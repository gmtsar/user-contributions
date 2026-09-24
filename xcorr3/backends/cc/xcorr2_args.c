#include "xcorr2.h"
#include "xcorr2_args.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>

/* Parse syntax before conversion; strtol alone accepts empty/whitespace input. */
static int option_integer(const char *name, const char *text, int minimum, int maximum) {
    const char *p = text;
    if (p && (*p == '+' || *p == '-')) ++p;
    if (!p || !*p) goto invalid;
    for (; *p; ++p) if (*p < '0' || *p > '9') goto invalid;
    errno = 0;
    char *end;
    long value = strtol(text, &end, 10);
    if (errno == ERANGE || *end || value < minimum || value > maximum) goto invalid;
    return (int)value;
invalid:
    fprintf(stderr, "xcorr_cc: -%s requires an integer in [%d, %d]\n", name, minimum, maximum);
    exit(EXIT_FAILURE);
}

/* The CUDA kernels/CUB counts still use int. Check before any int product,
 * grid allocation, GPU initialization or output creation. */
static const char *configuration_error(const struct st_xcorr *xc) {
    if (xc->nxl <= 0 || xc->nxl > INT_MAX - 3) return "nx is outside the sampling index range";
    if (xc->nyl <= 0 || xc->nyl > INT_MAX - 1) return "ny is outside the sampling index range";
    if (!TEST_2PWR(xc->xsearch) || xc->xsearch > INT_MAX / 6) return "xsearch must be a power of two within the index range";
    if (!TEST_2PWR(xc->ysearch) || xc->ysearch > INT_MAX / 6) return "ysearch must be a power of two within the index range";
    if (!TEST_2PWR(xc->ri)) return "range_interp must be a positive power of two";
    if (xc->interp_factor <= 0 || xc->interp_factor > INT_MAX / 8) return "interp is outside the interpolation index range";
    if (xc->interp_factor > 1 && (xc->xsearch < 8 || xc->ysearch < 8))
        return "xsearch and ysearch must be at least 8 with subpixel interpolation; use -nointerp for smaller windows";
    int64_t width = 4 * (int64_t)xc->xsearch, height = 4 * (int64_t)xc->ysearch;
    int64_t hi_side = 8 * (int64_t)xc->interp_factor;
    if (width * height > INT_MAX) return "xsearch/ysearch patch element count exceeds INT_MAX";
    if (width * xc->ri > INT_MAX || width * height > INT_MAX / xc->ri)
        return "range_interp element count exceeds INT_MAX";
    if (hi_side * hi_side > INT_MAX) return "interp element count exceeds INT_MAX";
    if ((int64_t)xc->m_nx * height > INT_MAX || (int64_t)xc->s_nx * height > INT_MAX)
        return "num_rng_bins/ysearch row-cache index exceeds INT_MAX";
    if ((xc->m_nx - 6 * (int64_t)xc->xsearch) / (xc->nxl + 3) <= 0)
        return "nx/xsearch produce a nonpositive range sampling increment";
    if ((xc->m_ny - 6 * (int64_t)xc->ysearch) / (xc->nyl + 1) <= 0)
        return "ny/ysearch produce a nonpositive azimuth sampling increment";
    return NULL;
}

void apply_args(const struct st_xcorr_args *args, struct st_xcorr *xc) {
    int num_patches, num_valid_az;
    double prf[2];
    const char *m_path = NULL, *s_path = NULL, *failure = NULL;
    struct prm_handler m_prm = {0}, s_prm = {0};
    memset(xc, 0, sizeof(*xc));
    if (!prm_open(&m_prm, args->m_prm) || !prm_open(&s_prm, args->s_prm)) goto fail;
    if (!prm_get_str(&m_prm, "SLC_file", &m_path) || !prm_get_str(&s_prm, "SLC_file", &s_path)) goto fail;
    if (!prm_get_int(&m_prm, "num_rng_bins", &xc->m_nx) ||
        !prm_get_int(&m_prm, "num_patches", &num_patches) ||
        !prm_get_int(&m_prm, "num_valid_az", &num_valid_az)) goto fail;
    if (xc->m_nx <= 0 || num_patches <= 0 || num_valid_az <= 0 || num_patches > INT_MAX / num_valid_az) {
        failure = "master PRM num_rng_bins/num_patches/num_valid_az must be positive and the line count must fit int"; goto fail;
    }
    xc->m_ny = num_patches * num_valid_az;
    if (!prm_get_int(&s_prm, "num_rng_bins", &xc->s_nx) ||
        !prm_get_int(&s_prm, "num_patches", &num_patches) ||
        !prm_get_int(&s_prm, "num_valid_az", &num_valid_az)) goto fail;
    if (xc->s_nx <= 0 || num_patches <= 0 || num_valid_az <= 0 || num_patches > INT_MAX / num_valid_az) {
        failure = "secondary PRM num_rng_bins/num_patches/num_valid_az must be positive and the line count must fit int"; goto fail;
    }
    xc->s_ny = num_patches * num_valid_az;
    if (!prm_get_f64(&m_prm, "PRF", &prf[0]) || !prm_get_f64(&s_prm, "PRF", &prf[1])) goto fail;
    if (prf[0] <= 0 || prf[1] <= 0) { failure = "PRF must be positive in both PRMs"; goto fail; }
    xc->astretcha = (prf[1] - prf[0]) / prf[0];
    if (!isfinite(xc->astretcha) || 1.0 + xc->astretcha <= 0) {
        failure = "PRF ratio produces an invalid azimuth stretch"; goto fail;
    }

    if (!args->noshift) {
        if (!prm_get_int(&s_prm, "rshift", &xc->x_offset) || !prm_get_int(&s_prm, "ashift", &xc->y_offset)) goto fail;
    } else
        xc->x_offset = xc->y_offset = 0;

    xc->xsearch = args->xsearch;
    xc->ysearch = args->ysearch;

    xc->nxl = args->nx;
    xc->nyl = args->ny;

    if (args->norange)
        xc->ri = 1;
    else
        xc->ri = args->range_interp;

    if (args->nointerp)
        xc->interp_factor = 1;
    else
        xc->interp_factor = args->interp;

    xc->n2x = xc->n2y = 8;

    /* Nyquist handling in FFT sub-pixel interpolation.
     * Default 0 = GMTSAR-compatible (fft_arrange_interpolate, no split) so
     * results reproduce the reference xcorr; -nyquist_split enables the
     * symmetric split. */
    xc->nyquist_split = args->nyquist_split ? 1 : 0;

    /* Integrated post-processing is explicitly requested, off by default. */
    xc->snr_thr = args->snr;
    xc->psnr_thr = args->psnr;
    /* Geocoding is opt-in (-geocode) so that a stray trans.dat cannot trigger
     * side effects when xcorr_cc stands in for xcorr inside align scripts.
     * -no_geocode is kept for backward compatibility and still forces it off. */
    xc->do_geocode = (args->geocode && !args->no_geocode) ? 1 : 0;
    xc->do_blockmedian = args->no_blockmedian ? 0 : 1;

    /* GPU throttle: rest briefly every batch of patches so a single GPU that
     * also drives the display stays responsive and gets thermal/power headroom.
     * Defaults leave a small headroom (~a few % throughput cost). Disable with
     * -throttle 0 when the display runs on another GPU (e.g. the iGPU). */
    xc->throttle_ms = args->throttle_ms;
    xc->throttle_every = args->throttle_every;

    failure = configuration_error(xc);
    if (failure) goto fail;
    xc->m_path = strdup(m_path);
    xc->s_path = strdup(s_path);
    if (!xc->m_path || !xc->s_path) { failure = "cannot allocate SLC paths"; goto fail; }

    prm_close(&m_prm);
    prm_close(&s_prm);
    return;
fail:
    fprintf(stderr, "xcorr_cc: %s\n", failure ? failure : (m_prm.error[0] ? m_prm.error : s_prm.error));
    prm_close(&m_prm);
    prm_close(&s_prm);
    free(xc->m_path); free(xc->s_path);
    xc->m_path = xc->s_path = NULL;
    exit(EXIT_FAILURE);
}

void parse_opts(struct st_xcorr_args *xa, int argc, char **argv) {
    enum {
        OPT_NX = 10,
        OPT_NY = 20,
        OPT_RANGE_INTERP = 30,
        OPT_XSEARCH = 40,
        OPT_YSEARCH = 50,
        OPT_INTERP = 60,
        OPT_NO_SHIFT = -10,
        OPT_NOINTERP = -20,
        OPT_NORANGE = -30,
        OPT_DEVICE = -40,
        OPT_NYQUIST_SPLIT = -50,
        OPT_NO_GEOCODE = -60,
        OPT_NO_BLOCKMEDIAN = -70,
        OPT_GEOCODE = -80,
        OPT_SNR = 70,
        OPT_PSNR = 80,
        OPT_THROTTLE = 90,
        OPT_THROTTLE_EVERY = 100,
        OPT_HELP = -100,
    };

    static const char *help = \
        "xcorr_cc - CUDA 加速的二维振幅互相关 / 像素偏移追踪 (POT)\n\n"
        "用法: xcorr_cc 主影像.PRM 辅影像.PRM [选项...]\n"
        "  主影像.PRM / 辅影像.PRM   已配准的 GMTSAR SLC 参数文件 (主=参考, 辅=搜索)\n"
        "\n"
        "采样与搜索:\n"
        "  -nx  n                距离向(x)采样点数\n"
        "  -ny  n                方位向(y)采样点数        (总计算量 = nx*ny)\n"
        "  -xsearch xs           距离向搜索半径, 2的幂[32/64/128/256], 默认 64\n"
        "  -ysearch ys           方位向搜索半径, 2的幂, 默认 64\n"
        "  -noshift              忽略 PRM 的 rshift/ashift(置0); 已配准影像测残余形变时用\n"
        "插值:\n"
        "  -range_interp ri      距离向 FFT 插值倍数(2的幂), 默认 2;  -norange 关闭\n"
        "  -interp factor        相关面亚像素插值倍数, 默认 16;      -nointerp 关闭\n"
        "  -nyquist_split        亚像素插值对称拆分 Nyquist(默认关 = 与原版 GMTSAR 逐点一致)\n"
        "抗异常值 + 集成后处理 (默认关闭, 需显式开启):\n"
        "  -geocode              开启集成后处理: 过滤/网格化/地理编码 (需当前目录有 trans.dat)\n"
        "  -snr  n               地理编码的相关阈值, 默认 10   (0 = 不按相关过滤)\n"
        "  -psnr n               峰显著性阈值, 剔除杂峰/极大异常值, 默认 5  (0 = 不过滤)\n"
        "  -no_blockmedian       网格化用原始 xyz2grd(保留逐点值), 而非默认 blockmedian + 3x3 中值去毛刺\n"
        "  -no_geocode           兼容旧脚本保留; 现在默认即不做地理编码\n"
        "GPU 负载冗余 (单卡既算又显时防止桌面卡死 + 留散热/供电余量):\n"
        "  -throttle ms          每算一批点让 GPU 歇 ms 毫秒, 默认 2 (0 = 跑满速, 显示器接核显时用)\n"
        "  -throttle_every n     每 n 个采样点休息一次, 默认 64 (n 越小越跟手但越慢)\n"
        "\n"
        "输出:\n"
        "  freq_xcorr.dat   6 列: x_像素  距离偏移  y_像素  方位偏移  相关(0-100)  peak_snr(峰显著性)\n"
        "  azi_offset_ll.grd / rng_offset_ll.grd   (仅 -geocode 且当前目录有 trans.dat 时生成; 地理坐标, 单位 m)\n"
        "\n"
        "示例:\n"
        "  xcorr_cc master.PRM aligned.PRM -nx 100 -ny 200 -xsearch 64 -ysearch 64 -noshift -psnr 6\n";

    static struct option long_options[] = {
        { "noshift", no_argument, NULL, OPT_NO_SHIFT },
        { "nx", required_argument, NULL, OPT_NX },
        { "ny", required_argument, NULL, OPT_NY },
        { "nointerp", no_argument, NULL, OPT_NOINTERP },
        { "norange", no_argument, NULL, OPT_NORANGE },
        { "range_interp", required_argument, NULL, OPT_RANGE_INTERP },
        { "xsearch", required_argument, NULL, OPT_XSEARCH },
        { "ysearch", required_argument, NULL, OPT_YSEARCH },
        { "interp", required_argument, NULL, OPT_INTERP },
        { "af", required_argument, NULL, OPT_DEVICE },
        { "nyquist_split", no_argument, NULL, OPT_NYQUIST_SPLIT },
        { "snr", required_argument, NULL, OPT_SNR },
        { "psnr", required_argument, NULL, OPT_PSNR },
        { "throttle", required_argument, NULL, OPT_THROTTLE },
        { "throttle_every", required_argument, NULL, OPT_THROTTLE_EVERY },
        { "geocode", no_argument, NULL, OPT_GEOCODE },
        { "no_geocode", no_argument, NULL, OPT_NO_GEOCODE },
        { "no_blockmedian", no_argument, NULL, OPT_NO_BLOCKMEDIAN },
        { "help", no_argument, NULL, OPT_HELP },
        { 0, 0, 0, 0 },
    };

    if (argc == 1) {
        fputs(help, stdout);
        exit(0);
    }

    memset(xa, 0, sizeof(struct st_xcorr_args));
    xa->nx = 16; xa->ny = 32;
    xa->xsearch = xa->ysearch = 64;
    xa->range_interp = 2; xa->interp = 16;
    xa->snr = 10; xa->psnr = 5;
    xa->throttle_ms = 2; xa->throttle_every = 64;
    optind = 0; /* Reinitialize getopt for host regression harnesses too. */
    opterr = 0;

    while (1) {
        int opt, long_index = -1, int_arg = 0;
        opt = getopt_long_only(argc, argv, "", long_options, &long_index);

        if (opt == -1) break;

        if (opt == '?' || opt == ':') {
            fprintf(stderr, "xcorr_cc: unknown option or missing/unexpected value near '%s'\n", argv[optind > 0 ? optind - 1 : 0]);
            exit(EXIT_FAILURE);
        }
        if (opt > 0 && long_index >= 0) {
            int minimum = (opt == OPT_SNR || opt == OPT_PSNR || opt == OPT_THROTTLE) ? 0 : 1;
            int maximum = opt == OPT_SNR ? 100 : INT_MAX;
            int_arg = option_integer(long_options[long_index].name, optarg, minimum, maximum);
        }

        switch (opt) {
            case OPT_NX:
                xa->nx = int_arg;
                break;
            case OPT_NY:
                xa->ny = int_arg;
                break;
            case OPT_XSEARCH:
                xa->xsearch = int_arg;
                break;
            case OPT_YSEARCH:
                xa->ysearch = int_arg;
                break;
            case OPT_INTERP:
                xa->interp = int_arg;
                break;
            case OPT_DEVICE:
                if (!strcmp(optarg, "cuda"))
                    xa->device = XCORR2_DEVICE_CUDA;
                else if (!strcmp(optarg, "opencl"))
                    xa->device = XCORR2_DEVICE_OPENCL;
                else if (!strcmp(optarg, "cpu"))
                    xa->device = XCORR2_DEVICE_CPU;
                else {
                    fprintf(stderr, "xcorr_cc: -af must be cuda, opencl or cpu\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case OPT_RANGE_INTERP:
                xa->range_interp = int_arg;
                break;
            case OPT_NO_SHIFT:
                xa->noshift = true;
                break;
            case OPT_NOINTERP:
                xa->nointerp = true;
                break;
            case OPT_NORANGE:
                xa->norange = true;
                break;
            case OPT_NYQUIST_SPLIT:
                xa->nyquist_split = true;
                break;
            case OPT_SNR:
                xa->snr = int_arg;
                break;
            case OPT_PSNR:
                xa->psnr = int_arg;
                break;
            case OPT_THROTTLE:
                xa->throttle_ms = int_arg;
                break;
            case OPT_THROTTLE_EVERY:
                xa->throttle_every = int_arg;
                break;
            case OPT_GEOCODE:
                xa->geocode = true;
                break;
            case OPT_NO_GEOCODE:
                xa->no_geocode = true;
                break;
            case OPT_NO_BLOCKMEDIAN:
                xa->no_blockmedian = true;
                break;
            case OPT_HELP:
                fputs(help, stdout);
                exit(0);
            default:
                fprintf(stderr, "xcorr_cc: unrecognized option\n");
                exit(EXIT_FAILURE);
        }
    }

    if (optind < argc)
        xa->m_prm = argv[optind++];
    else {
        fprintf(stderr, "错误: 未指定主影像 PRM 文件。\n\n");
        fputs(help, stderr);
        exit(EXIT_FAILURE);
    }

    if (optind < argc)
        xa->s_prm = argv[optind++];
    else {
        fprintf(stderr, "错误: 未指定辅影像 PRM 文件。\n\n");
        fputs(help, stderr);
        exit(EXIT_FAILURE);
    }
    if (optind != argc) {
        fprintf(stderr, "xcorr_cc: expected exactly two PRM inputs; unexpected '%s'\n", argv[optind]);
        exit(EXIT_FAILURE);
    }
}
