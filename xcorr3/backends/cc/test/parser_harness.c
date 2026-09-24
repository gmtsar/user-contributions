/* Host-only input regression harness: no CUDA runtime or image allocations. */
#include "xcorr2.h"
#include "xcorr2_args.h"
#include <stdlib.h>
#include <string.h>

static void json_string(const char *value) {
    const unsigned char *p = (const unsigned char *)value;
    putchar('"');
    for (; *p; ++p) {
        if (*p == '"' || *p == '\\') printf("\\%c", *p);
        else if (*p < 32) printf("\\u%04x", *p);
        else putchar(*p);
    }
    putchar('"');
}

int main(int argc, char **argv) {
    struct st_xcorr_args args;
    struct st_xcorr xc = {0};
    parse_opts(&args, argc, argv);
    apply_args(&args, &xc);
    printf("{\"m_nx\":%d,\"m_ny\":%d,\"s_nx\":%d,\"s_ny\":%d,"
           "\"x_offset\":%d,\"y_offset\":%d,\"xsearch\":%d,\"ysearch\":%d,"
           "\"nxl\":%d,\"nyl\":%d,\"astretcha\":%.17g,\"ri\":%d,"
           "\"interp_factor\":%d,\"nyquist_split\":%d,\"snr_thr\":%d,"
           "\"psnr_thr\":%d,\"do_geocode\":%d,\"do_blockmedian\":%d,"
           "\"throttle_ms\":%d,\"throttle_every\":%d,\"device\":%d,\"m_path\":",
           xc.m_nx, xc.m_ny, xc.s_nx, xc.s_ny, xc.x_offset, xc.y_offset,
           xc.xsearch, xc.ysearch, xc.nxl, xc.nyl, xc.astretcha, xc.ri,
           xc.interp_factor, xc.nyquist_split, xc.snr_thr, xc.psnr_thr,
           xc.do_geocode, xc.do_blockmedian, xc.throttle_ms, xc.throttle_every,
           args.device);
    json_string(xc.m_path);
    fputs(",\"s_path\":", stdout);
    json_string(xc.s_path);
    puts("}");
    free(xc.m_path);
    free(xc.s_path);
    return 0;
}
