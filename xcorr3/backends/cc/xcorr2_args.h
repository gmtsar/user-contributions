#ifndef __XCORR2_ARGS_H_INCLUDED__
#define __XCORR2_ARGS_H_INCLUDED__

struct st_xcorr {
    int m_nx, m_ny;
    int s_nx, s_ny;

    int x_offset, y_offset;
    int xsearch, ysearch;
    int nxl, nyl;
    double astretcha;

    int ri;
    int interp_factor;
    int n2x, n2y;  // high-res correlation window
    int nyquist_split;  // 0 = GMTSAR-compatible (no Nyquist split, default), 1 = split
    int snr_thr;        // correlation threshold for the integrated geocoding filter
    int psnr_thr;       // peak-significance SNR threshold for the integrated geocoding filter
    int do_geocode;     // 1 = auto post-process + geocode when trans.dat is present
    int do_blockmedian; // 1 = blockmedian + median despeckle grid; 0 = raw xyz2grd (manual-friendly)

    int throttle_ms;    // rest this many ms every throttle_every patches to leave GPU headroom (0 = full speed)
    int throttle_every; // patch batch size between GPU rests

    char *m_path;
    char *s_path;
};

enum enum_xcorr_af_device {
    XCORR2_DEVICE_DEFAULT = 0,
    XCORR2_DEVICE_CUDA,
    XCORR2_DEVICE_OPENCL,
    XCORR2_DEVICE_CPU
};

struct st_xcorr_args {
    const char *m_prm;
    const char *s_prm;

    int nx, ny;
    int xsearch, ysearch;
    int range_interp;
    int interp;

    bool noshift;
    bool nointerp;
    bool norange;
    bool nyquist_split;
    int snr;            // -snr  : correlation threshold for integrated geocoding (default 10)
    int psnr;           // -psnr : peak-significance SNR threshold (default 5)
    bool geocode;       // -geocode : enable the integrated post-processing + geocoding (needs trans.dat); off by default
    bool no_geocode;    // -no_geocode : kept for backward compatibility, forces geocoding off (now the default)
    bool no_blockmedian;// -no_blockmedian : grid with raw xyz2grd instead of blockmedian+median despeckle

    int throttle_ms;    // -throttle        : ms to rest each batch (default 2; 0 = full speed)
    int throttle_every; // -throttle_every  : positive patches per batch (default 64)

    enum enum_xcorr_af_device device;
};

void apply_args(const struct st_xcorr_args *args, struct st_xcorr *xc);
void parse_opts(struct st_xcorr_args *xa, int argc, char **argv);

#endif
