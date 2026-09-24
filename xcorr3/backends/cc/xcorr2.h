#ifndef __XCORR2_H_INCLUDED__
#define __XCORR2_H_INCLUDED__

#include <stdio.h>
#include <stdbool.h>
#include <complex.h>
#include <pthread.h>
#include <glib.h>

#define TEST_2PWR(n) ((n) > 0 && ((n) & ((n) - 1)) == 0)

// array_helper.c
void c64_array_print(const char *fmt, _Complex double *arr, int n, int m);
_Complex double *c64_array_slice(const _Complex double *mat, int n_cols,
                                int tl_y, int s_y, int tl_x, int s_x);
double *f64_array_slice(const double *mat, int n_cols,
                        int tl_y, int s_y, int tl_x, int s_x);

void f64_array_stats(const double *array, int ny, int nx,
                     double *average, double *max,
                     int *argmax_y, int *argmax_x);

// fft_helper.c
_Complex double *dft_interpolate_2d(_Complex double *in, int height, int width,
                                   int scale_h, int scale_w,
                                   pthread_mutex_t *fftw_lock);
double *rdft_interpolate_2d(double *in, int height, int width,
                            int scale_h, int scale_w,
                            pthread_mutex_t *fftw_lock);

// prm_helper.c
struct prm_handler {
    GHashTable *entry;
    char *filename;
    char error[512];
};

/* Open a fresh handler. On failure all parser resources are released and
 * error contains a contextual diagnostic. Close is safe after failed open. */
bool prm_open(struct prm_handler *handler, const char *fname);
void prm_close(struct prm_handler *handler);
/* Required fields return false and set error if malformed, absent or empty.
 * Failure leaves outputs untouched. Strings remain valid until prm_close(). */
bool prm_get_str(struct prm_handler *handler, const char *key, const char **out);
bool prm_get_int(struct prm_handler *handler, const char *key, int *out);
bool prm_get_f64(struct prm_handler *handler, const char *key, double *out);
#endif
