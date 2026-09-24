/*	$Id: xcorr.c 73 2013-04-19 17:59:45Z pwessel $	*/
/***************************************************************************/
/* xcorr does a 2-D cross correlation on complex or real images            */
/* either using a time convolution or wavenumber multiplication.           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * Creator:  Rob J. Mellors                                                *
 *           (San Diego State University)                                  *
 * Date   :  November 7, 2009                                              *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE				        	                           *
 *							                   *
 * 011810       Testing and very monor cosmetic modfifications DTS         *
 * 061520	Problem with sub-pixel interpolation RJM		   *
 * 		- fixed bug in 2D interpolationa			   *
 * 		- revised read_xcorr_data to read in all x position	   *
 * 		- reads directly into float rather than int	           *
 * 		- add range interpolation	                           *
 * 		- eliminated obsolete options and code	                   *
 * 		- renamed xcorr_utils.c print_results.c	                   *
 * 		- further testing....					   *
 ***************************************************************************/

/*-------------------------------------------------------*/
/* xcorr_mt: multi-process fork of the original GMTSAR xcorr.
 * Splits the (independent) correlation locations over N worker
 * processes and merges their output in the original location order.
 * Numerical compatibility is validated against the matching GMTSAR
 * build; byte identity is not guaranteed for every input or mode.
 * See VALIDATION_METHODS.md for inherited interpolation limitations
 * and the intentional correction of uninitialized -nointerp offsets.
 *
 * Derived from GMTSAR xcorr.c (R. J. Mellors et al., GMTSAR 6.6,
 * https://github.com/gmtsar/gmtsar), modified 2026-07 for
 * fork-based multi-process execution.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, version 3 of the
 * License (see LICENSE).                                        */
#define _POSIX_C_SOURCE 200809L
#include "gmtsar.h"
#include <errno.h>
#include <limits.h>
#include <signal.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

static void *checked_malloc(size_t size) {
    void *p = malloc(size);
    if (p == NULL) { die(" memory allocation failed\n", ""); exit(EXIT_FAILURE); }
    return p;
}

/* iloc half-open range [g_il0, g_il1) handled by this process;
 * the single-process path uses the full range [0, nlocs).      */
static int g_il0 = 0, g_il1 = 0;

char *USAGE = "xcorr [GMTSAR] - Compute 2-D cross-correlation of two images\n\n"
              "\nUsage: xcorr master.PRM aligned.PRM [-time] [-real] [-freq] [-nx n] [-ny "
              "n] [-xsearch xs] [-ysearch ys]\n"
              "master.PRM     	PRM file for reference image\n"
              "aligned.PRM     	 	PRM file of secondary image\n"
              "-time      		use time cross-correlation\n"
              "-freq      		use frequency cross-correlation (default)\n"
              "-real      		read float numbers instead of complex numbers\n"
              "-noshift  		ignore ashift and rshift in prm file (set to 0)\n"
              "-nx  nx    		number of locations in x (range) direction "
              "(int)\n"
              "-ny  ny    		number of locations in y (azimuth) direction "
              "(int)\n"
              "-nointerp     		do not interpolate correlation function\n"
              "-range_interp ri  	interpolate range by ri (power of two) [default: 2]\n"
              "-norange     		do not range interpolate \n"
              "-xsearch xs		search window size in x (range) direction (int "
              "power of 2 [32 64 128 256])\n"
              "-ysearch ys		search window size in y (azimuth) direction "
              "(int power of 2 [32 64 128 256])\n"
              "-interp  factor    	interpolate correlation function by factor "
              "(int) [default, 16]\n"
              "-v			verbose\n"
              "-nproc n  		number of worker processes (capped at online cores)\n"
              "           		[default: $OMP_NUM_THREADS, else cores-2]\n"
              "output: \n freq_xcorr.dat (default) \n time_xcorr.dat (if -time option))\n"
              "\nuse fitoffset.csh to convert output to PRM format\n"
              "\nExample:\n"
              "xcorr IMG-HH-ALPSRP075880660-H1.0__A.PRM "
              "IMG-HH-ALPSRP129560660-H1.0__A.PRM -nx 20 -ny 50 \n"
              "xcorr file1.grd file2.grd -nx 20 -ny 50 (takes grids with real numbers)\n";

/*-------------------------------------------------------------------------------*/
/* find "-nproc N" anywhere in argv, remove both tokens (so the original         */
/* parse_command_line never sees them) and return N, or -1 when absent           */
static int extract_nproc(int *argc, char **argv) {
	int i, pos = -1, nproc = -1;

	for (i = 1; i < *argc; i++) {
		if (!strcmp(argv[i], "-nproc")) {
			if (i + 1 >= *argc)
				die(" no option after -nproc!\n", "");
			char *end;
            long value;
            errno = 0;
            value = strtol(argv[i + 1], &end, 10);
            nproc = (!errno && *end == '\0' && value > 0 && value <= INT_MAX) ? (int)value : -1;
			pos = i;
			break;
		}
	}
	if (pos >= 0) {
		for (i = pos; i + 2 < *argc; i++)
			argv[i] = argv[i + 2];
		*argc -= 2;
		argv[*argc] = NULL;
		fprintf(stderr, " using %d worker processes\n", nproc);
	}
	return (nproc);
}
/*-------------------------------------------------------------------------------*/
/* precedence: -nproc N > OMP_NUM_THREADS > auto (online cores - 2, min 1);      */
/* invalid values fall through to the next level                                 */
static int resolve_nproc(int cli_nproc) {
	long n = -1, cores;
	char *e;

	cores = sysconf(_SC_NPROCESSORS_ONLN);
	if (cores < 1)
		cores = 1;

	if (cli_nproc > 0)
		n = cli_nproc;
	if (n < 1) {
		e = getenv("OMP_NUM_THREADS");
		if (e != NULL) {
            char *end;
            long value;
            errno = 0;
            value = strtol(e, &end, 10);
            if (!errno && *end == '\0' && value > 0 && value <= INT_MAX) n = value;
        }
	}
	if (n < 1) {
		/* default: leave 2 cores of headroom for the system */
		n = (cores > 2) ? cores - 2 : 1;
		fprintf(stderr, " auto: %ld online cores, keeping 2 in reserve -> %d worker processes\n", cores, (int)n);
	}
	return ((int)n);
}
/*-------------------------------------------------------------------------------*/
static double wall_seconds(struct timespec *t0, struct timespec *t1) {
	return ((double)(t1->tv_sec - t0->tv_sec) + (double)(t1->tv_nsec - t0->tv_nsec) / 1.0e9);
}
/*-------------------------------------------------------------------------------*/
int do_range_interpolate(void *API, struct FCOMPLEX *c, int nx, int ri, struct FCOMPLEX *work) {
	int i;

	/* interpolate c and put into work */
	fft_interpolate_1d(API, c, nx, work, ri);

	/* replace original with interpolated (only half) */
	for (i = 0; i < nx; i++) {
		c[i].r = work[i + nx / 2].r;
		c[i].i = work[i + nx / 2].i;
	}

	return (EXIT_SUCCESS);
}
/*-------------------------------------------------------------------------------*/
/* complex arrays used in fft correlation */
/* load complex arrays and mask out aligned */
/* c1 is master */
/* c2 is aligned */
/* c3 used in fft complex correlation */
/* c1, c2, and c3 are npy by npx */
/* d1, d2 are npy by nx (length of line in SLC) */
/*-------------------------------------------------------------------------------*/
void assign_values(void *API, struct xcorr *xc, int iloc) {
	int i, j, k, sx, mx;
	double mean1, mean2;

	/* master and aligned x offsets */
	mx = xc->loc[iloc].x - xc->npx / 2;
	sx = xc->loc[iloc].x + xc->x_offset - xc->npx / 2;

	for (i = 0; i < xc->npy; i++) {
		for (j = 0; j < xc->npx; j++) {
			k = i * xc->npx + j;

			xc->c3[k].i = xc->c3[k].r = 0.0f;

			xc->c1[k].r = xc->d1[i * xc->m_nx + mx + j].r;
			xc->c1[k].i = xc->d1[i * xc->m_nx + mx + j].i;

			xc->c2[k].r = xc->d2[i * xc->s_nx + sx + j].r;
			xc->c2[k].i = xc->d2[i * xc->s_nx + sx + j].i;
		}
	}

	/* range interpolate */
	if (xc->ri > 1) {
		for (i = 0; i < xc->npy; i++) {
			do_range_interpolate(API, &xc->c1[i * xc->npx], xc->npx, xc->ri, xc->ritmp);
			do_range_interpolate(API, &xc->c2[i * xc->npx], xc->npx, xc->ri, xc->ritmp);
		}
	}

	/* convert to amplitude and demean */
	mean1 = mean2 = 0.0;
	for (i = 0; i < xc->npy * xc->npx; i++) {
		xc->c1[i].r = Cabs(xc->c1[i]);
		xc->c1[i].i = 0.0f;

		xc->c2[i].r = Cabs(xc->c2[i]);
		xc->c2[i].i = 0.0f;

		mean1 += xc->c1[i].r;
		mean2 += xc->c2[i].r;
	}

	mean1 /= (double)(xc->npy * xc->npx);
	mean2 /= (double)(xc->npy * xc->npx);

	for (i = 0; i < xc->npy * xc->npx; i++) {
		xc->c1[i].r = xc->c1[i].r - (float)mean1;
		xc->c2[i].r = xc->c2[i].r - (float)mean2;
	}

	/* apply mask */
	for (i = 0; i < xc->npy * xc->npx; i++) {
		xc->c1[i].i = xc->c2[i].i = 0.0f;
		xc->c2[i].r = xc->c2[i].r * (float)xc->mask[i];

		xc->i1[i] = (int)(xc->c1[i].r);
		xc->i2[i] = (int)(xc->c2[i].r);
	}

	if (debug)
		fprintf(stderr, " mean %lf\n", mean1);
	if (debug)
		fprintf(stderr, " mean %lf\n", mean2);
}
/*-------------------------------------------------------------------------------*/
void do_correlation(void *API, struct xcorr *xc) {
	int i, j, iloc, istep;
	int row0, row1;

	/* correlation locations are independent: each worker process handles  */
	/* the iloc range [g_il0, g_il1) and writes its own ordered output     */
	istep = 1;

	/* allocate arrays   			*/
	allocate_arrays(xc);

	/* make mask 				*/
	make_mask(xc);

	iloc = 0;
	for (i = 0; i < xc->nyl; i += istep) {

		/* skip rows that hold no location of this process (no data read) */
		row0 = iloc;
		row1 = iloc + xc->nxl;
		if (row1 <= g_il0 || row0 >= g_il1) {
			iloc = row1;
			continue;
		}

		/* read in data for each row */
		read_xcorr_data(xc, iloc);

		for (j = 0; j < xc->nxl; j++) {

			if (iloc >= g_il0 && iloc < g_il1) {
				/* print_results always reads these, even with -nointerp. */
				xc->loc[iloc].xfrac = 0.0f;
				xc->loc[iloc].yfrac = 0.0f;

				if (debug)
					fprintf(stderr, " initial: iloc %d (%d,%d)\n", iloc, xc->loc[iloc].x, xc->loc[iloc].y);

				/* copy values from d1,d2 (real) to c1,c2 (complex) */
				assign_values(API, xc, iloc);

				if (debug)
					print_complex(xc->c1, xc->npy, xc->npx, 1);
				if (debug)
					print_complex(xc->c2, xc->npy, xc->npx, 1);

				/* correlate patch with data over offsets in time domain */
				if (xc->corr_flag < 2)
					do_time_corr(xc, iloc);

				/* correlate patch with data over offsets in freq domain */
				if (xc->corr_flag == 2)
					do_freq_corr(API, xc, iloc);

				/* oversample correlation surface  to obtain sub-pixel resolution */
				if (xc->interp_flag == 1)
					do_highres_corr(API, xc, iloc);

				/* write out results */
				print_results(xc, iloc);
			}

			iloc++;
		} /* end of x iloc loop */
	}     /* end of y iloc loop */
}
/*-------------------------------------------------------------------------------*/
/* want to avoid circular correlation so mask out most of b */
/* could adjust shape for different geometries */
/*-------------------------------------------------------------------------------*/
void make_mask(struct xcorr *xc) {
	int i, j, imask;
	imask = 0;

	for (i = 0; i < xc->npy; i++) {
		for (j = 0; j < xc->npx; j++) {
			xc->mask[i * xc->npx + j] = 1;
			if ((i < xc->ysearch) || (i >= (xc->npy - xc->ysearch))) {
				xc->mask[i * xc->npx + j] = imask;
			}
			if ((j < xc->xsearch) || (j >= (xc->npx - xc->xsearch))) {
				xc->mask[i * xc->npx + j] = imask;
			}
		}
	}
}
/*-------------------------------------------------------------------------------*/
void allocate_arrays(struct xcorr *xc) {
	int nx, ny, nx_exp, ny_exp;

	xc->d1 = (struct FCOMPLEX *)checked_malloc(xc->m_nx * xc->npy * sizeof(struct FCOMPLEX));
	xc->d2 = (struct FCOMPLEX *)checked_malloc(xc->s_nx * xc->npy * sizeof(struct FCOMPLEX));

	xc->i1 = (int *)checked_malloc(xc->npx * xc->npy * sizeof(int));
	xc->i2 = (int *)checked_malloc(xc->npx * xc->npy * sizeof(int));

	xc->c1 = (struct FCOMPLEX *)checked_malloc(xc->npx * xc->npy * sizeof(struct FCOMPLEX));
	xc->c2 = (struct FCOMPLEX *)checked_malloc(xc->npx * xc->npy * sizeof(struct FCOMPLEX));
	xc->c3 = (struct FCOMPLEX *)checked_malloc(xc->npx * xc->npy * sizeof(struct FCOMPLEX));

	xc->ritmp = (struct FCOMPLEX *)checked_malloc(xc->ri * xc->npx * sizeof(struct FCOMPLEX));
	xc->mask = (short *)checked_malloc(xc->npx * xc->npy * sizeof(short));

	/* this is size of correlation patch */
	xc->corr = (double *)checked_malloc(2 * xc->ri * (xc->nxc) * (xc->nyc) * sizeof(double));

	if (xc->interp_flag == 1) {
		nx = 2 * xc->n2x;
		ny = 2 * xc->n2y;
		nx_exp = nx * (xc->interp_factor);
		ny_exp = ny * (xc->interp_factor);
		xc->md = (struct FCOMPLEX *)checked_malloc(nx * ny * sizeof(struct FCOMPLEX));
		xc->cd_exp = (struct FCOMPLEX *)checked_malloc(nx_exp * ny_exp * sizeof(struct FCOMPLEX));
	}
}

/*-------------------------------------------------------*/
/* Private files prevent concurrent invocations from mixing worker results.
 * Publish only after every worker and every output write has succeeded. */
static char workdir[] = ".xcorr_mt-XXXXXX";
static int workdir_ready, cleanup_parts;

static void cleanup_output(void) {
	char path[128];
	int k;
	if (!workdir_ready) return;
	for (k = 0; k < cleanup_parts; k++) {
		snprintf(path, sizeof(path), "%s/part.%d", workdir, k);
		unlink(path);
	}
	snprintf(path, sizeof(path), "%s/result", workdir);
	unlink(path);
	rmdir(workdir);
}

static int close_output(FILE *f) {
	int bad = ferror(f);
	if (fclose(f) != 0) bad = 1;
	return bad;
}

int main(int argc, char **argv) {
	int input_flag = 0, nfiles = 2, cli_nproc, nproc, k, spawned = 0, failed = 0;
	struct xcorr *xc = checked_malloc(sizeof(*xc));
	struct timespec ts0, ts1;
	void *API;
	pid_t *pids;
	char result[128];

	verbose = debug = 0;
	/* Turn file-size-limit faults into checked stream errors, so cleanup runs. */
	if (signal(SIGXFSZ, SIG_IGN) == SIG_ERR) die(" cannot configure output error handling\n", "");
	xc->interp_flag = 0;
	xc->corr_flag = 2;
	xc->offset_flag = 0;
	cli_nproc = extract_nproc(&argc, argv);
	nproc = resolve_nproc(cli_nproc);
	if (argc < 3) die(USAGE, "");
	API = GMT_Create_Session(argv[0], 0U, 0U, NULL);
	if (API == NULL) return EXIT_FAILURE;
	set_defaults(xc);
	parse_command_line(argc, argv, xc, &nfiles, &input_flag, USAGE);
	if (input_flag == 0) handle_prm(API, argv, xc, nfiles);
	if (debug) print_params(xc);
	if (xc->corr_flag == 0) strcpy(xc->filename, "time_xcorr.dat");
	if (xc->corr_flag == 1) strcpy(xc->filename, "time_xcorr_Gatelli.dat");
	if (xc->corr_flag == 2) strcpy(xc->filename, "freq_xcorr.dat");
	if (xc->nxl < 1 || xc->nyl < 1 || xc->nxl > INT_MAX / xc->nyl)
		die(" invalid measurement grid\n", "");
	get_locations(xc);
	{
		long cores = sysconf(_SC_NPROCESSORS_ONLN);
		if (cores < 1) cores = 1;
		if (nproc > cores) nproc = (int)cores;
	}
	if (nproc > xc->nlocs) nproc = xc->nlocs;
	if (nproc < 1) die(" empty measurement grid\n", "");
	if (mkdtemp(workdir) == NULL) die(" cannot create private output directory: ", strerror(errno));
	workdir_ready = 1;
	cleanup_parts = nproc;
	if (atexit(cleanup_output) != 0) {
		cleanup_output();
		die(" cannot register output cleanup\n", "");
	}
	snprintf(result, sizeof(result), "%s/result", workdir);
	fprintf(stderr, " using %d worker processes\n", nproc);
	clock_gettime(CLOCK_MONOTONIC, &ts0);

	if (nproc == 1) {
		g_il0 = 0;
		g_il1 = xc->nlocs;
		xc->file = fopen(result, "w");
		if (xc->file == NULL) die(" cannot open output: ", strerror(errno));
		do_correlation(API, xc);
		if (close_output(xc->file)) die(" output write failed\n", "");
	} else {
		pids = checked_malloc((size_t)nproc * sizeof(*pids));
		fflush(NULL);
		for (k = 0; k < nproc; k++) {
			pid_t pid = fork();
			if (pid < 0) { failed = 1; break; }
			if (pid == 0) {
				char part[128];
				/* Children must never clean up another worker's files. */
				workdir_ready = 0;
				g_il0 = (int)((long)k * xc->nlocs / nproc);
				g_il1 = (int)((long)(k + 1) * xc->nlocs / nproc);
				if (xc->format == 0 || xc->format == 1) {
					fclose(xc->data1);
					fclose(xc->data2);
					xc->data1 = fopen(xc->data1_name, "r");
					xc->data2 = fopen(xc->data2_name, "r");
					if (!xc->data1 || !xc->data2) _exit(EXIT_FAILURE);
				}
				snprintf(part, sizeof(part), "%s/part.%d", workdir, k);
				xc->file = fopen(part, "w");
				if (!xc->file) _exit(EXIT_FAILURE);
				do_correlation(API, xc);
				if (close_output(xc->file)) _exit(EXIT_FAILURE);
				_exit(EXIT_SUCCESS);
			}
			pids[spawned++] = pid;
		}
		for (k = 0; k < spawned; k++) {
			int status;
			pid_t got;
			do { got = waitpid(pids[k], &status, 0); } while (got < 0 && errno == EINTR);
			if (got < 0 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) failed = 1;
		}
		free(pids);
		if (failed) die(" worker failed; existing output left unchanged\n", "");
		{
			FILE *out = fopen(result, "w");
			char buf[65536];
			if (!out) die(" cannot open merged output\n", "");
			for (k = 0; k < nproc; k++) {
				char part[128];
				FILE *in;
				size_t nr;
				snprintf(part, sizeof(part), "%s/part.%d", workdir, k);
				in = fopen(part, "r");
				if (!in) { fclose(out); die(" cannot open worker output\n", ""); }
				while ((nr = fread(buf, 1, sizeof(buf), in)) != 0) {
					if (fwrite(buf, 1, nr, out) != nr) { failed = 1; break; }
				}
				if (ferror(in)) failed = 1;
				if (fclose(in) != 0) failed = 1;
				if (failed) break;
			}
			if (close_output(out)) failed = 1;
			if (failed) die(" output merge failed\n", "");
		}
	}
	if (xc->format == 0 || xc->format == 1) { fclose(xc->data1); fclose(xc->data2); }
	if (GMT_Destroy_Session(API)) return EXIT_FAILURE;
	if (rename(result, xc->filename) != 0) die(" cannot publish output: ", strerror(errno));
	clock_gettime(CLOCK_MONOTONIC, &ts1);
	fprintf(stdout, " elapsed time: %lf \n", wall_seconds(&ts0, &ts1));
	return EXIT_SUCCESS;
}
