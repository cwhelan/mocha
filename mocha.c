/* The MIT License

   Copyright (C) 2015-2020 Giulio Genovese

   Author: Giulio Genovese <giulio.genovese@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <htslib/kseq.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/ksort.h>
#include "regidx.h"
#include "kmin.h"
#include "mocha.h"
#include "genome_rules.h"
#include "beta_binom.h"
#include "bcftools.h"

#define MOCHA_VERSION "2020-04-08"

/****************************************
 * CONSTANT DEFINITIONS                 *
 ****************************************/

// TODO replace SIGN with copysignf()
#define SIGN(x) (((x) > 0) - ((x) < 0))

#define BDEV_LRR_BAF_DFLT "-2.0,-4.0,-6.0,10.0,6.0,4.0"
#define BDEV_BAF_PHASE_DFLT "6.0,8.0,10.0,15.0,20.0,30.0,50.0,80.0,100.0,150.0,200.0"
#define MIN_DST_DFLT "400"
#define MEDIAN_BAF_ADJ_DFLT "5"
#define MEDIAN_LRR_ADJ_DFLT "5"
#define ORDER_LRR_GC_DFLT "2"
#define MAX_ORDER 5
#define XY_PRB_DFLT "1e-09"
#define ERR_PRB_DFLT "1e-04"
#define FLIP_PRB_DFLT "1e-02"
#define TEL_PRB_DFLT "1e-02"
#define CEN_PRB_DFLT "1e-04"
#define SHORT_ARM_CHRS_DFLT "13,14,15,21,22,chr13,chr14,chr15,chr21,chr22"
#define LRR_BIAS_DFLT "0.2"
#define LRR_HAP2DIP_DFLT                                                                       \
	"0.45" // https://www.illumina.com/documents/products/technotes/technote_cnv_algorithms.pdf

#define FLT_INCLUDE (1 << 0)
#define FLT_EXCLUDE (1 << 1)
#define WGS_DATA (1 << 2)
#define NO_LOG (1 << 3)
#define NO_ANNOT (1 << 4)
#define USE_SHORT_ARMS (1 << 5)
#define USE_CENTROMERES (1 << 6)
#define NO_BAF_FLIP (1 << 7)

#define LRR 0
#define BAF 1
#define AD0 0
#define AD1 1
#define LDEV 0
#define BDEV 1

#define LRR_BAF 0
#define BAF_PHASE 1

#define SEX_UNK 0
#define SEX_MAL 1
#define SEX_FEM 2

#define MOCHA_UNK 0
#define MOCHA_DEL 1
#define MOCHA_DUP 2
#define MOCHA_UPD 3
#define MOCHA_CNP_DEL 4
#define MOCHA_CNP_DUP 5
#define MOCHA_CNP_CNV 6

#define MOCHA_NOT 0
#define MOCHA_ARM 1
#define MOCHA_TEL 2

#define GT_NC 0
#define GT_AA 1
#define GT_AB 2
#define GT_BB 3

/****************************************
 * DATA STRUCTURES                      *
 ****************************************/

typedef struct {
	int pos;
	int allele_a;
	int allele_b;
	float adjust[4];
} locus_t;

typedef struct {
	float xy_log_prb;
	float err_log_prb;
	float flip_log_prb;
	float tel_log_prb;
	float cen_log_prb;
	// TODO make this an array
	float *bdev_lrr_baf, *bdev_baf_phase;
	int bdev_lrr_baf_n, bdev_baf_phase_n;
	int min_dst;
	float lrr_cutoff;
	float lrr_hap2dip;
	float lrr_auto2sex;
	float lrr_bias;
	int median_baf_adj;
	int median_lrr_adj;
	int order_lrr_gc;
	int flags;
	genome_rules_t *genome_rules;
	regidx_t *cnp_idx;
	regitr_t *cnp_itr;

	int rid;
	int n;
	locus_t *locus_arr;
	int m_locus;
	float *gc_arr;
	int m_gc;
	int n_flipped;
} model_t;

typedef struct {
	int sample_idx;
	int sex;
	int rid;
	int beg_pos;
	int end_pos;
	int length;
	int8_t p_arm;
	int8_t q_arm;
	int nsites;
	int nhets;
	int n50_hets;
	float bdev;
	float bdev_se;
	float ldev;
	float ldev_se;
	float lod_lrr_baf;
	float lod_baf_phase;
	int nflips;
	float baf_conc;
	float lod_baf_conc;
	int8_t type;
	float cf;
} mocha_t;

typedef struct {
	int n, m;
	mocha_t *a;
} mocha_table_t;

typedef struct {
	float lrr_median;
	float lrr_sd;
	float lrr_auto;
	float dispersion; // either rho(AD0, AD1) for WGS model or sd(BAF)
	float baf_conc;
	float baf_auto;
	float coeffs[MAX_ORDER + 1];
	float rel_ess;
} stats_t;

typedef struct {
	int idx;
	int sex;
	float adjlrr_sd;
	int nsites;
	int nhets;
	int x_nonpar_nhets;
	float x_nonpar_dispersion; // either rho(AD0, AD1) for WGS model or sd(BAF)
	float x_nonpar_lrr_median;
	float y_nonpar_lrr_median;
	float mt_lrr_median;
	stats_t stats;
	stats_t *stats_arr;
	int m_stats, n_stats;

	int n;
	int *vcf_imap_arr;
	int m_vcf_imap;
	int16_t *data_arr[2];
	int m_data[2];
	int8_t *phase_arr;
	int m_phase;
} sample_t;

/****************************************
 * INLINE FUNCTIONS AND CONSTANTS       *
 ****************************************/

// this macro from ksort.h defines the function
// void ks_introsort_int(size_t n, int a[]);
KSORT_INIT_GENERIC(int)

static inline float sqf(float x)
{
	return x * x;
}
static inline double sq(double x)
{
	return x * x;
}
// the x == y is necessary in case x == -INFINITY
static inline float log_mean_expf(float x, float y)
{
	return x == y ? x
		      : (x > y ? x + logf(1.0f + expf(y - x)) : y + logf(1.0f + expf(x - y)))
				- (float)M_LN2;
}

static beta_binom_t *beta_binom_null, *beta_binom_alt;
static inline float beta_binom_log_lkl(const beta_binom_t *self, int16_t ad0, int16_t ad1)
{
	return ad0 == bcf_int16_missing || ad1 == bcf_int16_missing
		       ? 0.0f
		       : beta_binom_log_unsafe(self, ad0, ad1);
}

/****************************************
 * CONVERT FLOAT TO INT16 AND VICEVERSA *
 ****************************************/

#define INT16_SCALE 1000 // BAF values from Illumina are scaled to 1000

static inline int16_t float_to_int16(float in)
{
	return isnan(in) ? bcf_int16_missing : (int16_t)roundf(INT16_SCALE * in);
}

static inline float int16_to_float(int16_t in)
{
	return in == bcf_int16_missing ? NAN : ((float)in) / INT16_SCALE;
}

/******************************************
 * LRR AND COVERAGE POLYNOMIAL REGRESSION *
 ******************************************/

// the following alternative code snippets were considered to perform GC regression:
// https://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
// https://github.com/natedomin/polyfit/blob/master/polyfit.c
// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/sgels_ex.c.htm
// this function needs to use doubles internally when dealing with WGS data
static int polyfit(const float *lrr, const float *gc, int n, const int *imap, int order,
		   float *coeffs)
{
	int m = order + 1;
	if (n < m || order > MAX_ORDER)
		return -1;
	double B[MAX_ORDER + 1] = {0.0};
	double P[((MAX_ORDER + 1) * 2) + 1] = {0.0};
	double A[(MAX_ORDER + 1) * 2 * (MAX_ORDER + 1)] = {0.0};

	// identify the column vector
	for (int i = 0; i < n; i++) {
		float x = imap ? gc[imap[i]] : gc[i];
		float y = lrr[i];
		if (isnan(x) || isnan(y))
			continue;
		float powx = 1.0f;

		for (int j = 0; j < m; j++) {
			B[j] += (double)(y * powx);
			powx *= x;
		}
	}

	// initialize the PowX array
	P[0] = (float)n;

	// compute the sum of the powers of X
	for (int i = 0; i < n; i++) {
		float x = imap ? gc[imap[i]] : gc[i];
		if (isnan(x))
			continue;
		float powx = x;

		for (int j = 1; j < ((2 * m) + 1); j++) {
			P[j] += (double)powx;
			powx *= x;
		}
	}

	// initialize the reduction matrix
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			A[(i * (2 * m)) + j] = P[i + j];
		}

		A[(i * (2 * m)) + (i + m)] = 1.0;
	}

	// move the identity matrix portion of the redux matrix to the
	// left side (find the inverse of the left side of the redux matrix)
	for (int i = 0; i < m; i++) {
		double x = A[(i * (2 * m)) + i];
		if (x != 0) {
			for (int k = 0; k < (2 * m); k++) {
				A[(i * (2 * m)) + k] /= x;
			}

			for (int j = 0; j < m; j++) {
				if (i != j) {
					double y = A[(j * (2 * m)) + i];
					for (int k = 0; k < (2 * m); k++) {
						A[(j * (2 * m)) + k] -=
							y * A[(i * (2 * m)) + k];
					}
				}
			}
		} else {
			// cannot work with singular matrices
			return -1;
		}
	}

	// calculate coefficients
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			double x = 0.0;
			for (int k = 0; k < m; k++) {
				x += (A[(i * (2 * m)) + (k + m)] * B[k]);
			}
			coeffs[i] = (float)x;
		}
	}

	return 0;
}

static void ad_to_lrr_baf(const int16_t *ad0, const int16_t *ad1, float *lrr, float *baf, int n)
{
	// this function keeps a list of logarithms of integers to minimize log calls
	static float *logf_arr = NULL;
	static int n_logf = 0, m_logf = 0;
	if (ad0 == NULL && ad1 == NULL && n == 0) {
		free(logf_arr);
		return;
	}

	for (int i = 0; i < n; i++) {
		if (ad0[i] == bcf_int16_missing && ad1[i] == bcf_int16_missing) {
			lrr[i] = NAN;
			baf[i] = NAN;
			continue;
		}
		int cov = (int)(ad0[i] == bcf_int16_missing ? 0 : ad0[i])
			  + (int)(ad1[i] == bcf_int16_missing ? 0 : ad1[i]);
		if (cov == 0) {
			lrr[i] = 0;
			baf[i] = NAN;
		} else {
			if (cov > n_logf) {
				hts_expand(float, cov, m_logf, logf_arr);
				for (int j = n_logf; j < cov; j++)
					logf_arr[j] = logf(j + 1);
				n_logf = cov;
			}
			lrr[i] = logf_arr[cov - 1];
			baf[i] = (ad0[i] == bcf_int16_missing || ad1[i] == bcf_int16_missing)
					 ? NAN
					 : (float)ad1[i] / (float)cov;
		}
	}
}

static void adjust_lrr(float *lrr, const float *gc, int n, const int *imap, const float *coeffs,
		       int order)
{
	for (int i = 0; i < n; i++) {
		float x = imap ? gc[imap[i]] : gc[i];
		float powx = 1.0f;
		for (int j = 0; j <= order; j++) {
			lrr[i] -= coeffs[j] * powx;
			powx *= x;
		}
	}
}

// computes total sum of squares
// this function needs to use doubles internally when dealing with WGS data
static float get_tss(const float *v, int n)
{
	double mean = 0.0;
	int j = 0;
	for (int i = 0; i < n; i++) {
		if (!isnan(v[i])) {
			mean += (double)v[i];
			j++;
		}
	}
	if (j <= 1)
		return NAN;
	mean /= (double)j;

	double tss = 0.0;
	for (int i = 0; i < n; i++) {
		if (!isnan(v[i]))
			tss += sq((double)v[i] - mean);
	}
	return (float)tss;
}

/*********************************
 * HMM AND OPTIMIZATION METHODS  *
 *********************************/

// compute Viterbi path from probabilities
static int8_t *retrace_viterbi(int T, int N, const float *log_prb, const int8_t *ptr)
{
	int i, t;
	int8_t *path = (int8_t *)malloc(T * sizeof(int8_t));

	// initialize last path state
	path[T - 1] = 0;
	for (i = 1; i < N; i++)
		if (log_prb[(int)path[T - 1]] < log_prb[i])
			path[T - 1] = (int8_t)i;

	// compute best path by tracing back the Markov chain
	for (t = T - 1; t > 0; t--)
		path[t - 1] = ptr[(t - 1) * N + (int)path[t]];

	return path;
}

// rescale Viterbi log probabilities to avoid underflow issues
static inline void rescale_log_prb(float *log_prb, int n)
{
	float max = -INFINITY;
	for (int i = 0; i < n; i++)
		max = max > log_prb[i] ? max : log_prb[i];
	for (int i = 0; i < n; i++)
		log_prb[i] -= max;
}

// compute the Viterbi path from BAF
// n is the length of the hidden Markov model
// m is the number of possible BAF deviations
static int8_t *log_viterbi_run(const float *emis_log_lkl, int T, int m, float xy_log_prb,
			       float flip_log_prb, float tel_log_prb, float cen_log_prb,
			       int last_p, int first_q)
{
	int t, i, j, changeidx;

	// determine the number of hidden states based on whether phase information is used
	int N = 1 + m + (isnan(flip_log_prb) ? 0 : m);

	// allocate memory necessary for running the algorithm
	float *log_prb = (float *)malloc(N * sizeof(float));
	float *new_log_prb = (float *)malloc(N * sizeof(float));
	int8_t *ptr = (int8_t *)malloc(N * (T - 1) * sizeof(int8_t));
	int8_t *path;

	// initialize and rescale the first state
	log_prb[0] = emis_log_lkl[0];
	for (i = 1; i < N; i++)
		log_prb[i] = xy_log_prb - (last_p == 0 ? cen_log_prb * 0.5f : tel_log_prb)
			     + emis_log_lkl[i];
	rescale_log_prb(log_prb, N);

	// compute best probabilities at each position
	for (t = 1; t < T; t++) {
		// this causes a penalty for mosaic chromosomal calls across the centromeres
		float exit_log_prb = t > last_p ? xy_log_prb + cen_log_prb * 0.5
						: xy_log_prb - cen_log_prb * 0.5;
		float enter_log_prb = t < first_q ? xy_log_prb + cen_log_prb * 0.5
						  : xy_log_prb - cen_log_prb * 0.5;

		for (i = 0; i < N; i++) {
			new_log_prb[i] = log_prb[i];
			ptr[(t - 1) * N + i] = (int8_t)i;
		}

		// compute whether a state switch should be considered for null state
		for (i = 1; i < N; i++) {
			if (new_log_prb[0] < log_prb[i] + exit_log_prb) {
				new_log_prb[0] = log_prb[i] + exit_log_prb;
				ptr[(t - 1) * N] = ptr[(t - 1) * N + i];
			}
			if (new_log_prb[i] < log_prb[0] + enter_log_prb) {
				new_log_prb[i] = log_prb[0] + enter_log_prb;
				ptr[(t - 1) * N + i] = ptr[(t - 1) * N];
			}
		}

		// compute whether a state switch should be considered for each other state
		// it will run twice if and only if phasing is used
		for (j = 0; j == 0 || (!isnan(flip_log_prb) && j == m); j += m) {
			float change_log_prb = log_prb[0] + enter_log_prb;
			changeidx = 0;
			for (i = 0; i < m; i++) {
				if (change_log_prb < log_prb[1 + j + i] + xy_log_prb * 1.5f) {
					change_log_prb = log_prb[1 + j + i] + xy_log_prb * 1.5f;
					changeidx = 1 + j + i;
				}
			}
			for (i = 0; i < m; i++) {
				if (new_log_prb[1 + j + i] < change_log_prb) {
					new_log_prb[1 + j + i] = change_log_prb;
					ptr[(t - 1) * N + 1 + j + i] =
						ptr[(t - 1) * N + changeidx];
				}
			}
		}

		// compute whether a phase flip should be considered for non-null states
		if (!isnan(flip_log_prb)) {
			for (i = 0; i < m; i++) {
				if (new_log_prb[1 + i]
				    < new_log_prb[1 + m + i] + flip_log_prb) {
					new_log_prb[1 + i] =
						new_log_prb[1 + m + i] + flip_log_prb;
					ptr[(t - 1) * N + 1 + i] = ptr[(t - 1) * N + 1 + m + i];
				}
				if (new_log_prb[1 + m + i]
				    < new_log_prb[1 + i] + flip_log_prb) {
					new_log_prb[1 + m + i] =
						new_log_prb[1 + i] + flip_log_prb;
					ptr[(t - 1) * N + 1 + m + i] = ptr[(t - 1) * N + 1 + i];
				}
			}
		}

		// update and rescale the current state
		new_log_prb[0] += emis_log_lkl[t * N];
		for (i = 0; i < m; i++) {
			new_log_prb[1 + i] += emis_log_lkl[t * N + 1 + i];
			if (!isnan(flip_log_prb))
				new_log_prb[1 + m + i] += emis_log_lkl[t * N + 1 + m + i];
		}
		for (i = 0; i < N; i++)
			log_prb[i] = new_log_prb[i];
		rescale_log_prb(log_prb, N);
	}

	// add closing cost to the last state
	for (i = 1; i < N; i++)
		log_prb[i] += xy_log_prb - (first_q == T ? cen_log_prb * 0.5 : tel_log_prb);
	rescale_log_prb(log_prb, N);

	path = retrace_viterbi(T, N, log_prb, ptr);

	// free memory
	free(log_prb);
	free(new_log_prb);
	free(ptr);

	// symmetrize the path
	if (!isnan(flip_log_prb))
		for (i = 0; i < T; i++)
			if (path[i] > m)
				path[i] = (int8_t)m - path[i];

	return path;
}

/*********************************
 * LRR AND BAF LIKELIHOODS       *
 *********************************/

// rescale emission probabilities to avoid problems with outliers
// TODO find a better name for this function
static void rescale_emis_log_lkl(float *log_prb, int n, float err_log_prb)
{
	float min_thr = -INFINITY;
	for (int i = 0; i < n; i++)
		if (min_thr < log_prb[i])
			min_thr = log_prb[i];
	min_thr += err_log_prb;
	for (int i = 0; i < n; i++)
		if (log_prb[i] < min_thr)
			log_prb[i] = min_thr;
}

static inline float norm_log_lkl(float x, float m, float s, float w)
{
	return isnan(x) ? 0.0f : -0.5f * sqf((x - m) / s) * w;
}

// lrr_bias is used in a different way from what done by Petr Danecek in bcftools/vcfcnv.c
static inline float lrr_baf_log_lkl(float lrr, float baf, float ldev, float bdev, float lrr_sd,
				    float baf_sd, float lrr_bias)
{
	return norm_log_lkl(lrr, ldev, lrr_sd, lrr_bias)
	       + log_mean_expf(norm_log_lkl(baf - 0.5f, bdev, baf_sd, 1.0f),
			       norm_log_lkl(baf - 0.5f, -bdev, baf_sd, 1.0f));
}

static inline float baf_phase_log_lkl(float baf, int8_t phase, float bdev, float baf_sd)
{
	return phase == 0 ? log_mean_expf(norm_log_lkl(baf - 0.5f, bdev, baf_sd, 1.0f),
					  norm_log_lkl(baf - 0.5f, -bdev, baf_sd, 1.0f))
			  : norm_log_lkl(baf - 0.5f, (float)SIGN(phase) * bdev, baf_sd, 1.0f);
}

// precomupute emission probabilities
static float *lrr_baf_emis_log_lkl(const float *lrr, const float *baf, int T, const int *imap,
				   float err_log_prb, float lrr_bias, float lrr_hap2dip,
				   float lrr_sd, float baf_sd, const float *bdev_lrr_baf_arr,
				   int m)
{
	float *ldev = (float *)malloc(m * sizeof(float));
	for (int i = 0; i < m; i++)
		ldev[i] = -logf(1.0f - 2.0f * bdev_lrr_baf_arr[i]) / (float)M_LN2 * lrr_hap2dip;
	int N = 1 + 2 * m;
	float *emis_log_lkl = (float *)malloc(N * T * sizeof(float));
	for (int t = 0; t < T; t++) {
		float x = imap ? lrr[imap[t]] : lrr[t];
		float y = imap ? baf[imap[t]] : baf[t];
		emis_log_lkl[t * N] =
			lrr_baf_log_lkl(x, y, 0.0f, 0.0f, lrr_sd, baf_sd, lrr_bias);
		for (int i = 0; i < m; i++) {
			emis_log_lkl[t * N + 1 + i] = lrr_baf_log_lkl(
				x, y, ldev[i], bdev_lrr_baf_arr[i], lrr_sd, baf_sd, lrr_bias);
		}
		// add states to distinguish LRR waves from true mosaic duplications/deletions
		for (int i = 0; i < m; i++) {
			if (bdev_lrr_baf_arr[i] < -1.0f / 6.0f
			    || bdev_lrr_baf_arr[i] >= 1.0f / 6.0f)
				emis_log_lkl[t * N + 1 + m + i] =
					emis_log_lkl[t * N] + err_log_prb;
			else
				emis_log_lkl[t * N + 1 + m + i] =
					lrr_baf_log_lkl(x, y, bdev_lrr_baf_arr[i], 0.0f, lrr_sd,
							baf_sd, lrr_bias);
		}
		rescale_emis_log_lkl(&emis_log_lkl[t * N], N, err_log_prb);
	}
	free(ldev);
	return emis_log_lkl;
}

// precomupute emission probabilities
static float *baf_phase_emis_log_lkl(const float *baf, const int8_t *gt_phase, int T,
				     const int *imap, float err_log_prb, float baf_sd,
				     const float *bdev, int m)
{
	int N = 1 + 2 * m;
	float *emis_log_lkl = (float *)malloc(N * T * sizeof(float));
	for (int t = 0; t < T; t++) {
		float x = imap ? baf[imap[t]] : baf[t];
		int8_t p = imap ? gt_phase[imap[t]] : gt_phase[t];
		emis_log_lkl[t * N] = baf_phase_log_lkl(x, (int8_t)1, 0.0f, baf_sd);
		for (int i = 0; i < m; i++) {
			emis_log_lkl[t * N + 1 + i] = baf_phase_log_lkl(x, p, bdev[i], baf_sd);
			if (p == 0)
				emis_log_lkl[t * N + 1 + m + i] = emis_log_lkl[t * N + 1 + i];
			else
				emis_log_lkl[t * N + 1 + m + i] =
					baf_phase_log_lkl(x, p, -bdev[i], baf_sd);
		}
		rescale_emis_log_lkl(&emis_log_lkl[t * N], N, err_log_prb);
	}
	return emis_log_lkl;
}

static int cnp_edge_is_not_cn2_lrr_baf(const float *lrr, const float *baf, int n, int a, int b,
				       float xy_log_prb, float err_log_prb, float lrr_bias,
				       float lrr_hap2dip, float lrr_sd, float baf_sd,
				       float ldev, float bdev)
{
	// test left edge
	float sum_log_lkl = 0.0f;
	for (int i = a - 1; i >= 0; i--) {
		float log_lkl =
			lrr_baf_log_lkl(lrr[i], baf[i], ldev, bdev, lrr_sd, baf_sd, lrr_bias)
			- lrr_baf_log_lkl(lrr[i], baf[i], 0.0f, 0.0f, lrr_sd, baf_sd, lrr_bias);
		if (log_lkl < err_log_prb)
			log_lkl = err_log_prb;
		else if (log_lkl > -err_log_prb)
			log_lkl = -err_log_prb;
		sum_log_lkl += log_lkl;
		if (sum_log_lkl > -xy_log_prb)
			return -1;
		if (sum_log_lkl < xy_log_prb)
			break;
	}

	// test right edge
	sum_log_lkl = 0.0f;
	for (int i = b + 1; i < n; i++) {
		float log_lkl =
			lrr_baf_log_lkl(lrr[i], baf[i], ldev, bdev, lrr_sd, baf_sd, lrr_bias)
			- lrr_baf_log_lkl(lrr[i], baf[i], 0, 0, lrr_sd, baf_sd, lrr_bias);
		if (log_lkl < err_log_prb)
			log_lkl = err_log_prb;
		else if (log_lkl > -err_log_prb)
			log_lkl = -err_log_prb;
		sum_log_lkl += log_lkl;
		if (sum_log_lkl > -xy_log_prb)
			return -1;
		if (sum_log_lkl < xy_log_prb)
			break;
	}

	return 0;
}

// return the LOD likelihood for a segment
static double lrr_baf_lod(const float *lrr_arr, const float *baf_arr, int n, const int *imap,
			  float err_log_prb, float lrr_bias, float lrr_hap2dip, float lrr_sd,
			  float baf_sd, double bdev_lrr_baf)
{
	if (n == 0 || bdev_lrr_baf < -0.5 || bdev_lrr_baf > 0.25)
		return -INFINITY; // kmin_brent does not handle NAN

	float ldev = -logf(1.0f - 2.0f * (float)bdev_lrr_baf) / (float)M_LN2 * lrr_hap2dip;
	float ret = 0.0f;
	for (int i = 0; i < n; i++) {
		float lrr = imap ? lrr_arr[imap[i]] : lrr_arr[i];
		float baf = imap ? baf_arr[imap[i]] : baf_arr[i];
		float log_lkl =
			lrr_baf_log_lkl(lrr, baf, ldev, (float)bdev_lrr_baf, lrr_sd, baf_sd,
					lrr_bias)
			- lrr_baf_log_lkl(lrr, baf, 0.0f, 0.0f, lrr_sd, baf_sd, lrr_bias);
		if (log_lkl < err_log_prb)
			log_lkl = err_log_prb;
		else if (log_lkl > -err_log_prb)
			log_lkl = -err_log_prb;
		ret += log_lkl;
	}
	return (double)ret * M_LOG10E;
}

// return the LOD likelihood for a segment
static double baf_lod(const float *baf_arr, int n, const int *imap, float err_log_prb,
		      float baf_sd, double bdev)
{
	if (n == 0 || bdev < 0.0 || bdev > 0.5)
		return -INFINITY; // kmin_brent does not handle NAN

	float ret = 0.0f;
	for (int i = 0; i < n; i++) {
		float baf = imap ? baf_arr[imap[i]] : baf_arr[i];
		if (isnan(baf))
			continue;
		float log_lkl =
			log_mean_expf(norm_log_lkl(baf - 0.5f, (float)bdev, baf_sd, 1.0f),
				      norm_log_lkl(baf - 0.5f, -(float)bdev, baf_sd, 1.0f))
			- norm_log_lkl(baf - 0.5f, 0.0f, baf_sd, 1.0f);
		if (log_lkl < err_log_prb)
			log_lkl = err_log_prb;
		else if (log_lkl > -err_log_prb)
			log_lkl = -err_log_prb;
		ret += log_lkl;
	}
	return (double)ret * M_LOG10E;
}

// return the LOD likelihood for a segment
static double baf_phase_lod(const float *baf_arr, const int8_t *gt_phase, int n,
			    const int *imap, const int8_t *bdev_phase, float err_log_prb,
			    float baf_sd, double bdev)
{
	if (n == 0 || bdev < 0.0 || bdev > 0.5)
		return -INFINITY; // kmin_brent does not handle NAN

	float ret = 0.0f;
	for (int i = 0; i < n; i++) {
		float baf = imap ? baf_arr[imap[i]] : baf_arr[i];
		int8_t p = imap ? gt_phase[imap[i]] : gt_phase[i];
		if (bdev_phase)
			p *= (int8_t)SIGN(bdev_phase[i]); // notice bdev_phase has no imap
		float log_lkl = baf_phase_log_lkl(baf, p, (float)bdev, baf_sd)
				- baf_phase_log_lkl(baf, 0, 0.0f, baf_sd);
		if (log_lkl < err_log_prb)
			log_lkl = err_log_prb;
		else if (log_lkl > -err_log_prb)
			log_lkl = -err_log_prb;
		ret += log_lkl;
	}
	return (double)ret * M_LOG10E;
}

// TODO find a better title for this function
static float compare_models(const float *baf, const int8_t *gt_phase, int n, const int *imap,
			    float xy_log_prb, float err_log_prb, float flip_log_prb,
			    float tel_log_prb, float baf_sd, const float *bdev, int m)
{
	if (n == 0)
		return NAN;
	float *emis_log_lkl =
		baf_phase_emis_log_lkl(baf, gt_phase, n, imap, err_log_prb, baf_sd, bdev, m);
	int8_t *path =
		log_viterbi_run(emis_log_lkl, n, m, xy_log_prb, flip_log_prb, tel_log_prb, 0.0f,
				0, 0); // TODO can I not pass these values instead of 0 0?
	free(emis_log_lkl);
	int nflips = 0;
	for (int i = 1; i < n; i++)
		if (path[i - 1] && path[i] && path[i - 1] != path[i])
			nflips++;
	double f(double x, void *data)
	{
		return -baf_phase_lod(baf, gt_phase, n, imap, path, err_log_prb, baf_sd, x);
	}
	double x, fx = kmin_brent(f, 0.1, 0.2, NULL, KMIN_EPS, &x);
	free(path);
	return -(float)fx + (float)nflips * flip_log_prb * (float)M_LOG10E;
}

void get_max_sum(const int16_t *ad0, const int16_t *ad1, int n, const int *imap, int *n1,
		 int *n2)
{
	*n1 = 0;
	*n2 = 0;
	for (int i = 0; i < n; i++) {
		int a = imap ? ad0[imap[i]] : ad0[i];
		int b = imap ? ad1[imap[i]] : ad1[i];
		if (a != bcf_int16_missing && b != bcf_int16_missing) {
			if (a > *n1)
				*n1 = a;
			if (b > *n1)
				*n1 = b;
			if (a + b > *n2)
				*n2 = a + b;
		}
	}
}

/*********************************
 * WGS AD LIKELIHOODS            *
 *********************************/

static inline float lrr_ad_log_lkl(float lrr, int16_t ad0, int16_t ad1, float ldev,
				   float lrr_sd, float lrr_bias, const beta_binom_t *beta_binom)
{
	return norm_log_lkl(lrr, ldev, lrr_sd, lrr_bias)
	       + log_mean_expf(beta_binom_log_lkl(beta_binom, ad0, ad1),
			       beta_binom_log_lkl(beta_binom, ad1, ad0));
}

static inline float ad_phase_log_lkl(int16_t ad0, int16_t ad1, int8_t phase,
				     const beta_binom_t *beta_binom)
{
	return phase == 0 ? log_mean_expf(beta_binom_log_lkl(beta_binom, ad0, ad1),
					  beta_binom_log_lkl(beta_binom, ad1, ad0))
			  : (phase > 0 ? beta_binom_log_lkl(beta_binom, ad0, ad1)
				       : beta_binom_log_lkl(beta_binom, ad1, ad0));
}

static float *lrr_ad_emis_log_lkl(const float *lrr, const int16_t *ad0, const int16_t *ad1,
				  int T, const int *imap, float err_log_prb, float lrr_bias,
				  float lrr_hap2dip, float lrr_sd, float ad_rho,
				  const float *bdev_lrr_baf_arr, int m)
{
	int N = 1 + 2 * m;
	int n1, n2;
	get_max_sum(ad0, ad1, T, NULL, &n1, &n2);
	float *emis_log_lkl = (float *)malloc(N * T * sizeof(float));
	for (int i = 0; i < 1 + m; i++) {
		float ldev = i == 0 ? 0.0f
				    : -logf(1.0f - 2.0f * bdev_lrr_baf_arr[i - 1])
					      / (float)M_LN2 * lrr_hap2dip;
		float bdev = i == 0 ? 0.0f : fabsf(bdev_lrr_baf_arr[i - 1]);
		beta_binom_t *beta_binom = i == 0 ? beta_binom_null : beta_binom_alt;
		beta_binom_update(beta_binom, 0.5f + bdev, ad_rho, n1, n2);

		for (int t = 0; t < T; t++) {
			float x = imap ? lrr[imap[t]] : lrr[t];
			int16_t a = imap ? ad0[imap[t]] : ad0[t];
			int16_t b = imap ? ad1[imap[t]] : ad1[t];
			emis_log_lkl[t * N + i] =
				lrr_ad_log_lkl(x, a, b, ldev, lrr_sd, lrr_bias, beta_binom);
		}
	}
	// generate states that should attract LRR waves with no BAF signal
	for (int i = 0; i < m; i++) {
		// do not make extreme states compete with LRR waves
		if (bdev_lrr_baf_arr[i] < -1.0f / 6.0f || bdev_lrr_baf_arr[i] >= 1.0f / 6.0f) {
			for (int t = 0; t < T; t++)
				emis_log_lkl[t * N + 1 + m + i] =
					emis_log_lkl[t * N] + err_log_prb;
		} else {
			float ldev = -logf(1.0f - 2.0f * bdev_lrr_baf_arr[i]) / (float)M_LN2
				     * lrr_hap2dip;
			for (int t = 0; t < T; t++) {
				float x = imap ? lrr[imap[t]] : lrr[t];
				int16_t a = imap ? ad0[imap[t]] : ad0[t];
				int16_t b = imap ? ad1[imap[t]] : ad1[t];
				emis_log_lkl[t * N + 1 + m + i] = lrr_ad_log_lkl(
					x, a, b, ldev, lrr_sd, lrr_bias, beta_binom_null);
			}
		}
	}
	for (int t = 0; t < T; t++)
		rescale_emis_log_lkl(&emis_log_lkl[t * N], N, err_log_prb);
	return emis_log_lkl;
}

static float *ad_phase_emis_log_lkl(const int16_t *ad0, const int16_t *ad1,
				    const int8_t *gt_phase, int T, const int *imap,
				    float err_log_prb, float ad_rho, const float *bdev_arr,
				    int m)
{
	int N = 1 + 2 * m;
	float *emis_log_lkl = (float *)malloc(N * T * sizeof(float));
	for (int i = 0; i < 1 + m; i++) {
		float bdev = i == 0 ? 0.0f : bdev_arr[i - 1];
		// TODO this function should come out of the loop
		int n1, n2;
		get_max_sum(ad0, ad1, T, imap, &n1, &n2);
		beta_binom_t *beta_binom = i == 0 ? beta_binom_null : beta_binom_alt;
		beta_binom_update(beta_binom, 0.5f + bdev, ad_rho, n1, n2);

		for (int t = 0; t < T; t++) {
			int16_t a = imap ? ad0[imap[t]] : ad0[t];
			int16_t b = imap ? ad1[imap[t]] : ad1[t];
			int8_t p = imap ? gt_phase[imap[t]] : gt_phase[t];
			emis_log_lkl[t * N + i] = ad_phase_log_lkl(a, b, p, beta_binom);
			if (i > 0)
				emis_log_lkl[t * N + m + i] =
					ad_phase_log_lkl(b, a, p, beta_binom);
		}
	}
	for (int t = 0; t < T; t++)
		rescale_emis_log_lkl(&emis_log_lkl[t * N], N, err_log_prb);
	return emis_log_lkl;
}

static int cnp_edge_is_not_cn2_lrr_ad(const float *lrr, int16_t *ad0, int16_t *ad1, int n,
				      int a, int b, float xy_log_prb, float err_log_prb,
				      float lrr_bias, float lrr_hap2dip, float lrr_sd,
				      float ad_rho, float ldev, float bdev)
{
	int n1, n2;
	get_max_sum(ad0, ad1, n, NULL, &n1, &n2);
	beta_binom_update(beta_binom_null, 0.5f, ad_rho, n1, n2);
	beta_binom_update(beta_binom_alt, 0.5f + bdev, ad_rho, n1, n2);

	// test left edge
	float sum_log_lkl = 0.0f;
	for (int i = a - 1; i >= 0; i--) {
		float log_lkl = lrr_ad_log_lkl(lrr[i], ad0[i], ad1[i], ldev, lrr_sd, lrr_bias,
					       beta_binom_alt)
				- lrr_ad_log_lkl(lrr[i], ad0[i], ad1[i], 0.0f, lrr_sd, lrr_bias,
						 beta_binom_null);

		if (log_lkl < err_log_prb)
			log_lkl = err_log_prb;
		else if (log_lkl > -err_log_prb)
			log_lkl = -err_log_prb;
		sum_log_lkl += log_lkl;
		if (sum_log_lkl > -xy_log_prb)
			return -1;
		if (sum_log_lkl < xy_log_prb)
			break;
	}

	// test right edge
	sum_log_lkl = 0.0f;
	for (int i = b + 1; i < n; i++) {
		float log_lkl = lrr_ad_log_lkl(lrr[i], ad0[i], ad1[i], ldev, lrr_sd, lrr_bias,
					       beta_binom_alt)
				- lrr_ad_log_lkl(lrr[i], ad0[i], ad1[i], 0.0f, lrr_sd, lrr_bias,
						 beta_binom_null);
		if (log_lkl < err_log_prb)
			log_lkl = err_log_prb;
		else if (log_lkl > -err_log_prb)
			log_lkl = -err_log_prb;
		sum_log_lkl += log_lkl;
		if (sum_log_lkl > -xy_log_prb)
			return -1;
		if (sum_log_lkl < xy_log_prb)
			break;
	}

	return 0;
}

// return the LOD likelihood for a segment
static double lrr_ad_lod(const float *lrr_arr, const int16_t *ad0_arr, const int16_t *ad1_arr,
			 int n, const int *imap, float err_log_prb, float lrr_bias,
			 float lrr_hap2dip, float lrr_sd, float ad_rho, double bdev_lrr_baf)
{
	if (n == 0 || bdev_lrr_baf < -0.5 || bdev_lrr_baf > 0.25)
		return -INFINITY; // kmin_brent does not handle NAN

	float ldev = -logf(1.0f - 2.0f * (float)bdev_lrr_baf) / (float)M_LN2 * lrr_hap2dip;
	int n1, n2;
	get_max_sum(ad0_arr, ad1_arr, n, imap, &n1, &n2);
	beta_binom_update(beta_binom_null, 0.5f, ad_rho, n1, n2);
	beta_binom_update(beta_binom_alt, 0.5f + (float)bdev_lrr_baf, ad_rho, n1, n2);
	float ret = 0.0f;
	for (int i = 0; i < n; i++) {
		float lrr = imap ? lrr_arr[imap[i]] : lrr_arr[i];
		int16_t ad0 = imap ? ad0_arr[imap[i]] : ad0_arr[i];
		int16_t ad1 = imap ? ad1_arr[imap[i]] : ad1_arr[i];
		float log_lkl =
			lrr_ad_log_lkl(lrr, ad0, ad1, ldev, lrr_sd, lrr_bias, beta_binom_alt)
			- lrr_ad_log_lkl(lrr, ad0, ad1, 0.0f, lrr_sd, lrr_bias,
					 beta_binom_null);
		if (log_lkl < err_log_prb)
			log_lkl = err_log_prb;
		else if (log_lkl > -err_log_prb)
			log_lkl = -err_log_prb;
		ret += log_lkl;
	}
	return (double)ret * M_LOG10E;
}

// return the LOD likelihood for a segment
static double ad_lod(const int16_t *ad0_arr, const int16_t *ad1_arr, int n, const int *imap,
		     float err_log_prb, float ad_rho, double bdev)
{
	if (n == 0 || bdev < 0.0 || bdev > 0.5)
		return -INFINITY; // kmin_brent does not handle NAN

	int n1, n2;
	get_max_sum(ad0_arr, ad1_arr, n, imap, &n1, &n2);
	beta_binom_update(beta_binom_null, 0.5f, ad_rho, n1, n2);
	beta_binom_update(beta_binom_alt, 0.5f + (float)bdev, ad_rho, n1, n2);
	float ret = 0.0f;
	for (int i = 0; i < n; i++) {
		int16_t ad0 = imap ? ad0_arr[imap[i]] : ad0_arr[i];
		int16_t ad1 = imap ? ad1_arr[imap[i]] : ad1_arr[i];
		float log_lkl = log_mean_expf(beta_binom_log_lkl(beta_binom_alt, ad0, ad1),
					      beta_binom_log_lkl(beta_binom_alt, ad1, ad0))
				- beta_binom_log_lkl(beta_binom_null, ad0, ad1);
		if (log_lkl < err_log_prb)
			log_lkl = err_log_prb;
		else if (log_lkl > -err_log_prb)
			log_lkl = -err_log_prb;
		ret += log_lkl;
	}
	return (double)ret * M_LOG10E;
}

// return the LOD likelihood for a segment
static double ad_phase_lod(const int16_t *ad0_arr, const int16_t *ad1_arr,
			   const int8_t *gt_phase, int n, const int *imap,
			   const int8_t *bdev_phase, float err_log_prb, float ad_rho,
			   double bdev)
{
	if (n == 0 || bdev < 0.0 || bdev > 0.5)
		return -INFINITY; // kmin_brent does not handle NAN

	int n1, n2;
	get_max_sum(ad0_arr, ad1_arr, n, imap, &n1, &n2);
	beta_binom_update(beta_binom_null, 0.5f, ad_rho, n1, n2);
	beta_binom_update(beta_binom_alt, 0.5f + (float)bdev, ad_rho, n1, n2);
	float ret = 0.0f;
	for (int i = 0; i < n; i++) {
		int16_t ad0 = imap ? ad0_arr[imap[i]] : ad0_arr[i];
		int16_t ad1 = imap ? ad1_arr[imap[i]] : ad1_arr[i];
		int8_t p = imap ? gt_phase[imap[i]] : gt_phase[i];
		if (bdev_phase)
			p *= (int8_t)SIGN(bdev_phase[i]); // notice bdev_phase has no imap
		float log_lkl = ad_phase_log_lkl(ad0, ad1, p, beta_binom_alt)
				- ad_phase_log_lkl(ad0, ad1, 0, beta_binom_null);
		if (log_lkl < err_log_prb)
			log_lkl = err_log_prb;
		else if (log_lkl > -err_log_prb)
			log_lkl = -err_log_prb;
		ret += log_lkl;
	}
	return (double)ret * M_LOG10E;
}

// TODO find a better title for this function
static float compare_wgs_models(const int16_t *ad0, const int16_t *ad1, const int8_t *gt_phase,
				int n, const int *imap, float xy_log_prb, float err_log_prb,
				float flip_log_prb, float tel_log_prb, float ad_rho,
				const float *bdev, int m)
{
	if (n == 0)
		return NAN;
	float *emis_log_lkl = ad_phase_emis_log_lkl(ad0, ad1, gt_phase, n, imap, err_log_prb,
						    ad_rho, bdev, m);
	int8_t *path =
		log_viterbi_run(emis_log_lkl, n, m, xy_log_prb, flip_log_prb, tel_log_prb, 0.0f,
				0, 0); // TODO can I not pass these values instead of 0 0?
	free(emis_log_lkl);
	int nflips = 0;
	for (int i = 1; i < n; i++)
		if (path[i - 1] && path[i] && path[i - 1] != path[i])
			nflips++;
	double f(double x, void *data)
	{
		return -ad_phase_lod(ad0, ad1, gt_phase, n, imap, path, err_log_prb, ad_rho, x);
	}
	double x, fx = kmin_brent(f, 0.1, 0.2, NULL, KMIN_EPS, &x);
	free(path);
	return -(float)fx + (float)nflips * flip_log_prb * (float)M_LOG10E;
}

// TODO change this or integrate with ad_lod
static double lod_lkl_beta_binomial(const int16_t *ad0_arr, const int16_t *ad1_arr, int n,
				    const int *imap, double ad_rho)
{
	if (n == 0 || ad_rho <= 0.0 || ad_rho >= 1.0)
		return -INFINITY;
	float ret = 0.0f;
	int n1, n2;
	get_max_sum(ad0_arr, ad1_arr, n, imap, &n1, &n2);
	beta_binom_update(beta_binom_null, 0.5f, ad_rho, n1, n2);
	for (int i = 0; i < n; i++) {
		int16_t ad0 = imap ? ad0_arr[imap[i]] : ad0_arr[i];
		int16_t ad1 = imap ? ad1_arr[imap[i]] : ad1_arr[i];
		ret += beta_binom_log_lkl(beta_binom_null, ad0, ad1);
	}
	return (double)ret * M_LOG10E;
}

/*********************************
 * BASIC STATISTICS FUNCTIONS    *
 *********************************/

// iterator of non-NaN values
static inline float next_not_missing(const float *v, const int *imap, int n, int *i)
{
	float x = NAN;
	while (*i < n) {
		x = imap ? v[imap[*i]] : v[*i];
		if (!isnan(x))
			break;
		(*i)++;
	}
	return x;
}

// compute BAF phase concordance for a float array with iterator
static void get_baf_conc(const float *baf, const int8_t *gt_phase, int n, const int *imap,
			 int *conc, int *disc)
{
	int i;
	float prev = NAN, next = NAN;
	*conc = 0, *disc = 0;
	for (i = 0; i < n; i++) {
		prev = ((imap ? baf[imap[i]] : baf[i]) - 0.5f)
		       * (imap ? gt_phase[imap[i]] : gt_phase[i]);
		if (!isnan(prev) && prev != 0.0f)
			break;
	}
	if (i == n)
		return;

	for (i++; i < n; i++) {
		next = ((imap ? baf[imap[i]] : baf[i]) - 0.5f)
		       * (imap ? gt_phase[imap[i]] : gt_phase[i]);
		if (!isnan(next) && next != 0.0f) {
			if (prev * next > 0.0f)
				(*conc)++;
			else if (prev * next < 0.0f)
				(*disc)++;
			prev = next;
		}
	}
}

// compute phased BAF autocorrelation for a float array with iterator
static float get_baf_auto_corr(const float *baf, const int8_t *gt_phase, int n, const int *imap)
{
	double var = 0.0, auto_corr = 0.0;
	float prev = NAN, next = NAN;
	for (int i = 0; i < n; i++) {
		next = ((imap ? baf[imap[i]] : baf[i]) - 0.5f)
		       * (imap ? gt_phase[imap[i]] : gt_phase[i]);
		if (!isnan(next)) {
			var += sq((double)next);
			if (!isnan(prev))
				auto_corr += prev * next;
			prev = next;
		}
	}
	auto_corr /= var;
	return auto_corr;
}

static float get_sample_sd(const float *v, int n, const int *imap);

// compute (adjusted) LRR autocorrelation for a float array with iterator
static float get_lrr_auto_corr(const float *lrr, int n, const int *imap)
{
	float value;
	double mean = 0.0;
	int i = 0, j = 0;
	for (value = next_not_missing(lrr, imap, n, &i); i < n;
	     i++, value = next_not_missing(lrr, imap, n, &i)) {
		mean += (double)value;
		j++;
	}
	if (j <= 1)
		return NAN;
	mean /= (double)j;

	double var = 0.0;
	i = 0;
	for (value = next_not_missing(lrr, imap, n, &i); i < n;
	     i++, value = next_not_missing(lrr, imap, n, &i)) {
		var += sq((double)value - mean);
	}

	double auto_corr = 0.0;
	i = 0;
	double prev = (double)next_not_missing(lrr, imap, n, &i) - mean, next;
	for (i++, value = next_not_missing(lrr, imap, n, &i); i < n;
	     i++, value = next_not_missing(lrr, imap, n, &i)) {
		next = (double)value - mean;
		auto_corr += prev * next;
		prev = next;
	}
	auto_corr /= var;
	return auto_corr;
}

// compute the n50 of a vector
static int get_n50(const int *v, int n, const int *imap)
{
	if (n <= 1)
		return -1;
	int i;
	int sum, sum2;
	int *w = (int *)malloc((n - 1) * sizeof(int));

	for (i = 0, sum = 0; i < n - 1; i++) {
		w[i] = imap ? v[imap[i + 1]] - v[imap[i]] : v[i + 1] - v[i];
		sum += w[i];
	}
	sum /= 2;

	ks_introsort_int((size_t)n - 1, w);

	for (i = 0, sum2 = 0; sum2 < sum && i < n - 1; i++)
		sum2 += w[i];
	int n50 = w[i - 1];
	free(w);
	return n50;
}

// compute sample standard deviation of a float array (with iterator)
// sqrt ( ( \sum x^2 - (\sum x)^2 / N ) / ( N - 1 ) )
// TODO this is not completely okay as the above way to compute the sd is error prone with
// floats
static float get_sample_sd(const float *v, int n, const int *imap)
{
	// float s = 0.0, s2 = 0.0;
	double mean = 0.0;
	int j = 0;
	for (int i = 0; i < n; i++) {
		double tmp = (double)(imap ? v[imap[i]] : v[i]);
		if (!isnan(tmp)) {
			mean += tmp;
			j++;
		}
	}
	if (j <= 1)
		return NAN;
	mean /= (double)j;

	double s2 = 0.0;
	for (int i = 0; i < n; i++) {
		double tmp = (double)(imap ? v[imap[i]] : v[i]);
		if (!isnan(tmp))
			s2 += sq(tmp - mean);
	}
	s2 /= (double)(j - 1);

	return (float)sqrt(s2);
}

// compute standard error of mean a float array (with iterator)
// sqrt ( ( \sum x^2 - (\sum x)^2 / N ) / ( N - 1 ) / N )
static float get_se_mean(const float *v, int n, const int *imap)
{
	int j = 0;
	for (int i = 0; i < n; i++) {
		float tmp = imap ? v[imap[i]] : v[i];
		if (!isnan(tmp))
			j++;
	}
	if (j <= 1)
		return NAN;

	return get_sample_sd(v, n, imap) / sqrtf(j);
}

/*********************************
 * SAMPLE METHODS                *
 *********************************/

static void mocha_print_ucsc(FILE *restrict stream, const mocha_t *mocha, int n,
			     const bcf_hdr_t *hdr)
{
	if (stream == NULL)
		return;
	const char *name[4];
	name[MOCHA_UNK] = "mCA_undetermined";
	name[MOCHA_DEL] = "mCA_loss";
	name[MOCHA_DUP] = "mCA_gain";
	name[MOCHA_UPD] = "mCA_neutral";
	const char *desc[4];
	desc[MOCHA_UNK] = "Undetermined";
	desc[MOCHA_DEL] = "Deletions";
	desc[MOCHA_DUP] = "Duplications";
	desc[MOCHA_UPD] = "CN-LOHs";
	kstring_t tmp = {0, 0, NULL};
	for (int i = 0; i < 4; i++) {
		fprintf(stream,
			"track name=%s description=\"%s\" visibility=4 priority=1 itemRgb=\"On\"\n",
			name[i], desc[i]);
		for (int j = 0; j < n; j++) {
			uint8_t red = 127, green = 127, blue = 127;
			switch (i) {
			case MOCHA_DEL:
				red -= (uint8_t)(127.0f * sqf(mocha[j].cf));
				green -= (uint8_t)(127.0f * sqf(mocha[j].cf));
				blue += (uint8_t)(128.0f * sqf(mocha[j].cf));
				break;
			case MOCHA_DUP:
				red += (uint8_t)(128.0f * sqf(mocha[j].cf));
				green -= (uint8_t)(127.0f * sqf(mocha[j].cf));
				blue -= (uint8_t)(127.0f * sqf(mocha[j].cf));
				break;
			case MOCHA_UPD:
				red += (uint8_t)(128.0f * mocha[j].cf);
				green += (uint8_t)(38.0f * mocha[j].cf);
				blue -= (uint8_t)(127.0f * mocha[j].cf);
				break;
			}
			if (i == mocha[j].type) {
				const char *sample_name =
					bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, mocha[j].sample_idx);
				tmp.l = 0;
				kputs(sample_name, &tmp);
				for (char *p = tmp.s; (p = strchr(p, ' ')); ++p)
					*p = '_';
				const char *seq_name = bcf_hdr_id2name(hdr, mocha[j].rid);
				if (strncmp(seq_name, "chr", 3) == 0)
					fprintf(stream,
						"%s\t%d\t%d\t%s\t%d\t.\t%d\t%d\t%d,%d,%d\n",
						seq_name, mocha[j].beg_pos, mocha[j].end_pos,
						tmp.s,
						isnan(mocha[j].cf) ? 0
								   : (int)(1e3 * mocha[j].cf),
						mocha[j].beg_pos, mocha[j].end_pos, red, green,
						blue);
				else
					fprintf(stream,
						"chr%s\t%d\t%d\t%s\t%d\t.\t%d\t%d\t%d,%d,%d\n",
						seq_name, mocha[j].beg_pos, mocha[j].end_pos,
						tmp.s,
						isnan(mocha[j].cf) ? 0
								   : (int)(1e3 * mocha[j].cf),
						mocha[j].beg_pos, mocha[j].end_pos, red, green,
						blue);
			}
		}
	}
	free(tmp.s);
	if (stream != stdout && stream != stderr)
		fclose(stream);
}

static void mocha_print(FILE *restrict stream, const mocha_t *mocha, int n,
			const bcf_hdr_t *hdr, int flags, char *genome, float lrr_hap2dip)
{
	if (stream == NULL)
		return;
	char sex[3];
	sex[SEX_UNK] = 'U';
	sex[SEX_MAL] = 'M';
	sex[SEX_FEM] = 'F';
	const char *type[6];
	type[MOCHA_UNK] = "Undetermined";
	type[MOCHA_DEL] = "Deletion";
	type[MOCHA_DUP] = "Duplication";
	type[MOCHA_UPD] = "CN-LOH";
	type[MOCHA_CNP_DEL] = "CNP Deletion";
	type[MOCHA_CNP_DUP] = "CNP Duplication";
	char arm_type[3];
	arm_type[MOCHA_NOT] = 'N';
	arm_type[MOCHA_ARM] = 'Y';
	arm_type[MOCHA_TEL] = 'T';
	fprintf(stream,
		"SAMPLE\tSEX\tCHROM\tBEG_%s\tEND_%s\tLENGTH\tP_ARM\tQ_ARM\tNSITES\tNHETS\tN50_HETS\tBDEV\tBDEV_SE\tREL_COV\tREL_COV_SE\tLOD_LRR_BAF\tLOD_BAF_PHASE\tNFLIPS\tBAF_CONC\tLOD_BAF_CONC\tTYPE\tCF\n",
		genome, genome);
	for (int i = 0; i < n; i++) {
		const char *sample_name = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, mocha->sample_idx);
		const char *seq_name = bcf_hdr_id2name(hdr, mocha->rid);
		fprintf(stream,
			"%s\t%c\t%s\t%d\t%d\t%d\t%c\t%c\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%d\t%.4f\t%.2f\t%s\t%.4f\n",
			sample_name, sex[mocha->sex], seq_name, mocha->beg_pos, mocha->end_pos,
			mocha->length, arm_type[mocha->p_arm], arm_type[mocha->q_arm],
			mocha->nsites, mocha->nhets, mocha->n50_hets, mocha->bdev,
			mocha->bdev_se, 2.0f * expf(mocha->ldev / lrr_hap2dip * (float)M_LN2),
			2.0f * expf(mocha->ldev / lrr_hap2dip * (float)M_LN2) * mocha->ldev_se
				/ lrr_hap2dip * (float)M_LN2,
			mocha->lod_lrr_baf, mocha->lod_baf_phase, mocha->nflips,
			mocha->baf_conc, mocha->lod_baf_conc, type[mocha->type], mocha->cf);
		mocha++;
	}
	if (stream != stdout && stream != stderr)
		fclose(stream);
}

// this function returns two values (a, b) such that:
// (i) a <= b (ii) pos[a-1] < beg (iii) beg <= pos[a] < end (iv) beg <= pos[b] < end (v)
// pos[b+1] >= end
static int get_cnp_edges(const int *pos, int n, int beg, int end, int *a, int *b)
{
	if (pos[0] >= end || pos[n - 1] < beg)
		return -1;

	int i = 0, j = n - 1, k;
	while (j - i > 1) {
		k = (i + j) / 2;
		if (pos[k] < beg)
			i = k;
		else
			j = k;
	}
	if (pos[j] >= end)
		return -1;

	*a = j;

	i = j;
	j = n - 1;
	while (j - i > 1) {
		k = (i + j) / 2;
		if (pos[k] < end)
			i = k;
		else
			j = k;
	}

	*b = i;
	return 0;
}

// classify mosaic chromosomal alteration type based on LRR and BAF
// LDEV = -log2( 1 - 2 x BDEV ) * LRR-hap2dip for duplications
// LDEV = -log2( 1 + 2 x BDEV ) * LRR-hap2dip for deletions
static int8_t mocha_type(float ldev, float ldev_se, float bdev, int nhets, float lrr_hap2dip,
			 int8_t p_arm, int8_t q_arm)
{
	// a LOD score can be computed from a chi-squared statistic by dividing by 2ln(10) ~ 4.6
	float z2_upd = sqf(ldev / ldev_se);
	// equivalent of a 2 LOD score bonus for ending in one but not two telomeres
	if ((p_arm != MOCHA_TEL && q_arm != MOCHA_TEL)
	    || (p_arm == MOCHA_TEL && q_arm == MOCHA_TEL))
		z2_upd += 4.0f * M_LN10;
	else
		z2_upd -= 4.0f * M_LN10;

	// if one model has 2 LOD scores point more than the other model, select the better
	// model
	if (ldev > 0) {
		if (nhets < 5 || isnan(bdev)) {
			return MOCHA_DUP;
		} else {
			float expected_ldev =
				-logf(1.0f - 2.0f * bdev > 2.0f / 3.0f ? 1.0f - 2.0f * bdev
								       : 2.0f / 3.0f)
				* M_LOG2E * lrr_hap2dip;
			float z2_dup = sqf((ldev - expected_ldev) / ldev_se);
			if (z2_upd > z2_dup + 4.0f * M_LN10)
				return MOCHA_DUP;
			if (z2_dup > z2_upd + 4.0f * M_LN10)
				return MOCHA_UPD;
		}
	} else {
		if (nhets < 5 || isnan(bdev)) {
			return MOCHA_DEL;
		} else {
			float expected_ldev = -logf(1.0f + 2.0f * bdev) * M_LOG2E * lrr_hap2dip;
			float z2_del = sqf((ldev - expected_ldev) / ldev_se);
			if (z2_upd > z2_del + 4.0f * M_LN10)
				return MOCHA_DEL;
			if (z2_del > z2_upd + 4.0f * M_LN10)
				return MOCHA_UPD;
		}
	}
	return MOCHA_UNK;
}

// best estimate for cell fraction using the following formula (if BDEV is available):
// BDEV = | 1 / 2 - 1 / CNF |
// CNF = 2 / ( 1 + 2 x BDEV ) for deletions
// CNF = 2 / ( 1 - 2 x BDEV ) for duplications
// CNF = 2 x 2^( LDEV / LRR-hap2dip )
// LDEV = - LRR-hap2dip / ln(2) * ln( 1 + 2 x BDEV ) for deletions
// LDEV = - LRR-hap2dip / ln(2) * ln( 1 - 2 x BDEV ) for duplications
static float mocha_cell_fraction(float ldev, float ldev_se, float bdev, int nhets, int8_t type,
				 float lrr_hap2dip)
{
	if (nhets < 5 || isnan(bdev)) {
		switch (type) {
		case MOCHA_DEL:
			return ldev < -lrr_hap2dip
				       ? 1.0f
				       : -2.0f
						 * (expf(ldev / lrr_hap2dip * (float)M_LN2)
						    - 1.0f);
		case MOCHA_DUP:
			return ldev * (float)M_LN2 > lrr_hap2dip * logf(1.5f)
				       ? 1.0f
				       : 2.0f
						 * (expf(ldev / lrr_hap2dip * (float)M_LN2)
						    - 1.0f);
		default:
			return NAN;
		}
	} else {
		switch (type) {
		case MOCHA_DEL:
			return 4.0f * bdev / (1.0f + 2.0f * bdev);
		case MOCHA_DUP:
			return bdev > 1.0f / 6.0f ? 1.0f : 4.0f * bdev / (1.0f - 2.0f * bdev);
		case MOCHA_UPD:
			return 2.0f * bdev;
		case MOCHA_UNK:
			return bdev < 0.05f ? 4.0f * bdev : NAN; // here it assumes it is either
								 // a deletion or a duplication
		default:
			return NAN;
		}
	}
}

static void get_mocha_stats(const int *pos, const float *lrr, const float *baf,
			    const int8_t *gt_phase, int n, int a, int b, int cen_beg,
			    int cen_end, int length, float baf_conc, mocha_t *mocha)
{
	mocha->nsites = b + 1 - a;

	if (a == 0)
		if (pos[a] < cen_beg)
			mocha->beg_pos = 0;
		else
			mocha->beg_pos = cen_end;
	else
		mocha->beg_pos = pos[a];

	if (b == n - 1)
		if (pos[b] >= cen_end)
			mocha->end_pos = length;
		else
			mocha->end_pos = cen_beg;
	else
		mocha->end_pos = pos[b];

	mocha->length = mocha->end_pos - mocha->beg_pos;

	if (mocha->beg_pos == 0)
		mocha->p_arm = MOCHA_TEL;
	else if (mocha->beg_pos < cen_beg)
		mocha->p_arm = MOCHA_ARM;
	else
		mocha->p_arm = MOCHA_NOT;

	if (mocha->end_pos == length)
		mocha->q_arm = MOCHA_TEL;
	else if (mocha->end_pos > cen_end)
		mocha->q_arm = MOCHA_ARM;
	else
		mocha->q_arm = MOCHA_NOT;

	mocha->ldev_se = get_se_mean(lrr + a, b + 1 - a, NULL);
	mocha->nhets = 0;
	for (int i = a; i <= b; i++)
		if (!isnan(baf[i]))
			mocha->nhets++;
	int conc, disc;
	get_baf_conc(baf + a, gt_phase + a, b + 1 - a, NULL, &conc, &disc);
	mocha->baf_conc = conc + disc > 0 ? (float)conc / (float)(conc + disc) : NAN;
	mocha->lod_baf_conc =
		((mocha->baf_conc > 0 ? (float)conc * logf(mocha->baf_conc / baf_conc) : 0)
		 + (mocha->baf_conc < 1
			    ? (float)disc * logf((1 - mocha->baf_conc) / (1 - baf_conc))
			    : 0))
		* (float)M_LOG10E;
	mocha->n50_hets = get_n50(pos + a, b + 1 - a, NULL);
	mocha->nflips = -1;
	mocha->bdev = NAN;
	mocha->bdev_se = NAN;
	mocha->lod_baf_phase = NAN;
}

// return segments called by the HMM or a suggestion of what state should be added to the HMM
static float get_path_segs(const int8_t *path, const float *hs_arr, int n, int hmm_model,
			   int middle, int **beg, int *m_beg, int **end, int *m_end, int *nseg)
{
	int a = 0, b = 0;
	*nseg = 0;
	for (b = 0; b < n; b++) {
		// check whether it is the end of a segment
		if (b != n - 1) {
			int x = abs(path[b]);
			int y = abs(path[b + 1]);
			if (x == y)
				continue;

			// swap the two elements
			if (x > y) {
				x = abs(path[b + 1]);
				y = abs(path[b]);
			}

			if (y - x == 1) {
				if (hmm_model == LRR_BAF && x > 0 && x != middle)
					return (hs_arr[x - 1] + hs_arr[y - 1]) * 0.5f;
				else if (hmm_model == BAF_PHASE)
					return x == 0 ? hs_arr[0] * 0.5f
						      : (hs_arr[x - 1] + hs_arr[y - 1]) * 0.5f;
			}
		}

		if (path[b]) {
			(*nseg)++;
			hts_expand(int, *nseg, *m_beg, *beg);
			(*beg)[(*nseg) - 1] = a;
			hts_expand(int, *nseg, *m_end, *end);
			(*end)[(*nseg) - 1] = b;
		}
		a = b + 1;
	}
	return 0;
}

// process one contig for one sample
static void sample_run(sample_t *self, mocha_table_t *mocha_table, const model_t *model)
{
	// do nothing if chromosome Y or MT are being tested
	if (model->rid == model->genome_rules->y_rid
	    || model->rid == model->genome_rules->mt_rid) {
		memset(self->data_arr[LDEV], 0, self->n * sizeof(int16_t));
		memset(self->data_arr[BDEV], 0, self->n * sizeof(int16_t));
		memset(self->phase_arr, 0, self->n * sizeof(int8_t));
		return;
	}

	mocha_t mocha;
	mocha.sample_idx = self->idx;
	mocha.sex = self->sex;
	mocha.rid = model->rid;

	int cen_beg = model->genome_rules->cen_beg[model->rid];
	int cen_end = model->genome_rules->cen_end[model->rid];
	int length = model->genome_rules->length[model->rid];

	// declutter code by copying these values onto the stack
	int n = self->n;
	int8_t *gt_phase = self->phase_arr;
	int16_t *ad0 = self->data_arr[AD0];
	int16_t *ad1 = self->data_arr[AD1];
	float *lrr = (float *)malloc(n * sizeof(float));
	float *baf = (float *)malloc(n * sizeof(float));
	if (model->flags & WGS_DATA) {
		ad_to_lrr_baf(ad0, ad1, lrr, baf, n);
	} else {
		for (int i = 0; i < n; i++) {
			lrr[i] = int16_to_float(self->data_arr[LRR][i]);
			baf[i] = int16_to_float(self->data_arr[BAF][i]);
		}
	}

	if (model->order_lrr_gc >= 0)
		adjust_lrr(lrr, model->gc_arr, n, self->vcf_imap_arr, self->stats.coeffs,
			   model->order_lrr_gc);

	int8_t *bdev_phase = (int8_t *)calloc(n, sizeof(int8_t));
	int *pos = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++)
		pos[i] = model->locus_arr[self->vcf_imap_arr[i]].pos;
	int *imap_arr = (int *)malloc(n * sizeof(int));
	int *hets_imap_arr = (int *)malloc(n * sizeof(int));
	float *pbaf_arr = (float *)malloc(n * sizeof(float));

	int16_t *ldev = (int16_t *)calloc(n, sizeof(int16_t));
	int16_t *bdev = (int16_t *)calloc(n, sizeof(int16_t));

	// TODO do I need special normalization for the X nonPAR region?
	if (model->rid == model->genome_rules->x_rid) {
		for (int i = 0; i < n; i++) {
			if (pos[i] > model->genome_rules->x_nonpar_beg
			    && pos[i] < model->genome_rules->x_nonpar_end) {
				lrr[i] = (self->sex == SEX_MAL) ? NAN
								: lrr[i] - model->lrr_auto2sex;
				baf[i] = (self->sex == SEX_MAL) ? NAN : baf[i];
			}
		}
	}

	if (model->cnp_itr) {
		while (regitr_overlap(model->cnp_itr)) {
			int a, b;
			if (get_cnp_edges(pos, n, model->cnp_itr->beg, model->cnp_itr->end, &a,
					  &b)
			    == 0) {
				int cnp_type = regitr_payload(model->cnp_itr, int);
				float exp_ldev = NAN;
				float exp_bdev = NAN;
				mocha.type = MOCHA_UNK;
				mocha.ldev = get_median(lrr + a, b + 1 - a, NULL);
				if (mocha.ldev > 0
				    && (cnp_type == MOCHA_CNP_DUP
					|| cnp_type == MOCHA_CNP_CNV)) {
					if (model->flags & WGS_DATA)
						mocha.lod_lrr_baf = lrr_ad_lod(
							lrr + a, ad0 + a, ad1 + a, b + 1 - a,
							NULL, model->err_log_prb,
							model->lrr_bias, model->lrr_hap2dip,
							self->adjlrr_sd, self->stats.dispersion,
							1.0f / 6.0f);
					else
						mocha.lod_lrr_baf = lrr_baf_lod(
							lrr + a, baf + a, b + 1 - a, NULL,
							model->err_log_prb, model->lrr_bias,
							model->lrr_hap2dip, self->adjlrr_sd,
							self->stats.dispersion, 1.0f / 6.0f);
					if (mocha.lod_lrr_baf
					    > -model->xy_log_prb * (float)M_LOG10E) {
						mocha.type = MOCHA_CNP_DUP;
						mocha.cf = NAN;
						exp_ldev = log2f(1.5f) * model->lrr_hap2dip;
						exp_bdev = 1.0f / 6.0f;
					}
				} else if (mocha.ldev <= 0
					   && (cnp_type == MOCHA_CNP_DEL
					       || cnp_type == MOCHA_CNP_CNV)) {
					if (model->flags & WGS_DATA)
						mocha.lod_lrr_baf = lrr_ad_lod(
							lrr + a, ad0 + a, ad1 + a, b + 1 - a,
							NULL, model->err_log_prb,
							model->lrr_bias, model->lrr_hap2dip,
							self->adjlrr_sd, self->stats.dispersion,
							-0.5f);
					else
						mocha.lod_lrr_baf = lrr_baf_lod(
							lrr + a, baf + a, b + 1 - a, NULL,
							model->err_log_prb, model->lrr_bias,
							model->lrr_hap2dip, self->adjlrr_sd,
							self->stats.dispersion, -0.5f);
					if (mocha.lod_lrr_baf
					    > -model->xy_log_prb * (float)M_LOG10E) {
						mocha.type = MOCHA_CNP_DEL;
						mocha.cf = NAN;
						exp_ldev = -model->lrr_hap2dip;
						exp_bdev = 0.5f;
					}
				}
				if (mocha.type == MOCHA_CNP_DUP
				    || mocha.type == MOCHA_CNP_DEL) {
					if (model->flags & WGS_DATA) {
						if (cnp_edge_is_not_cn2_lrr_ad(
							    lrr, ad0, ad1, n, a, b,
							    model->xy_log_prb,
							    model->err_log_prb, model->lrr_bias,
							    model->lrr_hap2dip, self->adjlrr_sd,
							    self->stats.dispersion, exp_ldev,
							    exp_bdev))
							continue;
					} else {
						if (cnp_edge_is_not_cn2_lrr_baf(
							    lrr, baf, n, a, b,
							    model->xy_log_prb,
							    model->err_log_prb, model->lrr_bias,
							    model->lrr_hap2dip, self->adjlrr_sd,
							    self->stats.dispersion, exp_ldev,
							    exp_bdev))
							continue;
					}
					get_mocha_stats(pos, lrr, baf, gt_phase, n, a, b,
							cen_beg, cen_end, length,
							self->stats.baf_conc, &mocha);
					// compute bdev, if possible
					if (mocha.nhets > 0) {
						double f(double x, void *data)
						{
							if (model->flags & WGS_DATA)
								return -ad_lod(
									ad0 + a, ad1 + a,
									b + 1 - a, NULL,
									model->err_log_prb,
									self->stats.dispersion,
									x);
							else
								return -baf_lod(
									baf + a, b + 1 - a,
									NULL,
									model->err_log_prb,
									self->stats.dispersion,
									x);
						}
						double x;
						kmin_brent(f, 0.1, 0.2, NULL, KMIN_EPS, &x);
						mocha.bdev = fabsf((float)x);
					} else
						mocha.bdev = NAN;
					mocha_table->n++;
					hts_expand(mocha_t, mocha_table->n, mocha_table->m,
						   mocha_table->a);
					mocha_table->a[mocha_table->n - 1] = mocha;
					for (int j = a; j <= b; j++) {
						// TODO add other stuff here, like setting ldev
						// and bdev
						lrr[j] = NAN; // do not use the data again
						baf[j] = NAN; // do not use the data again
					}
				}
			}
		}
	}

	float *hs_arr = NULL;
	int n_hs = 0, m_hs = 0;
	for (int hmm_model = 0; hmm_model < 2; hmm_model++) {
		// select data to use from the contig, depending on which HMM model is being
		// used
		int last_p = 0, first_q = 0;
		int n_imap = 0;
		for (int i = 0; i < n; i++)
			if ((hmm_model == LRR_BAF && !isnan(lrr[i]))
			    || (hmm_model == BAF_PHASE && !isnan(baf[i]))) {
				if (pos[i] < cen_beg)
					last_p++;
				if (pos[i] < cen_end)
					first_q++;
				n_imap++;
				imap_arr[n_imap - 1] = i;
			}
		if (n_imap == 0)
			continue;

		// compute emission probabilities and Viterbi path according to HMM model
		int middle = 0;
		// TODO eliminate this redundancy
		if (hmm_model == LRR_BAF) {
			n_hs = model->bdev_lrr_baf_n;
			hts_expand(float, n_hs, m_hs, hs_arr);
			for (int i = 0; i < model->bdev_lrr_baf_n; i++) {
				hs_arr[i] = model->bdev_lrr_baf[i];
				if (model->bdev_lrr_baf[i] < 0.0f)
					middle++;
			}
		} else if (hmm_model == BAF_PHASE) {
			n_hs = model->bdev_baf_phase_n;
			hts_expand(float, n_hs, m_hs, hs_arr);
			for (int i = 0; i < model->bdev_baf_phase_n; i++)
				hs_arr[i] = model->bdev_baf_phase[i];
		}
		int8_t *path;
		float ret;
		int *beg = NULL, m_beg = 0, *end = NULL, m_end = 0, nseg;
		do {
			if (n_hs + (hmm_model == LRR_BAF ? n_hs : 0) > 50)
				error("Too many states being tested for the HMM\n");

			float *emis_log_lkl;
			if (model->flags & WGS_DATA) {
				emis_log_lkl =
					hmm_model == LRR_BAF
						? lrr_ad_emis_log_lkl(
							lrr, ad0, ad1, n_imap, imap_arr,
							model->err_log_prb, model->lrr_bias,
							model->lrr_hap2dip, self->adjlrr_sd,
							self->stats.dispersion, hs_arr, n_hs)
						: ad_phase_emis_log_lkl(
							ad0, ad1, gt_phase, n_imap, imap_arr,
							model->err_log_prb,
							self->stats.dispersion, hs_arr, n_hs);
			} else {
				emis_log_lkl =
					hmm_model == LRR_BAF
						? lrr_baf_emis_log_lkl(
							lrr, baf, n_imap, imap_arr,
							model->err_log_prb, model->lrr_bias,
							model->lrr_hap2dip, self->adjlrr_sd,
							self->stats.dispersion, hs_arr, n_hs)
						: baf_phase_emis_log_lkl(
							baf, gt_phase, n_imap, imap_arr,
							model->err_log_prb,
							self->stats.dispersion, hs_arr, n_hs);
			}
			path = log_viterbi_run(
				emis_log_lkl, n_imap, n_hs + (hmm_model == LRR_BAF ? n_hs : 0),
				model->xy_log_prb,
				hmm_model == LRR_BAF ? NAN : model->flip_log_prb,
				model->tel_log_prb, model->cen_log_prb, last_p, first_q);
			free(emis_log_lkl);

			if (hmm_model == LRR_BAF)
				for (int i = 0; i < n_imap; i++)
					if (path[i] > n_hs)
						path[i] = 0;
			ret = get_path_segs(path, hs_arr, n_imap, hmm_model, middle, &beg,
					    &m_beg, &end, &m_end, &nseg);

			if (ret) // two consecutive hidden states were used, hinting that
				 // testing of a middle state might be necessary
			{
				free(path);
				n_hs++;
				if (middle && ret < 0.0f)
					middle++;
				hts_expand(float, n_hs, m_hs, hs_arr);
				hs_arr[n_hs - 1] = ret;
				ks_introsort_float(n_hs, hs_arr);
			}
		} while (ret);

		// loop through all the segments called by the Viterbi algorithm
		for (int i = 0; i < nseg; i++) {
			// compute edges of the call
			int a = imap_arr[beg[i]];
			if (beg[i] == 0)
				while (a > 0 && ldev[a - 1] == 0 && bdev[a - 1] == 0)
					a--; // extend call towards p telomere
			int b = imap_arr[end[i]];
			if (end[i] == n_imap - 1)
				while (b < n - 1 && ldev[b + 1] == 0 && bdev[b + 1] == 0)
					b++; // extend call towards q telomere
			mocha.ldev = get_median(lrr + a, b + 1 - a, NULL);
			get_mocha_stats(pos, lrr, baf, gt_phase, n, a, b, cen_beg, cen_end,
					length, self->stats.baf_conc, &mocha);

			double f(double x, void *data)
			{
				if (model->flags & WGS_DATA)
					return -lrr_ad_lod(lrr + a, ad0 + a, ad1 + a,
							   mocha.nsites, NULL,
							   model->err_log_prb, model->lrr_bias,
							   model->lrr_hap2dip, self->adjlrr_sd,
							   self->stats.dispersion, x);
				else
					return -lrr_baf_lod(lrr + a, baf + a, mocha.nsites,
							    NULL, model->err_log_prb,
							    model->lrr_bias, model->lrr_hap2dip,
							    self->adjlrr_sd,
							    self->stats.dispersion, x);
			}
			double x, fx = kmin_brent(f, -0.15, 0.15, NULL, KMIN_EPS, &x);
			mocha.lod_lrr_baf = -(float)fx;

			if (hmm_model == LRR_BAF) {
				// here you need to check whether the call would have been
				// better with the phased HMM model
				int n_hets_imap = 0;
				for (int j = beg[i]; j <= end[i]; j++)
					if (!isnan(baf[imap_arr[j]])) {
						n_hets_imap++;
						hets_imap_arr[n_hets_imap - 1] = imap_arr[j];
					}
				// TODO here it needs to pass information about the centromeres
				if (model->flags & WGS_DATA) {
					mocha.lod_baf_phase = compare_wgs_models(
						ad0, ad1, gt_phase, n_hets_imap, hets_imap_arr,
						model->xy_log_prb, model->err_log_prb,
						model->flip_log_prb, model->tel_log_prb,
						self->stats.dispersion, model->bdev_baf_phase,
						model->bdev_baf_phase_n);
				} else {
					mocha.lod_baf_phase = compare_models(
						baf, gt_phase, n_hets_imap, hets_imap_arr,
						model->xy_log_prb, model->err_log_prb,
						model->flip_log_prb, model->tel_log_prb,
						self->stats.dispersion, model->bdev_baf_phase,
						model->bdev_baf_phase_n);
				}
				if (mocha.lod_baf_phase > mocha.lod_lrr_baf)
					continue;

				// compute bdev, if possible
				if (n_hets_imap > 0) {
					double f(double x, void *data)
					{
						if (model->flags & WGS_DATA)
							return -ad_lod(ad0, ad1, n_hets_imap,
								       hets_imap_arr,
								       model->err_log_prb,
								       self->stats.dispersion,
								       x);
						else
							return -baf_lod(
								baf, n_hets_imap, hets_imap_arr,
								model->err_log_prb,
								self->stats.dispersion, x);
					}
					fx = kmin_brent(f, 0.1, 0.2, NULL, KMIN_EPS, &x);
					mocha.bdev = fabsf((float)x);
				} else
					mocha.bdev = NAN;
				mocha.bdev_se = NAN;
				for (int j = 0; j < n_hets_imap; j++)
					bdev_phase[hets_imap_arr[j]] =
						(int8_t)SIGN(baf[hets_imap_arr[j]] - 0.5f);
			} else {
				// penalizes the LOD by the number of phase flips
				mocha.nflips = 0;
				for (int j = beg[i]; j < end[i]; j++)
					if (path[j] != path[j + 1])
						mocha.nflips++;

				double f(double x, void *data)
				{
					if (model->flags & WGS_DATA)
						return -ad_phase_lod(
							ad0, ad1, gt_phase, mocha.nhets,
							imap_arr + beg[i], path + beg[i],
							model->err_log_prb,
							self->stats.dispersion, x);
					else
						return -baf_phase_lod(
							baf, gt_phase, mocha.nhets,
							imap_arr + beg[i], path + beg[i],
							model->err_log_prb,
							self->stats.dispersion, x);
				}
				double x, fx = kmin_brent(f, 0.1, 0.2, NULL, KMIN_EPS, &x);
				mocha.bdev = fabsf((float)x);
				mocha.lod_baf_phase = -(float)fx
						      + (float)mocha.nflips
								* model->flip_log_prb
								* (float)M_LOG10E;

				for (int j = 0; j < mocha.nhets; j++)
					pbaf_arr[j] = (baf[imap_arr[beg[i] + j]] - 0.5f)
						      * (float)SIGN(path[beg[i] + j]);
				mocha.bdev_se = get_se_mean(pbaf_arr, mocha.nhets, NULL);
				for (int j = beg[i]; j <= end[i]; j++)
					bdev_phase[imap_arr[j]] =
						(int8_t)SIGN(path[j]) * gt_phase[imap_arr[j]];
			}

			mocha.type =
				mocha_type(mocha.ldev, mocha.ldev_se, mocha.bdev, mocha.nhets,
					   model->lrr_hap2dip, mocha.p_arm, mocha.q_arm);
			mocha.cf = mocha_cell_fraction(mocha.ldev, mocha.ldev_se, mocha.bdev,
						       mocha.nhets, mocha.type,
						       model->lrr_hap2dip);
			mocha_table->n++;
			hts_expand(mocha_t, mocha_table->n, mocha_table->m, mocha_table->a);
			mocha_table->a[mocha_table->n - 1] = mocha;

			// update information that will be stored in the output VCF and make
			// remaining sites missing
			for (int j = a; j <= b; j++) {
				if (ldev[j] == 0)
					ldev[j] = float_to_int16(mocha.ldev);
				if (bdev[j] == 0)
					bdev[j] = float_to_int16(mocha.bdev);
				lrr[j] = NAN; // do not use the data again
				baf[j] = NAN; // do not use the data again
			}
		}
		free(path);
		free(beg);
		free(end);
	}

	// clean up
	free(hs_arr);
	free(imap_arr);
	free(hets_imap_arr);
	free(pbaf_arr);
	free(pos);
	free(lrr);
	free(baf);
	memcpy(self->data_arr[LDEV], ldev, n * sizeof(int16_t));
	memcpy(self->data_arr[BDEV], bdev, n * sizeof(int16_t));
	memcpy(self->phase_arr, bdev_phase, n * sizeof(int8_t));
	free(ldev);
	free(bdev);
	free(bdev_phase);
}

// computes the medoid contig for LRR regression
// TODO weight the coefficients appropriately
static int get_medoid(const float *coeffs, int n, int order)
{
	int medoid_idx = -1;
	float prev = INFINITY;
	for (int i = 0; i < n; i++) {
		float next = 0.0f;
		for (int j = 0; j < n; j++) {
			if (i != j) {
				for (int k = 0; k <= order; k++) {
					next += fabsf(coeffs[i * (order + 1) + k]
						      - coeffs[j * (order + 1) + k]);
				}
			}
		}
		if (next < prev) {
			prev = next;
			medoid_idx = i;
		}
	}
	return medoid_idx;
}

// groups numbers in two separate distributions
static float get_lrr_cutoff(const float *v, int n)
{
	if (n <= 1)
		return NAN;
	float *w = (float *)malloc(n * sizeof(float));
	int j = 0;
	for (int i = 0; i < n; i++) {
		if (!isnan(v[i]))
			w[j++] = v[i];
	}
	if (j <= 1) {
		free(w);
		return NAN;
	}
	ks_introsort_float((size_t)j, w);

	// identify a reasonable initial split allowing for some outliers
	int k = j / 2;
	int d = (int)sqrtf((float)j) - 1;
	while (k > 1 && w[k - 1] - w[d] > w[j - 1 - d] - w[k - 1])
		k--;
	while (k < j && w[k] - w[d] < w[j - 1 - d] - w[k])
		k++;

	// run k-means clustering EM
	while (k > 0
	       && w[k - 1] - (w[(k - 1) / 2] + w[k / 2]) * 0.5f
			  > (w[(j + k - 1) / 2] + w[(j + k) / 2]) * 0.5f - w[k - 1])
		k--;
	while (k < j
	       && w[k] - (w[(k - 1) / 2] + w[k / 2]) * 0.5f
			  < (w[(j + k - 1) / 2] + w[(j + k) / 2]) * 0.5f - w[k])
		k++;

	float cutoff = (k > 0 && k < j) ? (w[k - 1] + w[k]) * 0.5f : NAN;
	free(w);
	return cutoff;
}

// this function computes several contig stats and then clears the contig data from the sample
static void sample_stats(sample_t *self, const model_t *model)
{
	int n = self->n;
	if (n == 0)
		return;
	self->nsites += n;

	int16_t *ad0 = self->data_arr[AD0];
	int16_t *ad1 = self->data_arr[AD1];
	float *lrr = (float *)malloc(n * sizeof(float));
	float *baf = (float *)malloc(n * sizeof(float));
	if (model->flags & WGS_DATA) {
		ad_to_lrr_baf(ad0, ad1, lrr, baf, n);
	} else {
		for (int i = 0; i < n; i++) {
			lrr[i] = int16_to_float(self->data_arr[LRR][i]);
			baf[i] = int16_to_float(self->data_arr[BAF][i]);
		}
	}
	int *imap_arr = (int *)malloc(n * sizeof(int));

	if (model->rid == model->genome_rules->x_rid) {
		int n_imap = 0;
		for (int i = 0; i < n; i++) {
			if (!isnan(baf[i]))
				self->nhets++;
			int pos = model->locus_arr[self->vcf_imap_arr[i]].pos;
			if (pos > model->genome_rules->x_nonpar_beg
			    && pos < model->genome_rules->x_nonpar_end
			    && (pos < model->genome_rules->x_xtr_beg
				|| pos > model->genome_rules->x_xtr_end)) {
				if (!isnan(baf[i]))
					self->x_nonpar_nhets++;
				n_imap++;
				imap_arr[n_imap - 1] = i;
			}
		}
		self->x_nonpar_lrr_median = get_median(lrr, n_imap, imap_arr);

		if (model->flags & WGS_DATA) {
			double f(double x, void *data)
			{
				return -lod_lkl_beta_binomial(ad0, ad1, n_imap, imap_arr, x);
			}
			double x;
			kmin_brent(f, 0.1, 0.2, NULL, KMIN_EPS,
				   &x); // dispersions above 0.5 are not allowed
			self->x_nonpar_dispersion = (float)x;
		} else {
			self->x_nonpar_dispersion = get_sample_sd(baf, n_imap, imap_arr);
		}
	} else if (model->rid == model->genome_rules->y_rid) {
		int n_imap = 0;
		for (int i = 0; i < n; i++) {
			if (!isnan(baf[i]))
				self->nhets++;
			int pos = model->locus_arr[self->vcf_imap_arr[i]].pos;
			if (pos > model->genome_rules->y_nonpar_beg
			    && pos < model->genome_rules->y_nonpar_end
			    && (pos < model->genome_rules->y_xtr_beg
				|| pos > model->genome_rules->y_xtr_end)) {
				n_imap++;
				imap_arr[n_imap - 1] = i;
			}
		}
		self->y_nonpar_lrr_median = get_median(lrr, n_imap, imap_arr);
	} else if (model->rid == model->genome_rules->mt_rid) {
		self->mt_lrr_median = get_median(lrr, n, NULL);
	} else {
		// expand arrays if necessary
		self->n_stats++;
		hts_expand(stats_t, self->n_stats, self->m_stats, self->stats_arr);

		if (model->flags & WGS_DATA) {
			double f(double x, void *data)
			{
				return -lod_lkl_beta_binomial(ad0, ad1, n, NULL, x);
			}
			double x;
			kmin_brent(f, 0.1, 0.2, NULL, KMIN_EPS,
				   &x); // dispersions above 0.5 are not allowed
			self->stats_arr[self->n_stats - 1].dispersion = (float)x;
		} else {
			self->stats_arr[self->n_stats - 1].dispersion =
				get_sample_sd(baf, n, NULL);
		}
		for (int i = 0; i < n; i++)
			if (!isnan(baf[i]))
				self->nhets++;
		self->stats_arr[self->n_stats - 1].lrr_median = get_median(lrr, n, NULL);
		self->stats_arr[self->n_stats - 1].lrr_sd = get_sample_sd(lrr, n, NULL);

		int conc, disc;
		get_baf_conc(baf, self->phase_arr, n, NULL, &conc, &disc);
		self->stats_arr[self->n_stats - 1].baf_conc =
			(float)conc / (float)(conc + disc);
		self->stats_arr[self->n_stats - 1].baf_auto =
			get_baf_auto_corr(baf, self->phase_arr, n, NULL);
		if (model->order_lrr_gc == 0) {
			self->stats_arr[self->n_stats - 1].coeffs[0] = get_median(lrr, n, NULL);
		}
		// performs polynomial regression for LRR
		else if (model->order_lrr_gc > 0) {
			float tss = get_tss(lrr, n);
			polyfit(lrr, model->gc_arr, n, self->vcf_imap_arr, model->order_lrr_gc,
				self->stats_arr[self->n_stats - 1].coeffs);
			adjust_lrr(lrr, model->gc_arr, n, self->vcf_imap_arr,
				   self->stats_arr[self->n_stats - 1].coeffs,
				   model->order_lrr_gc);
			self->stats_arr[self->n_stats - 1].coeffs[0] +=
				get_median(lrr, n, NULL); // further adjusts by median
			float rss = get_tss(lrr, n);
			self->stats_arr[self->n_stats - 1].rel_ess = 1.0f - rss / tss;
		}
		// compute autocorrelation after GC correction
		self->stats_arr[self->n_stats - 1].lrr_auto = get_lrr_auto_corr(lrr, n, NULL);
	}

	free(lrr);
	free(baf);
	free(imap_arr);
}

// this function computes the median of contig stats
static void sample_summary(sample_t *self, int n, model_t *model)
{
	float *tmp_arr = (float *)malloc(n * sizeof(float));
	int m_tmp = n;

	for (int i = 0; i < n; i++) {
		hts_expand(float,
			   self[i].n_stats *(model->order_lrr_gc < 0 ? 1
								     : model->order_lrr_gc + 1),
			   m_tmp, tmp_arr);

		for (int j = 0; j < self[i].n_stats; j++)
			tmp_arr[j] = self[i].stats_arr[j].lrr_median;
		self[i].stats.lrr_median = get_median(tmp_arr, self[i].n_stats, NULL);
		for (int j = 0; j < self[i].n_stats; j++)
			tmp_arr[j] = self[i].stats_arr[j].lrr_sd;
		self[i].stats.lrr_sd = get_median(tmp_arr, self[i].n_stats, NULL);
		for (int j = 0; j < self[i].n_stats; j++)
			tmp_arr[j] = self[i].stats_arr[j].lrr_auto;
		self[i].stats.lrr_auto = get_median(tmp_arr, self[i].n_stats, NULL);
		for (int j = 0; j < self[i].n_stats; j++)
			tmp_arr[j] = self[i].stats_arr[j].dispersion;
		self[i].stats.dispersion = get_median(tmp_arr, self[i].n_stats, NULL);
		for (int j = 0; j < self[i].n_stats; j++)
			tmp_arr[j] = self[i].stats_arr[j].baf_conc;
		self[i].stats.baf_conc = get_median(tmp_arr, self[i].n_stats, NULL);
		for (int j = 0; j < self[i].n_stats; j++)
			tmp_arr[j] = self[i].stats_arr[j].baf_auto;
		self[i].stats.baf_auto = get_median(tmp_arr, self[i].n_stats, NULL);

		self[i].adjlrr_sd = self[i].stats.lrr_sd;
		if (model->order_lrr_gc == 0) {
			for (int j = 0; j < self[i].n_stats; j++)
				tmp_arr[j] = self[i].stats_arr[j].coeffs[0];
			self[i].stats.coeffs[0] = get_median(tmp_arr, self[i].n_stats, NULL);
		} else if (model->order_lrr_gc > 0 && self[i].n_stats > 0) {
			for (int j = 0; j < self[i].n_stats; j++)
				for (int k = 0; k <= model->order_lrr_gc; k++)
					tmp_arr[j * (model->order_lrr_gc + 1) + k] =
						self[i].stats_arr[j].coeffs[k];
			int medoid_idx =
				get_medoid(tmp_arr, self[i].n_stats, model->order_lrr_gc);
			for (int k = 0; k <= model->order_lrr_gc; k++)
				self[i].stats.coeffs[k] =
					tmp_arr[medoid_idx * (model->order_lrr_gc + 1) + k];
			for (int j = 0; j < self[i].n_stats; j++)
				tmp_arr[j] = self[i].stats_arr[j].rel_ess;
			self[i].stats.rel_ess = get_median(tmp_arr, self[i].n_stats, NULL);
			self[i].adjlrr_sd *= sqrtf(
				1.0f
				- self[i].stats.rel_ess); // not perfect, but good enough(?)
		}
		free(self[i].stats_arr);
	}

	if (model->flags & WGS_DATA) {
		if (isnan(model->lrr_cutoff))
			model->lrr_cutoff = -0.3f; // arbitrary cutoff between -M_LN2 and 0
		if (isnan(model->lrr_auto2sex))
			model->lrr_auto2sex = 0.0f;
	}

	// determine LRR cutoff between haploid and diploid
	if (isnan(model->lrr_cutoff)) {
		int j = 0;
		for (int i = 0; i < n; i++)
			if (!isnan(self[i].x_nonpar_lrr_median))
				tmp_arr[j++] = isnan(self[i].stats.lrr_median)
						       ? self[i].x_nonpar_lrr_median
						       : self[i].x_nonpar_lrr_median
								 - self[i].stats.lrr_median;
		model->lrr_cutoff = get_lrr_cutoff(tmp_arr, j);
	}

	// determine sex of samples
	for (int i = 0; i < n; i++) {
		float tmp = isnan(self[i].stats.lrr_median)
				    ? self[i].x_nonpar_lrr_median
				    : self[i].x_nonpar_lrr_median - self[i].stats.lrr_median;
		if (tmp < model->lrr_cutoff)
			self[i].sex = SEX_MAL;
		else if (tmp > model->lrr_cutoff)
			self[i].sex = SEX_FEM;
	}

	// determine LRR difference between autosomes and sex chromosomes
	if (isnan(model->lrr_auto2sex)) {
		int j = 0;
		for (int i = 0; i < n; i++)
			if (self[i].sex == SEX_FEM)
				tmp_arr[j++] = isnan(self[i].stats.lrr_median)
						       ? self[i].x_nonpar_lrr_median
						       : self[i].x_nonpar_lrr_median
								 - self[i].stats.lrr_median;
		if (j == 0) // if no females are present it doesn't matter
		{
			model->lrr_auto2sex = 0.0f;
		} else {
			float lrr_females = get_median(tmp_arr, j, NULL);
			model->lrr_auto2sex = lrr_females;
		}
	}
	free(tmp_arr);
}

static void sample_print(FILE *restrict stream, const sample_t *self, int n,
			 const bcf_hdr_t *hdr, int flags)
{
	if (stream == NULL)
		return;
	char sex[3];
	sex[SEX_UNK] = 'U';
	sex[SEX_MAL] = 'M';
	sex[SEX_FEM] = 'F';
	if (flags & WGS_DATA) {
		fprintf(stream,
			"SAMPLE\tCOV_MEDIAN\tCOV_SD\tCOV_AUTO\tBAF_CORR\tBAF_CONC\tBAF_AUTO\tNSITES\tNHETS\tX_NONPAR_NHETS\tX_NONPAR_BAF_CORR\tX_NONPAR_COV_MEDIAN\tY_NONPAR_COV_MEDIAN\tMT_COV_MEDIAN\tSEX\tREL_ESS\n");
		for (int i = 0; i < n; i++) {
			const char *sample_name =
				bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, self[i].idx);
			fprintf(stream,
				"%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%c\t%.4f\n",
				sample_name, expf(self[i].stats.lrr_median),
				expf(self[i].stats.lrr_median) * self[i].stats.lrr_sd,
				self[i].stats.lrr_auto, self[i].stats.dispersion,
				self[i].stats.baf_conc, self[i].stats.baf_auto, self[i].nsites,
				self[i].nhets, self[i].x_nonpar_nhets,
				self[i].x_nonpar_dispersion, expf(self[i].x_nonpar_lrr_median),
				expf(self[i].y_nonpar_lrr_median), expf(self[i].mt_lrr_median),
				sex[self[i].sex], self[i].stats.rel_ess);
		}
	} else {
		fprintf(stream,
			"SAMPLE\tLRR_MEDIAN\tLRR_SD\tLRR_AUTO\tBAF_SD\tBAF_CONC\tBAF_AUTO\tNSITES\tNHETS\tX_NONPAR_NHETS\tX_NONPAR_BAF_SD\tX_NONPAR_LRR_MEDIAN\tY_NONPAR_LRR_MEDIAN\tMT_LRR_MEDIAN\tSEX\tREL_ESS\n");
		for (int i = 0; i < n; i++) {
			const char *sample_name =
				bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, self[i].idx);
			fprintf(stream,
				"%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%c\t%.4f\n",
				sample_name, self[i].stats.lrr_median, self[i].stats.lrr_sd,
				self[i].stats.lrr_auto, self[i].stats.dispersion,
				self[i].stats.baf_conc, self[i].stats.baf_auto, self[i].nsites,
				self[i].nhets, self[i].x_nonpar_nhets,
				self[i].x_nonpar_dispersion, self[i].x_nonpar_lrr_median,
				self[i].y_nonpar_lrr_median, self[i].mt_lrr_median,
				sex[self[i].sex], self[i].stats.rel_ess);
		}
	}
	if (stream != stdout && stream != stderr)
		fclose(stream);
}

/*********************************
 * VCF READ AND WRITE METHODS    *
 *********************************/

// moves synced bcf reader to first useful line for reader0
static inline int bcf_sr_next_line_reader0(bcf_srs_t *sr)
{
	int nret = bcf_sr_next_line(sr);
	while (nret > 0 && !bcf_sr_has_line(sr, 0))
		nret = bcf_sr_next_line(sr);
	return nret;
}

// write one contig
static int put_contig(bcf_srs_t *sr, const sample_t *sample, const model_t *model,
		      htsFile *out_fh, bcf_hdr_t *out_hdr)
{
	int rid = model->rid;
	bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
	bcf_sr_seek(sr, bcf_hdr_id2name(hdr, rid), 0);
	int nsmpl = bcf_hdr_nsamples(out_hdr);

	int *synced_iter = (int *)calloc(nsmpl, sizeof(int));
	float *bdev = (float *)calloc(nsmpl, sizeof(float));
	float *ldev = (float *)calloc(nsmpl, sizeof(float));
	int *bdev_phase = (int *)calloc(nsmpl, sizeof(int));

	int i;
	for (i = 0; bcf_sr_next_line_reader0(sr); i++) {
		bcf1_t *line = bcf_sr_get_line(sr, 0);
		if (rid != line->rid)
			break;

		for (int j = 0; j < nsmpl; j++) {
			while (synced_iter[j] < sample[j].n - 1
			       && sample[j].vcf_imap_arr[synced_iter[j]] < i)
				synced_iter[j]++;
			if (sample[j].vcf_imap_arr[synced_iter[j]] == i) {
				if (sample[j].data_arr[LDEV])
					ldev[sample[j].idx] = int16_to_float(
						sample[j].data_arr[LDEV][synced_iter[j]]);
				if (sample[j].data_arr[BDEV])
					bdev[sample[j].idx] = int16_to_float(
						sample[j].data_arr[BDEV][synced_iter[j]]);
				if (sample[j].phase_arr)
					bdev_phase[sample[j].idx] =
						sample[j].phase_arr[synced_iter[j]];
			} else {
				// if no match variant found, match the end of the contig or
				// keep conservative
				if (i == 0 && sample[j].data_arr[BDEV])
					bdev[sample[j].idx] =
						int16_to_float(sample[j].data_arr[BDEV][0]);
				if (i == 0 && sample[j].data_arr[LDEV])
					ldev[sample[j].idx] =
						int16_to_float(sample[j].data_arr[LDEV][0]);
				if (sample[j].data_arr[BDEV]
				    && int16_to_float(sample[j].data_arr[BDEV][synced_iter[j]])
					       == 0.0f)
					bdev[sample[j].idx] = 0.0f;
				if (sample[j].data_arr[LDEV]
				    && int16_to_float(sample[j].data_arr[LDEV][synced_iter[j]])
					       == 0.0f)
					ldev[sample[j].idx] = 0.0f;
				if (sample[j].phase_arr)
					bdev_phase[sample[j].idx] = 0;
			}
		}
		if (!(model->flags & WGS_DATA)) {
			bcf_update_info_float(out_hdr, line, "ADJUST_BAF_LRR",
					      &model->locus_arr[i].adjust, 4);
		}
		if (!(model->flags & NO_ANNOT)) {
			bcf_update_format_float(out_hdr, line, "Ldev", ldev, (int)nsmpl);
			bcf_update_format_float(out_hdr, line, "Bdev", bdev, (int)nsmpl);
		}
		bcf_update_format_int32(out_hdr, line, "Bdev_Phase", bdev_phase, (int)nsmpl);

		if (bcf_write(out_fh, out_hdr, line) < 0)
			error("Unable to write to output VCF file\n");
	}

	free(synced_iter);
	free(ldev);
	free(bdev);
	free(bdev_phase);

	return i;
}

// write header
static bcf_hdr_t *print_hdr(htsFile *out_fh, bcf_hdr_t *hdr, int argc, char *argv[],
			    int record_cmd_line, int flags)
{
	bcf_hdr_t *out_hdr = bcf_hdr_dup(hdr);
	if (!(flags & WGS_DATA) && bcf_hdr_id2int(out_hdr, BCF_DT_ID, "ADJUST_BAF_LRR") < 0)
		bcf_hdr_append(
			out_hdr,
			"##INFO=<ID=ADJUST_BAF_LRR,Number=4,Type=Float,Description=\"Adjust values for BAF and LRR\">");
	if (!(flags & NO_ANNOT)) {
		if (bcf_hdr_id2int(out_hdr, BCF_DT_ID, "Ldev") < 0)
			bcf_hdr_append(
				out_hdr,
				"##FORMAT=<ID=Ldev,Number=1,Type=Float,Description=\"LRR deviation due to chromosomal alteration\">");
		if (bcf_hdr_id2int(out_hdr, BCF_DT_ID, "Bdev") < 0)
			bcf_hdr_append(
				out_hdr,
				"##FORMAT=<ID=Bdev,Number=1,Type=Float,Description=\"BAF deviation due to chromosomal alteration\">");
	}
	if (bcf_hdr_id2int(out_hdr, BCF_DT_ID, "Bdev_Phase") < 0)
		bcf_hdr_append(
			out_hdr,
			"##FORMAT=<ID=Bdev_Phase,Number=1,Type=Integer,Description=\"BAF deviation phase, if available\">");
	if (record_cmd_line)
		bcf_hdr_append_version(out_hdr, argc, argv, "bcftools_mocha");
	if (bcf_hdr_write(out_fh, out_hdr) < 0)
		error("Unable to write to output VCF file\n");
	return out_hdr;
}

// retrieve genotypes as NC, AA, AB, BB
// assumes little endian architecture
static int bcf_get_ab_genotypes(bcf_fmt_t *fmt, int8_t *gts, int nsmpl, int allele_a,
				int allele_b)
{
	if (!fmt || fmt->n != 2)
		return 0;

#define BRANCH(type_t, bcf_type_vector_end)                                                    \
	{                                                                                      \
		type_t *p = (type_t *)fmt->p;                                                  \
		for (int i = 0; i < nsmpl; i++, p += 2) {                                      \
			if (p[0] == bcf_type_vector_end || bcf_gt_is_missing(p[0])             \
			    || p[1] == bcf_type_vector_end || bcf_gt_is_missing(p[1])) {       \
				gts[i] = GT_NC;                                                \
			} else {                                                               \
				type_t allele_0 = bcf_gt_allele(p[0]);                         \
				type_t allele_1 = bcf_gt_allele(p[1]);                         \
				if (allele_0 == allele_a && allele_1 == allele_a)              \
					gts[i] = GT_AA;                                        \
				else if (allele_0 == allele_b && allele_1 == allele_b)         \
					gts[i] = GT_BB;                                        \
				else if ((allele_0 == allele_a && allele_1 == allele_b)        \
					 || (allele_0 == allele_b && allele_1 == allele_a))    \
					gts[i] = GT_AB;                                        \
				else                                                           \
					return -1;                                             \
			}                                                                      \
		}                                                                              \
	}
	switch (fmt->type) {
	case BCF_BT_INT8:
		BRANCH(int8_t, bcf_int8_vector_end);
		break;
	case BCF_BT_INT16:
		BRANCH(int16_t, bcf_int16_vector_end);
		break;
	case BCF_BT_INT32:
		BRANCH(int32_t, bcf_int32_vector_end);
		break;
	default:
		error("Unexpected type %d\n", fmt->type);
	}
#undef BRANCH

	return 1;
}

// check whether the BAF is flipped at homozygous sites
// if the BAF is completely flipped, complement the BAF values
// return -1 in case of error, 0 in case of an okay site, and 1 in case of a flipped site
static int bcf_check_baf_flipped(const int8_t *gts, const bcf_fmt_t *baf_fmt, int nsmpl)
{
	if (!baf_fmt || baf_fmt->n != 1 || baf_fmt->type != BCF_BT_FLOAT)
		return -1;

	int okay = 0, flipped = 0;

	float *p = (float *)baf_fmt->p;
	for (int i = 0; i < nsmpl; i++, p++) {
		if (gts[i] == GT_NC)
			continue;
		else if (gts[i] == GT_AA) {
			okay += p[0] < .5f;
			flipped += p[0] > .5f;
		} else if (gts[i] == GT_BB) {
			okay += p[0] > .5f;
			flipped += p[0] < .5f;
		}
	}

	if (okay > 0 && flipped > 0)
		return -1;
	else if (flipped > 0) {
		p = (float *)baf_fmt->p;
		for (int i = 0; i < nsmpl; i++, p++) {
			p[0] = 1.0f - p[0];
		}
		return 1;
	}

	return 0;
}

// read one contig
static void get_contig(bcf_srs_t *sr, sample_t *sample, model_t *model)
{
	int rid = model->rid;
	bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
	bcf_sr_seek(sr, bcf_hdr_id2name(hdr, rid), 0);

	bcf_fmt_t *baf_fmt = NULL, *lrr_fmt = NULL;
	bcf_info_t *info;
	int nsmpl = bcf_hdr_nsamples(hdr);

	int i;
	int allele_a_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ALLELE_A");
	if (allele_a_id >= 0 && bcf_hdr_id2type(hdr, BCF_HL_INFO, allele_a_id) != BCF_HT_INT)
		error("Error: input VCF file ALLELE_A info field is not of integer type\n");
	int allele_b_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ALLELE_B");
	if (allele_b_id >= 0 && bcf_hdr_id2type(hdr, BCF_HL_INFO, allele_b_id) != BCF_HT_INT)
		error("Error: input VCF file ALLELE_B info field is not of integer type\n");
	int gc_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GC");
	if (gc_id >= 0 && bcf_hdr_id2type(hdr, BCF_HL_INFO, gc_id) != BCF_HT_REAL)
		error("Error: input VCF file GC info field is not of float type\n");
	int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
	int baf_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "BAF");
	if (baf_id >= 0 && bcf_hdr_id2type(hdr, BCF_HL_FMT, baf_id) != BCF_HT_REAL)
		error("Error: input VCF file BAF format field is not of float type\n");
	int lrr_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "LRR");
	if (lrr_id >= 0 && bcf_hdr_id2type(hdr, BCF_HL_FMT, lrr_id) != BCF_HT_REAL)
		error("Error: input VCF file LRR format field is not of float type\n");
	int ad_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AD");
	if (ad_id >= 0 && bcf_hdr_id2type(hdr, BCF_HL_FMT, ad_id) != BCF_HT_INT)
		error("Error: input VCF file AD format field is not of integer type\n");

	// maybe this should be assert instead
	if (gt_id < 0)
		error("Error: input VCF file has no GT format field\n");
	if ((model->flags & WGS_DATA) && ad_id < 0)
		error("Error: input VCF file has no AD format field\n");
	if (!(model->flags & WGS_DATA) && (baf_id < 0 || lrr_id < 0))
		error("Error: input VCF file has no BAF or LRR format field\n");

	int8_t *gts = (int8_t *)malloc(nsmpl * sizeof(int8_t));
	int8_t *phase_arr = (int8_t *)malloc(nsmpl * sizeof(int8_t));
	float *float_arr = (float *)malloc(nsmpl * sizeof(float));
	int16_t *gt0 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
	int16_t *gt1 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
	int16_t *ad0 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
	int16_t *ad1 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
	int *last_het_pos = (int *)calloc(nsmpl, sizeof(int));
	int *last_pos = (int *)calloc(nsmpl, sizeof(int));

	model->n = 0;
	model->n_flipped = 0;
	for (int j = 0; j < nsmpl; j++)
		sample[j].n = 0;

	for (i = 0; bcf_sr_next_line_reader0(sr); i++) {
		bcf1_t *line = bcf_sr_get_line(sr, 0);
		if (rid != line->rid)
			break;
		int pos = line->pos + 1;

		hts_expand(locus_t, i + 1, model->m_locus, model->locus_arr);
		model->locus_arr[i].pos = pos;
		if (allele_a_id >= 0 && allele_b_id >= 0) {
			if ((info = bcf_get_info_id(line, allele_a_id)))
				model->locus_arr[i].allele_a = info->v1.i;
			else
				error("Error: ALLELE_A missing at position %s:%" PRId64 "\n",
				      bcf_hdr_id2name(hdr, line->rid), line->pos + 1);
			if ((info = bcf_get_info_id(line, allele_b_id)))
				model->locus_arr[i].allele_b = info->v1.i;
			else
				error("Error: ALLELE_B missing at position %s:%" PRId64 "\n",
				      bcf_hdr_id2name(hdr, line->rid), line->pos + 1);
		} else {
			model->locus_arr[i].allele_a = 0;
			model->locus_arr[i].allele_b = 1;
		}
		hts_expand(float, i + 1, model->m_gc, model->gc_arr);
		if (gc_id >= 0 && (info = bcf_get_info_id(line, gc_id)))
			model->gc_arr[i] = info->v1.f;
		else
			model->gc_arr[i] = NAN;

		// if failing inclusion/exclusion requirement, skip line
		if ((model->flags & FLT_EXCLUDE) && bcf_sr_get_line(sr, 1))
			continue;
		if ((model->flags & FLT_INCLUDE) && !bcf_sr_get_line(sr, 1))
			continue;

		// if site falls in short arm or centromere regions skip line
		if (!(model->flags & USE_SHORT_ARMS) && model->genome_rules->is_short_arm[rid]
		    && pos < model->genome_rules->cen_beg[rid])
			continue;
		if (!(model->flags & USE_CENTROMERES) && pos > model->genome_rules->cen_beg[rid]
		    && pos < model->genome_rules->cen_end[rid])
			continue;

		// if there are no genotypes, skip line
		bcf_fmt_t *gt_fmt = bcf_get_fmt_id(line, gt_id);
		if (bcf_get_ab_genotypes(gt_fmt, gts, nsmpl, model->locus_arr[i].allele_a,
					 model->locus_arr[i].allele_b)
		    < 0)
			error("Error: site %s:%" PRId64 " contains non-AB alleles\n",
			      bcf_hdr_id2name(hdr, line->rid), line->pos + 1);
		if (!bcf_get_genotype_phase(gt_fmt, phase_arr, nsmpl))
			continue;

		// if neither AD nor LRR and BAF formats are present, skip line
		if (model->flags & WGS_DATA) {
			if (!bcf_get_genotype_alleles(gt_fmt, gt0, gt1, nsmpl))
				continue;
			if (!bcf_get_allelic_depth(bcf_get_fmt_id(line, ad_id), gt0, gt1, ad0,
						   ad1, nsmpl))
				continue;
		} else {
			if (!(lrr_fmt = bcf_get_fmt_id(line, lrr_id))
			    || !(baf_fmt = bcf_get_fmt_id(line, baf_id)))
				continue;

			for (int j = 0; j < nsmpl; j++) {
				if (bcf_float_is_missing(
					    ((float *)(lrr_fmt->p
						       + lrr_fmt->size * sample[j].idx))[0]))
					((float *)(lrr_fmt->p
						   + lrr_fmt->size * sample[j].idx))[0] = NAN;
				if (bcf_float_is_missing(
					    ((float *)(baf_fmt->p
						       + baf_fmt->size * sample[j].idx))[0]))
					((float *)(baf_fmt->p
						   + baf_fmt->size * sample[j].idx))[0] = NAN;
			}

			// check whether you need to flip the BAF but only for bi-allelic sites
			if (!(model->flags & NO_BAF_FLIP)) {
				int ret = bcf_check_baf_flipped(gts, baf_fmt, nsmpl);
				if (ret < 0)
					error("Error: site %s:%" PRId64
					      " contains a mix of flipped and unflipped BAF values\n"
					      "Use bcftools query -f \"[%%CHROM\\t%%POS\\t%%SAMPLE\\t%%GT\\t%%BAF\\n]\" -r %s:%" PRId64
					      "-%" PRId64
					      " to investigate the issue\n"
					      "Use --no-BAF-flip to suppress BAF flipping\n",
					      bcf_hdr_id2name(hdr, line->rid), line->pos + 1,
					      bcf_hdr_id2name(hdr, line->rid), line->pos,
					      line->pos + 1);
				else
					model->n_flipped += ret;
			}

			// if allele A index is bigger than allele B index flip the BAF to make
			// sure it refers to the highest allele
			if (model->locus_arr[i].allele_a > model->locus_arr[i].allele_b) {
				for (int j = 0; j < nsmpl; j++)
					((float *)(baf_fmt->p
						   + baf_fmt->size * sample[j].idx))[0] =
						1.0f
						- ((float *)(baf_fmt->p
							     + baf_fmt->size
								       * sample[j].idx))[0];
			}

			// adjust BAF and LRR
			for (int j = 0; j < 4; j++)
				model->locus_arr[i].adjust[j] = 0.0f;

			int is_x_nonpar = rid == model->genome_rules->x_rid
					  && pos > model->genome_rules->x_nonpar_beg
					  && pos < model->genome_rules->x_nonpar_end
					  && (pos < model->genome_rules->x_xtr_beg
					      || pos > model->genome_rules->x_xtr_end);
			int is_y_or_mt = rid == model->genome_rules->y_rid
					 || rid == model->genome_rules->mt_rid;

			// compute BAF cluster center for heterozygous genotypes
			int k = 0;
			for (int j = 0; j < nsmpl; j++)
				if (gts[j] != GT_AB || (is_x_nonpar && sample[j].sex == SEX_MAL)
				    || is_y_or_mt)
					((float *)(baf_fmt->p
						   + baf_fmt->size * sample[j].idx))[0] = NAN;
				else
					float_arr[k++] =
						((float *)(baf_fmt->p
							   + baf_fmt->size * sample[j].idx))[0];
			if (model->median_baf_adj >= 0 && k >= model->median_baf_adj)
				model->locus_arr[i].adjust[0] =
					get_median(float_arr, k, NULL) - 0.5f;

			// compute LRR cluster center for each genotype
			if (model->rid != model->genome_rules->x_rid
			    && model->rid != model->genome_rules->y_rid) {
				for (int gt = GT_AA; gt <= GT_BB; gt++) {
					k = 0;
					for (int j = 0; j < nsmpl; j++)
						if (gts[j] == gt
						    && !(is_x_nonpar
							 && sample[j].sex == SEX_MAL)
						    && !is_y_or_mt)
							float_arr[k++] = ((
								float *)(lrr_fmt->p
									 + lrr_fmt->size
										   * sample[j]
											     .idx))
								[0];
					if (model->median_lrr_adj >= 0
					    && k >= model->median_lrr_adj)
						model->locus_arr[i].adjust[gt] =
							get_median(float_arr, k, NULL);
				}
			}

			for (int j = 0; j < nsmpl; j++) {
				((float *)(baf_fmt->p + baf_fmt->size * sample[j].idx))[0] -=
					model->locus_arr[i].adjust[0];
				if (gts[j])
					((float *)(lrr_fmt->p
						   + lrr_fmt->size * sample[j].idx))[0] -=
						model->locus_arr[i].adjust[gts[j]];
			}
		}

		// read line in memory
		model->n++;
		for (int j = 0; j < nsmpl; j++) {
			if (model->flags & WGS_DATA) {
				if (ad0[j] == bcf_int16_missing && ad1[j] == bcf_int16_missing)
					continue;

				// site too close to last het site or hom site too close to last
				// site
				if ((pos < last_het_pos[j] + model->min_dst)
				    || ((phase_arr[j] == bcf_int8_missing)
					&& (pos < last_pos[j] + model->min_dst)))
					continue;

				// substitute the last hom site with the current het site
				if (pos < last_pos[j] + model->min_dst)
					sample[j].n--;

				if (phase_arr[j] != bcf_int8_missing)
					last_het_pos[j] = pos;
				last_pos[j] = pos;

				sample[j].n++;
				hts_expand(int, sample[j].n, sample[j].m_vcf_imap,
					   sample[j].vcf_imap_arr);
				hts_expand(int8_t, sample[j].n, sample[j].m_phase,
					   sample[j].phase_arr);
				hts_expand(int16_t, sample[j].n, sample[j].m_data[AD0],
					   sample[j].data_arr[AD0]);
				hts_expand(int16_t, sample[j].n, sample[j].m_data[AD1],
					   sample[j].data_arr[AD1]);
				sample[j].vcf_imap_arr[sample[j].n - 1] = i;
				sample[j].phase_arr[sample[j].n - 1] = phase_arr[j];
				sample[j].data_arr[AD0][sample[j].n - 1] = ad0[j];
				sample[j].data_arr[AD1][sample[j].n - 1] = ad1[j];

			} else {
				float lrr = ((float *)(lrr_fmt->p
						       + lrr_fmt->size * sample[j].idx))[0];
				float baf = ((float *)(baf_fmt->p
						       + baf_fmt->size * sample[j].idx))[0];
				if (lrr == NAN && baf == NAN)
					continue;

				sample[j].n++;
				hts_expand(int, sample[j].n, sample[j].m_vcf_imap,
					   sample[j].vcf_imap_arr);
				hts_expand(int8_t, sample[j].n, sample[j].m_phase,
					   sample[j].phase_arr);
				hts_expand(int16_t, sample[j].n, sample[j].m_data[LRR],
					   sample[j].data_arr[LRR]);
				hts_expand(int16_t, sample[j].n, sample[j].m_data[BAF],
					   sample[j].data_arr[BAF]);
				sample[j].vcf_imap_arr[sample[j].n - 1] = i;
				sample[j].phase_arr[sample[j].n - 1] = phase_arr[j];
				sample[j].data_arr[LRR][sample[j].n - 1] = float_to_int16(lrr);
				sample[j].data_arr[BAF][sample[j].n - 1] = float_to_int16(baf);
			}
		}
	}
	free(gts);
	free(phase_arr);
	free(float_arr);
	free(gt0);
	free(gt1);
	free(ad0);
	free(ad1);
	free(last_het_pos);
	free(last_pos);
}

/*********************************
 * PLUGIN CODE                   *
 *********************************/

const char *about(void)
{
	return "Runs SNP-based method for detection of MOsaic CHromosomal Alterations (MoChA).\n";
}

static const char *usage_text(void)
{
	return "\n"
	       "About:   MOsaic CHromosomal Alterations caller, requires phased genotypes (GT)\n"
	       "         and either B-allele frequency (BAF) and Log R Ratio intensity (LRR)\n"
	       "         or allelic depth coverage (AD). (version " MOCHA_VERSION
	       " https://github.com/freeseek/mocha)\n"
	       "Usage:   bcftools +mocha [OPTIONS] <in.vcf>\n"
	       "\n"
	       "Required options:\n"
	       "    -r, --rules <assembly>[?]         predefined genome reference rules, 'list' to print available settings, append '?' for details\n"
	       "    -R, --rules-file <file>           genome reference rules, space/tab-delimited CHROM:FROM-TO,TYPE\n"
	       "\n"
	       "General Options:\n"
	       "    -s, --samples [^]<list>           comma separated list of samples to include (or exclude with \"^\" prefix)\n"
	       "    -S, --samples-file [^]<file>      file of samples to include (or exclude with \"^\" prefix)\n"
	       "        --force-samples               only warn about unknown subset samples\n"
	       "    -t, --targets [^]<region>         restrict to comma-separated list of regions. Exclude regions with \"^\" prefix\n"
	       "    -T, --targets-file [^]<file>      restrict to regions listed in a file. Exclude regions with \"^\" prefix\n"
	       "    -f, --apply-filters <list>        require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n"
	       "    -v, --variants [^]<file>          tabix-indexed [compressed] VCF/BCF file containing variants\n"
	       "                                      to include (or exclude with \"^\" prefix) in the analysis\n"
	       "        --threads <int>               number of extra output compression threads [0]\n"
	       "\n"
	       "Output Options:\n"
	       "    -o, --output <file>               write output to a file [no output]\n"
	       "    -O, --output-type <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n"
	       "        --no-version                  do not append version and command line to the header\n"
	       "    -a  --no-annotations              omit Ldev and Bdev FORMAT from output VCF (requires --output)\n"
	       "        --no-log                      suppress progress report on standard error\n"
	       "    -l  --log <file>                  write log to file [standard error]\n"
	       "    -m, --mosaic-calls <file>         write mosaic chromosomal alterations to a file [standard output]\n"
	       "    -g, --genome-stats <file>         write sample genome-wide statistics to a file [no output]\n"
	       "    -u, --ucsc-bed <file>             write UCSC bed track to a file [no output]\n"
	       "\n"
	       "HMM Options:\n"
	       "    -p  --cnp <file>                  list of regions to genotype in BED format\n"
	       "        --bdev-LRR-BAF <list>         comma separated list of inverse BAF deviations for LRR+BAF model [" BDEV_LRR_BAF_DFLT
	       "]\n"
	       "        --bdev-BAF-phase <list>       comma separated list of inverse BAF deviations for BAF+phase model\n"
	       "                                      [" BDEV_BAF_PHASE_DFLT
	       "]\n"
	       "    -d, --min-dist <int>              minimum base pair distance between consecutive sites for WGS data [" MIN_DST_DFLT
	       "]\n"
	       "        --no-BAF-flip                 do not correct BAF at flipped sites\n"
	       "        --median-BAF-adjust <int>     minimum number of heterozygous genotypes required to perform\n"
	       "                                      median BAF adjustment (-1 for no across samples BAF adjustment) [" MEDIAN_BAF_ADJ_DFLT
	       "]\n"
	       "        --median-LRR-adjust <int>     minimum number of same type genotypes required to perform\n"
	       "                                      median LRR adjustment (-1 for no across samples LRR adjustment) [" MEDIAN_LRR_ADJ_DFLT
	       "]\n"
	       "        --order-LRR-GC <int>          order of polynomial in local GC content to be used for polynomial\n"
	       "                                      regression of LRR (-1 for no LRR adjustment) [" ORDER_LRR_GC_DFLT
	       "]\n"
	       "        --xy-prob <float>             transition probability [" XY_PRB_DFLT
	       "]\n"
	       "        --err-prob <float>            uniform error probability [" ERR_PRB_DFLT
	       "]\n"
	       "        --flip-prob <float>           phase flip probability [" FLIP_PRB_DFLT
	       "]\n"
	       "        --telomere-advantage <float>  telomere advantage [" TEL_PRB_DFLT
	       "]\n"
	       "        --centromere-penalty <float>  centromere penalty [" CEN_PRB_DFLT
	       "]\n"
	       "        --short_arm_chrs <list>       list of chromosomes with short arms [" SHORT_ARM_CHRS_DFLT
	       "]\n"
	       "        --use_short_arms              use variants in short arms\n"
	       "        --use_centromeres             use variants in centromeres\n"
	       "        --LRR-cutoff <float>          LRR cutoff between haploid and diploid [estimated from X nonPAR]\n"
	       "        --LRR-hap2dip <float>         LRR difference between haploid and diploid [" LRR_HAP2DIP_DFLT
	       "]\n"
	       "        --LRR-auto2sex <float>        LRR difference between autosomes and diploid sex chromosomes [estimated from X nonPAR]\n"
	       "        --LRR-weight <float>          relative contribution from LRR for LRR+BAF model [" SHORT_ARM_CHRS_DFLT
	       "]\n"
	       "\n"
	       "Examples:\n"
	       "    bcftools +mocha -r GRCh37 input.bcf -v ^exclude.bcf -g stats.tsv -m mocha.tsv -p cnp.grch37.bed\n"
	       "    bcftools +mocha -r GRCh38 input.bcf -Ob -o output.bcf -g stats.tsv -m mocha.tsv -c 1.0 --LRR-weight 0.5\n"
	       "\n";
}

static float *read_list_invf(const char *str, int *n, float min, float max)
{
	char *tmp, **list = hts_readlist(str, 0, n);
	if (*n >= 128)
		error("Cannot handle list of 128 or more parameters: %s\n", str);
	float *ret = (float *)malloc(*n * sizeof(float));
	for (int i = 0; i < *n; i++) {
		ret[i] = 1.0f / strtof(list[i], &tmp);
		if (*tmp)
			error("Could not parse: %s\n", list[i]);
		if (ret[i] < min || ret[i] > max)
			error("Expected values from the interval [%f,%f], found 1/%s\n", min,
			      max, list[i]);
		free(list[i]);
	}
	free(list);
	ks_introsort_float((size_t)*n, ret);
	return ret;
}

static int cnp_parse(const char *line, char **chr_beg, char **chr_end, uint32_t *beg,
		     uint32_t *end, void *payload, void *usr)
{
	// Use the standard parser for CHROM,FROM,TO
	int ret = regidx_parse_bed(line, chr_beg, chr_end, beg, end, NULL, NULL);
	if (ret != 0)
		return ret;

	// Skip the fields that were parsed above
	char *ss = (char *)line;
	while (*ss && isspace(*ss))
		ss++;
	for (int i = 0; i < 3; i++) {
		while (*ss && !isspace(*ss))
			ss++;
		if (!*ss)
			return -2; // wrong number of fields
		while (*ss && isspace(*ss))
			ss++;
	}
	if (!*ss)
		return -2;

	// Parse the payload
	int *dat = (int *)payload;
	if (strncmp(ss, "DEL", 3) == 0)
		*dat = MOCHA_CNP_DEL;
	else if (strncmp(ss, "DUP", 3) == 0)
		*dat = MOCHA_CNP_DUP;
	else if (strncmp(ss, "CNV", 3) == 0)
		*dat = MOCHA_CNP_CNV;
	else
		*dat = MOCHA_UNK;
	return 0;
}

static FILE *get_file_handle(const char *str)
{
	FILE *ret;
	if (strcmp(str, "-") == 0)
		ret = stdout;
	else {
		ret = fopen(str, "w");
		if (!ret)
			error("Failed to open %s: %s\n", str, strerror(errno));
	}
	return ret;
}

int run(int argc, char *argv[])
{
	beta_binom_null = beta_binom_init();
	beta_binom_alt = beta_binom_init();

	// program options
	char *tmp = NULL;
	int rules_is_file = 0;
	int sample_is_file = 0;
	int targets_is_file = 0;
	int force_samples = 0;
	int output_type = FT_VCF;
	int n_threads = 0;
	int record_cmd_line = 1;
	char *sample_names = NULL;
	char *targets_list = NULL;
	char *output_fname = NULL;
	char *cnp_fname = NULL;
	char *filter_fname = NULL;
	char *rules = NULL;
	sample_t *sample = NULL;
	FILE *log_file = stderr;
	FILE *out_fm = stdout;
	FILE *out_fg = NULL;
	FILE *out_fu = NULL;
	bcf_hdr_t *hdr = NULL;
	bcf_hdr_t *out_hdr = NULL;
	htsFile *out_fh = NULL;
	mocha_table_t mocha_table = {0, 0, NULL};
	const char *short_arm_chrs = SHORT_ARM_CHRS_DFLT;
	const char *bdev_lrr_baf = BDEV_LRR_BAF_DFLT;
	const char *bdev_baf_phase = BDEV_BAF_PHASE_DFLT;

	// model parameters
	model_t model;
	memset(&model, 0, sizeof(model_t));
	model.xy_log_prb = logf(strtof(XY_PRB_DFLT, &tmp));
	model.err_log_prb = logf(strtof(ERR_PRB_DFLT, &tmp));
	model.flip_log_prb = logf(strtof(FLIP_PRB_DFLT, &tmp));
	model.tel_log_prb = logf(strtof(TEL_PRB_DFLT, &tmp));
	model.cen_log_prb = logf(strtof(CEN_PRB_DFLT, &tmp));
	model.min_dst = (int)strtol(MIN_DST_DFLT, &tmp, 0);
	model.lrr_bias = strtof(LRR_BIAS_DFLT, &tmp);
	model.lrr_cutoff = NAN;
	model.lrr_hap2dip = strtof(LRR_HAP2DIP_DFLT, &tmp);
	model.lrr_auto2sex = NAN;
	model.median_baf_adj = (int)strtol(MEDIAN_BAF_ADJ_DFLT, &tmp, 0);
	model.median_lrr_adj = (int)strtol(MEDIAN_LRR_ADJ_DFLT, &tmp, 0);
	model.order_lrr_gc = (int)strtol(ORDER_LRR_GC_DFLT, &tmp, 0);

	// create synced reader object
	bcf_srs_t *sr = bcf_sr_init();
	bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);

	static struct option loptions[] = {{"rules", required_argument, NULL, 'r'},
					   {"rules-file", required_argument, NULL, 'R'},
					   {"samples", required_argument, NULL, 's'},
					   {"samples-file", required_argument, NULL, 'S'},
					   {"force-samples", no_argument, NULL, 1},
					   {"targets", required_argument, NULL, 't'},
					   {"targets-file", required_argument, NULL, 'T'},
					   {"apply-filters", required_argument, NULL, 'f'},
					   {"variants", required_argument, NULL, 'v'},
					   {"threads", required_argument, NULL, 9},
					   {"output", required_argument, NULL, 'o'},
					   {"output-type", required_argument, NULL, 'O'},
					   {"no-version", no_argument, NULL, 8},
					   {"no-annotations", no_argument, NULL, 'a'},
					   {"mosaic-calls", required_argument, NULL, 'm'},
					   {"genome-stats", required_argument, NULL, 'g'},
					   {"ucsc-bed", required_argument, NULL, 'u'},
					   {"log", required_argument, NULL, 'l'},
					   {"no-log", no_argument, NULL, 10},
					   {"cnp", required_argument, NULL, 'p'},
					   {"bdev-LRR-BAF", required_argument, NULL, 11},
					   {"bdev-BAF-phase", required_argument, NULL, 12},
					   {"min-dist", required_argument, NULL, 'd'},
					   {"no-BAF-flip", no_argument, NULL, 13},
					   {"median-BAF-adjust", required_argument, NULL, 14},
					   {"median-LRR-adjust", required_argument, NULL, 15},
					   {"order-LRR-GC", required_argument, NULL, 16},
					   {"xy-prob", required_argument, NULL, 17},
					   {"err-prob", required_argument, NULL, 18},
					   {"flip-prob", required_argument, NULL, 19},
					   {"telomere-advantage", required_argument, NULL, 20},
					   {"centromere-penalty", required_argument, NULL, 21},
					   {"short_arm_chrs", required_argument, NULL, 22},
					   {"use_short_arms", no_argument, NULL, 23},
					   {"use_centromeres", no_argument, NULL, 24},
					   {"LRR-cutoff", required_argument, NULL, 25},
					   {"LRR-hap2dip", required_argument, NULL, 26},
					   {"LRR-auto2sex", required_argument, NULL, 27},
					   {"LRR-weight", required_argument, NULL, 28},
					   {NULL, 0, NULL, 0}};
	int c;
	while ((c = getopt_long(argc, argv, "h?r:R:s:S:t:T:f:v:o:O:am:g:u:l:p:c:b:d:", loptions,
				NULL))
	       >= 0) {
		switch (c) {
		case 'r':
			rules = optarg;
			break;
		case 'R':
			rules = optarg;
			rules_is_file = 1;
			break;
		case 's':
			sample_names = optarg;
			break;
		case 'S':
			sample_names = optarg;
			sample_is_file = 1;
			break;
		case 1:
			force_samples = 1;
			break;
		case 't':
			targets_list = optarg;
			break;
		case 'T':
			targets_list = optarg;
			targets_is_file = 1;
			break;
		case 'f':
			sr->apply_filters = optarg;
			break;
		case 'v':
			if (optarg[0] == '^') {
				filter_fname = optarg + 1;
				model.flags |= FLT_EXCLUDE;
			} else {
				filter_fname = optarg;
				model.flags |= FLT_INCLUDE;
			}
			break;
		case 9:
			n_threads = (int)strtol(optarg, &tmp, 0);
			if (*tmp)
				error("Could not parse: --threads %s\n", optarg);
			break;
		case 'o':
			output_fname = optarg;
			break;
		case 'O':
			switch (optarg[0]) {
			case 'b':
				output_type = FT_BCF_GZ;
				break;
			case 'u':
				output_type = FT_BCF;
				break;
			case 'z':
				output_type = FT_VCF_GZ;
				break;
			case 'v':
				output_type = FT_VCF;
				break;
			default:
				error("The output type \"%s\" not recognised\n", optarg);
			};
			break;
		case 8:
			record_cmd_line = 0;
			break;
		case 'a':
			model.flags |= NO_ANNOT;
			break;
		case 'l':
			log_file = get_file_handle(optarg);
			break;
		case 10:
			model.flags |= NO_LOG;
			break;
		case 'm':
			out_fm = get_file_handle(optarg);
			break;
		case 'g':
			out_fg = get_file_handle(optarg);
			break;
		case 'u':
			out_fu = get_file_handle(optarg);
			break;
		case 'p':
			cnp_fname = optarg;
			break;
		case 11:
			bdev_lrr_baf = optarg;
			break;
		case 12:
			bdev_baf_phase = optarg;
			break;
		case 'd':
			model.min_dst = (int)strtol(optarg, &tmp, 0);
			if (*tmp)
				error("Could not parse: --min-dist %s\n", optarg);
			break;
		case 13:
			model.flags |= NO_BAF_FLIP;
			break;
		case 14:
			model.median_baf_adj = (int)strtol(optarg, &tmp, 0);
			if (*tmp)
				error("Could not parse: --median-BAF-adjust %s\n", optarg);
			break;
		case 15:
			model.median_lrr_adj = (int)strtol(optarg, &tmp, 0);
			if (*tmp)
				error("Could not parse: --median-LRR-adjust %s\n", optarg);
			break;
		case 16:
			model.order_lrr_gc = (int)strtol(optarg, &tmp, 0);
			if (*tmp)
				error("Could not parse: --order-LRR-GC %s\n", optarg);
			break;
		case 17:
			model.xy_log_prb = logf(strtof(optarg, &tmp));
			if (*tmp)
				error("Could not parse: --xy-prob %s\n", optarg);
			break;
		case 18:
			model.err_log_prb = logf(strtof(optarg, &tmp));
			if (*tmp)
				error("Could not parse: --err-prob %s\n", optarg);
			break;
		case 19:
			model.flip_log_prb = logf(strtof(optarg, &tmp));
			if (*tmp)
				error("Could not parse: --flip-prob %s\n", optarg);
			break;
		case 20:
			model.tel_log_prb = logf(strtof(optarg, &tmp));
			if (*tmp)
				error("Could not parse: --telomere-advantage %s\n", optarg);
			break;
		case 21:
			model.cen_log_prb = logf(strtof(optarg, &tmp));
			if (*tmp)
				error("Could not parse: --centromere-penalty %s\n", optarg);
			break;
		case 22:
			short_arm_chrs = optarg;
			break;
		case 23:
			model.flags |= USE_SHORT_ARMS;
			break;
		case 24:
			model.flags |= USE_CENTROMERES;
			break;
		case 25:
			model.lrr_cutoff = strtof(optarg, &tmp);
			if (*tmp)
				error("Could not parse: --LRR-cutoff %s\n", optarg);
			break;
		case 26:
			model.lrr_hap2dip = strtof(optarg, &tmp);
			if (*tmp)
				error("Could not parse: --LRR-hap2dip %s\n", optarg);
			break;
		case 27:
			model.lrr_auto2sex = strtof(optarg, &tmp);
			if (*tmp)
				error("Could not parse: --LRR-auto2sex %s\n", optarg);
			break;
		case 28:
			model.lrr_bias = strtof(optarg, &tmp);
			if (*tmp)
				error("Could not parse: --LRR-weight %s\n", optarg);
			break;
		case 'h':
		case '?':
			error("%s", usage_text());
			break;
		default:
			error("Unknown argument: %s\n", optarg);
		}
	}

	if (!rules) {
		fprintf(log_file,
			"Genome reference assembly was not specified with --rules or --rules-file\n");
		error("%s", usage_text());
	}
	int len = strlen(rules);
	if (!rules_is_file && (strncmp(rules, "list", 4) == 0 || rules[len - 1] == '?'))
		genome_init_alias(log_file, rules, NULL);

	if (output_fname == NULL && (model.flags & NO_ANNOT)) {
		fprintf(log_file, "Option --no-annotations requires option --output\n");
		error("%s", usage_text());
	}

	if (model.order_lrr_gc > MAX_ORDER) {
		fprintf(log_file,
			"Polynomial order must not be greater than %d: --order-LRR-GC %d\n",
			MAX_ORDER, model.order_lrr_gc);
		error("%s", usage_text());
	}

	// parse parameters defining hidden states
	model.bdev_lrr_baf = read_list_invf(bdev_lrr_baf, &model.bdev_lrr_baf_n, -0.5f, 0.25f);
	model.bdev_baf_phase =
		read_list_invf(bdev_baf_phase, &model.bdev_baf_phase_n, 0.0f, 0.5f);

	// read list of regions to genotype
	if (cnp_fname) {
		model.cnp_idx = regidx_init(cnp_fname, cnp_parse, NULL, sizeof(int), NULL);
		if (!model.cnp_idx)
			error("Error: failed to initialize CNP regions: --cnp %s\n", cnp_fname);
		model.cnp_itr = regitr_init(model.cnp_idx);
	}

	// input VCF
	char *input_fname = NULL;
	if (optind >= argc) {
		if (!isatty(fileno((FILE *)stdin)))
			input_fname = "-";
	} else
		input_fname = argv[optind];
	if (!input_fname)
		error("%s", usage_text());

	// read in the regions from the command line
	if (targets_list) {
		if (bcf_sr_set_targets(sr, targets_list, targets_is_file, 0) < 0)
			error("Failed to read the targets: %s\n", targets_list);
	}

	if (bcf_sr_set_threads(sr, n_threads) < 0)
		error("Failed to create threads\n");
	if (!bcf_sr_add_reader(sr, input_fname))
		error("Failed to open %s: %s\n", input_fname, bcf_sr_strerror(sr->errnum));
	if (filter_fname)
		if (!bcf_sr_add_reader(sr, filter_fname))
			error("Failed to open %s: %s\n", filter_fname,
			      bcf_sr_strerror(sr->errnum));

	// check whether the necessary information has been included in the VCF
	hdr = bcf_sr_get_header(sr, 0);
	if (bcf_hdr_nsamples(hdr) == 0)
		error("Error: input VCF file has no samples\n");
	if (bcf_hdr_id2int(hdr, BCF_DT_ID, "GT") < 0)
		error("Error: input VCF file has no GT format field\n");
	if (!(bcf_hdr_id2int(hdr, BCF_DT_ID, "AD") < 0)) {
		model.flags |= WGS_DATA;
		model.lrr_hap2dip = (float)M_LN2;
	} else if ((bcf_hdr_id2int(hdr, BCF_DT_ID, "LRR") < 0
		    || bcf_hdr_id2int(hdr, BCF_DT_ID, "BAF") < 0))
		error("Error: input VCF file must contain either the AD format field or the LRR and BAF format fields\n");
	if (model.order_lrr_gc > 0 && (bcf_hdr_id2int(hdr, BCF_DT_ID, "GC") < 0))
		error("Error: input VCF has no GC info field: use \"--order-LRR-GC 0/-1\" to disable LRR adjustment through GC correction\n");

	// median adjustment is necessary with sequencing counts data
	if ((model.order_lrr_gc < 0) && (model.flags & WGS_DATA))
		model.order_lrr_gc = 0;

	fprintf(log_file, "Running MoChA version %s\n", MOCHA_VERSION);
	fprintf(log_file, "For updates: https://github.com/freeseek/mocha\n");

	// initialize genome parameters
	if (rules_is_file)
		model.genome_rules = genome_init_file(rules, hdr);
	else
		model.genome_rules = genome_init_alias(log_file, rules, hdr);
	if (!(model.flags & NO_LOG))
		fprintf(log_file, "Using genome assembly from %s\n", rules);
	readlist_short_arms(model.genome_rules, short_arm_chrs, hdr);

	// subset VCF file
	if (sample_names) {
		int ret = bcf_hdr_set_samples(hdr, sample_names, sample_is_file);
		if (ret < 0)
			error("Error parsing the sample list\n");
		else if (ret > 0) {
			if (force_samples)
				fprintf(log_file,
					"Warn: sample #%d not found in the header... skipping\n",
					ret);
			else
				error("Error: sample #%d not found in the header. Use \"--force-samples\" to ignore this error\n",
				      ret);
		}
		if (bcf_hdr_nsamples(hdr) == 0)
			error("Error: subsetting has removed all samples\n");
	}

	// output VCF
	if (output_fname) {
		out_fh = hts_open(output_fname, hts_bcf_wmode(output_type));
		if (out_fh == NULL)
			error("Cannot write to \"%s\": %s\n", output_fname, strerror(errno));
		if (n_threads)
			hts_set_opt(out_fh, HTS_OPT_THREAD_POOL, sr->p);
		out_hdr = print_hdr(out_fh, hdr, argc, argv, record_cmd_line, model.flags);
	}

	int nsmpl = bcf_hdr_nsamples(hdr);
	if (nsmpl == 0)
		error("Subsetting has removed all samples\n");
	if (!(model.flags & NO_LOG))
		fprintf(log_file, "Loading %d sample(s) from the VCF file\n", nsmpl);
	if (nsmpl < model.median_baf_adj && !(model.flags & WGS_DATA))
		error("Error: cannot perform median BAF adjustment with only %d sample(s): use \"--median-BAF-adjust -1\" to disable BAF adjustment\n",
		      nsmpl);
	if (nsmpl < model.median_lrr_adj && !(model.flags & WGS_DATA))
		error("Error: cannot perform median LRR adjustment with only %d sample(s): use \"--median-LRR-adjust -1\" to disable LRR adjustment\n",
		      nsmpl);

	sample = (sample_t *)calloc(nsmpl, sizeof(sample_t));
	for (int i = 0; i < nsmpl; i++) {
		sample[i].idx = i;
		sample[i].x_nonpar_lrr_median = NAN;
		sample[i].y_nonpar_lrr_median = NAN;
		sample[i].mt_lrr_median = NAN;
	}

	for (int rid = 0; rid < hdr->n[BCF_DT_CTG]; rid++) {
		model.rid = rid;
		get_contig(sr, sample, &model);
		if (model.n <= 0)
			continue;
		if (!(model.flags & NO_LOG)) {
			if (model.flags & NO_BAF_FLIP || model.n_flipped == 0)
				fprintf(log_file, "Read %d variants from contig %s\n", model.n,
					bcf_hdr_id2name(hdr, rid));
			else
				fprintf(log_file,
					"Read %d variants (%d BAF flipped) from contig %s\n",
					model.n, model.n_flipped, bcf_hdr_id2name(hdr, rid));
		}
		if (model.genome_rules->length[rid] < model.locus_arr[model.n - 1].pos)
			model.genome_rules->length[rid] = model.locus_arr[model.n - 1].pos;
		for (int j = 0; j < nsmpl; j++)
			sample_stats(sample + j, &model);
	}

	sample_summary(sample, nsmpl, &model);
	int cnt[3] = {0, 0, 0};
	for (int i = 0; i < nsmpl; i++)
		cnt[sample[i].sex]++;
	if (!(model.flags & NO_LOG))
		fprintf(log_file,
			"Estimated %d sample(s) of unknown sex, %d male(s) and %d female(s)\n",
			cnt[SEX_UNK], cnt[SEX_MAL], cnt[SEX_FEM]);
	sample_print(out_fg, sample, nsmpl, hdr, model.flags);

	if (isnan(model.lrr_cutoff))
		error("Error: Unable to estimate LRR-cutoff. Make sure "
		      "the X nonPAR region and both male and female samples are present in the VCF or specify the parameter\n");
	if (isnan(model.lrr_auto2sex))
		error("Error: Unable to estimate LRR-auto2sex. Make sure "
		      "both autosomes and the X nonPAR region are present in the VCF or specify the parameter\n");

	if (!(model.flags & NO_LOG))
		fprintf(log_file,
			"Model LRR parameters: LRR-cutoff=%.4f LRR-hap2dip=%.4f LRR-auto2sex=%.4f\n",
			model.lrr_cutoff, model.lrr_hap2dip, model.lrr_auto2sex);

	for (int rid = 0; rid < hdr->n[BCF_DT_CTG]; rid++) {
		model.rid = rid;
		get_contig(sr, sample, &model);
		if (model.n <= 0)
			continue;
		if (!(model.flags & NO_LOG)) {
			if (model.flags & NO_BAF_FLIP || model.n_flipped == 0)
				fprintf(log_file, "Read %d variants from contig %s\n", model.n,
					bcf_hdr_id2name(hdr, rid));
			else
				fprintf(log_file,
					"Read %d variants (%d BAF flipped) from contig %s\n",
					model.n, model.n_flipped, bcf_hdr_id2name(hdr, rid));
		}
		for (int j = 0; j < nsmpl; j++) {
			if (model.cnp_idx)
				regidx_overlap(model.cnp_idx, bcf_hdr_id2name(hdr, rid), 0,
					       model.genome_rules->length[rid], model.cnp_itr);
			sample_run(sample + j, &mocha_table, &model);
		}

		if (output_fname) {
			int nret = put_contig(sr, sample, &model, out_fh, out_hdr);
			if (!(model.flags & NO_LOG))
				fprintf(log_file, "Written %d variants for contig %s\n", nret,
					bcf_hdr_id2name(hdr, rid));
		}
	}

	// estimate LRR at common autosomal deletions and duplications
	if (!(model.flags & NO_LOG) && model.cnp_itr) {
		int n_cnp_del = 0, n_cnp_dup = 0, n_rare_dup = 0, n_rare_del = 0;
		float *cnp_ldev = (float *)malloc(mocha_table.n * sizeof(float));
		float *rare_ldev = (float *)malloc(mocha_table.n * sizeof(float));
		for (int i = 0; i < mocha_table.n; i++) {
			if (mocha_table.a[i].rid == model.genome_rules->x_rid)
				continue;
			if (mocha_table.a[i].type == MOCHA_CNP_DEL
			    && mocha_table.a[i].nhets == 0)
				cnp_ldev[n_cnp_del++] = mocha_table.a[i].ldev;
			else if (mocha_table.a[i].type == MOCHA_CNP_DUP
				 && mocha_table.a[i].nhets >= 5
				 && mocha_table.a[i].bdev >= 0.1f)
				cnp_ldev[mocha_table.n - (++n_cnp_dup)] = mocha_table.a[i].ldev;
			else if (mocha_table.a[i].type == MOCHA_DEL
				 && mocha_table.a[i].nhets == 0)
				rare_ldev[n_rare_del++] = mocha_table.a[i].ldev;
			else if (mocha_table.a[i].type == MOCHA_DUP
				 && mocha_table.a[i].nhets >= 5
				 && mocha_table.a[i].bdev >= 0.1f)
				rare_ldev[mocha_table.n - (++n_rare_dup)] =
					mocha_table.a[i].ldev;
		}
		float cnp_del_ldev = get_median(cnp_ldev, n_cnp_del, NULL);
		float cnp_dup_ldev =
			get_median(&cnp_ldev[mocha_table.n - n_cnp_dup], n_cnp_dup, NULL);
		float rare_del_ldev = get_median(rare_ldev, n_rare_del, NULL);
		float rare_dup_ldev =
			get_median(&rare_ldev[mocha_table.n - n_rare_dup], n_rare_dup, NULL);
		fprintf(log_file, "Adjusted LRR at CNP deletions=%.4f\n", cnp_del_ldev);
		fprintf(log_file, "Adjusted LRR at CNP duplications=%.4f\n", cnp_dup_ldev);
		fprintf(log_file, "Adjusted LRR at rare deletions=%.4f\n", rare_del_ldev);
		fprintf(log_file, "Adjusted LRR at rare duplications=%.4f\n", rare_dup_ldev);
		free(cnp_ldev);
		free(rare_ldev);
	}

	// clear sample data
	for (int j = 0; j < nsmpl; j++) {
		free(sample[j].vcf_imap_arr);
		free(sample[j].data_arr[BDEV]);
		free(sample[j].data_arr[LDEV]);
		free(sample[j].phase_arr);
	}

	// clear model data
	free(model.locus_arr);
	free(model.gc_arr);
	free(model.bdev_lrr_baf);
	free(model.bdev_baf_phase);
	genome_destroy(model.genome_rules);

	// free precomputed tables
	ad_to_lrr_baf(NULL, NULL, NULL, NULL, 0);

	// write table with mosaic chromosomal alterations (and UCSC bed track)
	mocha_print(out_fm, mocha_table.a, mocha_table.n, hdr, model.flags, rules,
		    model.lrr_hap2dip);
	mocha_print_ucsc(out_fu, mocha_table.a, mocha_table.n, hdr);
	free(mocha_table.a);

	// close output VCF
	if (output_fname) {
		bcf_hdr_destroy(out_hdr);
		hts_close(out_fh);
	}

	if (log_file != stdout && log_file != stderr)
		fclose(log_file);

	// clean up
	if (model.cnp_idx)
		regidx_destroy(model.cnp_idx);
	if (model.cnp_itr)
		regitr_destroy(model.cnp_itr);
	bcf_sr_destroy(sr);
	free(sample);

	beta_binom_destroy(beta_binom_null);
	beta_binom_destroy(beta_binom_alt);
	return 0;
}
