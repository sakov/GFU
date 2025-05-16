/******************************************************************************
 *
 * File:        ncminmax.c
 *
 * Created:     05/12/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: NCMINMAX finds minimum and maximum values of a variable.
 *             
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include "version.h"
#include "ncw.h"
#include "ncutils.h"
#include "utils.h"

#define PROGRAM_NAME "ncminmax"
#define PROGRAM_VERSION "0.12"
#define VERBOSE_DEF 0

#define MASKTYPE_NONE 0
#define MASKTYPE_BINARY 1
#define MASKTYPE_NLAYERS 2

int verbose = VERBOSE_DEF;
int strict = 0;
int doave = 0;

/**
 */
static void usage(int status)
{
    printf("  Usage: %s <file> <var> [-m <file> <var>] [-a] [-s] [-v {0*|1|2}]\n", PROGRAM_NAME);
    printf("         %s -v\n", PROGRAM_NAME);
    printf("  Options:\n");
    printf("    -m <file> <var> -- set mask (for 2D or 3D variables: either 2D with 0s and 1s;\n");
    printf("       or 2D with number of valid layers in a column\n");
    printf("    -a -- also report average\n");
    printf("    -s -- strict (no missing values allowed)\n");
    printf("    -v {0*|1|2} -- verbosity level | print version and exit\n");
    exit(status);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname, char** varname, char** mfname, char** mvarname)
{
    int i;

    if (argc == 1)
        usage(0);

    if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 'v') {
        printf("  %s v%s\n", PROGRAM_NAME, PROGRAM_VERSION);
        printf("  GFU v%s\n", VERSION);
        exit(0);
    }

    i = 1;
    while (i < argc) {
        if (argv[i][0] == '-') {
            if (argv[i][1] == 'v') {
                i++;
                if (i == argc || argv[i][0] == '-')
                    quit("no verbosity level specified after \"-v\"");
                str2int(argv[i], &verbose);
                i++;
            } else if (argv[i][1] == 's') {
                strict = 1;
                i++;
            } else if (argv[i][1] == 'a') {
                doave = 1;
                i++;
            } else if (argv[i][1] == 'm') {
                i++;
                if (i >= argc)
                    quit("no mask file name specified after \"-m\"\n");
                if (*mfname != NULL)
                    quit("-m: mask file name already specified\n");
                *mfname = argv[i];
                i++;
                if (i >= argc)
                    quit("no mask variable name specified after \"-m\"\n");
                *mvarname = argv[i];
                i++;
            } else
                quit("unknown option \"%s\"", argv[i]);
        } else if (*fname == NULL) {
            *fname = argv[i];
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no variable name specified");
            *varname = argv[i];
            i++;
        } else
            usage(1);
    }
}

/**
 */
int main(int argc, char* argv[])
{
    char* fname = NULL;
    char* varname = NULL;
    char* mfname = NULL;
    char* mvarname = NULL;

    int ncid = -1, varid = -1, ndims = -1;
    int dimids[NC_MAX_DIMS];
    size_t dimlens[NC_MAX_DIMS];
    size_t i, ij, size, slabsize;

    int mncid = -1, mvarid = -1;
    size_t msize;

    double* v = NULL;
    int* mask = NULL;
    int masktype = MASKTYPE_NONE;

    double min = DBL_MAX;
    double max = -DBL_MAX;
    double ave = 0.0;
    size_t n = 0;
    size_t imin = 0;
    size_t imax = 0;
    int nk, k;

    parse_commandline(argc, argv, &fname, &varname, &mfname, &mvarname);

    nk = ncu_getnfields(fname, varname);
    if (nk == 0)
        nk = 1;
    if (nk == 1 && verbose == 2)
        verbose = 1;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_var(ncid, varid, NULL, NULL, &ndims, dimids, NULL);
    for (i = 0; i < ndims; ++i)
        ncw_inq_dimlen(ncid, dimids[i], &dimlens[i]);
    size = ncw_get_varsize(ncid, varid);
    if (nk > 1)
        ncw_close(ncid);

    v = calloc(size / nk, sizeof(double));

    if (mfname != NULL) {
        ncw_open(mfname, NC_NOWRITE, &mncid);
        ncw_inq_varid(mncid, mvarname, &mvarid);
        msize = ncw_get_varsize(mncid, mvarid);
        if (msize != size / nk)
            quit("mask size %zu is not equal to layer size %zu", msize, size / nk);
        mask = calloc(msize, sizeof(int));
        ncw_get_var_int(mncid, mvarid, mask);
        ncw_close(mncid);

        masktype = MASKTYPE_BINARY;
        if (nk > 1) {
            for (i = 0; i < msize; ++i) {
                if (mask[i] > 1) {
                    masktype = MASKTYPE_NLAYERS;
                    break;
                }
            }
            if (masktype == MASKTYPE_BINARY)
                for (i = 0; i < msize; ++i)
                    if (mask[i] != 0)
                        mask[i] = nk;
        }
    }

    if (verbose > 1)
        printf("  %s:\n", fname);
    for (k = 0, ij = 0; k < nk; ++k) {
        double min_k = DBL_MAX;
        double max_k = -DBL_MAX;
        double ave_k = 0.0;
        size_t n_k = 0, imin_k = 0, imax_k = 0;

        if (nk > 1)
            ncu_readfield_double(fname, varname, k, -1, 1, nk, v);
        else {
            ncu_readvardouble(ncid, varid, size, v);
            ncw_close(ncid);
        }

        for (i = 0; i < size / nk; ++i, ++ij) {
            if (mask != NULL && mask[i] <= k)
                continue;
            if (isnan(v[i])) {
                if (!strict)
                    continue;
                else
                    quit("%s(%d) = missing", varname, i);
            }
            if (v[i] > max) {
                max = v[i];
                imax = ij;
            }
            if (v[i] < min) {
                min = v[i];
                imin = ij;
            }
            ave += v[i];
            n++;
            if (verbose > 1) {
                if (v[i] > max_k) {
                    max_k = v[i];
                    imax_k = i;
                }
                if (v[i] < min_k) {
                    min_k = v[i];
                    imin_k = i;
                }
                ave_k += v[i];
                n_k++;
            }
        }
        if (verbose == 1 && nk > 1) {
            printf(".");
            fflush(stdout);
        }
        if (verbose > 1 && n_k > 0) {
            if (doave) {
                ave_k /= (double) n_k;
                printf("    %s: %d: %.4g %.4g %.4g (", varname, k, min_k, ave_k, max_k);
            } else
                printf("    %s: %d: %.4g %.4g (", varname, k, min_k, max_k);
            for (i = ndims - 2, slabsize = size / nk; i < ndims; ++i) {
                size_t ii;

                slabsize /= dimlens[i];
                ii = imin_k / slabsize;
                printf((i == ndims - 2) ? "%lu" : ", %lu", ii);
                imin_k = imin_k - slabsize * ii;
            }
            printf(") (");
            for (i = ndims - 2, slabsize = size / nk; i < ndims; ++i) {
                size_t ii;

                slabsize /= dimlens[i];
                ii = imax_k / slabsize;
                printf((i == ndims - 2) ? "%lu" : ", %lu", ii);
                imax_k = imax_k - slabsize * ii;
            }
            printf(")\n");
        }
    }
    if (verbose == 1 && nk > 1)
        printf("\n");

    if (mask != NULL)
        free(mask);
    free(v);

    ave /= (double) n;
    if (verbose) {
        if (verbose == 1)
            printf("  %s:\n", fname);
        printf("    %s: min = %.4g\n", varname, min);
        printf("    %s: max = %.4g\n", varname, max);
        if (doave)
            printf("    %s: ave = %.4g\n", varname, ave);
        printf("    %s: size = ", varname);
        for (i = 0; i < ndims; ++i)
            printf((i > 0) ? " x %lu" : "%lu", dimlens[i]);
        printf("\n");
        printf("    %s: imin = %lu (", varname, imin);
        for (i = 0, slabsize = size; i < ndims; ++i) {
            size_t ii;

            slabsize /= dimlens[i];
            ii = imin / slabsize;
            printf((i == 0) ? "%lu" : ", %lu", ii);
            imin = imin - slabsize * ii;
        }
        printf(")\n");
        printf("    %s: imax = %lu (", varname, imax);
        slabsize = size;
        for (i = 0, slabsize = size; i < ndims; ++i) {
            size_t ii;

            slabsize /= dimlens[i];
            ii = imax / slabsize;
            printf((i == 0) ? "%lu" : ", %lu", ii);
            imax = imax - slabsize * ii;
        }
        printf(")\n");
        printf("    %s: %zu valid values (%.2f%%)\n", varname, n, (double) n / (double) size * 100.0);
    } else {
        if (doave)
            printf("  %.4g %.4g %.4g\n", min, ave, max);
        else
            printf("  %.4g %.4g\n", min, max);
    }
    return 0;
}
