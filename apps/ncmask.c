/******************************************************************************
 *
 * File:        ncmask.c
 *
 * Created:     02/05/2025
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: NCMASK writes 0, NaN, or fill value to areas defined by a layer
 *              mask from an external file.
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
#include <unistd.h>
#include "version.h"
#include "ncw.h"
#include "ncutils.h"
#include "utils.h"

#define PROGRAM_NAME "ncmask"
#define PROGRAM_VERSION "0.02"
#define VERBOSE_DEF 1

#define MASKTYPE_NONE 0
#define MASKTYPE_BINARY 1
#define MASKTYPE_NLAYERS 2

#define FILL_ZERO 0
#define FILL_NAN 1
#define FILL_FILLVALUE 2

int verbose = VERBOSE_DEF;

/**
 */
static void usage(int status)
{
    printf("  Usage: %s <file> <var> [0|{nan|fillvalue}] -m <file> <var> [-v]\n", PROGRAM_NAME);
    printf("         %s -v [0|1|2]\n", PROGRAM_NAME);
    printf("  Options:\n");
    printf("    <file> <var> [0|{nan|fillvalue}] - data file, variable and the fill value\n");
    printf("       (default = 0)\n");
    printf("    -m <file> <var> -- land mask (either 2D with 0s and 1s; or 2D with number\n");
    printf("       of valid layers in a column; or for 1D or 2D variables of the size of\n");
    printf("       the variable with 0s and 1s\n");
    printf("    -v [0|1|2] -- verbosity level (default = 1) | print version and exit\n");
    exit(status);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname, char** varname, char** mfname, char** mvarname, int* fill)
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
            if (i == argc && argv[i][0] == '-')
                quit("no variable name specified");
            *varname = argv[i];
            i++;
            if (i < argc && argv[i][0] != '-') {
                if (strcmp(argv[i], "0") == 0)
                    *fill = FILL_ZERO;
                else if (strcasecmp(argv[i], "nan") == 0)
                    *fill = FILL_NAN;
                else if (strcmp(argv[i], "fillvalue") == 0)
                    *fill = FILL_FILLVALUE;
                else
                    quit("could not understand fill value specification \"%s\"", argv[i]);
                i++;
            }
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
    int fill = FILL_ZERO;

    int ncid = -1, varid = -1;
    nc_type vtype;
    size_t i, size, slabsize;
    float fillvalue_f = NAN;
    double fillvalue_d = NAN;

    int mncid = -1, mvarid = -1;
    size_t msize;

    void* v = NULL;
    int* mask = NULL;
    int masktype = MASKTYPE_NONE;

    int layered;
    int nk, k;

    parse_commandline(argc, argv, &fname, &varname, &mfname, &mvarname, &fill);
    if (fname == NULL)
        quit("no data file specified");
    if (mfname == NULL)
        quit("no mask specified");

    layered = 1;
    nk = ncu_getnfields(fname, varname);
    if (nk == 0) {
        nk = 1;
        layered = 0;
    }

    if (verbose > 1) {
        printf("  data = %s\n", fname);
        printf("    variable = %s\n", varname);
        if (nk > 1)
            printf("    %d layers\n", nk);
    }

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vartype(ncid, varid, &vtype);
    if (verbose > 1) {
        int ndims;
        size_t dimlens[NC_MAX_DIMS];

        ncw_inq_vardims(ncid, varid, 4, &ndims, dimlens);
        printf("    size = ");
        for (i = 0; i < ndims; ++i)
            printf("%zu%s", dimlens[i], i < ndims - 1 ? " x " : "\n");
    }
    size = ncw_get_varsize(ncid, varid);
    slabsize = size / nk;
    if (fill == FILL_FILLVALUE) {
        if (vtype == NC_FLOAT)
            ncw_inq_var_fill(ncid, varid, NULL, &fillvalue_f);
        else
            ncw_inq_var_fill(ncid, varid, NULL, &fillvalue_d);
    }
    ncw_close(ncid);

    if (vtype == NC_FLOAT)
        v = calloc(slabsize, sizeof(float));
    else
        v = calloc(slabsize, sizeof(double));

    if (verbose > 1) {
        printf("  mask = %s\n", fname);
        printf("    variable = %s\n", mvarname);
    }
    ncw_open(mfname, NC_NOWRITE, &mncid);
    ncw_inq_varid(mncid, mvarname, &mvarid);
    if (verbose > 1) {
        int ndims;
        size_t dimlens[NC_MAX_DIMS];

        ncw_inq_vardims(mncid, mvarid, 4, &ndims, dimlens);
        printf("    size = ");
        for (i = 0; i < ndims; ++i)
            printf("%zu%s", dimlens[i], i < ndims - 1 ? " x " : "\n");
    }
    msize = ncw_get_varsize(mncid, mvarid);
    if (msize != slabsize)
        quit("mask size %zu is not equal layer size %zu", msize, slabsize);
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
    }
    if (verbose > 1) {
        printf("    type = %s\n", (masktype == MASKTYPE_NLAYERS) ? "no. of valid layers" : "binary");
    }

    if (verbose > 1) {
        printf("  applying the mask:");
        fflush(stdout);
    }

    for (k = 0; k < nk; ++k) {
        if (layered) {
            if (vtype == NC_FLOAT)
                ncu_readfield(fname, varname, k, -1, 1, nk, v);
            else
                ncu_readfield_double(fname, varname, k, -1, 1, nk, v);
        } else {
            if (vtype == NC_FLOAT)
                ncu_readvarfloat(ncid, varid, size, v);
            else if (vtype == NC_DOUBLE)
                ncu_readvardouble(ncid, varid, size, v);
        }
        if (masktype == MASKTYPE_BINARY) {
            if (vtype == NC_FLOAT) {
                float* vf = (float*) v;

                if (fill == FILL_ZERO) {
                    for (i = 0; i < slabsize; ++i)
                        if (mask[i] == 0)
                            vf[i] = 0.0;
                } else if (fill == FILL_NAN) {
                    for (i = 0; i < slabsize; ++i)
                        if (mask[i] == 0)
                            vf[i] = NAN;
                } else if (fill == FILL_FILLVALUE) {
                    for (i = 0; i < slabsize; ++i)
                        if (mask[i] == 0)
                            vf[i] = fillvalue_f;
                }
            } else {
                double* vd = (double*) v;

                if (fill == FILL_ZERO) {
                    for (i = 0; i < slabsize; ++i)
                        if (mask[i] == 0)
                            vd[i] = 0.0;
                } else if (fill == FILL_NAN) {
                    for (i = 0; i < slabsize; ++i)
                        if (mask[i] == 0)
                            vd[i] = NAN;
                } else if (fill == FILL_FILLVALUE) {
                    for (i = 0; i < slabsize; ++i)
                        if (mask[i] == 0)
                            vd[i] = fillvalue_d;
                }
            }
        } else if (masktype == MASKTYPE_NLAYERS) {
            if (vtype == NC_FLOAT) {
                float* vf = (float*) v;

                if (fill == FILL_ZERO) {
                    for (i = 0; i < slabsize; ++i)
                        if (mask[i] <= k)
                            vf[i] = 0.0;
                } else if (fill == FILL_NAN) {
                    for (i = 0; i < slabsize; ++i)
                        if (mask[i] <= k)
                            vf[i] = NAN;
                } else if (fill == FILL_FILLVALUE) {
                    for (i = 0; i < slabsize; ++i)
                        if (mask[i] <= k)
                            vf[i] = fillvalue_f;
                }
            } else {
                double* vd = (double*) v;

                if (fill == FILL_ZERO) {
                    for (i = 0; i < slabsize; ++i)
                        if (mask[i] <= k)
                            vd[i] = 0.0;
                } else if (fill == FILL_NAN) {
                    for (i = 0; i < slabsize; ++i)
                        if (mask[i] <= k)
                            vd[i] = NAN;
                } else if (fill == FILL_FILLVALUE) {
                    for (i = 0; i < slabsize; ++i)
                        if (mask[i] <= k)
                            vd[i] = fillvalue_d;
                }
            }
        }
        if (layered) {
            if (vtype == NC_FLOAT)
                ncu_writefield(fname, varname, k, -1, 1, nk, v);
            else
                ncu_writefield_double(fname, varname, k, -1, 1, nk, v);
        } else {
            ncw_open(fname, NC_WRITE, &ncid);
            if (vtype == NC_FLOAT)
                ncw_put_var_float(ncid, varid, v);
            else
                ncw_put_var_double(ncid, varid, v);
            ncw_close(ncid);
        }
        if (verbose) {
            printf(".");
            fflush(stdout);
        }
    }
    if (verbose)
        printf("\n");

    /*
     * copy command line and wdir to dst
     */
    {
        char* cmd = get_command(argc, argv);
        char attname[NC_MAX_NAME];
        char cwd[MAXSTRLEN];

        snprintf(attname, NC_MAX_NAME, "%s: command", PROGRAM_NAME);
        ncw_open(fname, NC_WRITE, &ncid);
        ncw_put_att_text(ncid, NC_GLOBAL, attname, cmd);
        free(cmd);
        if (getcwd(cwd, MAXSTRLEN) != NULL) {
            snprintf(attname, NC_MAX_NAME, "%s: wdir", PROGRAM_NAME);
            ncw_put_att_text(ncid, NC_GLOBAL, attname, cwd);
        }
        ncw_close(ncid);
    }

    free(v);
    if (mask != NULL)
        free(mask);
}
