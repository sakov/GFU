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
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include "version.h"
#include "ncw.h"
#include "utils.h"

#define PROGRAM_NAME "ncmask"
#define PROGRAM_VERSION "0.05"
#define VERBOSE_DEF 1

#define MASKTYPE_NONE 0
#define MASKTYPE_BINARY 1
#define MASKTYPE_NLAYERS 2

#define FILL_ZERO 0
#define FILL_NAN 1
#define FILL_FILLVALUE 2

#define MAXNDIMS 4

int verbose = VERBOSE_DEF;

/**
 */
static void usage(int status)
{
    printf("  Usage: %s <file> <var> [0*|nan|fillvalue] -m <file> <var> [-v {0|1*|2}]\n", PROGRAM_NAME);
    printf("         %s -v\n", PROGRAM_NAME);
    printf("  Options:\n");
    printf("    <file> <var> [0*|nan|fillvalue] - data file, variable and the fill value\n");
    printf("    -m <file> <var> -- land mask (for 2D or 3D variables: either 2D with 0s and 1s;\n");
    printf("       or 2D with number of valid layers in a column\n");
    printf("       of the variable with 0s and 1s\n");
    printf("    -v {0|1*|2} -- verbosity level | print version and exit\n");
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

    int ncid = -1, varid = -1, ndims = -1;
    nc_type vtype;
    size_t dimlens[NC_MAX_DIMS];
    int ndims_essential = -1;
    int dimids_essential[MAXNDIMS];

    void* fillvalue = NULL;
    int64_t zerovalue = 0;

    int mncid = -1, mvarid = -1;
    size_t msize;

    void* v = NULL;
    int* mask = NULL;
    int masktype = MASKTYPE_NONE;

    int ni = 1, nj = 1, nk = 1, nr = 1;
    int i, k, r;

    parse_commandline(argc, argv, &fname, &varname, &mfname, &mvarname, &fill);
    if (fname == NULL)
        quit("no data file specified");
    if (mfname == NULL)
        quit("no mask specified");

    /*
     * data
     */
    if (verbose > 1) {
        printf("  data = %s\n", fname);
        printf("    variable = %s\n", varname);
        if (nk > 1)
            printf("    %d layers\n", nk);
    }
    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vartype(ncid, varid, &vtype);
    if (fill == FILL_NAN && vtype != NC_FLOAT && vtype != NC_DOUBLE)
        quit("%s: fill value can not be set to \"nan\" for data of type \"%s\"", varname, ncw_nctype2str(vtype));
    if (ndims > MAXNDIMS)
        quit("%s: ndims = %d; should not exceed %d", varname, ndims, MAXNDIMS);
    ncw_inq_vardims(ncid, varid, MAXNDIMS, &ndims, dimlens);
    if (verbose > 1) {
        printf("    size = ");
        for (i = 0; i < ndims; ++i)
            printf("%zu%s", dimlens[i], i < ndims - 1 ? " x " : "\n");
    }
    for (i = 0, ndims_essential = 0; i < ndims; ++i)
        if (dimlens[i] > 1)
            dimids_essential[ndims_essential++] = i;
    if (ndims_essential == 4) {
        nr = dimlens[0];
        nk = dimlens[1];
        nj = dimlens[2];
        ni = dimlens[3];
        if (verbose > 1) {
            printf("    %d records\n", nr);
            printf("    %d layers\n", nk);
        }
    } else if (ndims_essential == 3) {
        nk = dimlens[dimids_essential[0]];
        nj = dimlens[dimids_essential[1]];
        ni = dimlens[dimids_essential[2]];
        if (verbose > 1)
            printf("    %d layers\n", nk);
    } else if (ndims_essential == 2) {
        nj = dimlens[dimids_essential[0]];
        ni = dimlens[dimids_essential[1]];
    } else if (ndims_essential == 1) {
        ni = dimlens[dimids_essential[0]];
    }

    /*
     * mask
     */
    if (verbose > 1) {
        printf("  mask = %s\n", mfname);
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
    if (msize != ni * nj)
        quit("mask size %zu is not equal layer size %zu", msize, ni * nj);
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
    if (verbose > 1)
        printf("    type = %s\n", (masktype == MASKTYPE_NLAYERS) ? "no. of valid layers" : "binary");

    if (verbose > 1) {
        printf("  applying:");
        fflush(stdout);
    }

    v = malloc(ni * nj * ncw_sizeof(vtype));

    if (fill == FILL_FILLVALUE) {
        fillvalue = malloc(ncw_sizeof(vtype));
        ncw_inq_var_fill(ncid, varid, NULL, fillvalue);
    } else if (fill == FILL_ZERO && vtype != NC_DOUBLE && vtype != NC_FLOAT) {
        if (ncw_att_exists(ncid, varid, "scale_factor") || ncw_att_exists(ncid, varid, "add_offset")) {
            double scale_factor = 1.0;
            double add_offset = 0.0;

            if (ncw_att_exists(ncid, varid, "scale_factor")) {
                ncw_check_attlen(ncid, varid, "scale_factor", 1);
                ncw_get_att_double(ncid, varid, "scale_factor", &scale_factor);
            }
            if (ncw_att_exists(ncid, varid, "add_offset")) {
                ncw_check_attlen(ncid, varid, "add_offset", 1);
                ncw_get_att_double(ncid, varid, "add_offset", &add_offset);
            }

            zerovalue = -(int64_t) (add_offset / scale_factor);
        }
    }

    for (r = 0; r < nr; ++r) {
        for (k = 0; k < nk; ++k) {
            size_t start[MAXNDIMS] = { 0, 0, 0, 0 };
            size_t count[MAXNDIMS] = { 1, 1, 1, 1 };

            if (ndims_essential == 4) {
                start[0] = k;
                start[1] = r;
                count[2] = nj;
                count[3] = ni;
            } else if (ndims_essential == 3) {
                start[dimids_essential[0]] = k;
                count[dimids_essential[1]] = nj;
                count[dimids_essential[2]] = ni;
            } else if (ndims_essential == 2) {
                count[dimids_essential[0]] = nj;
                count[dimids_essential[1]] = ni;
            } else if (ndims_essential == 1) {
                count[dimids_essential[0]] = ni;
            }

            ncw_get_vara(ncid, varid, start, count, v);
            if (vtype == NC_DOUBLE) {
                double* vv = (double*) v;

                if (fill == FILL_ZERO) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = 0.0;
                } else if (fill == FILL_NAN) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = NAN;
                } else if (fill == FILL_FILLVALUE) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = ((double*) fillvalue)[0];
                }
            } else if (vtype == NC_FLOAT) {
                float* vv = (float*) v;

                if (fill == FILL_ZERO) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = 0.0;
                } else if (fill == FILL_NAN) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = NAN;
                } else if (fill == FILL_FILLVALUE) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = ((float*) fillvalue)[0];
                }
            } else if (ncw_sizeof(vtype) == 1) {
                int8_t* vv = (int8_t*) v;

                if (fill == FILL_ZERO) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = zerovalue;
                } else if (fill == FILL_FILLVALUE) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = ((int8_t*) fillvalue)[0];
                }
            } else if (ncw_sizeof(vtype) == 2) {
                int16_t* vv = (int16_t*) v;

                if (fill == FILL_ZERO) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = zerovalue;
                } else if (fill == FILL_FILLVALUE) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = ((int16_t*) fillvalue)[0];
                }
            } else if (ncw_sizeof(vtype) == 4) {
                int32_t* vv = (int32_t*) v;

                if (fill == FILL_ZERO) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = zerovalue;
                } else if (fill == FILL_FILLVALUE) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = ((int32_t*) fillvalue)[0];
                }
            } else if (ncw_sizeof(vtype) == 8) {
                int64_t* vv = (int64_t*) v;

                if (fill == FILL_ZERO) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = zerovalue;
                } else if (fill == FILL_FILLVALUE) {
                    for (i = 0; i < ni * nj; ++i)
                        if (mask[i] <= k)
                            vv[i] = ((int64_t*) fillvalue)[0];
                }
            }
            ncw_put_vara(ncid, varid, start, count, v);
            if (verbose && nk > 1) {
                printf(".");
                fflush(stdout);
            }
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

    if (fillvalue != NULL)
        free(fillvalue);
    free(v);
    free(mask);
}
