/******************************************************************************
 *
 * File:        regrid_ll.c
 *
 * Created:     19/07/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: REGRID_LL is a utility for interpolating geophysical layered
 *              fields onto a different horizontal grid. It targets geographic
 *              grids with nodes defined in lat/lon and provides robust
 *              interpolation on sphere, including topologically complicated
 *              cases (e.g. tripolar grids, grids with discontinuities in
 *              node coordinates etc.). The source and destination horizontal
 *              grids can be either structured ([j][i]) or unstructured ([i]).
 *              The vertical layers are interpolated sequentially.
 *
 *              The data is assumed to be in NetCDF format.
 *             
 * Dependence:  For triangulation related matters the code uses nn library
 *              available from https://github.com/sakov/nn-c.       
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <nn.h>
#include <delaunay.h>
#include "version.h"
#include "ncw.h"
#include "ncutils.h"
#include "utils.h"

#define PROGRAM_NAME "regrid_ll"
#define PROGRAM_VERSION "0.00"

#define VERBOSE_DEF 1
#define DEG2RAD (M_PI / 180.0)

#define GRIDTYPE_UNDEF 0
#define GRIDTYPE_CURV 1
#define GRIDTYPE_RECT 2
#define GRIDTYPE_VECT 3

/**
 */
static void usage(int status)
{
    printf("  Usage: %s -i <src> -o <dst> -v <varname> -gi <src grid> <lon> <lat> [<numlayers>] -go <dst grid> <lon> <lat> [<numlayers>] [-n] [-s] [-V <verblevel>]\n", PROGRAM_NAME);
    printf("         %s -v\n", PROGRAM_NAME);
    printf("  Options:\n");
    printf("    -i <src> -- source file\n");
    printf("    -o <dst> -- destination file (clobbed)\n");
    printf("    -v <varname> -- variable to interpolate\n");
    printf("    -gi <src grid> <lon> <lat> [<numlayers>] -- source grid\n");
    printf("    -go <dst grid> <lon> <lat> [<numlayers>] -- destination grid\n");
    printf("    -n -- flag: do not use deepest valid values for filling the rest of\n");
    printf("          the column\n");
    printf("    -s -- flag: do not use the first and last columns of the source field\n");
    printf("          (e.g. with NEMO on ORCA grids)\n");
    printf("    -V <level> -- set verbosity to 0, 1, or 2 (default = 1)\n");
    printf("    -v -- print version and exit\n");
    exit(status);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname_src, char** fname_dst, char** varname, char** grdname_src, char** xname_src, char** yname_src, char** nkname_src, char** grdname_dst, char** xname_dst, char** yname_dst, char** nkname_dst, int* verbose, int* zerofill, int* skipfirstlast)
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
        if (argv[i][0] != '-') {
            printf("  error: argument \"%s\" does not follow usage\n", argv[i]);
            usage(1);
        }
        if (strcmp(&argv[i][1], "i") == 0) {
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no file name found after \"-i\"");
            *fname_src = argv[i];
            i++;
        } else if (strcmp(&argv[i][1], "o") == 0) {
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no file name found after \"-o\"");
            *fname_dst = argv[i];
            i++;
        } else if (strcmp(&argv[i][1], "v") == 0) {
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no variable name found after \"-v\"");
            *varname = argv[i];
            i++;
        } else if (strcmp(&argv[i][1], "gi") == 0) {
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no file name found after \"-gi\"");
            *grdname_src = argv[i];
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no X coordinate name found after \"-gi\"");
            *xname_src = argv[i];
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no Y coordinate name found after \"-gi\"");
            *yname_src = argv[i];
            i++;
            if (i < argc && argv[i][0] != '-') {
                *nkname_src = argv[i];
                i++;
            }
        } else if (strcmp(&argv[i][1], "go") == 0) {
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no file name found after \"-go\"");
            *grdname_dst = argv[i];
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no X coordinate name found after \"-go\"");
            *xname_dst = argv[i];
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no Y coordinate name found after \"-go\"");
            *yname_dst = argv[i];
            i++;
            if (i < argc && argv[i][0] != '-') {
                *nkname_dst = argv[i];
                i++;
            }
        } else if (strcmp(&argv[i][1], "V") == 0) {
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no verbosity level found after \"-V\"");
            *verbose = atoi(argv[i]);
            i++;
        } else if (strcmp(&argv[i][1], "n") == 0) {
            *zerofill = 1;
            i++;
        } else if (strcmp(&argv[i][1], "s") == 0) {
            *skipfirstlast = 1;
            i++;
        } else
            quit("unknown option \"%s\"", argv[i]);
    }
}

/** Converts from geographic to cartesian coordinates.
 * @param in Input: {lon, lat}
 * @param out Output: {x, y, z}
 */
static void ll2xyz(double in[2], double out[3])
{
    double lon = in[0] * DEG2RAD;
    double lat = in[1] * DEG2RAD;
    double coslat = cos(lat);

    out[0] = sin(lon) * coslat;
    out[1] = cos(lon) * coslat;
    out[2] = sin(lat);
}

/**
 */
int main(int argc, char* argv[])
{
    char* fname_src = NULL;
    char* fname_dst = NULL;
    char* varname = NULL;
    char* grdname_src = NULL;
    char* xname_src = NULL;
    char* yname_src = NULL;
    char* nkname_src = NULL;
    char* grdname_dst = NULL;
    char* xname_dst = NULL;
    char* yname_dst = NULL;
    char* nkname_dst = NULL;
    int verbose = VERBOSE_DEF;
    int gridtype = GRIDTYPE_UNDEF;
    int zerofill = 0;
    int skipfirstlast = 0;

    int npoint_filled_tot = 0;

    int ncid_src;
    int varid_src;
    int ndims_src;
    size_t dimlen_src[4];
    int dimids_src[4];
    int unlimdimid_src = -1;
    nc_type nctype;
    size_t ni_src = 0, nj_src = 0, nij_src = 0;

    float* xsrc = NULL;
    float* ysrc = NULL;
    int* nksrc = NULL;

    float* xsrc_south = NULL;
    float* ysrc_south = NULL;
    float* xsrc_north = NULL;
    float* ysrc_north = NULL;

    char fname_dst_tmp[MAXSTRLEN];

    float* xdst = NULL;
    float* ydst = NULL;
    int* nkdst = NULL;

    float* xdst_south = NULL;
    float* ydst_south = NULL;
    float* xdst_north = NULL;
    float* ydst_north = NULL;

    int ncid_dst;
    int varid_dst;
    size_t dimlen_dst[4];
    int dimids_dst[4];
    size_t ni_dst = 0, nj_dst = 0, nij_dst = 0;

    float* vsrc = NULL;
    float* vdst = NULL;
    float* vdst_last = NULL;

    point* points_south = NULL;
    point* points_north = NULL;

    size_t nk = 0;

    int i, k;

    parse_commandline(argc, argv, &fname_src, &fname_dst, &varname, &grdname_src, &xname_src, &yname_src, &nkname_src, &grdname_dst, &xname_dst, &yname_dst, &nkname_dst, &verbose, &zerofill, &skipfirstlast);

    if (fname_src == NULL)
        quit("no input file specified");
    if (fname_dst == NULL)
        quit("no output file specified");
    if (varname == NULL)
        quit("no variable name specified");
    if (grdname_src == NULL)
        quit("no input grid file specified");
    if (grdname_dst == NULL)
        quit("no output grid file specified");

    if (verbose) {
        printf("  src = \"%s\"\n", fname_src);
        printf("    varname =  \"%s\"\n", varname);
        fflush(stdout);
    }

    ncw_set_quitfn(quit);
    ncu_set_quitfn(quit);

    /*
     * src
     */
    ncw_open(fname_src, NC_NOWRITE, &ncid_src);
    ncw_inq_varid(ncid_src, varname, &varid_src);
    ncw_inq_vardims(ncid_src, varid_src, 4, &ndims_src, dimlen_src);
    ncw_inq_var(ncid_src, varid_src, NULL, &nctype, NULL, dimids_src, NULL);
    if (ncw_var_hasunlimdim(ncid_src, varid_src)) {
        ncw_inq_unlimdim(ncid_src, &unlimdimid_src);
        if (dimlen_src[0] != 1)
            quit("%s: %s: %s can not handle more than one record yet", fname_src, varname, PROGRAM_NAME);
    }
    if (verbose) {
        printf("    size = ");
        for (i = 0; i < ndims_src; ++i)
            printf("%zu%s", dimlen_src[i], i < ndims_src - 1 ? " x " : "\n");
        fflush(stdout);
    }

    /*
     * src grid
     */
    if (verbose) {
        printf("  src grid = \"%s\"\n", grdname_src);
        fflush(stdout);
    }
    {
        int ncid;
        int varid_x, varid_y;
        size_t dimlen[2];
        int ndims;

        ncw_open(grdname_src, NC_NOWRITE, &ncid);
        ncw_inq_varid(ncid, xname_src, &varid_x);
        ncw_inq_varid(ncid, yname_src, &varid_y);
        ncw_inq_varndims(ncid, varid_x, &ndims);
        if (ndims == 2) {
            gridtype = GRIDTYPE_CURV;
            ncw_inq_vardims(ncid, varid_x, 2, NULL, dimlen);
            if (dimlen[0] != dimlen_src[ndims_src - 2] || dimlen[1] != dimlen_src[ndims_src - 1])
                quit("%s,%s: dimensions of variable \"%s\" do not match grid dimensions of coordinate \"%s\"", fname_src, grdname_src, varname, xname_src);

            ni_src = dimlen[1];
            nj_src = dimlen[0];
            nij_src = ni_src * nj_src;
            xsrc = malloc(nij_src * sizeof(float));
            ysrc = malloc(nij_src * sizeof(float));
            ncu_readvarfloat(ncid, varid_x, nij_src, xsrc);
            ncu_readvarfloat(ncid, varid_y, nij_src, ysrc);
        } else if (ndims == 1) {
            ncw_inq_vardims(ncid, varid_x, 1, NULL, &dimlen[1]);
            ncw_inq_vardims(ncid, varid_y, 1, NULL, &dimlen[0]);
            if (dimlen[1] == dimlen_src[ndims_src - 1]
                && dimlen[0] == dimlen_src[ndims_src - 2]) {
                int ii, j;

                gridtype = GRIDTYPE_RECT;
                ni_src = dimlen[1];
                nj_src = dimlen[0];
                nij_src = ni_src * nj_src;
                xsrc = malloc(nij_src * sizeof(float));
                ysrc = malloc(nij_src * sizeof(float));
                ncu_readvarfloat(ncid, varid_x, ni_src, xsrc);
                ncu_readvarfloat(ncid, varid_y, nj_src, ysrc);

                for (ii = 0, j = 0; j < nj_src; ++j)
                    for (i = 0; i < ni_src; ++i, ++ii)
                        xsrc[ii] = xsrc[i];
                for (ii = 0, j = 0; j < nj_src; ++j, ii += ni_src)
                    ysrc[ii] = ysrc[j];
                for (ii = 0, j = 0; j < nj_src; ++j)
                    for (i = 0; i < ni_src; ++i, ++ii)
                        ysrc[ii] = ysrc[j * ni_src];
            } else if (dimlen[0] == dimlen_src[ndims_src - 1] && dimlen[1] == dimlen_src[ndims_src - 1]) {
                gridtype = GRIDTYPE_VECT;
                ni_src = dimlen[0];
                nj_src = 0;
                nij_src = dimlen[0];
                xsrc = malloc(nij_src * sizeof(float));
                ysrc = malloc(nij_src * sizeof(float));
                ncu_readvarfloat(ncid, varid_x, nij_src, xsrc);
                ncu_readvarfloat(ncid, varid_y, nij_src, ysrc);
            } else {
                quit("%s,%s: dimensions of variable \"%s\" do not match grid coordinate(s) \"%s\" or(and) \"%s\"", fname_src, grdname_src, varname, xname_src, yname_src);
            }
            if (nkname_src != NULL) {
                int varid;

                ncw_inq_varid(ncid, nkname_src, &varid);
                if (nj_src > 0) {
                    size_t dimlen[2] = { nj_src, ni_src };

                    ncw_check_vardims(ncid, varid, 2, dimlen);
                } else if (nj_src == 0) {
                    size_t dimlen = ni_src;

                    ncw_check_vardims(ncid, varid, 1, &dimlen);
                } else
                    quit("programming error");
                nksrc = malloc(nij_src * sizeof(int));
                ncw_get_var_int(ncid, varid, nksrc);
            }
        }
        ncw_close(ncid);
    }

    if (verbose) {
        printf("  dst = \"%s\"\n", fname_dst);
        fflush(stdout);
    }
    /*
     * dst grid
     */
    if (verbose) {
        printf("  dst grid = \"%s\"\n", grdname_dst);
        fflush(stdout);
    }
    {
        int ncid;
        int varid_x, varid_y;
        size_t dimlen[2];
        int ndims;

        ncw_open(grdname_dst, NC_NOWRITE, &ncid);
        ncw_inq_varid(ncid, xname_dst, &varid_x);
        ncw_inq_varid(ncid, yname_dst, &varid_y);
        ncw_inq_varndims(ncid, varid_x, &ndims);
        if (ndims == 2) {
            ncw_inq_vardims(ncid, varid_x, 2, NULL, dimlen);
            ni_dst = dimlen[1];
            nj_dst = dimlen[0];
            nij_dst = ni_dst * nj_dst;
            xdst = malloc(nij_dst * sizeof(float));
            ydst = malloc(nij_dst * sizeof(float));
            ncu_readvarfloat(ncid, varid_x, nij_dst, xdst);
            ncu_readvarfloat(ncid, varid_y, nij_dst, ydst);
        } else if (ndims == 1) {
            ncw_inq_vardims(ncid, varid_x, 1, NULL, &dimlen[1]);
            ncw_inq_vardims(ncid, varid_y, 1, NULL, &dimlen[0]);
            if (gridtype != GRIDTYPE_VECT) {
                int ii, j;

                ni_dst = dimlen[1];
                nj_dst = dimlen[0];
                nij_dst = ni_dst * nj_dst;
                xdst = malloc(nij_dst * sizeof(float));
                ydst = malloc(nij_dst * sizeof(float));
                ncu_readvarfloat(ncid, varid_x, ni_dst, xdst);
                ncu_readvarfloat(ncid, varid_y, nj_dst, ydst);

                for (ii = 0, j = 0; j < nj_dst; ++j)
                    for (i = 0; i < ni_dst; ++i, ++ii)
                        xdst[ii] = xdst[i];
                for (ii = 0, j = 0; j < nj_dst; ++j, ii += ni_dst)
                    ydst[ii] = ydst[j];
                for (ii = 0, j = 0; j < nj_dst; ++j)
                    for (i = 0; i < ni_dst; ++i, ++ii)
                        ydst[ii] = ydst[j * ni_dst];
            } else if (gridtype == GRIDTYPE_VECT) {
                if (dimlen[0] != dimlen[1])
                    quit("source grid is unstructured; destination grid is not (coordinates \"%s\" and \"%s\" of different length)", xname_dst, yname_dst);
                ni_dst = dimlen[0];
                nj_dst = 0;
                nij_dst = ni_dst;
                xdst = malloc(nij_dst * sizeof(float));
                ydst = malloc(nij_dst * sizeof(float));
                ncu_readvarfloat(ncid, varid_x, nij_dst, xdst);
                ncu_readvarfloat(ncid, varid_y, nij_dst, ydst);
            }
            if (nkname_dst != NULL) {
                int varid;

                ncw_inq_varid(ncid, nkname_dst, &varid);
                if (nj_dst > 0) {
                    size_t dimlen[2] = { nj_dst, ni_dst };

                    ncw_check_vardims(ncid, varid, 2, dimlen);
                } else if (nj_dst == 0) {
                    size_t dimlen = ni_dst;

                    ncw_check_vardims(ncid, varid, 1, &dimlen);
                } else
                    quit("programming error");
                nksrc = malloc(nij_dst * sizeof(int));
                ncw_get_var_int(ncid, varid, nksrc);
            }
        }
        ncw_close(ncid);
    }

    strcpy(fname_dst_tmp, fname_dst);
    strcat(fname_dst_tmp, ".tmp");
    ncw_create(fname_dst_tmp, NC_CLOBBER | NC_NETCDF4, &ncid_dst);
    /*
     * copy global attributes to destination
     */
    ncw_copy_atts(ncid_src, NC_GLOBAL, ncid_dst, NC_GLOBAL);
    for (i = 0; i < ndims_src; ++i) {
        char dimname[NC_MAX_NAME];

        ncw_inq_dimname(ncid_src, dimids_src[i], dimname);
        if (ndims_src - i == 1)
            dimlen_dst[i] = ni_dst;
        else if (ndims_src - i == 2)
            dimlen_dst[i] = nj_dst;
        else if (dimids_src[i] == unlimdimid_src)
            dimlen_dst[i] = 1;
        else {
            nk = dimlen_src[i];
            dimlen_dst[i] = dimlen_src[i];
        }
        ncw_def_dim(ncid_dst, dimname, dimlen_dst[i], &dimids_dst[i]);
    }
    if (verbose) {
        printf("    size = ");
        for (i = 0; i < ndims_src; ++i)
            printf("%zu%s", dimlen_dst[i], i < ndims_src - 1 ? " x " : "\n");
        fflush(stdout);
    }
    ncw_def_var(ncid_dst, varname, nctype, ndims_src, dimids_dst, &varid_dst);
    ncw_copy_atts(ncid_src, varid_src, ncid_dst, varid_dst);
    /*
     * copy command line and wdir to dst
     */
    {
        char* cmd = get_command(argc, argv);
        char attname[NC_MAX_NAME];
        char cwd[MAXSTRLEN];

        snprintf(attname, NC_MAX_NAME, "%s: command", PROGRAM_NAME);
        ncw_put_att_text(ncid_dst, NC_GLOBAL, attname, cmd);
        free(cmd);
        if (getcwd(cwd, MAXSTRLEN) != NULL) {
            snprintf(attname, NC_MAX_NAME, "%s: wdir", PROGRAM_NAME);
            ncw_put_att_text(ncid_dst, NC_GLOBAL, attname, cwd);
        }
    }

    ncw_close(ncid_dst);
    ncw_close(ncid_src);

    if (verbose) {
        printf("  converting src lon/lat to stereographic projections:");
        fflush(stdout);
    }
    xsrc_south = malloc(nij_src * sizeof(float));
    ysrc_south = malloc(nij_src * sizeof(float));
    xsrc_north = malloc(nij_src * sizeof(float));
    ysrc_north = malloc(nij_src * sizeof(float));
    for (i = 0; i < nij_src; ++i) {
        double ll[2] = { xsrc[i], -ysrc[i] };
        double xyz[3];

        ll2xyz(ll, xyz);
        xsrc_south[i] = xyz[0] / (1.0 - xyz[2]);
        ysrc_south[i] = xyz[1] / (1.0 - xyz[2]);

        ll[1] = ysrc[i];
        ll2xyz(ll, xyz);
        xsrc_north[i] = xyz[0] / (1.0 - xyz[2]);
        ysrc_north[i] = xyz[1] / (1.0 - xyz[2]);
    }
    free(xsrc);
    free(ysrc);
    if (verbose) {
        printf("\n");
        fflush(stdout);
    }

    if (verbose) {
        printf("  converting dst lon/lat to stereographic projections:");
        fflush(stdout);
    }
    xdst_south = malloc(nij_dst * sizeof(float));
    ydst_south = malloc(nij_dst * sizeof(float));
    xdst_north = malloc(nij_dst * sizeof(float));
    ydst_north = malloc(nij_dst * sizeof(float));
    for (i = 0; i < nij_dst; ++i) {
        double ll[2] = { xdst[i], -ydst[i] };
        double xyz[3];

        ll2xyz(ll, xyz);
        xdst_south[i] = xyz[0] / (1.0 - xyz[2]);
        ydst_south[i] = xyz[1] / (1.0 - xyz[2]);

        ll[1] = ydst[i];
        ll2xyz(ll, xyz);
        xdst_north[i] = xyz[0] / (1.0 - xyz[2]);
        ydst_north[i] = xyz[1] / (1.0 - xyz[2]);
    }
    free(xdst);
    if (verbose) {
        printf("\n");
        fflush(stdout);
    }

    if (verbose) {
        printf("  interpolating:");
        fflush(stdout);
    }

    vsrc = malloc(nij_src * sizeof(float));
    vdst = calloc(nij_dst, sizeof(float));
    if (nk > 1 && !zerofill) {
        vdst_last = malloc(nij_dst * sizeof(float));
        for (i = 0; i < nij_dst; ++i)
            vdst_last[i] = NAN;
    }

    if (nk == 0)
        nk = 1;

    points_south = malloc(nij_src * sizeof(point));
    points_north = malloc(nij_src * sizeof(point));
    for (k = 0; k < nk; ++k) {
        int npoint = 0;

        /*
         * stats
         */
        int npoint_dst = 0;
        int npoint_filled = 0;

        ncu_readfield(fname_src, varname, k, ni_src, nj_src, nk, vsrc);
        for (i = 0; i < nij_src; ++i) {
            point* p;
            int ii = i;

            /*
             * do not use the first and last columns
             */
            if (skipfirstlast && (i % ni_src == 0 || i % ni_src == ni_src - 1))
                continue;

            if (nksrc != NULL && k >= nksrc[i])
                continue;
            if (!isfinite(vsrc[ii]))
                continue;

            p = &points_south[npoint];
            p->x = xsrc_south[i];
            p->y = ysrc_south[i];
            p->z = vsrc[ii];

            p = &points_north[npoint];
            p->x = xsrc_north[i];
            p->y = ysrc_north[i];
            p->z = vsrc[ii];

            npoint++;
        }

        memset(vdst, 0, nij_dst * sizeof(float));
        if (npoint == 0)
            goto finalise_level;

        {
            delaunay* d_south = delaunay_build(npoint, points_south, 0, NULL, 0, NULL);
            delaunay* d_north = delaunay_build(npoint, points_north, 0, NULL, 0, NULL);
            void* interp_south = lpi_build(d_south);
            void* interp_north = lpi_build(d_north);

            for (i = 0; i < nij_dst; ++i) {
                point p;

                if (ydst[i] > 0.0) {
                    p.x = xdst_south[i];
                    p.y = ydst_south[i];
                } else {
                    p.x = xdst_north[i];
                    p.y = ydst_north[i];
                }

                if (nkdst == NULL || k < nkdst[i]) {
                    npoint_dst++;
                    if (ydst[i] > 0.0)
                        lpi_interpolate_point(interp_south, &p);
                    else
                        lpi_interpolate_point(interp_north, &p);

                    if (isfinite(p.z)) {
                        vdst[i] = (float) p.z;
                        if (vdst_last != NULL)
                            vdst_last[i] = (float) p.z;
                    } else {
                        npoint_filled_tot++;
                        npoint_filled++;
                        if (vdst_last != NULL && isfinite(vdst_last[i]))
                            vdst[i] = vdst_last[i];
                    }
                }
            }
            lpi_destroy(interp_south);
            delaunay_destroy(d_south);
            lpi_destroy(interp_north);
            delaunay_destroy(d_north);
        }

      finalise_level:

        ncu_writefield(fname_dst_tmp, varname, k, ni_dst, nj_dst, nk, vdst);
        if (verbose == 1) {
            printf("%c", (k + 1) % 10 ? '.' : '|');
            fflush(stdout);
        } else if (verbose > 1) {
            printf("\n    k = %d: %d in, %d out", k, npoint, npoint_dst);
            if (npoint_filled > 0)
                printf(" (%d filled)", npoint_filled);
            fflush(stdout);
        }
    }
    file_rename(fname_dst_tmp, fname_dst);
    if (verbose) {
        printf("\n");
        printf("  -> %s\n", fname_dst);
        if (verbose > 1)
            printf("  # cells filled = %d\n", npoint_filled_tot);
        fflush(stdout);
    }

    free(points_south);
    free(points_north);
    free(vsrc);
    free(vdst);
    free(ydst);
    free(xsrc_south);
    free(ysrc_south);
    free(xsrc_north);
    free(ysrc_north);
    free(xdst_south);
    free(ydst_south);
    free(xdst_north);
    free(ydst_north);
    if (vdst_last != NULL)
        free(vdst_last);

    return 0;
}
