/******************************************************************************
 *
 * File:        nccat.c        
 *
 * Created:     19/05/2023
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: NCCAT concatenates variables in NetCDF files over arbitrary
 *              dimensions.
 *
 * Revisions:   
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include "ncw.h"
#include "utils.h"
#include "version.h"

#define PROGRAM_NAME "nccat"
#define PROGRAM_VERSION "0.00"

#define NINC 10
#define VERBOSE 0

/**
 */
static void usage(int status)
{
    printf("  Usage: nccat [-v <var> [...]] [-d <dim> [...]] -i <src> [...] -o <dst> [-V <level>] \n");
    printf("         nccat -v\n");
    printf("  Options:\n");
    printf("    -v <var> [...] - variables to be concatenated (default: all)\n");
    printf("    -d <dim> [...] - dimensions to concatenate (default: those of different length)\n");
    printf("    -i <src> [...] - source files\n");
    printf("    -o <dst>       - destination file\n");
    printf("    -V <level>     - verbosity level (0 to 2)\n");
    printf("    -v             - print version and exit\n");
    exit(status);
}

/**
 */
static void parse_commandline(int argc, char* argv[], int* nvar, char*** vars, int* ndim, char*** dims, int* nsrc, char*** srcs, char** dst, int* verbose)
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
        if (argv[i][0] != '-')
            usage(1);
        if (argv[i][1] == 'v') {
            i++;
            while (i < argc && argv[i][0] != '-') {
                if (*nvar % NINC == 0)
                    *vars = realloc(*vars, (*nvar + NINC) * sizeof(void*));
                (*vars)[*nvar] = strdup(argv[i]);
                (*nvar)++;
                i++;
            }
        } else if (argv[i][1] == 'd') {
            i++;
            while (i < argc && argv[i][0] != '-') {
                if (*ndim % NINC == 0)
                    *dims = realloc(*dims, (*ndim + NINC) * sizeof(void*));
                (*dims)[*ndim] = argv[i];
                (*ndim)++;
                i++;
            }
        } else if (argv[i][1] == 'i') {
            i++;
            while (i < argc && argv[i][0] != '-') {
                if (*nsrc % NINC == 0)
                    *srcs = realloc(*srcs, (*nsrc + NINC) * sizeof(void*));
                (*srcs)[*nsrc] = argv[i];
                (*nsrc)++;
                i++;
            }
        } else if (argv[i][1] == 'o') {
            i++;
            *dst = argv[i];
            i++;
        } else if (argv[i][1] == 'V') {
            i++;
            if (!str2int(argv[i], verbose))
                quit("could not convert \"%s\" to int");
            i++;
        } else
            quit("can not interpret argument \"%s\"", argv[i]);
    }

    if (*nsrc == 0)
        quit("no input files specified");
    if (*nsrc == 1)
        quit("only one input file specified; nothing to do");
    if (*dst == NULL)
        quit("no output file specified");
}

/**
 */
static int time_units_changed(char* tunits0, int ncid, int varid)
{
    size_t attlen;
    char tunits[MAXSTRLEN];

    ncw_inq_attlen(ncid, varid, "units", &attlen);
    assert(attlen < MAXSTRLEN);
    ncw_get_att_text(ncid, varid, "units", tunits);

    return strcmp(tunits0, tunits);
}

/**
 */
static void adjusttime(char* tunits0, int ncid, int varid, size_t vlen, void* v)
{
    size_t attlen;
    char tunits[MAXSTRLEN];
    double multiple, offset;
    nc_type nctype;
    size_t i;

    ncw_inq_attlen(ncid, varid, "units", &attlen);
    assert(attlen < MAXSTRLEN);
    ncw_get_att_text(ncid, varid, "units", tunits);
    tunits_convert(tunits0, tunits, &multiple, &offset);
    ncw_inq_vartype(ncid, varid, &nctype);
    if (nctype == NC_FLOAT) {
        for (i = 0; i < vlen; ++i)
            ((float*) v)[i] = ((float*) v)[i] * multiple + offset;
    } else if (nctype == NC_DOUBLE) {
        for (i = 0; i < vlen; ++i)
            ((double*) v)[i] = ((double*) v)[i] * multiple + offset;
    } else
        quit("time variable can not be adjusted for data types other than NC_FLOAT or NC_DOUBLE");
}

/**
 */
int main(int argc, char** argv)
{
    int nvar = 0;
    char** vars = NULL;
    int ndim_force = 0;
    char** dims_force = NULL;
    int nsrc = 0;
    char** srcs = NULL;
    char* dst = NULL;
    int verbose = VERBOSE;

    int* ncids_src = NULL;
    int* varids_src = NULL;
    int sid, vid;

    char* tmpdst = NULL;
    int ncid_dst;
    int i;

    parse_commandline(argc, argv, &nvar, &vars, &ndim_force, &dims_force, &nsrc, &srcs, &dst, &verbose);

    if (verbose > 1) {
        printf("  %s v%s\n", PROGRAM_NAME, PROGRAM_VERSION);
        printf("  GFU v%s\n", VERSION);
    }

    ncids_src = malloc(nsrc * sizeof(int));
    for (sid = 0; sid < nsrc; ++sid)
        ncw_open(srcs[sid], NC_NOWRITE, &ncids_src[sid]);
    for (i = 0; i < ndim_force; ++i)
        if (!ncw_dim_exists(ncids_src[0], dims_force[i]))
            quit("%s: no dimension \"%s\"", srcs[i], dims_force[i]);

    if (nvar == 0) {
        ncw_inq_nvars(ncids_src[0], &nvar);
        if (nvar == 0)
            quit("%s: no variables found", ncids_src[0]);
        vars = malloc(nvar * sizeof(void*));
        for (vid = 0; vid < nvar; ++vid) {
            char varname[NC_MAX_NAME];

            ncw_inq_varname(ncids_src[0], vid, varname);
            vars[vid] = strdup(varname);
        }
    }

    tmpdst = malloc(strlen(dst) + 32);
    snprintf(tmpdst, strlen(dst) + 32, "%s.pid%zu.nccat.tmp", dst, (size_t) getpid());
    ncw_create(tmpdst, NC_CLOBBER | NC_NETCDF4, &ncid_dst);
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

    /*
     * main cycle
     */
    varids_src = malloc(nsrc * sizeof(int));
    for (vid = 0; vid < nvar; ++vid) {
        int ndim;
        size_t** dimlens_src = NULL;
        size_t dimlens_dst[NC_MAX_DIMS];
        int dimids_src[NC_MAX_DIMS];
        int dimids_dst[NC_MAX_DIMS];
        nc_type nctype;
        int varid_dst;
        int did_merge;
        int did;
        void* v_dst = NULL;
        void** v_src = NULL;
        size_t typesize, vlen_dst, slicesize;
        size_t offset_dst;
        size_t* offsets_src = NULL;
        int istime = 0;
        char tunits0[MAXSTRLEN];
        int i, ni;

        for (sid = 0; sid < nsrc; ++sid) {
            ncw_inq_varid(ncids_src[sid], vars[vid], &varids_src[sid]);
            if (sid == 0)
                ncw_inq_varndims(ncids_src[sid], varids_src[sid], &ndim);
            else
                ncw_check_varndims(ncids_src[sid], varids_src[sid], ndim);
        }

        dimlens_src = alloc2d(nsrc, ndim, sizeof(size_t));
        for (sid = 0; sid < nsrc; ++sid)
            ncw_inq_vardims(ncids_src[sid], varids_src[sid], ndim, NULL, dimlens_src[sid]);
        did_merge = -1;
        for (did = 0; did < ndim; ++did) {
            for (sid = 1; sid < nsrc; ++sid) {
                if (dimlens_src[sid][did] != dimlens_src[0][did]) {
                    if (did_merge < 0)
                        did_merge = did;
                    else if (did != did_merge)
                        quit("can not concatenate variable \"%s\": dimension sizes in source files are different for more than one dimension", vars[vid]);
                }
            }
        }

        ncw_inq_vardimid(ncids_src[0], varids_src[0], dimids_src);
        for (i = 0; i < ndim_force; ++i) {
            int dimid_force;

            ncw_inq_dimid(ncids_src[0], dims_force[i], &dimid_force);
            for (did = 0; did < ndim; ++did) {
                if (dimid_force == dimids_src[did]) {
                    if (did_merge < 0)
                        did_merge = did;
                    else {
                        if (did != did_merge) {
                            char dimname[NC_MAX_NAME];

                            ncw_inq_dimname(ncids_src[0], dimids_src[did_merge], dimname);
                            quit("can not merge variable \"%s\" on more than one dimension (\"%s\" and \"%s\")", dims_force[i], dimname);
                        }
                    }
                }
            }
        }

        /*
         * set dimensions
         */
        for (did = 0; did < ndim; ++did)
            dimlens_dst[did] = dimlens_src[0][did];
        if (did_merge >= 0)
            for (sid = 1; sid < nsrc; ++sid)
                dimlens_dst[did_merge] += dimlens_src[sid][did_merge];
        for (did = 0; did < ndim; ++did) {
            char dimname[NC_MAX_NAME];

            ncw_inq_dimname(ncids_src[0], dimids_src[did], dimname);
            if (ncw_dim_exists(ncid_dst, dimname)) {
                ncw_check_dimlen(ncid_dst, dimname, dimlens_dst[did]);
                ncw_inq_dimid(ncid_dst, dimname, &dimids_dst[did]);
            } else
                ncw_def_dim(ncid_dst, dimname, dimlens_dst[did], &dimids_dst[did]);
        }

        if (verbose > 0) {
            printf("  %s%s", (verbose > 1) ? "  " : "", vars[vid]);
            for (did = 0; did < ndim; ++did) {
                char dimname[NC_MAX_NAME];

                ncw_inq_dimname(ncids_src[0], dimids_src[did], dimname);
                printf("%s%s%s", (did == 0) ? "(" : ", ", dimname, (did == ndim - 1) ? ")" : "");
            }
            if (did_merge >= 0) {
                char dimname[NC_MAX_NAME];

                ncw_inq_dimname(ncids_src[0], dimids_src[did_merge], dimname);
                printf(" - merged by \"%s\"", dimname);
                if (verbose > 1) {
                    for (sid = 0; sid < nsrc; ++sid)
                        printf("%s%zu%s", (sid == 0) ? " (" : " + ", dimlens_src[sid][did_merge], (sid == nsrc - 1) ? ")" : "");
                }
                printf("\n");
            } else
                printf(" - just copied\n");
        }

        /*
         * (1) copy variable definition
         */
        ncw_inq_vartype(ncids_src[0], varids_src[0], &nctype);
        ncw_def_var(ncid_dst, vars[vid], nctype, ndim, dimids_dst, &varid_dst);
        ncw_copy_atts(ncids_src[0], varids_src[0], ncid_dst, varid_dst);
        ncw_enddef(ncid_dst);

        /*
         * (2) copy variable data
         */
        typesize = ncw_sizeof(nctype);
        vlen_dst = dimlens_dst[0];
        for (did = 1; did < ndim; ++did)
            vlen_dst *= dimlens_dst[did];

        v_dst = malloc(vlen_dst * typesize);
        /*
         * check if merging is needed for this variable
         */
        if (did_merge < 0) {
            /*
             * no merging  -- just copy the data from the first source
             */
            ncw_get_var(ncids_src[0], varids_src[0], v_dst);
            ncw_put_var(ncid_dst, varid_dst, v_dst);
            ncw_redef(ncid_dst);

            goto nextvar;
        }

        istime = varistime(ncids_src[0], vid);
        if (istime) {
            size_t attlen;

            ncw_inq_attlen(ncids_src[0], varids_src[0], "units", &attlen);
            assert(attlen < MAXSTRLEN);
            ncw_get_att_text(ncids_src[0], varids_src[0], "units", tunits0);
        }

        /*
         * read data
         */
        v_src = malloc(nsrc * sizeof(void*));
        for (sid = 0; sid < nsrc; ++sid) {
            size_t vlen_src = vlen_dst / dimlens_dst[did_merge] * dimlens_src[sid][did_merge];

            v_src[sid] = malloc(vlen_src * typesize);
            ncw_get_var(ncids_src[sid], varids_src[sid], v_src[sid]);
            if (sid > 0 && istime && time_units_changed(tunits0, ncids_src[sid], varids_src[sid]))
                adjusttime(tunits0, ncids_src[sid], varids_src[sid], vlen_src, v_src[sid]);
        }

        slicesize = 1;
        for (did = did_merge + 1; did < ndim; ++did)
            slicesize *= dimlens_dst[did];

        /*
         * merge source data into destination array 
         */
        ni = 1;
        for (did = 0; did < did_merge; ++did)
            ni *= dimlens_dst[did];
        offset_dst = 0;
        offsets_src = calloc(nsrc, sizeof(size_t));
        for (i = 0; i < ni; ++i) {
            for (sid = 0; sid < nsrc; ++sid) {
                size_t size = slicesize * dimlens_src[sid][did] * typesize;

                memcpy(&((char*) v_dst)[offset_dst], &((char*) v_src[sid])[offsets_src[sid]], size);
                offset_dst += size;
                offsets_src[sid] += size;
            }
        }

        /*
         * write
         */
        ncw_put_var(ncid_dst, varid_dst, v_dst);
        ncw_redef(ncid_dst);

        /*
         * free memory
         */
        free(offsets_src);
        for (sid = 0; sid < nsrc; ++sid)
            free(v_src[sid]);
        free(v_src);
      nextvar:
        free(v_dst);
        free(dimlens_src);
        free(vars[vid]);
    }

    ncw_enddef(ncid_dst);
    ncw_close(ncid_dst);
    file_rename(tmpdst, dst);
    free(tmpdst);

    free(varids_src);
    for (sid = 0; sid < nsrc; ++sid)
        ncw_close(ncids_src[sid]);
    free(ncids_src);
    free(srcs);
    free(vars);
    if (dims_force != NULL)
        free(dims_force);

    return 0;
}
