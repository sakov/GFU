/******************************************************************************
 *
 * File:        ncd2f.c
 *
 * Created:     15/08/2024
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: NCD2F reads a variable (or all variables) of type double from a
 *              NetCDF file, casts it (them) to float and saves to another file.
 *
 *              The data is assumed to be in NetCDF format.
 *             
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include "version.h"
#include "ncw.h"
#include "ncutils.h"
#include "utils.h"

#define PROGRAM_NAME "ncd2f"
#define PROGRAM_VERSION "0.01"

#define VERBOSE_DEF 1

#define NDIM_MIN_DEF 2
#define DIMNAME_NTRIES 10
#define NVAR_INC 100
#define MAXSIZE (1024 * 1024 * 256)
#define TEMPVARSUF "_d2f_tmp"

int verbose = VERBOSE_DEF;

/**
 */
static void usage(int status)
{
    printf("  Usage: %s -i <src> [{-v <var> [...] | -d <N>}] -o <dst>\n", PROGRAM_NAME);
    printf("         %s -v\n", PROGRAM_NAME);
    printf("  Options:\n");
    printf("    -i <src>       -- source file\n");
    printf("    -o <dst>       -- destination file\n");
    printf("    -v <var> [...] -- variable to cast (default: all variables of type double\n");
    printf("                      that have 2 or more dimensions)\n");
    printf("    -d <N>         -- minimal number of dimensions for a variable to be casted\n");
    printf("    -v             -- print version and exit\n");
    exit(status);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname_src, char** fname_dst, int* nvar, char*** vars, int* ndim_min)
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
        if (strcmp(argv[i], "-i") == 0) {
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no file name found after \"-i\"");
            *fname_src = argv[i];
            i++;
        } else if (strcmp(argv[i], "-o") == 0) {
            i++;
            if (i == argc || argv[i][0] == '-')
                quit("no file name found after \"-o\"");
            *fname_dst = argv[i];
            i++;
        } else if (strcmp(argv[i], "-v") == 0) {
            i++;
            if (i >= argc)
                quit("no variable specified after \"-v\"\n");
            if (*nvar % NVAR_INC == 0)
                *vars = realloc(*vars, (*nvar + NVAR_INC) * sizeof(void*));
            (*vars)[*nvar] = strdup(argv[i]);
            (*nvar)++;
            i++;
        } else if (argv[i][1] == 'd') {
            if (*vars != NULL)
                quit("can not use both \"-v\" and \"-d\"");
            i++;
            if (!str2int(argv[i], ndim_min))
                quit("could not convert \"%s\" to int", argv[i]);
            i++;
        } else
            quit("unknown option \"%s\"", argv[i]);
    }
}

/**
 */
static int copy_vardef_newtype(int ncid_src, int varid_src, int ncid_dst, char* varname_dst, nc_type newtype)
{
    int unlimdimid_src = -1;
    char varname[NC_MAX_NAME];
    int varid_dst;
    nc_type type;
    int ndims;
    int dimids_src[NC_MAX_DIMS], dimids_dst[NC_MAX_DIMS];
    int natts;
    int status;
    int i, ii;

    status = nc_redef(ncid_dst);

    ncw_inq_unlimdim(ncid_src, &unlimdimid_src);
    ncw_inq_varname(ncid_src, varid_src, varname);
    ncw_inq_var(ncid_src, varid_src, NULL, &type, &ndims, dimids_src, &natts);

    for (i = 0; i < ndims; ++i) {
        char dimname[NC_MAX_NAME];
        size_t len;

        dimids_dst[i] = -1;
        ncw_inq_dim(ncid_src, dimids_src[i], dimname, &len);

#if defined(NCW_SKIPSINGLE)
        /*
         * skip "normal" (not unlimited) "inner" dimensions of length 1
         */
        if (len == 1 && i < ndims - 1 && dimids_src[i] != unlimdimid_src)
            continue;
#endif

        if (!ncw_dim_exists(ncid_dst, dimname))
            ncw_def_dim(ncid_dst, dimname, (dimids_src[i] == unlimdimid_src) ? NC_UNLIMITED : len, &dimids_dst[i]);
        else {
            int dimid_dst;
            size_t len_dst;
            char dimname_dst[MAXSTRLEN];
            int format;
            int j;

            ncw_inq_dimid(ncid_dst, dimname, &dimid_dst);
            ncw_inq_dimlen(ncid_dst, dimid_dst, &len_dst);
            if (len == len_dst) {       /* (all good) */
                dimids_dst[i] = dimid_dst;
                continue;
            }

            ncw_inq_format(ncid_dst, &format);
            if (format == NC_FORMAT_NETCDF4) {
                int nunlimdims = -1;
                int unlimdimids[NC_MAX_DIMS];
                int d;

                (void) nc_inq_unlimdims(ncid_dst, &nunlimdims, unlimdimids);
                for (d = 0; d < nunlimdims; ++d)
                    if (dimid_dst == unlimdimids[d])
                        break;
                if (d < nunlimdims) {
                    dimids_dst[i] = dimid_dst;
                    continue;
                }
            } else {
                int unlimdimid_dst = -1;

                ncw_inq_unlimdim(ncid_dst, &unlimdimid_dst);
                if (dimid_dst == unlimdimid_dst) {
                    dimids_dst[i] = dimid_dst;
                    continue;
                }
            }

            /*
             * So... there is this dimension in the destination file, but it has
             * the wrong length. We will define and use another dimension then.
             */
            for (j = 0; j < DIMNAME_NTRIES; ++j) {
                {
                    char buf[MAXSTRLEN];

                    snprintf(buf, MAXSTRLEN, "%s%d", dimname, j);
                    if (strlen(buf) < NC_MAX_NAME)
                        strcpy(dimname_dst, buf);
                    else {
                        int diff = strlen(buf) - NC_MAX_NAME + 1;

                        dimname[NC_MAX_NAME - diff] = 0;
                        sprintf(dimname_dst, "%s%d", dimname, j);
                    }
                }
                if (ncw_dim_exists(ncid_dst, dimname_dst)) {
                    ncw_inq_dimid(ncid_dst, dimname_dst, &dimid_dst);
                    ncw_inq_dimlen(ncid_dst, dimid_dst, &len_dst);
                    if (len_dst == len) {       /* (all good) */
                        dimids_dst[i] = dimid_dst;
                        break;
                    }
                } else
                    break;      /* this name is unique, dimids_dst[i] remains 
                                 * -1 */
            }
            if (j == DIMNAME_NTRIES) {  /* (error) */
                char fname_dst[MAXSTRLEN];

                strncpy(fname_dst, ncw_get_path(ncid_dst), MAXSTRLEN - 1);
                ncw_close(ncid_dst);    /* (to be able to examine the file) */
                quit("\"%s\": ncw_copy_vardef(): technical problem copying \"%s\" from \"%s\"\n", fname_dst, varname, ncw_get_path(ncid_src));
            }
            if (dimids_dst[i] >= 0)     /* (all good) */
                continue;
            ncw_def_dim(ncid_dst, dimname_dst, (dimids_src[i] == unlimdimid_src) ? NC_UNLIMITED : len, &dimids_dst[i]);
        }
    }

    for (i = 0, ii = 0; i < ndims; ++i) {
        if (dimids_dst[i] < 0)
            continue;
        dimids_dst[ii] = dimids_dst[i];
        ii++;
    }
    ndims = ii;

    ncw_def_var(ncid_dst, varname_dst, newtype, ndims, dimids_dst, &varid_dst);
    ncw_copy_atts(ncid_src, varid_src, ncid_dst, varid_dst);

    ncw_def_var_deflate(ncid_dst, varid_dst, 0, 1, 1);
    
    if (status == NC_NOERR)
        nc_enddef(ncid_dst);

    return varid_dst;
}

/**
 */
int main(int argc, char* argv[])
{
    char* fname_src = NULL;
    char* fname_dst = NULL;
    char** varnames_src = NULL;
    int nvar = 0;
    int ndim_min = NDIM_MIN_DEF;

    char** varnames_dst = NULL;
    char* fname_dst_tmp = NULL;
    int ncid_src, ncid_dst;
    int vid;

    parse_commandline(argc, argv, &fname_src, &fname_dst, &nvar, &varnames_src, &ndim_min);

    if (fname_src == NULL)
        quit("no input file specified");
    if (fname_dst == NULL)
        quit("no output file specified");
    
    ncw_set_quitfn(quit);
    ncu_set_quitfn(quit);

    ncw_open(fname_src, NC_NOWRITE, &ncid_src);
    
    /*
     * get variables to convert
     */
    if (varnames_src == NULL) {
        int nvar_all, vid;
        
        nvar = 0;
        ncw_inq_nvars(ncid_src, &nvar_all);
        for (vid = 0; vid < nvar_all; ++vid) {
            char name[NC_MAX_NAME];
            nc_type type;
            int ndim;
            
            ncw_inq_var(ncid_src, vid, name, &type, &ndim, NULL, NULL);
            if (type != NC_DOUBLE || ndim < ndim_min)
                continue;
            if (nvar % NVAR_INC == 0)
                varnames_src = realloc(varnames_src, (nvar + NVAR_INC) * sizeof(void*));
            varnames_src[nvar] = strdup(name);
            nvar++;
        }
    }
    
    /*
     * open destination
     */
    if (verbose)
        printf("  %s:\n", fname_dst);
    if (file_exists(fname_dst)) {
        ncw_open(fname_dst, NC_WRITE, &ncid_dst);
        for (vid = 0; vid < nvar; ++vid)
            if (ncw_var_exists(ncid_dst, varnames_src[vid]))
                quit("%s: variable \"%s\" already exists", fname_dst, varnames_src[vid]);
    } else {
        fname_dst_tmp = malloc(MAXSTRLEN);
        strncpy(fname_dst_tmp, fname_dst, MAXSTRLEN - 1);
        if (strlen(fname_dst_tmp) < MAXSTRLEN - 4)
            strcat(fname_dst_tmp, ".tmp");
        else
            quit("destination too long");
        ncw_create(fname_dst_tmp, NC_CLOBBER | NC_NETCDF4, &ncid_dst);
    }
    
    /*
     * if appending to an existing file -- write to a temporary variable name
     * first
     */
    if (fname_dst_tmp == NULL) {
        int vid;
        
        varnames_dst = malloc(nvar * sizeof(void*));
        for (vid = 0; vid < nvar; ++vid) {
            varnames_dst[vid] = malloc(strlen(varnames_src[vid]) + strlen(TEMPVARSUF) + 1);
            strcpy(varnames_dst[vid], varnames_src[vid]);
            strcat(varnames_dst[vid], TEMPVARSUF);
        }
    } else
        varnames_dst = varnames_src;

    /*
     * copy global attributes to dst
     */
    ncw_copy_atts(ncid_src, NC_GLOBAL, ncid_dst, NC_GLOBAL);
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

    for (vid = 0; vid < nvar; ++vid) {
        int varid_src, varid_dst;
        size_t size;
        float* v;
 
        ncw_open(fname_src, NC_NOWRITE, &ncid_src);
        if (ncw_var_exists(ncid_dst, varnames_src[vid]))
            quit("%s: variable \"%s\" already exists", fname_dst, varnames_src[vid]);
        if (verbose) {
            printf("    %s:", varnames_src[vid]);
            fflush(stdout);
        }
        ncw_inq_varid(ncid_src, varnames_src[vid], &varid_src);
        if (!ncw_var_exists(ncid_dst, varnames_dst[vid]))
            varid_dst = copy_vardef_newtype(ncid_src, varid_src, ncid_dst, varnames_dst[vid], NC_FLOAT);
        else
            ncw_inq_varid(ncid_dst, varnames_dst[vid], &varid_dst);
#if 0        
        {
            int natt, i;

            ncw_inq_varnatts(ncid_dst, varid_dst, &natt);
            for (i = 0; i < natt; ++i) {
                char attname[NC_MAX_NAME];

                ncw_inq_attname(ncid_dst, varid_dst, i, attname);
                if (strcmp(attname, "_FillValue") == 0 || strcmp(attname, "valid_min") == 0 || strcmp(attname, "valid_max") == 0 || strcmp(attname, "valid_range") == 0 || strcmp(attname, "scale_factor") == 0 || strcmp(attname, "add_offset") == 0 || strcmp(attname, "missing_value") == 0) {
                    ncw_del_att(ncid_dst, varid_dst, attname);
                    i--;
                    natt--;
                }
            }
        }
#endif

        size = ncw_get_varsize(ncid_src, varid_src);
        if (size <= MAXSIZE) {
            int use_putvara = 0;
            
            v = malloc(size * ncw_sizeof(NC_FLOAT));
            ncu_readvarfloat(ncid_src, varid_src, size, v);
            
            if (ncw_var_hasunlimdim(ncid_dst, varid_dst)) {
                int nr_src = ncw_inq_nrecords(ncid_src);
                int nr_dst = ncw_inq_nrecords(ncid_dst);

                if (nr_src > nr_dst)
                    use_putvara = 1;
            }

            if (!use_putvara)
                ncw_put_var_float(ncid_dst, varid_dst, v);
            else {
                int ndims;
                size_t dimlen[NC_MAX_VAR_DIMS];
                size_t start[NC_MAX_VAR_DIMS];
                int i;

                ncw_inq_vardims(ncid_src, varid_src, NC_MAX_DIMS, &ndims, dimlen);
                for (i = 0; i < ndims; ++i)
                    start[i] = 0;
                ncw_put_vara_float(ncid_dst, varid_dst, start, dimlen, v);
            }
            free(v);
        } else {
            int nk = ncu_getnfields(fname_src, varnames_src[vid]);
            int k;

            assert(nk > 0);
            v = malloc(size / nk * ncw_sizeof(NC_FLOAT));
            for (k = 0; k < nk; ++k) {
                ncu_readfield(fname_src, varnames_src[vid], k, -1, -1, nk, v);
                ncu_writefield((fname_dst_tmp == NULL) ? fname_dst : fname_dst_tmp, varnames_dst[vid], k, -1, -1, nk, v);
                if (verbose) {
                    printf(".");
                    fflush(stdout);
                }
            }
            free(v);
        }
        ncw_sync(ncid_dst);
        ncw_redef(ncid_dst);
        if (verbose)
            printf("\n");
    }
    if (fname_dst_tmp == NULL) {
        for (vid = 0; vid < nvar; ++vid)
            ncw_rename_var(ncid_dst, varnames_dst[vid], varnames_src[vid]);
     }
    ncw_close(ncid_dst);

    if (fname_dst_tmp != NULL)
        file_rename(fname_dst_tmp, fname_dst);
    
    if (varnames_dst != varnames_src) {
        for (vid = 0; vid < nvar; ++vid)
            free(varnames_dst[vid]);
        free(varnames_dst);
    }
    for (vid = 0; vid < nvar; ++vid)
        free(varnames_src[vid]);
    free(varnames_src);

    return 0;
}
