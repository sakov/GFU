/******************************************************************************
 *
 * File:        ncave.c        
 *
 * Created:     28/09/2022
 *
 * Authors:     Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: NCAVE targets averaging very large ensemble dumps. Compared to
 *              NCEA/NCRA it (1) conserves memory by averaging on layer-by-layer
 *              basis, and (2) when run on multiple CPUs processes layers in
 *              parallel. Note that for optimal efficiency the input dumps must
 *              be chunked by layers.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ftw.h>
#include <unistd.h>
#include <libgen.h>
#if defined(MPI)
#include  <mpi.h>
#endif
#include "ncw.h"
#include "ncutils.h"
#include "version.h"
#include "distribute.h"
#include "utils.h"

#define PROGRAM_NAME "ncave"
#define PROGRAM_VERSION "0.01"

#define ALIGN __attribute__((aligned(32)))

#define NSRC_INC 10
#define NVAR_INC 10
#define NFIELD_INC 100
#define BIGNUM 1.0e+20
#define TMPVARNAME "ncave.tmpname"

int verbose = 0;
int force = 0;

int nprocesses = 1;
int rank = 0;

char tmpdirname[MAXSTRLEN] = ".";

typedef struct {
    int fid;
    char* varname;
    int ni, nj, nk;
    int k;
    size_t start[4];
    size_t count[4];
    size_t n;
} field;

/**
 */
static void printlog(const char* format, ...)
{
    va_list args;

    if (rank == 0) {
        va_start(args, format);
        (void) vprintf(format, args);
        va_end(args);
        fflush(stdout);
    }
}

/**
 */
static int dir_exists(char dirname[])
{
    DIR* d;

    d = opendir(dirname);
    if (d == NULL)
        return 0;
    closedir(d);
    return 1;
}

/**
 */
static int dir_createifabsent(char dirname[])
{
    int status;

    if (dir_exists(dirname))
        return 1;
    status = mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (status == 0)
        return 1;
    else {
        int errno_saved = errno;

        quit("could not create directory \"%s\": %s", dirname, strerror(errno_saved));
    }
    return 0;
}

/**
 */
static int rmentry(const char* entry, const struct stat* sb, int typeflag, struct FTW* ftwbuf)
{
    return remove(entry);
}

/**
 */
static void dir_rmallifexists(char dirname[])
{
    int status;

    if (!dir_exists(dirname))
        return;
    status = nftw(dirname, rmentry, 64, FTW_DEPTH | FTW_PHYS);
    if (status != 0) {
        int errno_saved = errno;

        quit("dir_rmallifexists(): \"%s\": %s", dirname, strerror(errno_saved));
    }
}

/**
 */
static void usage(int exitstatus)
{
    printf("  Usage: ncave [-v <var>] [...] [-c <var>] [...] [-V] [-f] {<src> [...] <dst>}\n");
    printf("         ncave -v\n");
    printf("  Parameters:\n");
    printf("    -v <var>            -- variable to be averaged over all input files\n");
    printf("                           (default: all variables with 2 or more dimensions)\n");
    printf("    -c <var>            -- variable to be copied from the first input file\n");
    printf("    {<src> [...] <dst>} -- list of input files followed by the output  file\n");
    printf("    -f                  -- overwrite destination if exists\n");
    printf("    -V                  -- verbose\n");
    printf("    -v                  -- print version and exit\n");
    exit(exitstatus);
}

/**
 */
static void parse_commandline(int argc, char* argv[], int* nsrc, char*** srcs, int* nvar, char*** vars, int* ncvar, char*** cvars)
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
            if (argv[i][1] == 'V') {
                verbose = 1;
                i++;
            } else if (argv[i][1] == 'f') {
                force = 1;
                i++;
            } else if (argv[i][1] == 'v') {
                i++;
                if (i >= argc)
                    quit("no variable specified after \"-v\"\n");
                if (*nvar % NVAR_INC == 0)
                    *vars = realloc(*vars, (*nvar + NVAR_INC) * sizeof(void*));
                (*vars)[*nvar] = strdup(argv[i]);
                (*nvar)++;
                i++;
            } else if (argv[i][1] == 'c') {
                i++;
                if (i >= argc)
                    quit("no variable specified after \"-c\"\n");
                if (*ncvar % NVAR_INC == 0)
                    *cvars = realloc(*cvars, (*ncvar + NVAR_INC) * sizeof(void*));
                (*cvars)[*ncvar] = strdup(argv[i]);
                (*ncvar)++;
                i++;
            } else {
                printf("  ncave: ERROR: unknown option \"%s\"\n", argv[i]);
                usage(1);
            }
        } else {
            if (*nsrc != 0) {
                printf("  ncave: ERROR: input and output files need to be specified in a continuous sequence\n");
                usage(1);
            }
            while (i < argc && argv[i][0] != '-') {
                if (*nsrc % NSRC_INC == 0)
                    *srcs = realloc(*srcs, (*nsrc + NSRC_INC) * sizeof(void*));
                (*srcs)[*nsrc] = strdup(argv[i]);
                (*nsrc)++;
                i++;
            }
        }
    }
    if (*nsrc == 0)
        quit("no input specified");
    if (*nsrc == 1)
        quit("no output specified");
    if (*nvar > 0 && *ncvar > 0) {
        for (i = 0; i < *nvar; ++i) {
            int j;

            for (j = 0; j < *ncvar; ++j)
                if (strcmp((*vars)[i], (*cvars)[j]) == 0)
                    quit("variable \"%s\" is specified to be both averaged and copied", (*vars)[i]);
        }
    }
}

/**
 */
static void getvars(char* fname, int* nvar, char*** vars)
{
    int ncid;
    int nvartotal, vid;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_nvars(ncid, &nvartotal);
    assert(nvartotal > 0);
    *vars = malloc(nvartotal * sizeof(char*));
    for (vid = 0; vid < nvartotal; ++vid) {
        int ndims, nd;
        size_t dimlen[4];
        char varname[NC_MAX_NAME];
        int i, i1;

        ncw_inq_varname(ncid, vid, varname);
        ncw_inq_vardims(ncid, vid, 4, &ndims, dimlen);
        i1 = 0;
        if (ndims > 4)
            quit("%s: %s: do not know how to treat a %d-dimensional variable", fname, varname, ndims);
        if (ncw_var_hasunlimdim(ncid, vid)) {
            if (dimlen[0] != 1)
                quit("%s: %s: unlimited dimension length is allowed to be 1 only", fname, varname);
            i1 = 1;
        }
        for (i = i1; i < ndims; ++i)
            if (dimlen[i] > 1)
                break;
        nd = ndims - i;

        if (nd < 2)
            continue;

        (*vars)[*nvar] = strdup(varname);
        (*nvar)++;
    }
    ncw_close(ncid);
}

/**
 */
static void getfields(char* fname, int nvar, char** vars, int* nfield, field** fields)
{
    int ncid;
    int v;

    ncw_open(fname, NC_NOWRITE, &ncid);
    for (v = 0; v < nvar; ++v) {
        int varid;
        int ndim, nd;
        size_t dimlen[4];
        int i, i1;

        ncw_inq_varid(ncid, vars[v], &varid);
        ncw_inq_vardims(ncid, varid, 4, &ndim, dimlen);
        i1 = 0;
        if (ncw_var_hasunlimdim(ncid, varid))
            i1 = 1;
        for (i = i1; i < ndim; ++i)
            if (dimlen[i] > 1)
                break;
        nd = ndim - i;

        if (nd <= 2) {
            field* f;

            if (*nfield % NFIELD_INC == 0)
                *fields = realloc(*fields, (*nfield + NFIELD_INC) * sizeof(field));
            f = &(*fields)[*nfield];
            f->fid = *nfield;
            f->varname = vars[v];
            if (nd < 2) {
                f->nj = -1;
                f->ni = -1;
                f->nk = -1;
                f->k = -1;
                ncw_inq_varsize(ncid, varid, &f->n);
            } else if (nd == 2) {
                f->nj = dimlen[i];
                f->ni = dimlen[i + 1];
                f->nk = 0;
                f->k = 0;
                f->n = f->ni * f->nj;
            }
            for (i = 0; i < ndim; ++i)
                f->count[i] = dimlen[i];
            (*nfield)++;
        } else {
            field* f;
            int k;

            for (k = 0; k < dimlen[i]; ++k) {
                if (*nfield % NFIELD_INC == 0)
                    *fields = realloc(*fields, (*nfield + NFIELD_INC) * sizeof(field));
                f = &(*fields)[*nfield];
                f->fid = *nfield;
                f->varname = vars[v];
                f->nk = dimlen[i];
                f->nj = dimlen[i + 1];
                f->ni = dimlen[i + 2];
                f->k = k;
                f->n = f->ni * f->nj;
                if (i == 1)
                    f->count[0] = 1;
                f->count[i] = dimlen[i];
                f->start[i] = k;
                f->count[i] = 1;
                f->count[i + 1] = f->nj;
                f->count[i + 2] = f->ni;
                (*nfield)++;
            }
        }
    }
    ncw_close(ncid);
}

/**
 */
static void gettilename(field* f, char* dst, char** tilename)
{
    char path[MAXSTRLEN];
    char* bname;
    size_t len;

    strncpy(path, dst, MAXSTRLEN - 1);
    path[MAXSTRLEN - 1] = 0;
    bname = basename(path);
    len = strlen(tmpdirname) + strlen(bname) + strlen(f->varname) + 10 + 1;
    *tilename = malloc(len);
    snprintf(*tilename, len, "%s/%s-%s-%03d.tmp", tmpdirname, bname, f->varname, (f->k >= 0) ? f->k : 0);
}

/**
 */
int main(int argc, char* argv[])
{
    int nsrc = 0;
    char** srcs = NULL;
    char* dst;

    /*
     * variables to average
     */
    int nvar = 0;
    char** vars = NULL;

    /*
     * variables to copy
     */
    int ncvar = 0;
    char** cvars = NULL;
    int nfield = 0;
    field* fields = NULL;
    int i;

    parse_commandline(argc, argv, &nsrc, &srcs, &nvar, &vars, &ncvar, &cvars);

    dst = srcs[nsrc - 1];
    nsrc--;
    if (file_exists(dst) && !force)
        quit("destination \"%s\" exists", dst);

#if defined(MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    if (verbose) {
        printlog("  ncave v%s\n", VERSION);
        printlog("  MPI: initialised %d process(es)\n", nprocesses);
    }
#if defined(DEBUG)
    if (verbose)
        printlog("  master PID = %3d\n", getpid());
#endif

    if (verbose) {
        printlog("  input = %s\n", srcs[0]);
        for (i = 1; i < nsrc; ++i)
            printlog("          %s\n", srcs[i]);
        printlog("  output = %s\n", srcs[nsrc]);
    }

    if (nvar == 0) {
        if (ncvar == 0) {
            getvars(srcs[0], &nvar, &vars);
            if (verbose)
                printlog("  found %d variable(s) to average:", nvar);
        }
    } else if (verbose)
        printlog("  averaging %d variable(s):", nvar);
    if (nvar == 0 && ncvar == 0)
        quit("no variables to average");
    if (verbose) {
        for (i = 0; i < nvar; ++i)
            printlog(" %s", vars[i]);
        printlog("\n");
    }
    if (ncvar > 0 && verbose) {
        printlog("  copying %d variable(s):", ncvar);
        for (i = 0; i < ncvar; ++i)
            printlog(" %s", cvars[i]);
        printlog("\n");
    }

    getfields(srcs[0], nvar, vars, &nfield, &fields);
    if (verbose)
        printlog("  %d field(s)\n", nfield);

    /*
     * create temporary directory for tiles
     */
    {
        char path[MAXSTRLEN];
        char* bname;

        strncpy(path, dst, MAXSTRLEN - 1);
        path[MAXSTRLEN - 1] = 0;
        bname = basename(path);
        if (snprintf(tmpdirname, MAXSTRLEN, ".%s.ncave.tmp", bname) < 0)
            quit("temporary directory name too long (> 2048)");
        if (rank == 0)
            dir_createifabsent(tmpdirname);
#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    /*
     * calculate average fields and write them to tiles
     */
    if (nfield > 0) {
        distribute_iterations(0, nfield - 1, nprocesses, rank);
        if (verbose)
            printlog("  writing tiles:");
        for (i = my_first_iteration; i <= my_last_iteration; ++i) {
            field* f = &fields[i];
            float* vin = malloc(f->n * sizeof(float));
            float* vout = calloc(f->n, sizeof(float));
            int j, ij;

            for (j = 0; j < nsrc; ++j) {
                ncu_readfield(srcs[j], f->varname, f->k, f->ni, f->nj, f->nk, vin);
                for (ij = 0; ij < f->n; ++ij)
                    vout[ij] += vin[ij];
            }
            for (ij = 0; ij < f->n; ++ij)
                vout[ij] /= (float) nsrc;

            {
                char* tilename = NULL;
                int ncid, dimid, varid;

                gettilename(f, dst, &tilename);
                ncw_create(tilename, NC_NOWRITE, &ncid);
                ncw_def_dim(ncid, "n", f->n, &dimid);
                ncw_def_var(ncid, f->varname, NC_FLOAT, 1, &dimid, &varid);
                ncw_enddef(ncid);
                ncw_put_var_float(ncid, varid, vout);
                ncw_close(ncid);
                if (verbose > 0) {
                    printf(".");
                    fflush(stdout);
                }
                free(tilename);
            }
            free(vin);
            free(vout);
        }
#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (verbose)
            printlog("\n");
    }

    /*
     * assemble tiles
     */
    if (verbose)
        printlog("  assembling:");
    if (rank == 0) {
        char dst_tmp[MAXSTRLEN];
        int ncid_src, ncid_dst;

        /*
         * set the temporary destination
         */
        strncpy(dst_tmp, dst, MAXSTRLEN - 1);
        dst_tmp[MAXSTRLEN - 1] = 0;
        if (strlen(dst_tmp) < MAXSTRLEN - 4)
            strcat(dst_tmp, ".tmp");
        else
            quit("destination too long");

        ncw_open(srcs[0], NC_NOWRITE, &ncid_src);
        ncw_create(dst_tmp, NC_CLOBBER | NC_NETCDF4, &ncid_dst);

        for (i = 0; i < nvar; ++i) {
            int varid_src;

            ncw_inq_varid(ncid_src, vars[i], &varid_src);
            ncw_copy_vardef(ncid_src, varid_src, ncid_dst);
        }
        for (i = 0; i < ncvar; ++i) {
            int varid_src;

            ncw_inq_varid(ncid_src, cvars[i], &varid_src);
            ncw_copy_vardef(ncid_src, varid_src, ncid_dst);
        }
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
        ncw_enddef(ncid_dst);

        for (i = 0; i < ncvar; ++i) {
            int varid_src;

            ncw_inq_varid(ncid_src, cvars[i], &varid_src);
            ncw_copy_vardata(ncid_src, varid_src, ncid_dst);
        }

        ncw_close(ncid_src);
        ncw_close(ncid_dst);

        for (i = 0; i < nfield; ++i) {
            field* f = &fields[i];
            float* v = malloc(f->n * sizeof(float));
            char* tilename = NULL;
            int ncid_tile;

            gettilename(f, dst, &tilename);
            ncw_open(tilename, NC_NOWRITE, &ncid_tile);
            ncw_get_var_float(ncid_tile, 0, v);
            ncw_close(ncid_tile);

            ncu_writefield(dst_tmp, f->varname, f->k, f->ni, f->nj, f->nk, v);
            if (verbose)
                printlog(".");

            free(v);
            free(tilename);
        }

        file_rename(dst_tmp, dst);
        dir_rmallifexists(tmpdirname);

        if (verbose)
            printlog("\n");
    }

    /*
     * clean up
     */
    if (nfield > 0)
        free(fields);
    nsrc++;
    for (i = 0; i < nsrc; ++i)
        free(srcs[i]);
    free(srcs);
    if (nvar > 0) {
        for (i = 0; i < nvar; ++i)
            free(vars[i]);
        free(vars);
    }
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (verbose)
        printlog("  finished\n");
#if defined(MPI)
    MPI_Finalize();
#endif

    return 0;
}
