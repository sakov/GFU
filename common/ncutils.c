/******************************************************************************
 *
 * File:        ncutils.c        
 *
 * Created:     6/9/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Intermediate level NetCDF r/w procedures. Among other things
 *              these procedures are supposed to handle the following
 *              attributes: _FillValue, missing_value, valid_range, valid_min,
 *              valid_max, add_offset, and scale_factor.
 *
 *              This file is an excerpt from ncutils.c in EnKF-C package.
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include "ncw.h"
#include "ncutils.h"

static void ncu_quit_fn_def(char* format, ...);
static ncu_quit_fn quit = ncu_quit_fn_def;

/**
 */
static void ncu_quit_fn_def(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "\n\n  error: ncu: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n\n");
    exit(1);
}

/**
 */
void ncu_set_quitfn(ncu_quit_fn quitfn)
{
    quit = quitfn;
}

/**
 */
void ncu_readvarfloat(int ncid, int varid, size_t n, float v[])
{
    char attnames[2][20] = { "_FillValue", "missing_value" };
    nc_type vartype = -1;
    int typesize = 0;
    void* vv = NULL;
    size_t i, s;

    ncw_check_varsize(ncid, varid, n);
    ncw_get_var_float(ncid, varid, v);

    ncw_inq_vartype(ncid, varid, &vartype);
    typesize = ncw_sizeof(vartype);

    for (s = 0; s < 2; ++s) {
        char* attname = attnames[s];

        if (ncw_att_exists2(ncid, varid, attname)) {
            ncw_check_attlen(ncid, varid, attname, 1);

            if (vv == NULL) {
                if (vartype != NC_FLOAT) {
                    vv = malloc(n * typesize);
                    ncw_get_var(ncid, varid, vv);
                } else
                    vv = v;
            }

            if (typesize == 1) {
                int8_t value;

                ncw_get_att(ncid, varid, attname, &value);
                for (i = 0; i < n; ++i)
                    if (((int8_t*) vv)[i] == value)
                        v[i] = NAN;
            } else if (typesize == 2) {
                int16_t value;

                ncw_get_att(ncid, varid, attname, &value);
                for (i = 0; i < n; ++i)
                    if (((int16_t*) vv)[i] == value)
                        v[i] = NAN;
            } else if (typesize == 4) {
                int32_t value;

                ncw_get_att(ncid, varid, attname, &value);
                for (i = 0; i < n; ++i)
                    if (((int32_t*) vv)[i] == value)
                        v[i] = NAN;
            } else if (typesize == 8) {
                int64_t value;

                ncw_get_att(ncid, varid, attname, &value);
                for (i = 0; i < n; ++i)
                    if (((int64_t*) vv)[i] == value)
                        v[i] = NAN;
            } else
                quit("programming error");
        }
    }
    if (!ncw_att_exists2(ncid, varid, "_FillValue")) {
        if (vv == NULL) {
            if (typesize != sizeof(float)) {
                vv = malloc(n * typesize);
                ncw_get_var(ncid, varid, vv);
            } else
                vv = v;
        }
        if (typesize == 1) {
            int nofill;
            int8_t value;

            ncw_inq_var_fill(ncid, varid, &nofill, &value);
            if (nofill == 0)
                for (i = 0; i < n; ++i)
                    if (((int8_t*) vv)[i] == value)
                        v[i] = NAN;
        } else if (typesize == 2) {
            int nofill;
            int16_t value;

            ncw_inq_var_fill(ncid, varid, &nofill, &value);
            if (nofill == 0)
                for (i = 0; i < n; ++i)
                    if (((int16_t*) vv)[i] == value)
                        v[i] = NAN;
        } else if (typesize == 4) {
            int nofill;
            int32_t value;

            ncw_inq_var_fill(ncid, varid, &nofill, &value);
            if (nofill == 0)
                for (i = 0; i < n; ++i)
                    if (((int32_t*) vv)[i] == value)
                        v[i] = NAN;
        } else if (typesize == 8) {
            int nofill;
            int64_t value;

            ncw_inq_var_fill(ncid, varid, &nofill, &value);
            if (nofill == 0)
                for (i = 0; i < n; ++i)
                    if (((int64_t*) vv)[i] == value)
                        v[i] = NAN;
        } else
            quit("programming error");
    }
    if (ncw_att_exists2(ncid, varid, "valid_min")) {
        ncw_check_attlen(ncid, varid, "valid_min", 1);

        if (vv == NULL) {
            if (typesize != sizeof(float)) {
                vv = malloc(n * typesize);
                ncw_get_var(ncid, varid, vv);
            } else
                vv = v;
        }

        if (vartype == NC_BYTE || vartype == NC_CHAR) {
            char value;

            ncw_get_att(ncid, varid, "valid_min", &value);
            for (i = 0; i < n; ++i)
                if (((char*) vv)[i] < value)
                    v[i] = NAN;
        } else if (vartype == NC_UBYTE) {
            char value;

            ncw_get_att(ncid, varid, "valid_min", &value);
            for (i = 0; i < n; ++i)
                if (((unsigned char*) vv)[i] < value)
                    v[i] = NAN;
        } else if (vartype == NC_SHORT) {
            int16_t value;

            ncw_get_att(ncid, varid, "valid_min", &value);
            for (i = 0; i < n; ++i)
                if (((int16_t*) vv)[i] < value)
                    v[i] = NAN;
        } else if (vartype == NC_USHORT) {
            uint16_t value;

            ncw_get_att(ncid, varid, "valid_min", &value);
            for (i = 0; i < n; ++i)
                if (((uint16_t *) vv)[i] < value)
                    v[i] = NAN;
        } else if (vartype == NC_INT || vartype == NC_LONG) {
            int32_t value;

            ncw_get_att(ncid, varid, "valid_min", &value);
            for (i = 0; i < n; ++i)
                if (((int32_t*) vv)[i] < value)
                    v[i] = NAN;
        } else if (vartype == NC_UINT) {
            uint32_t value;

            ncw_get_att(ncid, varid, "valid_min", &value);
            for (i = 0; i < n; ++i)
                if (((uint32_t *) vv)[i] < value)
                    v[i] = NAN;
        } else if (vartype == NC_INT64) {
            int64_t value;

            ncw_get_att(ncid, varid, "valid_min", &value);
            for (i = 0; i < n; ++i)
                if (((int64_t*) vv)[i] < value)
                    v[i] = NAN;
        } else if (vartype == NC_UINT64) {
            uint64_t value;

            ncw_get_att(ncid, varid, "valid_min", &value);
            for (i = 0; i < n; ++i)
                if (((uint64_t *) vv)[i] < value)
                    v[i] = NAN;
        } else if (vartype == NC_FLOAT) {
            float value;

            assert(sizeof(float) == 4);
            ncw_get_att(ncid, varid, "valid_min", &value);
            for (i = 0; i < n; ++i)
                if (((float*) vv)[i] < value)
                    v[i] = NAN;
        } else if (vartype == NC_DOUBLE) {
            double value;

            assert(sizeof(double) == 8);
            ncw_get_att(ncid, varid, "valid_min", &value);
            for (i = 0; i < n; ++i)
                if (((double*) vv)[i] < value)
                    v[i] = NAN;
        } else
            quit("programming error");
    }
    if (ncw_att_exists2(ncid, varid, "valid_max")) {
        ncw_check_attlen(ncid, varid, "valid_max", 1);

        if (vv == NULL) {
            if (typesize != sizeof(float)) {
                vv = malloc(n * typesize);
                ncw_get_var(ncid, varid, vv);
            } else
                vv = v;
        }

        if (vartype == NC_BYTE || vartype == NC_CHAR) {
            char value;

            ncw_get_att(ncid, varid, "valid_max", &value);
            for (i = 0; i < n; ++i)
                if (((char*) vv)[i] > value)
                    v[i] = NAN;
        } else if (vartype == NC_UBYTE) {
            char value;

            ncw_get_att(ncid, varid, "valid_max", &value);
            for (i = 0; i < n; ++i)
                if (((unsigned char*) vv)[i] > value)
                    v[i] = NAN;
        } else if (vartype == NC_SHORT) {
            int16_t value;

            ncw_get_att(ncid, varid, "valid_max", &value);
            for (i = 0; i < n; ++i)
                if (((int16_t*) vv)[i] > value)
                    v[i] = NAN;
        } else if (vartype == NC_USHORT) {
            uint16_t value;

            ncw_get_att(ncid, varid, "valid_max", &value);
            for (i = 0; i < n; ++i)
                if (((uint16_t *) vv)[i] > value)
                    v[i] = NAN;
        } else if (vartype == NC_INT || vartype == NC_LONG) {
            int32_t value;

            ncw_get_att(ncid, varid, "valid_max", &value);
            for (i = 0; i < n; ++i)
                if (((int32_t*) vv)[i] > value)
                    v[i] = NAN;
        } else if (vartype == NC_UINT) {
            uint32_t value;

            ncw_get_att(ncid, varid, "valid_max", &value);
            for (i = 0; i < n; ++i)
                if (((uint32_t *) vv)[i] > value)
                    v[i] = NAN;
        } else if (vartype == NC_INT64) {
            int64_t value;

            ncw_get_att(ncid, varid, "valid_max", &value);
            for (i = 0; i < n; ++i)
                if (((int64_t*) vv)[i] > value)
                    v[i] = NAN;
        } else if (vartype == NC_UINT64) {
            uint64_t value;

            ncw_get_att(ncid, varid, "valid_max", &value);
            for (i = 0; i < n; ++i)
                if (((uint64_t *) vv)[i] > value)
                    v[i] = NAN;
        } else if (vartype == NC_FLOAT) {
            float value;

            assert(sizeof(float) == 4);
            ncw_get_att(ncid, varid, "valid_max", &value);
            for (i = 0; i < n; ++i)
                if (((float*) vv)[i] > value)
                    v[i] = NAN;
        } else if (vartype == NC_DOUBLE) {
            double value;

            assert(sizeof(double) == 8);
            ncw_get_att(ncid, varid, "valid_max", &value);
            for (i = 0; i < n; ++i)
                if (((double*) vv)[i] > value)
                    v[i] = NAN;
        } else
            quit("programming error");
    }
    if (ncw_att_exists2(ncid, varid, "valid_range")) {
        ncw_check_attlen(ncid, varid, "valid_range", 2);

        if (vv == NULL) {
            if (typesize != sizeof(float)) {
                vv = malloc(n * typesize);
                ncw_get_var(ncid, varid, vv);
            } else
                vv = v;
        }

        if (vartype == NC_BYTE || vartype == NC_CHAR) {
            char value[2];

            ncw_get_att(ncid, varid, "valid_range", value);
            for (i = 0; i < n; ++i)
                if (((char*) vv)[i] < value[0] || ((char*) vv)[i] > value[1])
                    v[i] = NAN;
        } else if (vartype == NC_UBYTE) {
            char value[2];

            ncw_get_att(ncid, varid, "valid_range", value);
            for (i = 0; i < n; ++i)
                if (((unsigned char*) vv)[i] < value[0]
                    || ((unsigned char*) vv)[i] > value[1])
                    v[i] = NAN;
        } else if (vartype == NC_SHORT) {
            int16_t value[2];

            ncw_get_att(ncid, varid, "valid_range", value);
            for (i = 0; i < n; ++i)
                if (((int16_t*) vv)[i] < value[0]
                    || ((int16_t*) vv)[i] > value[1])
                    v[i] = NAN;
        } else if (vartype == NC_USHORT) {
            uint16_t value[2];

            ncw_get_att(ncid, varid, "valid_range", value);
            for (i = 0; i < n; ++i)
                if (((uint16_t *) vv)[i] < value[0]
                    || ((uint16_t *) vv)[i] > value[1])
                    v[i] = NAN;
        } else if (vartype == NC_INT || vartype == NC_LONG) {
            int32_t value[2];

            ncw_get_att(ncid, varid, "valid_range", value);
            for (i = 0; i < n; ++i)
                if (((int32_t*) vv)[i] < value[0]
                    || ((int32_t*) vv)[i] > value[1])
                    v[i] = NAN;
        } else if (vartype == NC_UINT) {
            uint32_t value[2];

            ncw_get_att(ncid, varid, "valid_range", value);
            for (i = 0; i < n; ++i)
                if (((uint32_t *) vv)[i] < value[0]
                    || ((uint32_t *) vv)[i] > value[1])
                    v[i] = NAN;
        } else if (vartype == NC_INT64) {
            int64_t value[2];

            ncw_get_att(ncid, varid, "valid_range", value);
            for (i = 0; i < n; ++i)
                if (((int64_t*) vv)[i] < value[0]
                    || ((int64_t*) vv)[i] > value[1])
                    v[i] = NAN;
        } else if (vartype == NC_UINT64) {
            uint64_t value[2];

            ncw_get_att(ncid, varid, "valid_range", value);
            for (i = 0; i < n; ++i)
                if (((uint64_t *) vv)[i] < value[0]
                    || ((uint64_t *) vv)[i] > value[1])
                    v[i] = NAN;
        } else if (vartype == NC_FLOAT) {
            float value[2];

            assert(sizeof(float) == 4);
            ncw_get_att(ncid, varid, "valid_range", value);
            for (i = 0; i < n; ++i)
                if (((float*) vv)[i] < value[0] || ((float*) vv)[i] > value[1])
                    v[i] = NAN;
        } else if (vartype == NC_DOUBLE) {
            double value[2];

            assert(sizeof(double) == 8);
            ncw_get_att(ncid, varid, "valid_range", value);
            for (i = 0; i < n; ++i)
                if (((double*) vv)[i] < value[0]
                    || ((double*) vv)[i] > value[1])
                    v[i] = NAN;
        } else
            quit("programming error");
    }
    if (vv != NULL && vartype != NC_FLOAT)
        free(vv);

    if (ncw_att_exists(ncid, varid, "scale_factor")) {
        float scale_factor;

        ncw_check_attlen(ncid, varid, "scale_factor", 1);
        ncw_get_att_float(ncid, varid, "scale_factor", &scale_factor);
        for (i = 0; i < n; ++i)
            v[i] *= scale_factor;
    }
    if (ncw_att_exists(ncid, varid, "add_offset")) {
        float add_offset;

        ncw_check_attlen(ncid, varid, "scale_factor", 1);
        ncw_get_att_float(ncid, varid, "add_offset", &add_offset);
        for (i = 0; i < n; ++i)
            v[i] += add_offset;
    }
}

/** Reads one horizontal field (layer) for a variable from a NetCDF file.
 ** Verifies that the field dimensions are ni x nj.
 */
void ncu_readfield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v)
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    size_t start[4], count[4];
    size_t i, n;
    int hasrecorddim;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);
    if (hasrecorddim)
        if (dimlen[0] == 0)
            quit("ncu_readfield(): \"%s\": %s: empty record dimension", fname, varname);

    if (ndims == 4) {
        if (nj == 0) {
            quit("ncu_readfield(): \"%s\": %s: expected positive \"j\" dimension for a 4-dimensional variable\n", fname, varname);
        }
        if (!hasrecorddim && dimlen[0] != 1)
            quit("ncu_readfield(): \"%s\": %s: for a 4-dimensional variable expected the first dimension to be either unlimited or of length 1\n", fname, varname);
        start[0] = dimlen[0] - 1;
        if (dimlen[1] != nk) {
            if (dimlen[1] != 1)
                quit("ncu_readfield(): \"%s\": vertical dimension of variable \"%s\" (nk = %d) does not match grid dimension (nk = %d)", fname, varname, dimlen[1], nk);
            else
                k = 0;          /* ignore k */
        }
        if (dimlen[3] != ni || dimlen[2] != nj)
            quit("ncu_readfield(): \"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[3], dimlen[2], ni, nj);
        start[1] = k;
        start[2] = 0;
        start[3] = 0;
        count[0] = 1;
        count[1] = 1;
        count[2] = dimlen[2];
        count[3] = dimlen[3];
    } else if (ndims == 3) {
        if (nj > 0) {
            if (!hasrecorddim) {
                if (dimlen[0] != nk && !(dimlen[0] == 1 && (k == 0 || k == nk - 1)))
                    quit("ncu_readfield(): \"%s\": vertical dimension of variable \"%s\" (nk = %d) does not match grid dimension (nk = %d)", fname, varname, dimlen[0], nk);
                start[0] = (dimlen[0] == 1) ? 0 : k;
                start[1] = 0;
                start[2] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
                count[2] = dimlen[2];
            } else {
                /*
                 * ignore k in this case
                 */
                start[0] = dimlen[0] - 1;
                start[1] = 0;
                start[2] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
                count[2] = dimlen[2];
            }
            if (dimlen[2] != ni || dimlen[1] != nj)
                quit("ncu_readfield(): \"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[2], dimlen[1], ni, nj);
        } else {
            if (!hasrecorddim && dimlen[0] != 1)
                quit("ncu_readfield(): \"%s\": %s: for a 3-dimensional variable on unstructured horizontal grid expected the first dimension to be either unlimited or of length 1\n", fname, varname);
            start[0] = dimlen[0] - 1;
            if (dimlen[1] != nk) {
                if (dimlen[1] != 1)
                    quit("ncu_readfield(): \"%s\": vertical dimension of variable \"%s\" (nk = %d) does not match grid dimension (nk = %d)", fname, varname, dimlen[1], nk);
                else
                    k = 0;      /* ignore k */
            }
            if (dimlen[2] != ni)
                quit("ncu_readfield(): \"%s\": horizontal dimension of variable \"%s\" (ni = %d) does not match grid dimension (ni = %d)", fname, varname, dimlen[2], ni);
            start[1] = k;
            start[2] = 0;
            count[0] = 1;
            count[1] = 1;
            count[2] = dimlen[2];
        }
    } else if (ndims == 2) {
        if (nj > 0) {
            if (hasrecorddim)
                quit("ncu_readfield(): %s: can not read layer from a 1D variable \"%s\"", fname, varname);
            /*
             * ignore k
             */
            start[0] = 0;
            start[1] = 0;
            count[0] = dimlen[0];
            count[1] = dimlen[1];
            if (dimlen[1] != ni || dimlen[0] != nj)
                quit("ncu_readfield(): \"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[1], dimlen[0], ni, nj);
        } else {
            if (!hasrecorddim) {
                if (dimlen[0] != nk && !(dimlen[0] == 1 && (k == 0 || k == nk - 1)))
                    quit("ncu_readfield(): \"%s\": vertical dimension of variable \"%s\" (nk = %d) does not match grid dimension (nk = %d)", fname, varname, dimlen[0], nk);
                start[0] = k;
                start[1] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
            } else {
                /*
                 * ignore k in this case
                 */
                start[0] = dimlen[0] - 1;
                start[1] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
            }
            if (dimlen[1] != ni)
                quit("ncu_readfield(): \"%s\": horizontal dimension of variable \"%s\" (ni = %d) does not match grid dimension (ni = %d)", fname, varname, dimlen[1], ni);
        }
    } else if (ndims == 1) {
        if (nj > 0) {
            quit("ncu_readfield(): %s: can not read 2D field for \"%s\": # of dimensions = %d", fname, varname, ndims);
        } else {
            if (hasrecorddim)
                quit("ncu_readfield(): %s: can not read layer from a 0D variable \"%s\"", fname, varname);
            /*
             * ignore k in this case
             */
            start[0] = 0;
            count[0] = dimlen[0];
        }
    } else
        quit("ncu_readfield(): %s: can not read 2D field for \"%s\": # of dimensions = %d", fname, varname, ndims);

    ncw_get_vara_float(ncid, varid, start, count, v);

    n = 1;
    for (i = 0; i < ndims; ++i)
        n *= count[i];

    {
        void* vv = NULL;
        nc_type vartype = -1;
        int typesize = 0;
        char buf[128];
        void* attval = buf;

        ncw_inq_vartype(ncid, varid, &vartype);
        typesize = ncw_sizeof(vartype);
        if (vartype != NC_FLOAT) {
            vv = malloc(n * typesize);
            ncw_get_vara(ncid, varid, start, count, vv);
        } else
            vv = v;

        if (ncw_att_exists2(ncid, varid, "_FillValue")) {
            ncw_check_attlen(ncid, varid, "_FillValue", 1);
            ncw_get_att(ncid, varid, "_FillValue", attval);
        } else {
            int nofill;

            ncw_inq_var_fill(ncid, varid, &nofill, attval);
            if (nofill)
                goto skip;
        }
        if (typesize == 1) {
            for (i = 0; i < n; ++i)
                if (((int8_t*) vv)[i] == ((int8_t*) attval)[0])
                    v[i] = NAN;
        } else if (typesize == 2) {
            for (i = 0; i < n; ++i)
                if (((int16_t*) vv)[i] == ((int16_t*) attval)[0])
                    v[i] = NAN;
        } else if (typesize == 4) {
            for (i = 0; i < n; ++i)
                if (((int32_t*) vv)[i] == ((int32_t*) attval)[0])
                    v[i] = NAN;
        } else if (typesize == 8) {
            for (i = 0; i < n; ++i)
                if (((int64_t*) vv)[i] == ((int64_t*) attval)[0])
                    v[i] = NAN;
        } else
            quit("programming error");

      skip:
        if (ncw_att_exists2(ncid, varid, "missing_value")) {
            ncw_check_attlen(ncid, varid, "missing_value", 1);
            ncw_get_att(ncid, varid, "missing_value", attval);
            if (typesize == 1) {
                for (i = 0; i < n; ++i)
                    if (((int8_t*) vv)[i] == ((int8_t*) attval)[0])
                        v[i] = NAN;
            } else if (typesize == 2) {
                for (i = 0; i < n; ++i)
                    if (((int16_t*) vv)[i] == ((int16_t*) attval)[0])
                        v[i] = NAN;
            } else if (typesize == 4) {
                for (i = 0; i < n; ++i)
                    if (((int32_t*) vv)[i] == ((int32_t*) attval)[0])
                        v[i] = NAN;
            } else if (typesize == 8) {
                for (i = 0; i < n; ++i)
                    if (((int64_t*) vv)[i] == ((int64_t*) attval)[0])
                        v[i] = NAN;
            } else
                quit("programming error");
        }

        if (ncw_att_exists2(ncid, varid, "valid_min")) {
            ncw_check_attlen(ncid, varid, "valid_min", 1);
            ncw_get_att(ncid, varid, "valid_min", attval);
            if (vartype == NC_BYTE || vartype == NC_CHAR) {
                for (i = 0; i < n; ++i)
                    if (((char*) vv)[i] < ((char*) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_UBYTE) {
                for (i = 0; i < n; ++i)
                    if (((unsigned char*) vv)[i] < ((unsigned char*) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_SHORT) {
                for (i = 0; i < n; ++i)
                    if (((int16_t*) vv)[i] < ((int16_t*) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_USHORT) {
                for (i = 0; i < n; ++i)
                    if (((uint16_t *) vv)[i] < ((uint16_t *) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_INT || vartype == NC_LONG) {
                for (i = 0; i < n; ++i)
                    if (((int32_t*) vv)[i] < ((int32_t*) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_UINT) {
                for (i = 0; i < n; ++i)
                    if (((uint32_t *) vv)[i] < ((uint32_t *) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_INT64) {
                for (i = 0; i < n; ++i)
                    if (((int64_t*) vv)[i] < ((int64_t*) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_UINT64) {
                for (i = 0; i < n; ++i)
                    if (((uint64_t *) vv)[i] < ((uint64_t *) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_FLOAT) {
                assert(sizeof(float) == 4);
                for (i = 0; i < n; ++i)
                    if (((float*) vv)[i] < ((float*) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_DOUBLE) {
                for (i = 0; i < n; ++i)
                    if (((double*) vv)[i] < ((double*) attval)[0])
                        v[i] = NAN;
            } else
                quit("programming error");
        }

        if (ncw_att_exists2(ncid, varid, "valid_max")) {
            ncw_check_attlen(ncid, varid, "valid_max", 1);
            ncw_get_att(ncid, varid, "valid_max", attval);
            if (vartype == NC_BYTE || vartype == NC_CHAR) {
                for (i = 0; i < n; ++i)
                    if (((char*) vv)[i] > ((char*) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_UBYTE) {
                for (i = 0; i < n; ++i)
                    if (((unsigned char*) vv)[i] > ((unsigned char*) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_SHORT) {
                for (i = 0; i < n; ++i)
                    if (((int16_t*) vv)[i] > ((int16_t*) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_USHORT) {
                for (i = 0; i < n; ++i)
                    if (((uint16_t *) vv)[i] > ((uint16_t *) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_INT || vartype == NC_LONG) {
                for (i = 0; i < n; ++i)
                    if (((int32_t*) vv)[i] > ((int32_t*) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_UINT) {
                for (i = 0; i < n; ++i)
                    if (((uint32_t *) vv)[i] > ((uint32_t *) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_INT64) {
                for (i = 0; i < n; ++i)
                    if (((int64_t*) vv)[i] > ((int64_t*) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_UINT64) {
                for (i = 0; i < n; ++i)
                    if (((uint64_t *) vv)[i] > ((uint64_t *) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_FLOAT) {
                assert(sizeof(float) == 4);
                for (i = 0; i < n; ++i)
                    if (((float*) vv)[i] > ((float*) attval)[0])
                        v[i] = NAN;
            } else if (vartype == NC_DOUBLE) {
                for (i = 0; i < n; ++i)
                    if (((double*) vv)[i] > ((double*) attval)[0])
                        v[i] = NAN;
            } else
                quit("programming error");
        }

        if (ncw_att_exists2(ncid, varid, "valid_range")) {
            ncw_check_attlen(ncid, varid, "valid_range", 2);
            ncw_get_att(ncid, varid, "valid_range", attval);
            if (vartype == NC_BYTE || vartype == NC_CHAR) {
                for (i = 0; i < n; ++i)
                    if (((char*) vv)[i] < ((char*) attval)[0]
                        || ((char*) vv)[i] > ((char*) attval)[1])
                        v[i] = NAN;
            } else if (vartype == NC_UBYTE) {
                for (i = 0; i < n; ++i)
                    if (((unsigned char*) vv)[i] < ((unsigned char*) attval)[0]
                        || ((unsigned char*) vv)[i] > ((unsigned char*) attval)[1])
                        v[i] = NAN;
            } else if (vartype == NC_SHORT) {
                for (i = 0; i < n; ++i)
                    if (((int16_t*) vv)[i] < ((int16_t*) attval)[0]
                        || ((int16_t*) vv)[i] > ((int16_t*) attval)[1])
                        v[i] = NAN;
            } else if (vartype == NC_USHORT) {
                for (i = 0; i < n; ++i)
                    if (((uint16_t *) vv)[i] < ((uint16_t *) attval)[0]
                        || ((uint16_t *) vv)[i] > ((uint16_t *) attval)[1])
                        v[i] = NAN;
            } else if (vartype == NC_INT || vartype == NC_LONG) {
                for (i = 0; i < n; ++i)
                    if (((int32_t*) vv)[i] < ((int32_t*) attval)[0]
                        || ((int32_t*) vv)[i] > ((int32_t*) attval)[1])
                        v[i] = NAN;
            } else if (vartype == NC_UINT) {
                for (i = 0; i < n; ++i)
                    if (((uint32_t *) vv)[i] < ((uint32_t *) attval)[0]
                        || ((uint32_t *) vv)[i] > ((uint32_t *) attval)[1])
                        v[i] = NAN;
            } else if (vartype == NC_INT64) {
                for (i = 0; i < n; ++i)
                    if (((int64_t*) vv)[i] < ((int64_t*) attval)[0]
                        || ((int64_t*) vv)[i] > ((int64_t*) attval)[1])
                        v[i] = NAN;
            } else if (vartype == NC_UINT64) {
                for (i = 0; i < n; ++i)
                    if (((uint64_t *) vv)[i] < ((uint64_t *) attval)[0]
                        || ((uint64_t *) vv)[i] > ((uint64_t *) attval)[1])
                        v[i] = NAN;
            } else if (vartype == NC_FLOAT) {
                assert(sizeof(float) == 4);
                for (i = 0; i < n; ++i)
                    if (((float*) vv)[i] < ((float*) attval)[0]
                        || ((float*) vv)[i] > ((float*) attval)[1])
                        v[i] = NAN;
            } else if (vartype == NC_DOUBLE) {
                for (i = 0; i < n; ++i)
                    if (((double*) vv)[i] < ((double*) attval)[0]
                        || ((double*) vv)[i] > ((double*) attval)[1])
                        v[i] = NAN;
            } else
                quit("programming error");
        }

        if (vartype != NC_FLOAT)
            free(vv);
    }

    if (ncw_att_exists(ncid, varid, "scale_factor")) {
        float scale_factor;

        ncw_check_attlen(ncid, varid, "scale_factor", 1);
        ncw_get_att_float(ncid, varid, "scale_factor", &scale_factor);
        for (i = 0; i < n; ++i)
            v[i] *= scale_factor;
    }

    if (ncw_att_exists(ncid, varid, "add_offset")) {
        float add_offset;

        ncw_check_attlen(ncid, varid, "add_offset", 1);
        ncw_get_att_float(ncid, varid, "add_offset", &add_offset);
        for (i = 0; i < n; ++i)
            v[i] += add_offset;
    }

    ncw_close(ncid);
}

/** Writes one horizontal field (layer) for a variable to a NetCDF file.
 */
void ncu_writefield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v)
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    size_t start[4], count[4];
    size_t i, n;
    int hasrecorddim;

    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);

    if (ndims == 4) {
        if (nj == 0)
            quit("ncu_writefield(): \"%s\": %s: expected positive \"j\" dimension for a 4-dimensional variable\n", fname, varname);
        if (!hasrecorddim && dimlen[0] != 1)
            quit("ncu_writefield(): \"%s\": %s: for a 4-dimensional variable expected the first dimension to be either unlimited or of length 1\n", fname, varname);
        start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
        if (dimlen[1] != nk)
            quit("ncu_writefield(): \"%s\": vertical dimension of variable \"%s\" (nk = %d) does not match grid dimension (nk = %d)", fname, varname, dimlen[1], nk);
        start[1] = k;
        start[2] = 0;
        start[3] = 0;
        count[0] = 1;
        count[1] = 1;
        count[2] = dimlen[2];
        count[3] = dimlen[3];
        if (dimlen[3] != ni || dimlen[2] != nj)
            quit("ncu_writefield(): \"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[3], dimlen[2], ni, nj);
    } else if (ndims == 3) {
        if (nj > 0) {
            if (!hasrecorddim) {
                if (dimlen[0] != nk && !(dimlen[0] == 1 && (k == 0 || k == nk - 1)))
                    quit("ncu_writefield(): \"%s\": vertical dimension of variable \"%s\" (nk = %d) does not match grid dimension (nk = %d)", fname, varname, dimlen[0], nk);
                start[0] = k;
                start[1] = 0;
                start[2] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
                count[2] = dimlen[2];
            } else {
                /*
                 * 2D variable, ignore k
                 */
                start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
                start[1] = 0;
                start[2] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
                count[2] = dimlen[2];
            }
            if (dimlen[2] != ni || dimlen[1] != nj)
                quit("ncu_writefield(): \"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[2], dimlen[1], ni, nj);
        } else {
            if (!hasrecorddim && dimlen[0] != 1)
                quit("ncu_writefield(): \"%s\": %s: for a 3-dimensional variable on unstructured horizontal grid expected the first dimension to be either unlimited or of length 1\n", fname, varname);
            start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
            if (dimlen[1] != nk) {
                if (dimlen[1] != 1)
                    quit("ncu_writefield(): \"%s\": vertical dimension of variable \"%s\" (nk = %d) does not match grid dimension (nk = %d)", fname, varname, dimlen[1], nk);
                else
                    k = 0;      /* ignore k */
            }
            if (dimlen[2] != ni)
                quit("ncu_writefield(): \"%s\": horizontal dimension of variable \"%s\" (ni = %d) does not match grid dimension (ni = %d)", fname, varname, dimlen[2], ni);
            start[1] = k;
            start[2] = 0;
            count[0] = 1;
            count[1] = 1;
            count[2] = dimlen[2];
        }
    } else if (ndims == 2) {
        if (nj > 0) {
            if (hasrecorddim)
                quit("ncu_writefield(): \"%s\": can not write a layer to a 1D variable \"%s\"", fname, varname);
            /*
             * ignore k
             */
            start[0] = 0;
            start[1] = 0;
            count[0] = dimlen[0];
            count[1] = dimlen[1];
            if (dimlen[1] != ni || dimlen[0] != nj)
                quit("ncu_writefield(): \"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[1], dimlen[0], ni, nj);
        } else {
            if (!hasrecorddim) {
                if (dimlen[0] != nk && !(dimlen[0] == 1 && (k == 0 || k == nk - 1)))
                    quit("ncu_writefield(): \"%s\": vertical dimension of variable \"%s\" (nk = %d) does not match grid dimension (nk = %d)", fname, varname, dimlen[0], nk);
                start[0] = k;
                start[1] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
            } else {
                /*
                 * ignore k in this case
                 */
                start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
                start[1] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
            }
            if (dimlen[1] != ni)
                quit("ncu_writefield(): \"%s\": horizontal dimension of variable \"%s\" (ni = %d) does not match grid dimension (ni = %d)", fname, varname, dimlen[1], ni);
        }
    } else if (ndims == 1) {
        if (nj > 0) {
            quit("ncu_writefield(): %s: can not write 2D field for \"%s\": # of dimensions = %d", fname, varname, ndims);
        } else {
            if (hasrecorddim)
                quit("ncu_writefield(): %s: can not write layer to a 0D variable \"%s\"", fname, varname);
            /*
             * ignore k in this case
             */
            start[0] = 0;
            count[0] = dimlen[0];
        }
    } else
        quit("ncu_writefield(): \"%s\": can not write 2D field for \"%s\": # of dimensions = %d", fname, varname, ndims);

    n = 1;
    for (i = 0; i < ndims; ++i)
        n *= count[i];

    if (ncw_att_exists(ncid, varid, "add_offset")) {
        float add_offset;

        ncw_get_att_float(ncid, varid, "add_offset", &add_offset);
        for (i = 0; i < n; ++i)
            v[i] -= add_offset;
    }

    if (ncw_att_exists(ncid, varid, "scale_factor")) {
        float scale_factor;

        ncw_get_att_float(ncid, varid, "scale_factor", &scale_factor);
        for (i = 0; i < n; ++i)
            v[i] /= scale_factor;
    }

    {
        float attval[2];

        if (ncw_att_exists2(ncid, varid, "valid_min")) {
            ncw_check_attlen(ncid, varid, "valid_min", 1);
            ncw_get_att_float(ncid, varid, "valid_min", attval);
            for (i = 0; i < n; ++i)
                if (v[i] < attval[0])
                    v[i] = attval[0];
        }
        if (ncw_att_exists2(ncid, varid, "valid_max")) {
            ncw_check_attlen(ncid, varid, "valid_max", 1);
            ncw_get_att_float(ncid, varid, "valid_max", attval);
            for (i = 0; i < n; ++i)
                if (v[i] > attval[0])
                    v[i] = attval[0];
        }
        if (ncw_att_exists2(ncid, varid, "valid_range")) {
            ncw_check_attlen(ncid, varid, "valid_range", 2);
            ncw_get_att_float(ncid, varid, "valid_range", attval);
            for (i = 0; i < n; ++i)
                if (v[i] < attval[0])
                    v[i] = attval[0];
                else if (v[i] > attval[1])
                    v[i] = attval[1];
        }
        /*
         * This section does not always work as intended. E.g. for "int" dst
         * variable type _FillValue = -2147483647 is converted to -2.14748365e8f
         * and then in nc_put_vara_float() to -2147483648. Removing the code
         * below makes in that case no difference because NaNf is also converted
         * to -2147483648 in nc_put_vara_float().
         *
         * Yet sometimes the section below turns useful, e.g. in cases when src
         * and dst have the same data type but different fill values.
         */
        {
            int nofill = 0;
            
            if (ncw_att_exists2(ncid, varid, "_FillValue")) {
                ncw_check_attlen(ncid, varid, "_FillValue", 1);
                ncw_get_att_float(ncid, varid, "_FillValue", attval);
            } else
                ncw_inq_var_fill(ncid, varid, &nofill, attval);
            if (!nofill)
                for (i = 0; i < n; ++i)
                    if (isnan(v[i]))
                        v[i] = attval[0];
        }
        if (ncw_att_exists2(ncid, varid, "missing_value")) {
            ncw_check_attlen(ncid, varid, "missing_value", 1);
            ncw_get_att_float(ncid, varid, "missing_value", attval);
            for (i = 0; i < n; ++i)
                if (isnan(v[i]))
                    v[i] = attval[0];
        }
    }

    ncw_put_vara_float(ncid, varid, start, count, v);
    ncw_close(ncid);
}
