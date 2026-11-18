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
 *              This file initially was (and still mainly is) an excerpt from
 *              ncutils.c in EnKF-C package.
 *
 * Revisions:   PS 20261116 Consolidated typed read/write attribute handling
 *              and 2D slice logic.
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

/** Converts a value of the specified NetCDF type to double.
 */
static double ncu_todouble(nc_type type, const void* p)
{
    switch (type) {
    case NC_BYTE:
        return (double) ((const int8_t*) p)[0];
    case NC_CHAR:
        return (double) ((const uint8_t*) p)[0];
    case NC_UBYTE:
        return (double) ((const uint8_t*) p)[0];
    case NC_SHORT:
        return (double) ((const int16_t*) p)[0];
    case NC_USHORT:
        return (double) ((const uint16_t*) p)[0];
    case NC_INT:
        return (double) ((const int32_t*) p)[0];
    case NC_UINT:
        return (double) ((const uint32_t *) p)[0];
    case NC_INT64:
        return (double) ((const int64_t*) p)[0];
    case NC_UINT64:
        return (double) ((const uint64_t*) p)[0];
    case NC_FLOAT:
        return (double) ((const float*) p)[0];
    case NC_DOUBLE:
        return (double) ((const double*) p)[0];
    default:
        quit("ncu_todouble(): unexpected NetCDF type %d", (int) type);
    }
    return NAN;
}

typedef union {
    int8_t i8;
    uint8_t u8;
    int16_t i16;
    uint16_t u16;
    int32_t i32;
    uint32_t u32;
    int64_t i64;
    uint64_t u64;
    float f32;
    double f64;
} ncu_val;

/** Applies _FillValue, missing_value, valid_min, valid_max, and valid_range
 ** masks from raw typed buffer `vv` onto either a float[] or double[] array.
 */
static void ncu_apply_masks(int ncid, int varid, nc_type vartype, size_t n, const void* vv, void* out, int is_double)
{
    int has_fill = 0, has_miss = 0, has_vmin = 0, has_vmax = 0, has_vrange = 0;
    double fill = NAN, miss = NAN, vmin = NAN, vmax = NAN;
    double vrange[2] = { NAN, NAN };

    if (ncw_att_exists2(ncid, varid, "_FillValue")) {
        ncw_check_attlen(ncid, varid, "_FillValue", 1);
        ncw_get_att_double(ncid, varid, "_FillValue", &fill);
        has_fill = 1;
    } else {
        int nofill = 1;
        ncu_val fill_buf = { 0 };

        ncw_inq_var_fill(ncid, varid, &nofill, &fill_buf);
        if (!nofill) {
            fill = ncu_todouble(vartype, &fill_buf);
            has_fill = 1;
        }
    }

    if (ncw_att_exists2(ncid, varid, "missing_value")) {
        ncw_check_attlen(ncid, varid, "missing_value", 1);
        ncw_get_att_double(ncid, varid, "missing_value", &miss);
        has_miss = 1;
    }
    if (ncw_att_exists2(ncid, varid, "valid_min")) {
        ncw_check_attlen(ncid, varid, "valid_min", 1);
        ncw_get_att_double(ncid, varid, "valid_min", &vmin);
        has_vmin = 1;
    }
    if (ncw_att_exists2(ncid, varid, "valid_max")) {
        ncw_check_attlen(ncid, varid, "valid_max", 1);
        ncw_get_att_double(ncid, varid, "valid_max", &vmax);
        has_vmax = 1;
    }
    if (ncw_att_exists2(ncid, varid, "valid_range")) {
        ncw_check_attlen(ncid, varid, "valid_range", 2);
        ncw_get_att_double(ncid, varid, "valid_range", vrange);
        has_vrange = 1;
    }

    if (!has_fill && !has_miss && !has_vmin && !has_vmax && !has_vrange)
        return;

#define APPLY_MASK_LOOP(CTYPE) do {                                            \
        const CTYPE* raw = (const CTYPE*) vv;                                  \
        double vr_min = vrange[0];                                             \
        double vr_max = vrange[1];                                             \
        for (size_t i = 0; i < n; ++i) {                                       \
            double val = (double) raw[i];                                      \
            int invalid = (has_fill   && val == fill) ||                       \
                          (has_miss   && val == miss) ||                       \
                          (has_vmin   && val < vmin)  ||                       \
                          (has_vmax   && val > vmax)  ||                       \
                          (has_vrange && (val < vr_min || val > vr_max));      \
            if (invalid) {                                                     \
                if (is_double) ((double*) out)[i] = NAN;                       \
                else           ((float*)  out)[i] = NAN;                       \
            }                                                                  \
        }                                                                      \
    } while (0)

    switch (vartype) {
    case NC_BYTE:
        APPLY_MASK_LOOP(int8_t);

        break;
    case NC_CHAR:
    case NC_UBYTE:
        APPLY_MASK_LOOP(uint8_t);

        break;
    case NC_SHORT:
        APPLY_MASK_LOOP(int16_t);

        break;
    case NC_USHORT:
        APPLY_MASK_LOOP(uint16_t);

        break;
    case NC_INT:
        APPLY_MASK_LOOP(int32_t);

        break;
    case NC_UINT:
        APPLY_MASK_LOOP(uint32_t);
        break;
    case NC_INT64:
        APPLY_MASK_LOOP(int64_t);

        break;
    case NC_UINT64:
        APPLY_MASK_LOOP(uint64_t);

        break;
    case NC_FLOAT:
        APPLY_MASK_LOOP(float);

        break;
    case NC_DOUBLE:
        APPLY_MASK_LOOP(double);

        break;
    default:
        quit("ncu_apply_masks(): unexpected NetCDF type %d", (int) vartype);
    }

#undef APPLY_MASK_LOOP
}

/** Allocates raw buffer if needed, reads raw values, applies masks, and unpacks
 ** scale_factor / add_offset onto an already-populated float[] or double[]
 ** array.
 */
static void ncu_postprocess_read(int ncid, int varid, const size_t start[], const size_t count[], size_t n, void* out, int is_double, int is_full_var)
{
    nc_type vartype;
    nc_type target_type = is_double ? NC_DOUBLE : NC_FLOAT;
    int needs_raw = ncw_att_exists2(ncid, varid, "_FillValue") || ncw_att_exists2(ncid, varid, "missing_value") || ncw_att_exists2(ncid, varid, "valid_range") || ncw_att_exists2(ncid, varid, "valid_min") || ncw_att_exists2(ncid, varid, "valid_max");

    ncw_inq_vartype(ncid, varid, &vartype);

    /*
     * Check implicit NetCDF fill value even if _FillValue attribute is
     * absent 
     */
    if (!needs_raw && !ncw_att_exists2(ncid, varid, "_FillValue")) {
        int nofill = 1;
        ncu_val dummy[2];

        ncw_inq_var_fill(ncid, varid, &nofill, dummy);
        if (!nofill)
            needs_raw = 1;
    }

    if (needs_raw) {
        void* vv = out;

        if (vartype != target_type) {
            int typesize = ncw_sizeof(vartype);

            vv = malloc(n * (size_t) typesize);
            if (vv == NULL)
                quit("ncutils: could not allocate memory (%zu bytes)", n * (size_t) typesize);
            if (is_full_var)
                ncw_get_var(ncid, varid, vv);
            else
                ncw_get_vara(ncid, varid, start, count, vv);
        }

        ncu_apply_masks(ncid, varid, vartype, n, vv, out, is_double);

        if (vv != out)
            free(vv);
    }

    if (ncw_att_exists(ncid, varid, "scale_factor")) {
        ncw_check_attlen(ncid, varid, "scale_factor", 1);
        if (is_double) {
            double sf;

            ncw_get_att_double(ncid, varid, "scale_factor", &sf);
            for (size_t i = 0; i < n; ++i)
                ((double*) out)[i] *= sf;
        } else {
            float sf;

            ncw_get_att_float(ncid, varid, "scale_factor", &sf);
            for (size_t i = 0; i < n; ++i)
                ((float*) out)[i] *= sf;
        }
    }

    if (ncw_att_exists(ncid, varid, "add_offset")) {
        ncw_check_attlen(ncid, varid, "add_offset", 1);
        if (is_double) {
            double ao;

            ncw_get_att_double(ncid, varid, "add_offset", &ao);
            for (size_t i = 0; i < n; ++i)
                ((double*) out)[i] += ao;
        } else {
            float ao;

            ncw_get_att_float(ncid, varid, "add_offset", &ao);
            for (size_t i = 0; i < n; ++i)
                ((float*) out)[i] += ao;
        }
    }
}

/** Rounds values to be written to a NetCDF variable of an integer type
 ** (otherwise the NetCDF library truncates them).
 */
static void ncu_roundforint(nc_type type, size_t n, void* v, int is_double)
{
    size_t i;

    if (type == NC_FLOAT || type == NC_DOUBLE || type == NC_CHAR)
        return;
    if (is_double) {
        double* dv = (double*) v;

        for (i = 0; i < n; ++i)
            dv[i] = round(dv[i]);
    } else {
        float* fv = (float*) v;

        for (i = 0; i < n; ++i)
            fv[i] = roundf(fv[i]);
    }
}

/** Applies inverse add_offset/scale_factor packing, rounds for integer types,
 ** and clamps to valid_min/valid_max/valid_range before writing.
 ** Leaves NaNs intact so ncu_put_vara() can write exact fill values.
 */
static void ncu_preprocess_write(int ncid, int varid, size_t n, void* v, int is_double)
{
    nc_type vartype;
    size_t i;

    if (ncw_att_exists(ncid, varid, "add_offset")) {
        ncw_check_attlen(ncid, varid, "add_offset", 1);
        if (is_double) {
            double ao;

            ncw_get_att_double(ncid, varid, "add_offset", &ao);
            for (i = 0; i < n; ++i)
                ((double*) v)[i] -= ao;
        } else {
            float ao;

            ncw_get_att_float(ncid, varid, "add_offset", &ao);
            for (i = 0; i < n; ++i)
                ((float*) v)[i] -= ao;
        }
    }

    if (ncw_att_exists(ncid, varid, "scale_factor")) {
        ncw_check_attlen(ncid, varid, "scale_factor", 1);
        if (is_double) {
            double sf;

            ncw_get_att_double(ncid, varid, "scale_factor", &sf);
            for (i = 0; i < n; ++i)
                ((double*) v)[i] /= sf;
        } else {
            float sf;

            ncw_get_att_float(ncid, varid, "scale_factor", &sf);
            for (i = 0; i < n; ++i)
                ((float*) v)[i] /= sf;
        }
    }

    ncw_inq_vartype(ncid, varid, &vartype);
    ncu_roundforint(vartype, n, v, is_double);

#define CLAMP_WRITE_ATTS(CTYPE, GET_ATT_FN) do {                               \
        CTYPE* arr = (CTYPE*) v;                                               \
        CTYPE attval[2];                                                       \
        if (ncw_att_exists2(ncid, varid, "valid_min")) {                       \
            ncw_check_attlen(ncid, varid, "valid_min", 1);                     \
            GET_ATT_FN(ncid, varid, "valid_min", attval);                      \
            for (i = 0; i < n; ++i)                                            \
                if (arr[i] < attval[0]) arr[i] = attval[0];                    \
        }                                                                      \
        if (ncw_att_exists2(ncid, varid, "valid_max")) {                       \
            ncw_check_attlen(ncid, varid, "valid_max", 1);                     \
            GET_ATT_FN(ncid, varid, "valid_max", attval);                      \
            for (i = 0; i < n; ++i)                                            \
                if (arr[i] > attval[0]) arr[i] = attval[0];                    \
        }                                                                      \
        if (ncw_att_exists2(ncid, varid, "valid_range")) {                     \
            ncw_check_attlen(ncid, varid, "valid_range", 2);                   \
            GET_ATT_FN(ncid, varid, "valid_range", attval);                    \
            for (i = 0; i < n; ++i)                                            \
                if (arr[i] < attval[0]) arr[i] = attval[0];                    \
                else if (arr[i] > attval[1]) arr[i] = attval[1];               \
        }                                                                      \
    } while (0)

    if (is_double)
        CLAMP_WRITE_ATTS(double, ncw_get_att_double);

    else
        CLAMP_WRITE_ATTS(float, ncw_get_att_float);

#undef CLAMP_WRITE_ATTS
}

/** Resolves the fill value for a variable (checking _FillValue, implicit
 ** NetCDF fill, and missing_value in order). Returns NAN if none is set.
 */
static double ncu_get_write_fill(int ncid, int varid, nc_type vartype)
{
    if (ncw_att_exists2(ncid, varid, "_FillValue")) {
        double fill = NAN;

        ncw_check_attlen(ncid, varid, "_FillValue", 1);
        ncw_get_att_double(ncid, varid, "_FillValue", &fill);
        return fill;
    } else {
        int nofill = 1;
        ncu_val buf[2];

        ncw_inq_var_fill(ncid, varid, &nofill, buf);
        if (!nofill)
            return ncu_todouble(vartype, buf);
    }

    if (ncw_att_exists2(ncid, varid, "missing_value")) {
        double miss = NAN;

        ncw_check_attlen(ncid, varid, "missing_value", 1);
        ncw_get_att_double(ncid, varid, "missing_value", &miss);
        return miss;
    }

    return NAN;
}

/** Writes a float array to a NetCDF variable, replacing NaNs with the exact
 ** fill value. Uses double for 32/64-bit integer types so large integer fill
 ** values are never rounded by float precision.
 */
static void ncu_put_vara(int ncid, int varid, const size_t start[], const size_t count[], size_t n, void* v, int is_double)
{
    nc_type vartype;
    double fill;
    size_t i;

    ncw_inq_vartype(ncid, varid, &vartype);
    fill = ncu_get_write_fill(ncid, varid, vartype);

    if (is_double) {
        double* dv = (double*) v;

        if (!isnan(fill)) {
            for (i = 0; i < n; ++i)
                if (isnan(dv[i]))
                    dv[i] = fill;
        }
        ncw_put_vara_double(ncid, varid, start, count, dv);
    } else if (vartype == NC_INT || vartype == NC_UINT || vartype == NC_INT64 || vartype == NC_UINT64) {
        float* fv = (float*) v;
        double* dv = malloc(n * sizeof(double));

        if (dv == NULL)
            quit("ncu_put_vara(): could not allocate memory (%zu bytes)", n * sizeof(double));
        for (i = 0; i < n; ++i)
            dv[i] = isnan(fv[i]) ? fill : (double) fv[i];
        ncw_put_vara_double(ncid, varid, start, count, dv);
        free(dv);
    } else {
        float* fv = (float*) v;

        if (!isnan(fill)) {
            float ffill = (float) fill;

            for (i = 0; i < n; ++i)
                if (isnan(fv[i]))
                    fv[i] = ffill;
        }
        ncw_put_vara_float(ncid, varid, start, count, fv);
    }
}

/** Validates variable dimensions against grid dimensions (ni, nj, nk) and
 ** computes the NetCDF hyperslab start[], count[], and total element count *n
 ** for reading or writing a single horizontal layer k.
 */
static void ncu_get_field_slice(const char* caller, const char* fname, const char* varname, int ncid, int varid, int k, int ni, int nj, int nk, int is_write, size_t start[4], size_t count[4], size_t* n)
{
    int ndims;
    size_t dimlen[4];
    int hasrecorddim;
    int i;

    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);

    if (!is_write && hasrecorddim && dimlen[0] == 0)
        quit("%s(): \"%s\": %s: empty record dimension", caller, fname, varname);

    if (ndims == 4) {
        if (nj == 0)
            quit("%s(): \"%s\": %s: expected positive \"j\" dimension for a 4-dimensional variable\n", caller, fname, varname);
        if (!hasrecorddim && dimlen[0] != 1)
            quit("%s(): \"%s\": %s: for a 4-dimensional variable expected the first dimension to be either unlimited or of length 1\n", caller, fname, varname);
        start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
        if (dimlen[1] != (size_t) nk) {
            if (dimlen[1] != 1)
                quit("%s(): \"%s\": vertical dimension of variable \"%s\" (nk = %zu) does not match grid dimension (nk = %d)", caller, fname, varname, dimlen[1], nk);
            else
                k = 0;          /* ignore k */
        }
        if ((ni >= 0 && nj >= 0) && (dimlen[3] != (size_t) ni || dimlen[2] != (size_t) nj))
            quit("%s(): \"%s\": horizontal dimensions of variable \"%s\" (ni = %zu, nj = %zu) do not match grid dimensions (ni = %d, nj = %d)", caller, fname, varname, dimlen[3], dimlen[2], ni, nj);
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
                if (nk >= 0 && dimlen[0] != (size_t) nk && !(dimlen[0] == 1 && (k == 0 || k == nk - 1)))
                    quit("%s(): \"%s\": vertical dimension of variable \"%s\" (nk = %zu) does not match grid dimension (nk = %d)", caller, fname, varname, dimlen[0], nk);
                start[0] = (dimlen[0] == 1) ? 0 : k;
                start[1] = 0;
                start[2] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
                count[2] = dimlen[2];
            } else {
                /*
                 * 2D variable with record dimension, ignore k
                 */
                start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
                start[1] = 0;
                start[2] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
                count[2] = dimlen[2];
            }
            if ((ni >= 0 && nj >= 0) && (dimlen[2] != (size_t) ni || dimlen[1] != (size_t) nj))
                quit("%s(): \"%s\": horizontal dimensions of variable \"%s\" (ni = %zu, nj = %zu) do not match grid dimensions (ni = %d, nj = %d)", caller, fname, varname, dimlen[2], dimlen[1], ni, nj);
        } else {
            if (!hasrecorddim && dimlen[0] != 1)
                quit("%s(): \"%s\": %s: for a 3-dimensional variable on unstructured horizontal grid expected the first dimension to be either unlimited or of length 1\n", caller, fname, varname);
            start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
            if (nk >= 0 && dimlen[1] != (size_t) nk) {
                if (dimlen[1] != 1)
                    quit("%s(): \"%s\": vertical dimension of variable \"%s\" (nk = %zu) does not match grid dimension (nk = %d)", caller, fname, varname, dimlen[1], nk);
                else
                    k = 0;      /* ignore k */
            }
            if (ni >= 0 && dimlen[2] != (size_t) ni)
                quit("%s(): \"%s\": horizontal dimension of variable \"%s\" (ni = %zu) does not match grid dimension (ni = %d)", caller, fname, varname, dimlen[2], ni);
            start[1] = k;
            start[2] = 0;
            count[0] = 1;
            count[1] = dimlen[1];
            count[2] = dimlen[2];
        }
    } else if (ndims == 2) {
        if (nj > 0) {
            if (hasrecorddim)
                quit("%s(): \"%s\": can not %s a layer %s a 1D variable \"%s\"", caller, fname, is_write ? "write" : "read", is_write ? "to" : "from", varname);
            /*
             * ignore k
             */
            start[0] = 0;
            start[1] = 0;
            count[0] = dimlen[0];
            count[1] = dimlen[1];
            if ((ni >= 0 && nj >= 0) && (dimlen[1] != (size_t) ni || dimlen[0] != (size_t) nj))
                quit("%s(): \"%s\": horizontal dimensions of variable \"%s\" (ni = %zu, nj = %zu) do not match grid dimensions (ni = %d, nj = %d)", caller, fname, varname, dimlen[1], dimlen[0], ni, nj);
        } else {
            if (!hasrecorddim) {
                if (dimlen[0] != (size_t) nk && !(dimlen[0] == 1 && (k == 0 || k == nk - 1)))
                    quit("%s(): \"%s\": vertical dimension of variable \"%s\" (nk = %zu) does not match grid dimension (nk = %d)", caller, fname, varname, dimlen[0], nk);
                start[0] = (dimlen[0] == 1) ? 0 : k;
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
            if (ni >= 0 && dimlen[1] != (size_t) ni)
                quit("%s(): \"%s\": horizontal dimension of variable \"%s\" (ni = %zu) does not match grid dimension (ni = %d)", caller, fname, varname, dimlen[1], ni);
        }
    } else if (ndims == 1) {
        if (nj > 0) {
            quit("%s(): %s: can not %s 2D field for \"%s\": # of dimensions = %d", caller, fname, is_write ? "write" : "read", varname, ndims);
        } else {
            if (hasrecorddim)
                quit("%s(): %s: can not %s layer %s a 0D variable \"%s\"", caller, fname, is_write ? "write" : "read", is_write ? "to" : "from", varname);
            /*
             * ignore k in this case
             */
            start[0] = 0;
            count[0] = dimlen[0];
        }
    } else
        quit("%s(): \"%s\": can not %s 2D field for \"%s\": # of dimensions = %d", caller, fname, is_write ? "write" : "read", varname, ndims);

    *n = 1;
    for (i = 0; i < ndims; ++i)
        *n *= count[i];
}

/**
 */
int ncu_getnfields(char fname[], char varname[])
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    int hasrecorddim;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = (ncw_var_hasunlimdim(ncid, varid) || dimlen[0] == 1);
    ncw_close(ncid);

    if (ndims > 4)
        quit("%s: %s: do not know how to handle more than 4-dimensional variables\n", fname, varname);
    if (ndims == 4) {
        if (!hasrecorddim)
            quit("%s: %s: expect an unlimited dimension to be present for a 4-dimensional variable\n", fname, varname);
        return (int) dimlen[1];
    }
    if (ndims == 3)
        return (hasrecorddim) ? 1 : (int) dimlen[0];
    if (ndims == 2)
        return (hasrecorddim) ? 0 : 1;

    return 0;
}

/**
 */
void ncu_readvarfloat(int ncid, int varid, size_t n, float v[])
{
    ncw_check_varsize(ncid, varid, n);
    ncw_get_var_float_fixerange(ncid, varid, v);
    ncu_postprocess_read(ncid, varid, NULL, NULL, n, v, 0, 1);
}

/**
 */
void ncu_readvardouble(int ncid, int varid, size_t n, double v[])
{
    ncw_check_varsize(ncid, varid, n);
    ncw_get_var_double(ncid, varid, v);
    ncu_postprocess_read(ncid, varid, NULL, NULL, n, v, 1, 1);
}

/** Reads one horizontal field (layer) for a variable from a NetCDF file.
 ** Verifies that the field dimensions are ni x nj.
 */
void ncu_readfield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v)
{
    int ncid, varid;
    size_t start[4], count[4], n;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncu_get_field_slice("ncu_readfield", fname, varname, ncid, varid, k, ni, nj, nk, 0, start, count, &n);
    ncw_get_vara_float_fixerange(ncid, varid, start, count, v);
    ncu_postprocess_read(ncid, varid, start, count, n, v, 0, 0);
    ncw_close(ncid);
}

/** Writes one horizontal field (layer) for a variable to a NetCDF file.
 */
void ncu_writefield_double(char fname[], char varname[], int k, int ni, int nj, int nk, double* v)
{
    int ncid, varid;
    size_t start[4], count[4], n;

    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncu_get_field_slice("ncu_writefield_double", fname, varname, ncid, varid, k, ni, nj, nk, 1, start, count, &n);
    ncu_preprocess_write(ncid, varid, n, v, 1);
    ncu_put_vara(ncid, varid, start, count, n, v, 1);
    ncw_close(ncid);
}

/** Writes one horizontal field (layer) for a variable to a NetCDF file.
 */
void ncu_writefield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v)
{
    int ncid, varid;
    size_t start[4], count[4], n;

    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncu_get_field_slice("ncu_writefield", fname, varname, ncid, varid, k, ni, nj, nk, 1, start, count, &n);
    ncu_preprocess_write(ncid, varid, n, v, 0);
    ncu_put_vara(ncid, varid, start, count, n, v, 0);
    ncw_close(ncid);
}

/**
 */
void ncu_readfield_double(char fname[], char varname[], int k, int ni, int nj, int nk, double* v)
{
    int ncid, varid;
    size_t start[4], count[4], n;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncu_get_field_slice("ncu_readfield_double", fname, varname, ncid, varid, k, ni, nj, nk, 0, start, count, &n);
    ncw_get_vara_double(ncid, varid, start, count, v);
    ncu_postprocess_read(ncid, varid, start, count, n, v, 1, 0);
    ncw_close(ncid);
}
