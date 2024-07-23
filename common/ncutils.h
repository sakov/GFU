/******************************************************************************
 *
 * File:        ncutils.h        
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
 *              This file is an excerpt from ncutils.h in EnKF-C package.
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_NCUTILS_H)

typedef void (*ncu_quit_fn) (char* format, ...);

void ncu_set_quitfn(ncu_quit_fn quit_fn);

/*
 * generic read procedures
 */
void ncu_readvarfloat(int ncid, int varid, size_t n, float v[]);

/*
 * model r/w procedures
 */
void ncu_readfield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v);
void ncu_writefield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v);

#define _NCUTILS_H
#endif
