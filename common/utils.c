#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <stdint.h>
#include <errno.h>
#include <assert.h>
#if defined(MPI)
#include <mpi.h>
#endif
#include "ncw.h"
#include "ncutils.h"
#include "utils.h"

#define BASEYEAR 1970
#define BASEMONTH 1
#define BASEDAY 1

/**
 */
void quit(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "\n  error: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n\n");
#if defined(MPI)
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif

    exit(1);
}

/**
 */
int str2double(char* token, double* value)
{
    char* end = NULL;

    if (token == NULL) {
        *value = NAN;
        return 0;
    }

    *value = strtod(token, &end);

    if (end == token) {
        *value = NAN;
        return 0;
    }

    return 1;
}

/**
 */
int str2int(char* token, int* value)
{
    long int tmp;
    char* end = NULL;

    if (token == NULL) {
        *value = INT_MAX;
        return 0;
    }

    tmp = strtol(token, &end, 10);

    if (end == token || tmp > INT_MAX || tmp < INT_MIN) {
        *value = INT_MAX;
        return 0;
    }

    *value = (int) tmp;
    return 1;
}

/**
 */
void print_command(int argc, char* argv[])
{
    int i;
    char cwd[MAXSTRLEN];

    printf("    command = \"%s", argv[0]);
    for (i = 1; i < argc; ++i)
        printf(" %s", argv[i]);
    printf("\"\n");
    if (getcwd(cwd, MAXSTRLEN))
        printf("    dir = \"%s\"\n", cwd);
}

/**
 */
void print_time(const char offset[])
{
    time_t t;
    struct tm tm;

    t = time(NULL);
    tm = *localtime(&t);

    printf("%s%d-%02d-%02d %02d:%02d:%02d\n", offset, tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
}

/**
 */
char* get_command(int argc, char* argv[])
{
    char* cmd = NULL;
    int len = 0;
    int i;

    len = strlen(argv[0]);
    for (i = 1; i < argc; ++i)
        len += strlen(argv[i]) + 1;
    len++;

    cmd = malloc(len);
    strcpy(cmd, argv[0]);
    len = strlen(argv[0]);
    for (i = 1; i < argc; ++i) {
        cmd[len] = ' ';
        len++;
        strcpy(&cmd[len], argv[i]);
        len += strlen(argv[i]);
    }
    cmd[len] = 0;

    return cmd;
}

/**
 */
int file_exists(char* fname)
{
    FILE* f;

    f = fopen(fname, "r");
    if (f == NULL)
        return 0;
    fclose(f);
    return 1;
}

/**
 */
void file_rename(char oldname[], char newname[])
{
    int status = -1;

    status = rename(oldname, newname);
    if (status != 0) {
        int errno_saved = errno;

        quit("file_rename(): could not rename \"%s\" to \"%s\": %s", oldname, newname, strerror(errno_saved));
    }
}

/**
 */
void* alloc2d(size_t nj, size_t ni, size_t unitsize)
{
    size_t size;
    void* p;
    void** pp;
    int i;

    if (ni == 0 || nj == 0 || unitsize == 0)
        quit("alloc2d(): invalid size (nj = %zt, ni = %zt, unitsize = %zt)", nj, ni, unitsize);
    if (SIZE_MAX / nj / sizeof(void*) == 0 || SIZE_MAX / nj / ni / unitsize == 0 || nj * sizeof(void*) > SIZE_MAX - nj * ni * unitsize)
        quit("alloc2d: allocation size overflow: nj = %zt, ni = %zt, unitsize = %zt", nj, ni, unitsize);

    size = nj * sizeof(void*) + nj * ni * unitsize;
    if ((p = malloc(size)) == NULL) {
        int errno_saved = errno;

        quit("alloc2d(): %s", strerror(errno_saved));
    }
    memset(p, 0, size);

    pp = p;
    p = &((size_t*) p)[nj];
    for (i = 0; i < nj; i++)
        pp[i] = &((char*) p)[i * ni * unitsize];

    return pp;
}

/*
** scalar date routines    --    public domain by Ray Gardner
** Numerically, these will work over the range 1/01/01 thru 14699/12/31.
** Practically, these only work from the beginning of the Gregorian 
** calendar thru 14699/12/31.  The Gregorian calendar took effect in
** much of Europe in about 1582, some parts of Germany in about 1700, in
** England and the colonies in about 1752ff, and in Russia in 1918.
*/
static int isleap(unsigned yr)
{
    return yr % 400 == 0 || (yr % 4 == 0 && yr % 100 != 0);
}

static unsigned months_to_days(unsigned month)
{
    return (month * 3057 - 3007) / 100;
}

static long years_to_days(unsigned yr)
{
    return yr * 365L + yr / 4 - yr / 100 + yr / 400;
}

static long ymd_to_scalar(unsigned yr, unsigned mo, unsigned day)
{
    long scalar;

    scalar = day + months_to_days(mo);
    if (mo > 2)                 /* adjust if past February */
        scalar -= isleap(yr) ? 1 : 2;
    yr--;
    scalar += years_to_days(yr);
    return scalar;
}

/**
 */
static long int daydiff(unsigned int y1, unsigned int m1, unsigned int d1, unsigned int y2, unsigned int m2, unsigned int d2)
{
    long int dn1, dn2;

    dn1 = ymd_to_scalar(y1, m1, d1);
    dn2 = ymd_to_scalar(y2, m2, d2);

    return dn1 - dn2;
}

#define NTNAMES 4
static char* TNAMES[] = { "t",
    "time",
    "Time",
    "TIME"
};

#define NTUNITS 3
static char* TUNITS[] = { "seconds",
    "hours",
    "days"
};

/** Determine whether the variable is time variable.
 ** A variable is consedered to be "time" variable if
 **   (1) It has one of predefined names
 **   (2) It has dimension of 1 (or 0 perhaps)
 **   (3) It has atribute "units" that contains one of {"seconds", "hours",
 **       "days"} and "since".
 ** Note: we only need to give special treatment to "time" variable in case
 **   when the initial time is different across the merged files.
 */
int varistime(int ncid, int varid)
{
    char varname[NC_MAX_NAME];
    int ndims;
    size_t attlen;
    char attname[MAXSTRLEN];
    int i;

    ncw_inq_varndims(ncid, varid, &ndims);
    if (ndims > 1)
        return 0;

    ncw_inq_varname(ncid, varid, varname);
    for (i = 0; i < NTNAMES; ++i)
        if (strcmp(varname, TNAMES[i]) == 0)
            break;
    if (i == NTNAMES)
        return 0;

    if (!ncw_att_exists(ncid, varid, "units"))
        return 0;

    ncw_inq_attlen(ncid, varid, "units", &attlen);
    assert(attlen < MAXSTRLEN);
    ncw_get_att_text(ncid, varid, "units", attname);
    if (strstr(attname, "since") == NULL)
        return 0;

    for (i = 0; i < NTUNITS; ++i)
        if (strstr(attname, TUNITS[i]) != NULL)
            break;
    return (i < NTUNITS);
}

/**
 */
void tunits_convert(char* tunits1, char* tunits2, double* tunits_multiple, double* tunits_offset)
{
    char tunits1_copy[MAXSTRLEN];
    char* tunits[2] = { tunits1_copy, tunits2 };
    double multiple[2], offset[2];
    int i;

    strcpy(tunits1_copy, tunits1);
    for (i = 0; i < 2; ++i) {
        char* startdate;
        int year, month, day, h, m, s;
        char* token;
        char seps_date[] = " -\n";
        char seps_time[] = " :\n";

        if (strncasecmp(tunits[i], "fraction of a ", 14) == 0)
            tunits[i] += 14;

        if (strncasecmp(tunits[i], "sec", 3) == 0)
            multiple[i] = 86400.0;
        else if (strncasecmp(tunits[i], "hou", 3) == 0)
            multiple[i] = 24.0;
        else if (strncasecmp(tunits[i], "day", 3) == 0)
            multiple[i] = 1.0;
        else
            quit("can not interpret time units \"%s\"", tunits[i]);

        startdate = strstr(tunits[i], "since");
        if (startdate == NULL)
            quit("can not interpret time units \"%s\"", tunits[i]);
        startdate += strlen("since ");
        if ((token = strtok(startdate, seps_date)) == NULL)
            quit("can not interpret time units \"%s\"", tunits[i]);
        if (!str2int(token, &year))
            quit("could not convert \"%s\" to time units", tunits[i]);
        if ((token = strtok(NULL, seps_date)) == NULL)
            quit("can not interpret time units \"%s\"", tunits[i]);
        if (!str2int(token, &month))
            quit("could not convert \"%s\" to time units", tunits[i]);
        if ((token = strtok(NULL, seps_date)) == NULL)
            quit("can not interpret time units \"%s\"", tunits[i]);
        if (!str2int(token, &day))
            quit("could not convert \"%s\" to time units", tunits[i]);
        h = 0;
        m = 0;
        s = 0;
        if ((token = strtok(NULL, seps_time)) != NULL)
            if (!str2int(token, &h))
                quit("could not convert \"%s\" to time units", tunits[i]);
        if ((token = strtok(NULL, seps_time)) != NULL)
            if (!str2int(token, &m))
                quit("could not convert \"%s\" to time units", tunits[i]);
        if ((token = strtok(NULL, seps_time)) != NULL)
            if (!str2int(token, &s))
                quit("could not convert \"%s\" to time units", tunits[i]);

        offset[i] = (double) daydiff(year, month, day, BASEYEAR, BASEMONTH, BASEDAY) + (double) h / 24.0 + (double) m / 1440.0 + (double) s / 86400.0;
    }

    *tunits_multiple = multiple[0] / multiple[1];
    *tunits_offset = (offset[1] - offset[0]) * multiple[1];
}
