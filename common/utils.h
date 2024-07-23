#if !defined(_UTILS_H)

#define MAXSTRLEN 4096

void quit(char* format, ...);
int str2double(char* token, double* value);
int str2int(char* token, int* value);
void print_command(int argc, char* argv[]);
char* get_command(int argc, char* argv[]);
void print_time(const char offset[]);
int file_exists(char* fname);
void file_rename(char oldname[], char newname[]);
void* alloc2d(size_t nj, size_t ni, size_t unitsize);
int varistime(int ncid, int varid);
void tunits_convert(char* tunits1, char* tunits2, double* tunits_multiple, double* tunits_offset);

#define UTILS_H
#endif
