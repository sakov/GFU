# Set MPI status to MPI to compile with MPI, or leave it empty
MPISTATUS_NCAVE = MPI

CC = gcc
CFLAGS = -g  -Wall -pedantic -std=c99 -D_GNU_SOURCE -O2

INCS = -I $(HOME)/local/include -I common -I apps
LIBNC = -lnetcdf -lhdf5 -lhdf5_hl
# on Ubuntu the above may look as follows:
# LIBNC = -lnetcdf -lhdf5_serial -lhdf5_serial_hl
LIBS = $(LIBNC) -lm
LIBNN = -L $(HOME)/local/lib -lnn

# do not edit below

VERSION := $(shell sed 's/[^"]*"\([^"]*\)".*/\1/' common/version.h)
CCMPI = OMPI_MPICC=$(CC) mpicc
CFLAGSMPI = $(CFLAGS) -DMPI

PROGRAMS =\
bin/regrid_ll\
bin/nccat\
bin/ncave\
bin/ncd2f\
bin/ncminmax\
bin/ncmask

SRC_REGRID_LL =\
apps/regrid_ll.c\
common/ncutils.c\
common/utils.c\
common/ncw.c

HDR_REGRID_LL =\
common/ncw.h\
common/ncutils.h\
common/utils.h\
common/version.h

SRC_NCCAT =\
apps/nccat.c\
common/utils.c\
common/ncw.c

HDR_NCCAT =\
common/utils.h\
common/ncw.h\
common/version.h

SRC_NCAVE =\
apps/ncave.c\
common/utils.c\
common/distribute.c\
common/ncutils.c\
common/ncw.c

HDR_NCAVE =\
common/distribute.h\
common/ncw.h\
common/ncutils.h\
common/version.h\
common/utils.h

SRC_NCD2F =\
apps/ncd2f.c\
common/ncutils.c\
common/utils.c\
common/ncw.c\
common/stringtable.c

HDR_NCD2F =\
common/ncw.h\
common/ncutils.h\
common/utils.h\
common/version.h\
common/stringtable.h

SRC_NCMINMAX =\
apps/ncminmax.c\
common/ncutils.c\
common/utils.c\
common/ncw.c

HDR_NCMINMAX =\
common/ncw.h\
common/ncutils.h\
common/utils.h\
common/version.h

SRC_NCMASK =\
apps/ncmask.c\
common/ncutils.c\
common/utils.c\
common/ncw.c

HDR_NCMASK =\
common/ncw.h\
common/ncutils.h\
common/utils.h\
common/version.h

default: bin $(PROGRAMS)

bin:
	mkdir -p bin

bin/regrid_ll: Makefile $(SRC_REGRID_LL) $(HDR_REGRID_LL)
	$(CC) $(CFLAGS) $(INCS) -o $@ $(SRC_REGRID_LL) $(LIBNN) $(LIBS)

bin/nccat: Makefile $(SRC_NCCAT) $(HDR_NCCAT)
	$(CC) $(CFLAGS) $(INCS) -o $@ $(SRC_NCCAT) $(LIBS)

bin/ncave: Makefile $(SRC_NCAVE) $(HDR_NCAVE)
	$(CC$(MPISTATUS_NCAVE)) $(CFLAGS$(MPISTATUS_NCAVE)) $(INCS) -o $@ $(SRC_NCAVE) $(LIBS)

bin/ncd2f: Makefile $(SRC_NCD2F) $(HDR_NCD2F)
	$(CC) $(CFLAGS) $(INCS) -o $@ $(SRC_NCD2F) $(LIBS)

bin/ncminmax: Makefile $(SRC_NCMINMAX) $(HDR_NCMINMAX)
	$(CC) $(CFLAGS) $(INCS) -o $@ $(SRC_NCMINMAX) $(LIBS)

bin/ncmask: Makefile $(SRC_NCMASK) $(HDR_NCMASK)
	$(CC) $(CFLAGS) $(INCS) -o $@ $(SRC_NCMASK) $(LIBS)

clean:
	rm -f bin/*

tar:
	make clean; cd ..; tar -czvf gfu-v$(VERSION).tar.gz gfu; echo "  ->../gfu-v$(VERSION).tar.gz"

indent:
	indent -T delaunay -T nc_type -T nctype2str -T field -T int8_t -T int16_t -T int32_t -T int64_t -T uint16_t -T uint32_t -T uint64_t -T size_t -T stringtable */*.[ch]; rm -f */*.[ch]~
