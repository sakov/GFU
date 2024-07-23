# Set MPI status to MPI to compile with MPI, or leave it empty
MPISTATUS_NCAVE = MPI

CC = gcc
CFLAGS = -g  -Wall -pedantic -std=c99 -D_GNU_SOURCE -O2

INCS = -I $(HOME)/local/include -I common -I apps
LIBS = -lnetcdf -lhdf5 -lhdf5_hl -lm
# on Ubuntu the above may look as follows:
# LIB_NC = -lnetcdf -lhdf5_serial -lhdf5_serial_hl
LIBNN = -L $(HOME)/local/lib -lnn

# do not edit below

VERSION := $(shell sed 's/[^"]*"\([^"]*\)".*/\1/' common/version.h)
CCMPI = OMPI_MPICC=$(CC) mpicc
CFLAGSMPI = $(CFLAGS) -DMPI

PROGRAMS =\
bin/regrid_ll\
bin/nccat\
bin/ncave

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

default: bin $(PROGRAMS)

bin:
	mkdir -p bin

bin/regrid_ll: Makefile $(SRC_REGRID_LL) $(HDR_REGRID_LL)
	$(CC) $(CFLAGS) $(INCS) -o $@ $(SRC_REGRID_LL) $(LIBNN) $(LIBS)

bin/nccat: Makefile $(SRC_NCCAT) $(HDR_NCCAT)
	$(CC) $(CFLAGS) $(INCS) -o $@ $(SRC_NCCAT) $(LIBS)

bin/ncave: Makefile $(SRC_NCAVE) $(HDR_NCAVE)
	$(CC$(MPISTATUS_NCAVE)) $(CFLAGS$(MPISTATUS_NCAVE)) $(INCS) -o $@ $(SRC_NCAVE) $(LIBS)

clean:
	rm -f bin/*

tar:
	make clean; cd ..; tar -czvf gfu-v$(VERSION).tar.gz gfu; echo "  ->../gfu-v$(VERSION).tar.gz"

indent:
	indent -T size_t -T nc_type -T field -T int8_t -T int16_t -T int32_t -T int64_t -T delaunay */*.[ch]; rm -f */*.[ch]~
