.SUFFIXES:
.SUFFIXES: .cpp .o

.cpp.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $<

CC=icc

CFLAGS=-std=c++11 -Wall

# enable debug output
CFLAGS+= -D__DEBUG

# armadillo settings
CFLAGS+= -DARMA_ALLOW_FAKE_GCC

# you may define a larger MAXNATOM for larger unit cell
#CFLAGS+= -DMAXNATOM=1024 # the default is 1024

INCLUDES=-I. -I$(ARMAINCLUDE) -I$(SPGINCLUDE)
LIBS=$(ARMALIBS) $(SPGLIBS)


ARMAINCLUDE=/opt/software/armadillo-10.8.0/include
ARMALIBS=-L/opt/software/armadillo-10.8.0 -larmadillo

SPGINCLUDE=/opt/software/spglib/1.16.5/include
SPGLIBS=-L/opt/software/spglib/1.16.5/lib -lsymspg
