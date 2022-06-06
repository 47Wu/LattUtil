include ../make.inc

objs = utils.o input.o regstruct.o supercell.o

all: lattutil

lattutil: lattutil.cpp $(objs)
	$(CC) $(CFLAGS) -o LattUtil.x lattutil.cpp $(INCLUDES) $(objs) $(LIBS)
	cp LattUtil.x ../bin/

clean:
	rm -rf *.o *.mod *.x
	rm -f ../bin/*.x

remake: clean
	./cp_src.sh
	make all
