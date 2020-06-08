
IDIR =../include
NDIR=/usr/local/Cellar/netcdf/4.7.4
ODIR=../obj
LDIR =../lib

CC=gcc
CFLAGS=-I$(IDIR) -I$(NDIR)/include

LIBS=-lnetcdf -lm

_DEPS = main.h locate.h mesh.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o locate.o mesh.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 