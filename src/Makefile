
IDIR =../include
NDIR=/usr/local/Cellar/netcdf/4.7.4
ODIR=../obj
LDIR =../lib
XDIR=../exec

CC=gcc
CFLAGS=-I$(IDIR) -I$(NDIR)/include

LIBS=-lnetcdf -lm

_DEPS = main.h mesh.h interpolate.h particle.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o mesh.o interpolate.o particle.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(XDIR)/main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 