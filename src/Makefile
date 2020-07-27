
IDIR =../include
NDIR=/usr/local/Cellar/netcdf/4.7.4
CAIRO_DIR=/usr/local/Cellar/cairo/1.16.0_3
FFMPEG_DIR=/usr/local/Cellar/ffmpeg/4.3_1
ODIR=../obj
LDIR =../lib
XDIR=../exec

CC=gcc
CFLAGS=-I$(IDIR) -I$(NDIR)/include -I$(CAIRO_DIR)/include/cairo -I$(FFMPEG_DIR)/include -DNEWFFMPEG=1

LIBS=-lnetcdf -lcairo -lavformat -lavcodec -lswresample -lswscale -lavutil -lm

_DEPS = main.h mesh.h interpolate.h particle.h visualize.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o mesh.o interpolate.o particle.o visualize.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(XDIR)/main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 