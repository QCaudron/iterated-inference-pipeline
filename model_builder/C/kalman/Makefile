CC=gcc
CFLAGS=-fopenmp -std=gnu99 -Wall -O3 -I ../core/ -DGSL_RANGE_CHECK_OFF
LDFLAGS=-lplomtpl -lplom -lm -lgsl -lgslcblas -ljansson -lzmq
EXEC=kalman ksimplex kmcmc
ALL_SRC= $(wildcard *.c)
ALL_OBJ= $(ALL_SRC:.c=.o)

SRC=$(filter-out main.c simplex_kalman.c kmcmc_main.c kmcmc.c, $(ALL_SRC))
OBJ= $(SRC:.c=.o)
BINDIR=../..

all: $(EXEC)

obj: $(ALL_OBJ)

kalman: $(OBJ) main.o
	$(CC) $(CFLAGS) -L../lib -o $@ $^ $(LDFLAGS)

ksimplex: $(OBJ) simplex_kalman.o
	$(CC) $(CFLAGS) -L../lib -o $@ $^ $(LDFLAGS)

kmcmc: $(OBJ) kmcmc_main.o kmcmc.o
	$(CC) $(CFLAGS) -L../lib -o $@ $^ $(LDFLAGS)

%.o: %.c kalman.h
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY: clean mrproper install

install:
	test -d $(BINDIR) || mkdir $(BINDIR); mv $(EXEC) $(BINDIR)/

clean:
	rm *.o

mrproper: clean
	rm $(BINDIR)/$(EXEC)
