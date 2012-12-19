CC=gcc
CFLAGS=-fopenmp -std=gnu99 -Wall -O3 -I ../core/ -I ../simulation/ -I ../kalman/ -I ../pmcmc/ -DGSL_RANGE_CHECK_OFF
LDFLAGS=-lm -lgsl -lgslcblas -ljansson -lzmq
LIB=libplomtpl.a
SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)
LIBDIR=../lib

all: $(EXEC) $(LIB)

$(LIB): $(OBJ)
	ar -rcs $@ $(OBJ)

%.o: %.c plom.h simulation.h kalman.h pmcmc.h
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY: clean mrproper install

install:
	test -d $(LIBDIR) || mkdir $(LIBDIR); mv $(LIB) $(LIBDIR)/

clean:
	rm *.o

mrproper: clean
	rm $(LIBDIR)/$(LIB)
