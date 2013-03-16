CC=gcc
CFLAGS=-fopenmp -std=gnu99 -Wall -O3 -DGSL_RANGE_CHECK_OFF
LDFLAGS=-lplom -lplomtpl -lm -lgsl -lgslcblas -ljansson -lzmq
LIB=libplomtpl.a
SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)
LIBDIR=../lib

EXEC=simul simplex smc mif pmcmc kalman ksimplex kmcmc worker

all: $(LIB)

install: $(EXEC)
	mv $(EXEC) ../../

$(LIB): $(OBJ)
	ar -rcs $@ $(OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

simul: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@  $(LDFLAGS) -lplom_simul 

simplex: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@  $(LDFLAGS) -lplom_simplex

smc: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@  $(LDFLAGS) -lplom_smc 

mif: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@  $(LDFLAGS) -lplom_mif 

pmcmc: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@  $(LDFLAGS) -lplom_pmcmc 

kalman: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@  $(LDFLAGS) -lplom_kalman

ksimplex: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@  $(LDFLAGS) -lplom_ksimplex

kmcmc: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@  $(LDFLAGS) -lplom_kmcmc

worker: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@  $(LDFLAGS) -lplom_worker 

.PHONY: clean

clean:
	rm *.o $(LIB)

uninstall: clean
	rm ../../{simul,simplex,smc,mif,pmcmc,kalman,ksimplex,kmcmc,worker}
