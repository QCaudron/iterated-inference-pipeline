CC=gcc
CFLAGS=-fopenmp -std=gnu99 -Wall -O3 -DGSL_RANGE_CHECK_OFF -I kalman -I pmcmc -I simulation -I mif -I simplex -I core
LIB=libplom.a libplom_smc.a libplom_simplex.a libplom_mif.a libplom_pmcmc.a libplom_kalman.a libplom_ksimplex.a libplom_kmcmc.a libplom_simul.a libplom_worker.a
ALL_SRC= $(wildcard */*.c)
ALL_SRC_NO_TEMPLATE=$(filter-out templates/core_template.c templates/kalman_template.c templates/simulation_template.c, $(ALL_SRC))
SRC=$(filter-out smc/main.c simplex/main.c mif/main.c worker/main.c pmcmc/mainc.c kalman/main.c kalman/kmcmc_main.c kalman/simplex_kalmanc.c simulation/main.c, $(ALL_SRC_NO_TEMPLATE))
INCLUDES=$(wildcard */*.h)
OBJ= $(SRC:.c=.o)
ALL_OBJ= $(ALL_SRC_NO_TEMPLATE:.c=.o)

PREFIX=/usr/local

all: $(LIB)

libplom.a: $(OBJ)
	ar -rcs $@ $(OBJ)

%.o: %.c $(INCLUDES)
	$(CC) $(CFLAGS) -o $@ -c $<

libplom_smc.a: smc/main.o
	ar -rcs $@ $^

libplom_simplex.a: simplex/main.o
	ar -rcs $@ $^

libplom_mif.a: mif/main.o
	ar -rcs $@ $^

libplom_pmcmc.a: pmcmc/main.o
	ar -rcs $@ $^

libplom_kalman.a: kalman/main.o
	ar -rcs $@ $^

libplom_ksimplex.a: kalman/simplex_kalman.o
	ar -rcs $@ $^

libplom_kmcmc.a: kalman/kmcmc_main.o
	ar -rcs $@ $^

libplom_simul.a: simulation/main.o
	ar -rcs $@ $^

libplom_worker.a: worker/main.o
	ar -rcs $@ $^

.PHONY: clean uninstall install

install:
	test -d $(PREFIX)/lib || mkdir $(PREFIX)/lib; mv $(LIB) $(PREFIX)/lib/; \
	test -d $(PREFIX)/include || mkdir $(PREFIX)/include; cp $(INCLUDES) $(PREFIX)/include/

clean:
	rm $(ALL_OBJ)

uninstall:
	rm $(PREFIX)/{lib/libplom.a,lib/libplom_smc.a,lib/libplom_simplex.a,lib/libplom_mif.a,lib/libplom_pmcmc.a,lib/libplom_kalman.a,lib/libplom_ksimplex.a,lib/libplom_kmcmc.a,lib/libplom_simul.a,lib/libplom_worker.a,include/plom.h,include/mcmc_util.h,include/kalman.h,include/simplex.h,include/simulation.h,include/mif.h,include/pmcmc.h}
