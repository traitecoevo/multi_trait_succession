SHELL = /bin/sh

CC=g++
CFLAGS=-O3 -Wall #-O3#-pg  #-g for debugging -pg for profiling

INCS = -Iutils -I/usr/local/include
UTILS= utils/arg_parser.cpp utils/MatrixCPP.cpp utils/Utils.cpp utils/Collect.cpp utils/RungeKutta_CashKarp_45.cpp utils/gsl_root_class.cpp  utils/gsl_random_numbers.cpp
SOURCES= evolve.cpp Strategy.cpp AdaptiveDynamics_2D.cpp AdaptiveDynamics_nD.cpp PopulationDynamics.cpp FitnessSolver.cpp  EBT_Metapopulation.cpp EBT_Base.cpp  site.cpp
OBJECTS= $(UTILS:.cpp=.o) $(SOURCES:.cpp=.o)
LIBS = -L/usr/local/lib -lgsl -lgslcblas
EXECUTABLE=evolve

all: $(SOURCES) $(EXECUTABLE)
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) -c $(INCS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
