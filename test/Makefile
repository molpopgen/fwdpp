CXX=c++
DEBUG=-DNDEBUG
#DEBUG=-g
#DEBUG=
OPT=-O2  -funroll-all-loops
#OPT=-O2 -pg
#OPT=
CXXFLAGS=$(OPT) -Wall -W -I..  $(DEBUG)

TARGETS=phiab strobeck_morgan pmono gsl_rng_test gsl_rng_test2 gsl_rng_test_std gsl_rng_test2_std
all: gsl_rng_test.o gsl_rng_test2.o gsl_rng_test_std.o gsl_rng_test2_std.o phiab.o strobeck_morgan.o pmono.o
	$(CXX) $(CXXFLAGS) -o gsl_rng_test gsl_rng_test.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o gsl_rng_test2 gsl_rng_test2.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o gsl_rng_test_std gsl_rng_test_std.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o gsl_rng_test2_std gsl_rng_test2_std.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o strobeck_morgan strobeck_morgan.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CC) $(CCFLAGS) -o phiab phiab.o -lm
	$(CXX) $(CXXFLAGS) -o pmono pmono.o -lsequence

clean:
	rm -f *.o $(TARGETS)
