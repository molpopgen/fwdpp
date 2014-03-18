CXX=c++
DEBUG=-DNDEBUG
#DEBUG=-g
#DEBUG=
OPT=-O2  -funroll-all-loops
#OPT=-O2 -pg
#OPT=
CXXFLAGS=$(OPT) -Wall -W -I..  $(DEBUG)

TARGETS=diploid diploid_ind diploid_fixed_sh diploid_twopop_mig diploid_binaryIO RHH bneck_selection migsel_ind diploid_binaryIO_ind diploid_fixed_sh_ind gsl_rng_test gsl_rng_test2 TFL2013
all: diploid.o diploid_ind.o diploid_binaryIO.o diploid_binaryIO_ind.o diploid_fixed_sh.o diploid_twopop_mig.o RHH.o diploid_fixed_sh_ind.o bneck_selection.o migsel_ind.o gsl_rng_test.o gsl_rng_test2.o TFL2013.o
	$(CXX) $(CXXFLAGS) -o diploid diploid.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_ind diploid_ind.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o RHH RHH.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_binaryIO diploid_binaryIO.o -lsequence -lz -lgsl -lboost_system -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_binaryIO_ind diploid_binaryIO_ind.o -lsequence -lz -lgsl -lboost_system -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_fixed_sh diploid_fixed_sh.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_fixed_sh_ind diploid_fixed_sh_ind.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_twopop_mig diploid_twopop_mig.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o bneck_selection bneck_selection.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o migsel_ind migsel_ind.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o TFL2013 TFL2013.o -lsequence -lz -lgsl -lgslcblas -lboost_system -lboost_iostreams $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o gsl_rng_test gsl_rng_test.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o gsl_rng_test2 gsl_rng_test2.o -lgsl -lgslcblas $(LDFLAGS)

clean:
	rm -f *.o $(TARGETS)
