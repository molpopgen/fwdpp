CXX=c++
DEBUG=-DNDEBUG
#DEBUG=-g
#DEBUG=
OPT=-O2  -funroll-all-loops
#OPT=-O2 -pg
#OPT=
CXXFLAGS=$(OPT) -Wall -W -I..  $(DEBUG) -std=c++11

TARGETS=diploid diploid_ind diploid_fixed_sh diploid_twopop_mig diploid_binaryIO RHH bneck_selection migsel_ind diploid_binaryIO_ind diploid_fixed_sh_ind TFL2013 diploid_ind_2locus
all: diploid.o diploid_ind.o diploid_binaryIO.o diploid_binaryIO_ind.o diploid_fixed_sh.o diploid_twopop_mig.o RHH.o diploid_fixed_sh_ind.o bneck_selection.o migsel_ind.o TFL2013.o diploid_ind_2locus.o
TARGETS=diploid diploid_ind diploid_fixed_sh diploid_twopop_mig diploid_binaryIO diploid_gzbinaryIO RHH bneck_selection migsel_ind diploid_binaryIO_ind diploid_gzbinaryIO_ind diploid_fixed_sh_ind TFL2013 pfix
all: diploid.o diploid_ind.o diploid_binaryIO.o diploid_gzbinaryIO.o diploid_binaryIO_ind.o diploid_gzbinaryIO_ind.o diploid_fixed_sh.o diploid_twopop_mig.o RHH.o diploid_fixed_sh_ind.o bneck_selection.o migsel_ind.o TFL2013.o pfix.o
	$(CXX) $(CXXFLAGS) -o diploid diploid.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_ind diploid_ind.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_ind_2locus diploid_ind_2locus.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o RHH RHH.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_binaryIO diploid_binaryIO.o -lsequence -lz -lgsl -lboost_system -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_gzbinaryIO diploid_gzbinaryIO.o -lsequence -lz -lgsl -lboost_system -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_binaryIO_ind diploid_binaryIO_ind.o -lsequence -lz -lgsl -lboost_system -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_gzbinaryIO_ind diploid_gzbinaryIO_ind.o -lsequence -lz -lgsl -lboost_system -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_fixed_sh diploid_fixed_sh.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_fixed_sh_ind diploid_fixed_sh_ind.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o diploid_twopop_mig diploid_twopop_mig.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o bneck_selection bneck_selection.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o migsel_ind migsel_ind.o -lsequence -lz -lgsl -lgslcblas -lboost_system $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o TFL2013 TFL2013.o -lsequence -lz -lgsl -lgslcblas -lboost_system -lboost_iostreams $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o pfix pfix.o -lgsl -lgslcblas $(LDFLAGS)

clean:
	rm -f *.o $(TARGETS)
