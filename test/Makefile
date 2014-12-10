CXX=c++
#DEBUG=-DNDEBUG
#DEBUG=-g
#DEBUG=
OPT=-I.. -O2  -funroll-all-loops
#OPT=-O2 -pg
#OPT=
CXXFLAGS=$(OPT) -std=c++11 -Wall -W -I..  $(DEBUG)

TARGETS=phiab strobeck_morgan pmono procsm
all:  phiab.o strobeck_morgan.o pmono.o procsm.o
	$(CXX) $(CXXFLAGS) -o strobeck_morgan strobeck_morgan.o -lsequence  $(LDFLAGS) -lz -lgsl -lgslcblas -lboost_system 
	$(CC) $(CCFLAGS) -o phiab phiab.o  $(LDFLAGS) -lm 
	$(CXX) $(CXXFLAGS) -o pmono pmono.o  $(LDFLAGS) -lsequence -lz 
	$(CXX) $(CXXFLAGS) -o procsm procsm.o  $(LDFLAGS) -lsequence -lz

clean:
	rm -f *.o $(TARGETS)
