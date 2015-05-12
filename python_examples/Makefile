CXXFLAGS=-fPIC -Wall -W -O3 -I.. -I. `python-config --includes` -std=c++11 
LIBS=`python-config --ldflags` -lboost_python -lboost_system -lgsl -lgslcblas

all: fwdpp_boost_python.o social_evol.o
	$(CXX) -shared -o fwdpp_boost_python.so fwdpp_boost_python.o $(LIBS)
	$(CXX) -shared -o social_evol.so social_evol.o $(LIBS)

clean:
	rm -f *.o *.so
