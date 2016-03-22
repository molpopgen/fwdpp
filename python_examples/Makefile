CXXFLAGS=-fPIC -Wall -W  -I.. -I. `python-config --cflags --includes` -std=c++11 
LIBS=`python-config --ldflags` -lboost_python -lboost_system -lgsl -lgslcblas

all: fwdpp_boost_python.o social_evol.o fwdpp_pybind11.o
	$(CXX) -shared -o fwdpp_boost_python.so fwdpp_boost_python.o $(LIBS)
	$(CXX) -shared -o fwdpp_pybind11.so fwdpp_pybind11.o $(LIBS)
	$(CXX) -shared -o social_evol.so social_evol.o $(LIBS)

clean:
	rm -f *.o *.so
