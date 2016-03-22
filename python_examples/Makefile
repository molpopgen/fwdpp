CXXFLAGS=-fPIC -Wall -W  -I.. -I. `python-config --cflags --includes` -std=c++11 
LIBS_BOOST_PYTHON=`python-config --ldflags` -lboost_python -lboost_system -lgsl -lgslcblas
LIBS_PYBIND11=`python-config --ldflags` -lgsl -lgslcblas

all: fwdpp_boost_python.o social_evol.o fwdpp_pybind11.o
	$(CXX) -shared -o fwdpp_boost_python.so fwdpp_boost_python.o $(LIBS_BOOST_PYTHON)
	$(CXX) -shared -o fwdpp_pybind11.so fwdpp_pybind11.o $(LIBS_PYBIND11)
	$(CXX) -shared -o social_evol.so social_evol.o $(LIBS)

clean:
	rm -f *.o *.so
