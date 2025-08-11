CXX=g++
CFLAGS=-O3 -Wall -std=c++20 -shared -undefined dynamic_lookup -Xclang -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp
# CFLAGS=-O3 -Wall -std=c++20 -shared -fPIC -fopenmp

simple_add: simple_add.cpp
	$(CXX) $(CFLAGS) $$(python -m pybind11 --includes) simple_add.cpp -o simple_add$$(python -m pybind11 --extension-suffix)

dist_py: dist.o dist_py.cpp
	$(CXX) $(CFLAGS) $$(python -m pybind11 --includes) dist_py.cpp dist.o -o dist_py$$(python -m pybind11 --extension-suffix)

dist.o: dist.cpp dist.h
	$(CXX) -c $(CFLAGS) dist.cpp -o dist.o $$(python -m pybind11 --includes)

dist: dist.cpp dist.h
	$(CXX) $(CFLAGS) dist.cpp -o dist $$(python -m pybind11 --includes)

testopenmp: testopenmp.cpp
	$(CXX) $(CFLAGS) testopenmp.cpp -o testopenmp
