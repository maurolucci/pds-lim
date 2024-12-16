# Makefile

CC = g++ -std=c++20

SRCPATH = src
INCPATH = include
DEPPATH = _deps

LIBS   = -L/opt/gurobi1103/linux64/lib -lgurobi_c++ -lgurobi110
DEPS = $(DEPPATH)/tinyxml2/CMakeFiles/tinyxml2.dir/tinyxml2.cpp.o $(DEPPATH)/fmt/CMakeFiles/fmt.dir/src/format.cc.o
INCS = -I$(INCPATH) -I$(DEPPATH)/fmt/include -I$(DEPPATH)/tinyxml2 -I/opt/gurobi1103/linux64/include

all: deps pds-lim

pds-lim: main.cpp pds.o graphio.o gurobi_solve.o $(DEPS)
	$(CC) -o $@.out $^ $(INCS) $(LIBS)

deps:
	cd $(DEPPATH)/fmt && make
	cd $(DEPPATH)/tinyxml2 && make

pds.o: $(SRCPATH)/pds.cpp $(INCPATH)/pds.hpp
	$(CC) -g -c $< $(INCS)

graphio.o: $(SRCPATH)/graphio.cpp $(INCPATH)/graphio.hpp $(INCPATH)/pds.hpp
	$(CC) -g -c $< $(INCS)

gurobi_solve.o: $(SRCPATH)/gurobi_solve.cpp $(INCPATH)/gurobi_solve.hpp $(INCPATH)/gurobi_common.hpp
	$(CC) -g -c $< $(INCS)