# Makefile

CC = g++ -std=c++20

SRCPATH = src
INCPATH = include
DEPPATH = _deps

LIBS   = -L/opt/gurobi1203/linux64/lib -lgurobi_c++ -lgurobi120 -L/usr/lib/x86_64-linux-gnu -lboost_program_options
DEPS = $(DEPPATH)/tinyxml2/CMakeFiles/tinyxml2.dir/tinyxml2.cpp.o $(DEPPATH)/fmt/CMakeFiles/fmt.dir/src/format.cc.o
INCS = -I$(INCPATH) -I$(DEPPATH)/fmt/include -I$(DEPPATH)/tinyxml2 -I/opt/gurobi1203/linux64/include

FLAGS = -Wall -Wextra -O3

all: deps pds-lim

pds-lim: main.cpp pds.o graphio.o gurobi_solve.o efps_solve.o fps_solve.o fort_solve.o $(DEPS)
	$(CC) -o $@ $^ $(INCS) $(LIBS) $(FLAGS)

deps:
	cd $(DEPPATH)/fmt && make
	cd $(DEPPATH)/tinyxml2 && make

pds.o: $(SRCPATH)/pds.cpp $(INCPATH)/pds.hpp
	$(CC) -g -c $< $(INCS) $(FLAGS)

graphio.o: $(SRCPATH)/graphio.cpp $(INCPATH)/graphio.hpp $(INCPATH)/pds.hpp
	$(CC) -g -c $< $(INCS) $(FLAGS)

gurobi_solve.o: $(SRCPATH)/gurobi_solve.cpp $(INCPATH)/gurobi_solve.hpp $(INCPATH)/gurobi_common.hpp $(INCPATH)/pds.hpp
	$(CC) -g -c $< $(INCS) $(FLAGS)

efps_solve.o: $(SRCPATH)/efps_solve.cpp $(INCPATH)/efps_solve.hpp $(INCPATH)/gurobi_common.hpp $(INCPATH)/pds.hpp
	$(CC) -g -c $< $(INCS) $(FLAGS)

fps_solve.o: $(SRCPATH)/fps_solve.cpp $(INCPATH)/fps_solve.hpp $(INCPATH)/gurobi_common.hpp $(INCPATH)/pds.hpp
	$(CC) -g -c $< $(INCS) $(FLAGS)

fort_solve.o: $(SRCPATH)/fort_solve.cpp $(INCPATH)/fort_solve.hpp $(INCPATH)/gurobi_common.hpp $(INCPATH)/pds.hpp
	$(CC) -g -c $< $(INCS) $(FLAGS)

