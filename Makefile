SHELL = /bin/bash
SHELLFLAGS="-O extglob -c"

CXX = g++
INCLUDE = include/
BIN = bin/
VPATH = src/
CXXFLAG = -std=c++11 -lboost_program_options -I ${INCLUDE}

all: src/Main.cpp $(objects)
	${CXX} ${CXXFLAG} -o ${BIN}/superimposer.bin $< ${objects}

src/Align_cliques.o: src/Align_cliques.cpp
	$(CXX) $(CXXFLAG) -c -o $@ $<
src/AlignmentFunctions.o: src/AlignmentFunctions.cpp 
	$(CXX) $(CXXFLAG) -c -o $@ $<
src/CliqParser.o: src/CliqParser.cpp
	$(CXX) $(CXXFLAG) -c -o $@ $<
src/jacobi_eigenvalue.o: src/jacobi_eigenvalue.cpp
	$(CXX) $(CXXFLAG) -c -o $@ $<

objects : src/Align_cliques.o  src/AlignmentFunctions.o  src/CliqParser.o  src/jacobi_eigenvalue.o
objects = src/Align_cliques.o  src/AlignmentFunctions.o  src/CliqParser.o  src/jacobi_eigenvalue.o
