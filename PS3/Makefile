# Makefile for Problem Set 3: Extended Huckel Program

# Compiler
CPP = g++

# Compiler flags
CPPFLAGS = -O2 -std=c++11 

factorial: factorial.h
	$(CPP) $(CPPFLAGS) -c factorial.cpp

molecule: molecule.h
	$(CPP) $(CPPFLAGS) -c molecule.cpp 

eht_matrices: eht_matrices.h
	$(CPP) $(CPPFLAGS) -c eht_matrices.cpp

test: factorial molecule eht_matrices
	$(CPP) $(CPPFLAGS) -o test test.cpp factorial.o molecule.o eht_matrices.o -larmadillo

energy_diff: factorial molecule eht_matrices
	$(CPP) $(CPPFLAGS) -o energy_diff energy_diff.cpp factorial.o molecule.o eht_matrices.o -larmadillo

all: test energy_diff

clean:
	rm -f *.o test energy_diff

	