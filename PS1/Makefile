# Makefile for Problem Set 1 

# Compiler
CPP = g++

# Compiler flags
CPPFLAGS = -O2 -std=c++11 


gold_lj: gold_lj.h
	$(CPP) $(CPPFLAGS) -c gold_lj.cpp -larmadillo

test: gold_lj
	$(CPP) $(CPPFLAGS) -o test test.cpp gold_lj.o -larmadillo

all: gold_lj test

clean:
	rm -f *.o test