# Makefile for Problem Set 2: Analytical and Numerical Integration for Overlap Integrals

# Compiler
CPP = g++

# Compiler flags
CPPFLAGS = -O2 -std=c++11 

factorial: factorial.h
	$(CPP) $(CPPFLAGS) -c factorial.cpp

num_integration: num_integration.h
	$(CPP) $(CPPFLAGS) -c num_integration.cpp 

analytical_integration: analytical_integration.h
	$(CPP) $(CPPFLAGS) -c analytical_integration.cpp

test_numerical: num_integration
	$(CPP) $(CPPFLAGS) -o test_numerical test_numerical.cpp num_integration.o -larmadillo

test_analytical: analytical_integration
	$(CPP) $(CPPFLAGS) -o test_analytical test_analytical.cpp analytical_integration.o factorial.o -larmadillo

all: factorial num_integration analytical_integration test_numerical test_analytical

clean:
	rm -f *.o test_numerical test_analytical


# Citation:
# https://kcwong-joe.medium.com/passing-a-function-as-a-parameter-in-c-a132e69669f6
# "std::bind can be used if we would like to pass a function from one object to 
# another object in the std::function parameter"