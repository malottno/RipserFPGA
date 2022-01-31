build: ripser


all: output ripser ripser-coeff ripser-debug

output: main.o ripser.o
	c++ -std=c++11 -Wall main.o ripser.o -o output
	
main.o: main.cpp
	c++ -c main.cpp
	
ripser.o: ripser.cpp ripser.hpp
	c++ -c ripser.cpp


clean:
	rm -f ripser ripser-coeff ripser-debug
