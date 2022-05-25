all: output

output: main.o
	c++ -std=c++11 -Wall main.cpp -o ripserMain

clean:
	rm -f main.o
