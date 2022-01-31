all: output

output: main.o
	c++ -std=c++11 -Wall main.cpp ripser.cpp -o output

clean:
	rm -f main.o
