all: output

output: main.o
	c++ -std=c++11 -Wall main.cpp -o output ripser.o

ripser:
	c++ -c -std=c++11 -Wall ripser.cpp -o ripser

clean:
	rm -f main.o
