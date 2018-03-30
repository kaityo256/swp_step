CC=g++
CPPFLAGS=-O3 -std=c++11 -march=native

all: a.out

a.out: force.cpp conf.hpp
	$(CC) $(CPPFLAGS) $< -o $@

test: a.out
	./a.out > test.txt
	diff test.txt orig.txt

clean:
	rm -f a.out

clear: clean
	rm -f pair.dat
