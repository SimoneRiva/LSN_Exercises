CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex05.exe : random.o error.o ex05.o
	$(CC) -o $@ $^

ex05.o : ex05.cpp ex05.h
	$(CC) -c ex05.cpp $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp $(CFLAGS)

error.o : error.cpp error.h
	$(CC) -c error.cpp $(CFLAGS)

clean :
	rm *.o ex05.exe
