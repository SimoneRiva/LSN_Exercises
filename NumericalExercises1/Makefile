CC = g++
CFLAGS = -Wall -O3 --std=c++11

compile: ex01.1.exe ex01.2.exe ex01.3.exe

ex01.1.exe : random.o error.o ex01.1.o
	$(CC) -o $@ $^

ex01.2.exe : random.o ex01.2.o
	$(CC) -o $@ $^

ex01.3.exe : random.o error.o ex01.3.o
	$(CC) -o $@ $^

ex01.1.o : ex01.1.cpp funzioneBase.h
	$(CC) -c ex01.1.cpp $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp $(CFLAGS)

error.o : error.cpp error.h
	$(CC) -c error.cpp $(CFLAGS)

ex01.2.o : ex01.2.cpp
	$(CC) -c ex01.2.cpp $(CFLAGS)

ex01.3.o : ex01.3.cpp
	$(CC) -c ex01.3.cpp $(CFLAGS)
