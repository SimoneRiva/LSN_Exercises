CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex03.exe : random.o error.o ex03.o
	$(CC) -o $@ $^

ex03.o : ex03.cpp funzioneBase.h
	$(CC) -c ex03.cpp $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp $(CFLAGS)

error.o : error.cpp error.h
	$(CC) -c error.cpp $(CFLAGS)
