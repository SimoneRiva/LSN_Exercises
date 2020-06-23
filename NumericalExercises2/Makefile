CC = g++
CFLAGS = -Wall -O3 --std=c++11

compile : ex02.1.exe ex02.2.exe

ex02.1.exe : random.o error.o ex02.1.o
	$(CC) -o $@ $^

ex02.2.exe : random.o ranWalk.o ex02.2.o
	$(CC) -o $@ $^

ex02.1.o : ex02.1.cpp funzioneBase.h
	$(CC) -c ex02.1.cpp $(CFLAGS)

ex02.2.o : ex02.2.cpp
	$(CC) -c ex02.2.cpp $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp $(CFLAGS)

error.o : error.cpp error.h
	$(CC) -c error.cpp $(CFLAGS)

ranWalk.o : ranWalk.cpp ranWalk.h
	$(CC) -c ranWalk.cpp $(CFLAGS)
