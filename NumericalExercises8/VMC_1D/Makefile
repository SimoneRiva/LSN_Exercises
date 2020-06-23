CC = g++
CFLAGS = -Wall -O3 --std=c++11

VMC_1D.exe : random.o error.o VMC_1D.o
	$(CC) -o $@ $^

VMC_1D.o : VMC_1D.cpp VMC_1D.h funzioneBase.h
	$(CC) -c VMC_1D.cpp $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp $(CFLAGS)

error.o : error.cpp error.h
	$(CC) -c error.cpp $(CFLAGS)

clean :
	rm *.o VMC_1D.exe 
