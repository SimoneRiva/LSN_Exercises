LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

qmc1d: qmc1d.o
	g++ -O3 -Wall -o $@ $^ ${LIBS}
 
%.o : %.cpp
	g++ -Wall -c $< ${INCS}

clean:
	rm *.o qmc1d potential.dat kinetic.dat probability.dat
