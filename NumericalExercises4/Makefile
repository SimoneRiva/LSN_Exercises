CC = g++
CFLAGS = -std=c++0x

MolDyn_NVE.exe : MolDyn_NVE.o error.o
	$(CC) $(CFLAGS) -o $@ $^

MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) $(CFLAGS) -c MolDyn_NVE.cpp

error.o : error.cpp error.h
	$(CC) $(CFLAGS) -c error.cpp

copy: copy_config copy_old

copy_config: config.final
	cp config.final config.0

copy_old: old.final
	cp old.final old.0
