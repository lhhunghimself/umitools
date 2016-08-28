CC=gcc
CPP=g++
LD=g++ 
ifdef DEBUG
CFLAGS= -I./ -I./optparse -ggdb -I -Wall -std=c++11
else
CFLAGS= -I./ -I./optparse -Ofast -std=c++11 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-functions -m64 -mavx -msse3 -msse4 -frerun-loop-opt -static $(HDRS) $(DEFINES)
endif
LIBS= -lm -lc -lgcc -lrt -lz

.PHONY: clean all
all: speedTest umiappend

OBJ_speedTest = speedTest.o optparse/optparse.o
OBJ_umiappend = umiappend.o optparse/optparse.o

umiappend : $(OBJ_umiappend)
	$(LD) -o umiappend $(LFLAGS) $(OBJ_umiappend) $(LIBS)
umiappend.o : umiappend.cpp 
	$(CPP) $(CFLAGS) -c umiappend.cpp
speedTest : $(OBJ_speedTest)
	$(LD) -o speedTest $(LFLAGS) $(OBJ_speedTest) $(LIBS)
speedTest.o : speedTest.cpp 
	$(CPP) $(CFLAGS) -c speedTest.cpp

optparse.o :
	cd optparse && make

clean:
	-@rm -rf *.o 2>/dev/null || true
	-@rm -rf core.* 2>/dev/null || true
default: umiappend
