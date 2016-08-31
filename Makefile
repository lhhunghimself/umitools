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
all: umimerge speedTest umiappend umisample

OBJ_umimerge = umimerge.o
OBJ_speedTest = speedTest.o optparse/optparse.o
OBJ_umiappend = umiappend.o optparse/optparse.o
OBJ_umisample = umisample.o optparse/optparse.o

umimerge : $(OBJ_umimerge)
	$(LD) -o umimerge $(LFLAGS) $(OBJ_umimerge) $(LIBS)
umimerge.o : umimerge.cpp 
	$(CPP) $(CFLAGS) -c umimerge.cpp	
		
umiappend : $(OBJ_umiappend)
	$(LD) -o umiappend $(LFLAGS) $(OBJ_umiappend) $(LIBS)
umisample : $(OBJ_umisample)
	$(LD) -o umisample $(LFLAGS) $(OBJ_umisample) $(LIBS)
umiappend.o : umiappend.cpp 
	$(CPP) $(CFLAGS) -c umiappend.cpp	
umisample.o : umisample.cpp 
	$(CPP) $(CFLAGS) -c umisample.cpp
speedTest : $(OBJ_speedTest)
	$(LD) -o speedTest $(LFLAGS) $(OBJ_speedTest) $(LIBS)
speedTest.o : speedTest.cpp 
	$(CPP) $(CFLAGS) -c speedTest.cpp

optparse.o :
	cd optparse && make

clean:
	-@rm -rf *.o 2>/dev/null || true
	-@rm -rf core.* 2>/dev/null || true
default: all
