CC=gcc
CPP=g++
LD=g++ -L /usr/local/lib
ifdef DEBUG
CFLAGS= -I./ -I./optparse -ggdb -I -Wall -fopenmp -std=c++11
else
CFLAGS= -I./ -I./optparse -Ofast -std=c++11 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-functions -fopenmp -m64 -mavx -msse3 -msse4 -frerun-loop-opt -static $(HDRS) $(DEFINES)
endif
LIBS= -lboost_filesystem -lboost_system -lgomp -lm -lc -lgcc -lrt -lz

.PHONY: clean all
all: umimerge_parallel umimerge speedTest umiappend umisample umisplit

OBJ_umisplit = umisplit.o optparse/optparse.o
OBJ_umimerge = umimerge.o
OBJ_speedTest = speedTest.o optparse/optparse.o
OBJ_umiappend = umiappend.o optparse/optparse.o
OBJ_umisample = umisample.o optparse/optparse.o
OBJ_umimerge_parallel = umimerge_parallel.o


umimerge_parallel : $(OBJ_umimerge_parallel)
	$(LD) -o umimerge_parallel $(LFLAGS) $(OBJ_umimerge_parallel) $(LIBS)
umimerge_parallel.o : umimerge_parallel.cpp 
	$(CPP) $(CFLAGS) -c umimerge_parallel.cpp	
umimerge : $(OBJ_umimerge)
	$(LD) -o umimerge $(LFLAGS) $(OBJ_umimerge) $(LIBS)
umimerge.o : umimerge.cpp 
	$(CPP) $(CFLAGS) -c umimerge.cpp	
		
umiappend : $(OBJ_umiappend)
	$(LD) -o umiappend $(LFLAGS) $(OBJ_umiappend) $(LIBS)
umiappend.o : umiappend.cpp 
	$(CPP) $(CFLAGS) -c umiappend.cpp		
umisplit : $(OBJ_umisplit)
	$(LD) -o umisplit $(LFLAGS) $(OBJ_umisplit) $(LIBS)
umisplit.o : umisplit.cpp 
	$(CPP) $(CFLAGS) -c umisplit.cpp			
umisample : $(OBJ_umisample)
	$(LD) -o umisample $(LFLAGS) $(OBJ_umisample) $(LIBS)
umisample.o : umisample.cpp 
	$(CPP) $(CFLAGS) -c umisample.cpp
speedTest : $(OBJ_speedTest)
	$(LD) -o speedTest $(LFLAGS) $(OBJ_speedTest) $(LIBS)
speedTest.o : speedTest.cpp 
	$(CPP) $(CFLAGS) -c speedTest.cpp

optparse.o :
	cd optparse && make

clean:
	-@rm  umimerge_parallel umimerge speedTest umiappend umisample umisplit
	-@rm -rf *.o 2>/dev/null || true
	-@rm -rf core.* 2>/dev/null || true
default: all
