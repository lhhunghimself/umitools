CC=gcc
LD=gcc 
CFLAGS= -I./ -I./optparse -Ofast -fexpensive-optimizations -ffast-math -finline-functions -frerun-loop-opt -static $(HDRS) $(DEFINES)
LIBS= -lm -lc -lgcc -lrt -lz

.PHONY: clean all

OBJ_umiappend = umiappend.o optparse/optparse.o
OBJ_umiextract = umiextract.o optparse/optparse.o
#make optparse
umiextract : $(OBJ_umiextract)
	$(LD) -o umiextract $(LFLAGS) $(OBJ_umiextract) $(LIBS)
umiextract.o : umiextract.c 
	$(CC) $(CFLAGS) -c umiextract.c
umiappend : $(OBJ_umiappend)
	$(LD) -o umiappend $(LFLAGS) $(OBJ_umiappend) $(LIBS)
umiappend.o : umiappend.c 
	$(CC) $(CFLAGS) -c umiappend.c
optparse.o :
	cd optparse && make
clean:
	-@rm -rf *.o 2>/dev/null || true
	-@rm -rf core.* 2>/dev/null || true
	-@rm -rf 'umiappend' 2>/dev/null || true
default: append
