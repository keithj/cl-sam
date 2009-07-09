
CC = gcc
CFLAGS = -g -Wall -fPIC

LD = gcc
LDFLAGS = -shared -Wl

RM = /bin/rm -f

all: libbgzf.so

libbgzf.so: bgzf.o
	$(LD) $(LDFLAGS) -o libbgzf.so bgzf.o -lc -lz

bgzf.o: bgzf.c
	$(CC) $(CFLAGS) -c bgzf.c

clean:
	$(RM) *.o
