-include ../Makefile.inc

all: libhmatlab.so

libhmatlab.so: hmatlab.c hmatlab.h h2matlab.c
	$(CC) -o $@ $(CFLAGS) -shared -fPIC hmatlab.c h2matlab.c $(LDFLAGS) $(LIBS)

clean:
	rm -f libhmatlab.so
