-include ../Makefile.inc

all: libhmatlab.so

libhmatlab.so: hmatlab.c hmatlab.h
	$(CC) -o $@ $(CFLAGS) -shared -fPIC $< $(LDFLAGS) $(LIBS)

clean:
	rm -f libhmatlab.so
