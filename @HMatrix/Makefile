-include ../Makefile.inc

MEX_FILES = \
    size$(MEX_EXT) 

all: $(MEX_FILES) _private

_private:
	$(MAKE) -C private

%$(MEX_EXT): %.c
	$(MEX) -o $@ $(CFLAGS) $(MATLAB_CFLAGS) $< $(LDFLAGS) $(LIBS) -lhmatlab > /dev/null

clean:
	rm -f *$(MEX_EXT)
	$(MAKE) -C private clean
