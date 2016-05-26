-include ../Makefile.inc

MEX_FILES = \
    matrix_size$(MEX_EXT) \
    hmatrix_rank$(MEX_EXT) \
    hmatrix_trace$(MEX_EXT) \
    norm$(MEX_EXT)

all: $(MEX_FILES) _private

_private:
	$(MAKE) -C private

%$(MEX_EXT): %.c
	$(MEX) -o $@ $(CFLAGS) $(MATLAB_CFLAGS) $< $(LDFLAGS) $(LIBS) $(MEX_LIBS) > /dev/null

clean:
	rm -f *$(MEX_EXT)
	$(MAKE) -C private clean
