-include ../Makefile.inc

MEX_FILES = \
	matrix_size$(MEX_EXT) \
	h2matrix_rank$(MEX_EXT)

all: $(MEX_FILES) _private

_private:
	$(MAKE) -C private

%$(MEX_EXT): %.c
	$(MEX) $< $(LIBS) $(MEX_LIBS) CFLAGS="$(CFLAGS) $(MATLAB_CFLAGS)" LDFLAGS="$(LDFLAGS)" > /dev/null

clean:
	rm -f *$(MEX_EXT)
	$(MAKE) -C private clean
