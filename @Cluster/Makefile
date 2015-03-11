-include ../Makefile.inc

MEX_FILES = \
	cluster_size$(MEX_EXT)


all: $(MEX_FILES) private-dir

%$(MEX_EXT): %.c
	$(MEX) -o $@ $(CFLAGS) $(MATLAB_CFLAGS) $< $(LDFLAGS) $(LIBS) $(MEX_LIBS) > /dev/null

private-dir: 
	$(MAKE) -C private

clean:
	rm -f *$(MEX_EXT)
	$(MAKE) -C private clean
