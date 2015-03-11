-include ../../Makefile.inc

MEX_FILES = \
	delete_cluster$(MEX_EXT) \
	std_subdivision_scheme$(MEX_EXT)

all: $(MEX_FILES)

%$(MEX_EXT): %.c
	$(MEX) -o $@ $(CFLAGS) $(MATLAB_CFLAGS) $< $(LDFLAGS) $(LIBS) $(MEX_LIBS) > /dev/null

clean:
	rm -f *$(MEX_EXT)
