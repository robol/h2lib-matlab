-include ../../Makefile.inc

MEX_FILES = \
    delete_h2matrix$(MEX_EXT) 		\
    mvm_h2matrix_avector$(MEX_EXT)   	\
    tridiag$(MEX_EXT)       		\
    h2matrix_sum$(MEX_EXT)		\
    h2matrix_minus$(MEX_EXT) 		\
    h2matrix_prod$(MEX_EXT)             \
    h2matrix_mldivide$(MEX_EXT)

all: $(MEX_FILES)

%$(MEX_EXT): %.c
	$(MEX) $< CFLAGS="$(CFLAGS) $(MATLAB_CFLAGS)" LDFLAGS="$(LDFLAGS) $(LIBS) $(MEX_LIBS)" > /dev/null

clean:
	rm -f *$(MEX_EXT)
