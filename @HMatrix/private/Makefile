-include ../../Makefile.inc

MEX_FILES = \
    allocate_hmatrix$(MEX_EXT) \
    hmatrix_create$(MEX_EXT) 		\
    mvm_hmatrix_avector$(MEX_EXT)   	\
    tridiag$(MEX_EXT)       	\
    hmatrix_sum$(MEX_EXT)	       	\
    hmatrix_prod$(MEX_EXT)	        \
    hmatrix_inv$(MEX_EXT)		\
    delete_hmatrix$(MEX_EXT) 		\
    hmatrix_minus$(MEX_EXT)	        \
    hmatrix_uminus$(MEX_EXT)	       \
    hmatrix_mldivide$(MEX_EXT)          \
    hmatrix_mldivide_full$(MEX_EXT)    \
    band$(MEX_EXT)                      \
    generators$(MEX_EXT)

all: $(MEX_FILES)

%$(MEX_EXT): %.c
	$(MEX) -o $@ $(CFLAGS) $(MATLAB_CFLAGS) $< $(LDFLAGS) $(LIBS) $(MEX_LIBS) > /dev/null

clean:
	rm -f *$(MEX_EXT)
