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
    hmatrix_mldivide$(MEX_EXT)         \
    hmatrix_mldivide_full$(MEX_EXT)    \
    hmatrix_scale$(MEX_EXT)            \
    hmatrix_shift$(MEX_EXT)            \
    band$(MEX_EXT)                     \
    generators$(MEX_EXT) \
    hmatrix_from_full$(MEX_EXT)

all: $(MEX_FILES)

%$(MEX_EXT): %.c
	$(MEX) $< $(LIBS) $(MEX_LIBS) CFLAGS="$(CFLAGS) $(MATLAB_CFLAGS)" LDFLAGS="$(LDFLAGS)" 

clean:
	rm -f *$(MEX_EXT)
