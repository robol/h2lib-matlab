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
    hmatrix_rank$(MEX_EXT)              \
    hmatrix_minus$(MEX_EXT)

all: $(MEX_FILES)

%$(MEX_EXT): %.c
	$(MEX) -o $@ $(CFLAGS) $(MATLAB_CFLAGS) $< $(LDFLAGS) $(LIBS) -lhmatlab > /dev/null

clean:
	rm -f *$(MEX_EXT)
