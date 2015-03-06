-include Makefile.inc

CFLAGS= -g -IH2Lib/Library 
LDFLAGS= -L. -LH2Lib -Wl,-RH2Lib -Wl,-R.
LIBS= -lh2

MATLAB_MEX=hmatrix_create.mexa64 \
    mvm_hmatrix_avector.mexa64   \
    hmatrix_tridiag.mexa64       \
    hmatrix_sum.mexa64	         \
    hmatrix_prod.mexa64	         \
    hmatrix_inv.mexa64		\
    hmatrix_size.mexa64

all: $(MATLAB_MEX)

%.mexa64: %.c libhmatlab.so
	$(MEX) -o $@ $(CFLAGS) -D__MATLAB_TRICK $< $(LDFLAGS) $(LIBS) -lhmatlab > /dev/null

libhmatlab.so: hmatlab.h hmatlab.c
	$(CC) -o libhmatlab.so -shared -fPIC $(CFLAGS) hmatlab.c $(LDFLAGS) $(LIBS)

clean:
	rm -f *.mex *.mexa64 libhmatlab.so
