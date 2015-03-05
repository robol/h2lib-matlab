-include Makefile.inc

CFLAGS= -g -IH2Lib/Library -D__MATLAB_TRICK
LDFLAGS= -LH2Lib
LIBS= -lh2

MATLAB_MEX=hmatrix_create.mexa64 mvm_hmatrix_avector.mexa64

all: $(MATLAB_MEX)

%.mexa64: %.c
	$(MEX) -o $@ $(CFLAGS) $< $(LDFLAGS) $(LIBS) 

clean:
	rm -f *.mex *.mexa64
