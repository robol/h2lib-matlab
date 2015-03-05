-include Makefile.inc

CFLAGS= -IH2Lib
LDFLAGS= -LH2Lib
LIBS= -lh2

MATLAB_MEX=hmatrix_create.mexa64

all: $(MATLAB_MEX)

%.mexa64: %.c
	$(MEX) -o $@ $(CFLAGS) $< $(LDFLAGS) $(LIBS) 

clean:
	rm -f *.mex *.mexa64
