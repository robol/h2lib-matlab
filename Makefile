-include Makefile.inc

ifdef USE_COMPLEX
  CFLAGS += -DUSE_COMPLEX
endif

all: Makefile.inc hmatlab hmatrix cluster h2matrix

Makefile.inc:
	$(error "Makefile.inc not present, please run ./configure before make")

hmatlab:
	$(MAKE) -C lib

hmatrix: hmatlab
	$(MAKE) -C @HMatrix

h2matrix: hmatlab
	$(MAKE) -C @H2Matrix

cluster: hmatlab
	$(MAKE) -C @Cluster

check:
	$(MAKE) -C tests check

clean:
	$(MAKE) -C lib clean
	$(MAKE) -C @HMatrix clean
	$(MAKE) -C @Cluster clean

distclean: clean
	rm -f Makefile.inc
