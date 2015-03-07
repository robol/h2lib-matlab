-include Makefile.inc

all: Makefile.inc hmatlab hmatrix cluster

Makefile.inc:
	$(error "Makefile.inc not present, please run ./configure before make")

hmatlab:
	$(MAKE) -C lib

hmatrix: hmatlab
	$(MAKE) -C @HMatrix

cluster: hmatlab
	$(MAKE) -C @Cluster

check:
	$(MAKE) -C tests check

clean:
	$(MAKE) -C lib clean
	$(MAKE) -C @HMatrix clean
	$(MAKE) -C @Cluster clean

