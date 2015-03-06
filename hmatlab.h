#include <avector.h>
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>

#ifndef _HMATLAB_H
#define _HMATLAB_H

pcluster create_cluster (int * a, int n);

phmatrix create_tridiag_hmatrix (double * a, double * b, double * c, pccluster rc, pccluster cc);

#endif
