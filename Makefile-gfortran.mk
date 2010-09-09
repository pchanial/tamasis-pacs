FFLAGS_MOD = -J
FFLAGS_DEBUG = -g -fbacktrace -O3 -fcheck=all -ffree-form -fopenmp -Wall -fPIC -cpp -DGFORTRAN
FFLAGS_RELEASE = -fbacktrace -O3 -ffree-form -fopenmp -Wall -fPIC -cpp -DGFORTRAN
FFLAGS_F95 = -O3 -fPIC -std=f95
LDFLAGS = -lgomp $(shell pkg-config --libs cfitsio) $(shell pkg-config --libs wcslib) $(shell pkg-config --libs fftw3) -latlas -L/usr/lib/atlas
FCOMPILER=gnu95
