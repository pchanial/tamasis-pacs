FFLAGS_MOD = -module
FFLAGS_DEBUG = -debug -fpp -O0 -static -fPIC -free -openmp -ftz -traceback -DIFORT -check all -ftrapuv -fp-model precise
FFLAGS_RELEASE = -fpp -fast -fPIC -free -openmp -ftz -DIFORT -fp-model precise
FFLAGS_F95 = -fast -fPIC -ftz -fp-model precise
# for static linking: use -static -static-intel
LDFLAGS = $(filter-out -lm,-liomp5 $(shell pkg-config --libs cfitsio) $(shell pkg-config --libs wcslib) $(shell pkg-config --libs fftw3) -latlas -L/usr/lib/atlas)
FCOMPILER = intelem

ifeq ($(PROF_GEN),1)
    FFLAGS_PROF = -prof_gen -prof_dir/home/pchanial/profiles
endif

ifeq ($(PROF_USE),1)
    FFLAGS_PROF = -prof_use -prof_dir/home/pchanial/profiles
endif

FCVERSION = $(lastword $(shell ifort -v))
SUPPORT_QUAD = 1