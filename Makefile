# Makefile for project TAMASIS
# Author: P. Chanial


# set up default compiler
ifeq "$(origin FC)" "default"
    ifneq ($(shell which ifort),)
        FC=ifort
    else ifneq ($(shell which gfortran),)
        FC=gfortran
    endif
endif

# set up compiler flags
ifeq "$(FC)" "gfortran"
    FFLAGS_DEBUG = -g -fbacktrace -Warray-temporaries -O0 -fcheck=all -ffree-form -fopenmp -Wall -fPIC -cpp -DGFORTRAN
    FFLAGS_RELEASE = -fbacktrace -O3 -ffree-form -fopenmp -Wall -fPIC -cpp -DGFORTRAN
    LDFLAGS = -lgomp $(shell pkg-config --libs cfitsio) $(shell pkg-config --libs wcslib) $(shell pkg-config --libs fftw3)
    FCOMPILER=gnu95
else ifeq ($(FC),ifort)
    FFLAGS_DEBUG = -debug -fpp -O0 -static -fPIC -free -openmp -ftz -traceback -DIFORT -check all -ftrapuv
    FFLAGS_RELEASE = -fpp -fast -fPIC -free -openmp -ftz -DIFORT
    # for static linking: use -static -static-intel
    LDFLAGS = -liomp5 $(shell pkg-config --libs cfitsio) $(shell pkg-config --libs wcslib) $(shell pkg-config --libs fftw3)
    FCOMPILER = intelem
else
    $(error Unsupported compiler '$(FC)'.)
endif

ifeq ($(PROF_GEN),1)
    DEBUG = 1
endif

ifeq ($(DEBUG),1)
    FFLAGS = $(FFLAGS_DEBUG) -DDEBUG
else
    FFLAGS = $(FFLAGS_RELEASE)
endif

ifeq ($(PROF_GEN),1)
    FFLAGS += -prof_gen -prof_dir/home/pchanial/profiles
endif

ifeq ($(PROF_USE),1)
    FFLAGS += -prof_use -prof_dir/home/pchanial/profiles
endif

INCLUDES = -Iinclude -Iinclude/wcslib-4.4.4-Fortran90

MODULES = $(wildcard module_*.f)
SOURCES = $(wildcard test_*.f) pacs_photproject.f
EXECS = $(SOURCES:.f=)
FORTRANTESTS = $(wildcard test_*.f)
PYTHONTESTS = $(wildcard tamasis/tests/test_*.py)

# apply a function to each element of a list
map = $(foreach a,$(2),$(call $(1),$(a)))

# recursively find dependencies
finddeps = $(1).o $(if $($(1)),$(call map,finddeps,$($(1))))

# define module dependencies
module_cfitsio = module_stdio
module_deglitching = module_math module_pointingmatrix module_precision module_string
module_filtering = module_precision
module_fitstools = module_cfitsio module_precision module_string
module_instrument = module_precision module_string
module_madcap = module_filtering module_pointingmatrix module_precision module_string
module_math = module_precision
module_observation = module_fitstools module_math module_precision module_string
module_optionparser = module_string
module_pacsinstrument = module_fitstools module_math module_pacsobservation module_pointingmatrix module_projection module_tamasis module_wcs 
module_pacsobservation = module_fitstools module_math module_observation module_precision module_string
module_pointingmatrix = module_math module_precision module_projection
module_preprocessor = module_math module_sort
module_projection = module_precision module_sort module_stack
module_sort = module_math
module_wcs = module_fitstools module_math module_string module_wcslib

# define executable dependencies
pacs_photproject = module_fitstools module_deglitching module_optionparser module_pacsinstrument module_pacsobservation module_pointingmatrix module_preprocessor module_string module_tamasis
test_cfitsio = module_cfitsio
test_compression = module_compression module_math module_precision
test_deglitching = module_deglitching module_math module_pointingmatrix module_precision
test_fitstools = module_fitstools
test_madcap = module_filtering module_fitstools module_madcap
test_math = module_math module_precision
test_ngc6946_bpj = module_fitstools module_math module_pacsinstrument module_pacsobservation module_pointingmatrix module_preprocessor module_projection module_tamasis
test_optionparser = module_optionparser
test_pacsinstrument = module_pacsinstrument module_pacsobservation module_tamasis
test_pacsobservation = module_pacsinstrument  module_pacsobservation module_string module_tamasis
test_pacspointing = module_math module_pacsobservation module_precision
test_pointingmatrix = module_pointingmatrix
test_preprocessor = module_fitstools module_math module_preprocessor
test_projection = module_projection module_sort
test_read_config = module_instrument
test_sort = module_fitstools module_math module_sort
test_stack = module_stack
test_string = module_string
test_module_string = module_string
test_stdio = module_cfitsio module_stdio
test_wcs = module_fitstools module_math module_wcs
test_wcslib1 = module_cfitsio module_wcslib 
test_wcslib2 = module_fitstools module_math module_wcslib
test_wcslibc = module_cfitsio module_wcslibc

.PHONY : all tests
all : $(EXECS) tamasis/tamasisfortran.so
#all : $(EXECS)

# if %.mod doesn't exist, make %.o. It will create %.mod with the same 
# timestamp. If it does, do nothing
%.mod : %.o
	@if [ ! -f $@ ]; then \
	    rm $< ;\
	    $(MAKE) $< ;\
	fi

%.o : %.f
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $<

%: %.o
	$(FC) -o $@ $^ $(LDFLAGS)

.SECONDEXPANSION:
$(MODULES:.f=.o) $(SOURCES:.f=.o):%.o: $$(addsuffix .mod,$$($$*))
$(EXECS):%:$$(sort $$(call finddeps,$$*))

tamasis/tamasisfortran.so: tamasisfortran.f90 $(MODULES:.f=.o)
	unset LDFLAGS ; \
	f2py --fcompiler=${FCOMPILER} --f90exec=$(FC) --f90flags="$(FFLAGS)" -DF2PY_REPORT_ON_ARRAY_COPY=1 -c $^ -m tamasisfortran $(LDFLAGS) ; \
	mv -f tamasisfortran.so tamasis/
# /opt/core-3.1-amd64/ifc/11.1/lib/intel64/libiomp5.a /opt/core-3.1-amd64/ifc/11.1/lib/intel64/libifport.a /opt/core-3.1-amd64/ifc/11.1/lib/intel64/libifcore.a /opt/core-3.1-amd64/ifc/11.1/lib/intel64/libimf.a /opt/core-3.1-amd64/ifc/11.1/lib/intel64/libsvml.a

clean:
	@rm -f *.o *.mod *.so *~ $(EXECS);\
	find . -name "*pyc" -exec rm '{}' ';'

tests: tests_fortran tests_python

tests_fortran: $$(filter-out test_wcslib%,$(FORTRANTESTS:.f=))
	@for test in $^; do \
	echo;\
	echo "Running Fortran test: "$$test"...";\
	echo "=====================";\
	./$$test; \
	done

tests_python:
	@for test in $(PYTHONTESTS); do \
	echo;\
	echo "Running Python test: "$$test"...";\
	echo "====================";\
	python $$test; \
	done

