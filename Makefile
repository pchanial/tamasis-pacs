# Makefile for project TAMASIS
# Author: P. Chanial


ifeq "$(origin FC)" "default"
    ifneq ($(shell which ifort),)
        FC=ifort
    else ifneq ($(shell which gfortran),)
        FC=gfortran
    endif
endif

ifeq "$(FC)" "gfortran"
    FFLAGS_DEBUG = -g -fbacktrace -Warray-temporaries -O3 -fcheck=all -fopenmp -Wall -fPIC
    FFLAGS_RELEASE = -fbacktrace -O3 -fopenmp -Wall -fPIC
    LDFLAGS  = -lgomp $(shell pkg-config --libs cfitsio) $(shell pkg-config --libs wcslib)
    FCOMPILER=gnu95
else ifeq ($(FC),ifort)
    FFLAGS_DEBUG = -fpp -O2 -static -fPIC -openmp -traceback
    FFLAGS_RELEASE = -fpp -fast -openmp -ftz -ip  -ipo
    LDFLAGS  = -liomp5 $(shell pkg-config --libs cfitsio) $(shell pkg-config --libs wcslib)
    FCOMPILER = intelem
else
    $(error Unsupported compiler '$(FC)'.)
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

INCLUDES = wcslib-4.4.4-Fortran90

MODULES = $(wildcard module_*.f90)
SOURCES = $(wildcard test_*.f90) pacs_photproject.f90
EXECS = $(SOURCES:.f90=)

# apply a function to each element of a list
map = $(foreach a,$(2),$(call $(1),$(a)))

# recursively find dependencies
finddeps = $(1).o $(if $($(1)),$(call map,finddeps,$($(1))))

# define module dependencies
module_cfitsio = module_stdio
module_deglitching = module_math module_pointingmatrix module_precision
module_fitstools = module_cfitsio module_precision module_string
module_instrument = module_precision module_string
module_math = module_precision
module_optionparser = module_string
module_pacsinstrument = module_fitstools module_math module_pacsobservation module_pacspointing module_pointingmatrix module_projection module_wcs 
module_pacsobservation = module_fitstools module_precision module_string
module_pacspointing = module_fitstools module_pacsobservation module_precision module_string
module_pointingmatrix = module_math module_precision module_projection
module_preprocessor = module_math
module_projection = module_precision module_sort module_stack
module_wcs = module_fitstools module_math module_string module_wcslib

# define executable dependencies
pacs_photproject = module_fitstools module_deglitching module_optionparser module_pacsinstrument module_pacsobservation module_pacspointing module_pointingmatrix module_preprocessor
test_cfitsio = module_cfitsio
test_deglitching = module_deglitching module_math module_pointingmatrix module_precision
test_fitstools = module_fitstools
test_math = module_math module_precision
test_ngc6946_bpj = module_fitstools module_math module_pacsinstrument module_pacsobservation module_pacspointing module_pointingmatrix module_preprocessor module_projection
test_optionparser = module_optionparser
test_pacsinstrument = module_pacsinstrument module_pacspointing
test_pacsobservation = module_pacsinstrument  module_pacsobservation module_string
test_pacspointing = module_math module_pacsobservation module_pacspointing
test_pointingmatrix = module_pointingmatrix
test_projection = module_projection module_sort
test_read_config = module_instrument
test_sort = module_sort
test_stack = module_stack
test_string = module_string
test_module_string = module_string
test_stdio = module_cfitsio module_stdio
test_wcs = module_fitstools module_math module_wcs
test_wcslib1 = module_cfitsio module_wcslib 
test_wcslib2 = module_fitstools module_math module_wcslib
test_wcslibc = module_cfitsio module_wcslibc

.PHONY : all tests
all : $(EXECS) tamasisfortran.so
#all : $(EXECS)

# if %.mod doesn't exist, make %.o. It will create %.mod with the same 
# timestamp. If it does, do nothing
%.mod : %.o
	@if [ ! -f $@ ]; then \
	    rm $< ;\
	    $(MAKE) $< ;\
	fi

%.o : %.f90
	$(FC) $(FFLAGS) -I$(INCLUDES) -c -o $@ $<

%: %.o
	$(FC) -o $@ $^ $(LDFLAGS)

.SECONDEXPANSION:
$(MODULES:.f90=.o) $(SOURCES:.f90=.o):%.o: $$(addsuffix .mod,$$($$*))
$(EXECS):%:$$(sort $$(call finddeps,$$*))

tamasisfortran.so: tamasisfortran.f90 $(MODULES:.f90=.o)
	unset LDFLAGS ; \
	f2py --fcompiler=${FCOMPILER} --f90exec=$(FC) --f90flags="$(FFLAGS)" -DF2PY_REPORT_ON_ARRAY_COPY=1 -c $^ -m tamasisfortran $(LDFLAGS)

clean:
	rm -f *.o *.mod *.so *~ $(EXECS)

tests: $$(filter-out test_wcslib%,$$(filter test_%,$$(EXECS)))
	@for test in $^; do \
	echo;\
	echo "Running test: "$$test"...";\
	echo "=============";\
	./$$test; \
	done
