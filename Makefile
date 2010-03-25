# Makefile for project TAMASIS
# Author: P. Chanial

FC=ifort
FFLAGS_DEBUG = -g -fast -warn all -openmp
FFLAGS_RELEASE = -fast -openmp -ftz -ip  -ipo -march=native 
LDFLAGS  = -liomp5 $(shell pkg-config --libs cfitsio) $(shell pkg-config --libs wcslib)
FCOMPILER = intelem

FC=gfortran	
FFLAGS_DEBUG = -g -O3 -fcheck=all -fopenmp -Wall -fPIC
FFLAGS_RELEASE = -O3 -fopenmp -Wall -fPIC
LDFLAGS  = -lgomp $(shell pkg-config --libs cfitsio) $(shell pkg-config --libs wcslib)
FCOMPILER=gnu95

INCLUDES = wcslib-4.4.4-Fortran90
FFLAGS   = $(FFLAGS_DEBUG)

MODULES = precision.f90 string.f90 $(wildcard module_*.f90)
SOURCES = $(wildcard test_*.f90) pacs_photproject.f90
EXECS = $(SOURCES:.f90=)

# apply a function to each element of a list
map = $(foreach a,$(2),$(call $(1),$(a)))

# recursively find dependencies
finddeps = $(1).o $(if $($(1)),$(call map,finddeps,$($(1))))

# define module dependencies
module_cfitsio = module_stdio
module_deglitching = module_math module_pointingmatrix precision
module_fitstools = string module_cfitsio
module_instrument = precision string
module_math = precision
module_optionparser = string
module_pacsinstrument = string module_fitstools module_math module_pacspointing module_pointingmatrix module_projection module_wcs 
module_pacspointing = module_math module_fitstools
module_pointingmatrix = module_pointingelement precision
module_preprocessor = module_math
module_projection = precision module_sort module_stack
module_wcs = module_fitstools module_wcslib string

# define executable dependencies
pacs_photproject = module_fitstools module_deglitching module_optionparser module_pacsinstrument module_pacspointing module_pointingmatrix module_preprocessor
test_cfitsio = module_cfitsio
test_deglitching = precision module_deglitching module_math module_pacsinstrument module_pointingmatrix
test_fitstools = module_fitstools
test_math = module_math precision
test_ngc6946_bpj = module_fitstools module_pacsinstrument module_pacspointing module_pointingmatrix module_preprocessor module_projection module_math
test_optionparser = module_optionparser
test_pacs = module_pacsinstrument module_pacspointing module_fitstools
test_pacsinstrument = module_pacsinstrument module_pacspointing
test_pacspointing = module_math module_pacspointing
test_projection = module_projection module_sort
test_read_config = module_instrument
test_sort = module_sort
test_stack = module_stack
test_stdio = module_stdio module_cfitsio
test_wcs = module_wcs module_fitstools module_math
test_wcslib1 = module_wcslib module_cfitsio 
test_wcslib2 = module_wcslib module_fitstools module_math
test_wcslibc = module_wcslibc module_cfitsio

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
	rm -f *.o *.mod *.so $(EXECS)

tests: $$(filter-out test_wcslib%,$$(filter test_%,$$(EXECS)))
	@for test in $^; do \
	echo;\
	echo "Running test: "$$test"...";\
	echo "=============";\
	./$$test; \
	done
