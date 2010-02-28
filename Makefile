# Makefile for project TAMASIS
# Author: P. Chanial

FC=gfortran
FFLAGS_DEBUG   = -g -O3 -fcheck=all -fopenmp -Wall -fPIC
FFLAGS_RELEASE = -O3 -fopenmp -fPIC
FFLAGS   = $(FFLAGS_DEBUG)
LDFLAGS  = -lgomp $(shell pkg-config --libs cfitsio) $(shell pkg-config --libs wcslib)
INCLUDES = wcslib-4.4.4-Fortran90

MODULES = precision.f90 string.f90 $(wildcard module_*.f90)
SOURCES = $(wildcard test_*.f90)
EXECS = $(SOURCES:.f90=)

# apply a function for each element of a list
map = $(foreach a,$(2),$(call $(1),$(a)))

# recursively find dependencies
finddeps = $(1).o $(if $($(1)),$(call map,finddeps,$($(1))))

# define module dependencies
module_fitstools = string module_cfitsio module_wcslib
module_instrument = precision string
module_pacsinstrument = string module_fitstools module_pacspointing module_pointingmatrix module_projection module_wcs module_wcslib 
module_pacspointing = precision module_fitstools
module_pointingmatrix = module_pointingelement
module_projection = precision module_sort module_stack
module_wcs = module_fitstools module_wcslib

# define executable dependencies
test_cfitsio = module_cfitsio module_stdio
test_fitstools = module_fitstools module_wcslib
test_ngc6946_bpj = module_fitstools module_pacsinstrument module_pacspointing module_pointingmatrix module_preprocessor module_projection module_wcslib precision
test_pacs = module_pacsinstrument module_pacspointing module_fitstools module_wcslib
test_pointing = module_pacspointing
test_projection = module_projection module_sort
test_read_config = module_instrument
test_sort = module_sort
test_stack = module_stack
test_stdio = module_stdio module_cfitsio
test_wcslib1 = module_wcslib module_cfitsio 
test_wcslib2 = module_wcslib module_fitstools
test_wcslibc = module_wcslibc module_cfitsio

.PHONY : all
all : $(EXECS) tamasisfortran.so

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
	f2py --fcompiler=gnu95 --f90exec=$(FC) --f90flags="-fopenmp" -DF2PY_REPORT_ON_ARRAY_COPY=1 -c $^ -m tamasisfortran $(LDFLAGS)

clean:
	rm -f *.o *.mod *.so $(EXECS)

tests:
	@for test in $(EXECS); do \
	echo;\
	echo "Running test: "$$test"...";\
	echo "=============";\
	./$$test; \
	done
