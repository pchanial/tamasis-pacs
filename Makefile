FC=gfortran
FFLAGS_DEBUG   = -g -O3 -fcheck=all -fopenmp -Wall -fPIC
FFLAGS_RELEASE = -O3 -fopenmp -fPIC
FFLAGS   = $(FFLAGS_DEBUG)
LDFLAGS  = -lgomp $(shell pkg-config --libs cfitsio) $(shell pkg-config --libs wcslib)
INCLUDES = wcslib-4.4.4-Fortran90

MODULES = precision.f90 string.f90 $(wildcard module_*.f90)
SOURCES = $(wildcard test_*.f90)
EXECS = $(SOURCES:.f90=)

.PHONY : all

all : $(EXECS) tamasisfortran.so

%.mod : %.o
	@if [ ! -f $@ ]; then \
	    rm $< ;\
	    $(MAKE) $< ;\
	fi

%.o : %.f90
	$(FC) $(FFLAGS) -I$(INCLUDES) -c -o $@ $<

module_fitstools.o : module_cfitsio.mod module_wcslib.mod
module_pacsinstrument.o : string.mod module_fitstools.mod module_pacspointing.mod module_pointingmatrix.mod module_projection.mod module_wcslib.mod 
module_pacspointing.o : precision.mod module_fitstools.mod
module_pointingmatrix.o : module_pointingelement.mod
module_projection.o : precision.mod module_sort.mod module_stack.mod

test_cfitsio.o : module_cfitsio.mod module_stdio.mod
test_fitstools.o : module_fitstools.mod module_wcslib.mod
test_ngc6946_bpj.o : module_fitstools.mod module_pacsinstrument.mod module_pacspointing.mod module_pointingmatrix.mod module_preprocessor.mod module_projection.mod module_wcslib.mod precision.mod
test_pacs.o : module_pacsinstrument.mod module_pacspointing.mod module_fitstools.mod module_wcslib.mod
test_pointing.o : module_pacspointing.mod
test_projection.o : module_projection.mod module_sort.mod
test_read_config.o : module_instrument.mod
test_sort.o : module_sort.mod
test_stack.o : module_stack.mod
test_stdio.o : module_stdio.mod module_cfitsio.mod
test_wcslib1.o : module_wcslib.mod module_cfitsio.mod 
test_wcslib2.o : module_wcslib.mod module_fitstools.mod
test_wcslibc.o : module_wcslibc.mod module_cfitsio.mod

test_cfitsio : module_cfitsio.o module_stdio.o
test_fitstools : module_cfitsio.o module_fitstools.o module_wcslib.o
test_ngc6946_bpj : module_cfitsio.o module_fitstools.o module_pacsinstrument.o module_pacsinstrument.o module_pacspointing.o module_pointingelement.o module_pointingmatrix.o module_projection.o module_preprocessor.o module_sort.o module_stack.o module_wcslib.o precision.o string.o
test_pacs : module_pacsinstrument.o string.o module_fitstools.o module_cfitsio.o module_stack.o module_sort.o module_projection.o module_pacspointing.o module_pointingelement.o module_pointingmatrix.o module_wcslib.o
test_pointing : module_pacspointing.o precision.o module_fitstools.o module_cfitsio.o
test_projection : module_projection.o module_sort.o module_stack.o
test_read_config :  module_instrument.o precision.o string.o
test_sort : module_sort.o
test_stack : module_stack.o
test_stdio : module_stdio.o module_cfitsio.o
test_wcslib1 : module_wcslib.o module_cfitsio.o
test_wcslib2 : module_wcslib.o module_cfitsio.o module_fitstools.o
test_wcslibc : module_wcslibc.o module_cfitsio.o module_fitstools.o string.o

%: %.o
	$(FC) -o $@ $^ $(LDFLAGS)

tamasisfortran.so: tamasisfortran.f90 $(MODULES:.f90=.o)
	unset LDFLAGS ; \
	#$(FC) $(FFLAGS) -I$(INCLUDES) -c -o $@ $< ; \
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
