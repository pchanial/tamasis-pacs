# Makefile for project TAMASIS
# Author: P. Chanial

# set up default compiler
ifeq "$(origin FC)" "default"
    ifneq (,$(shell which ifort))
        FC=ifort
    else ifneq ($(shell which gfortran),)
        FC=gfortran
    endif
endif

# set up compiler flags
ifeq "$(FC)" "gfortran"
    include Makefile-gfortran.mk
else ifeq ($(FC),ifort)
    include Makefile-ifort.mk
else
    $(error Unsupported compiler '$(FC)')
endif

# add libraries built by this Makefile
LDFLAGS += -L$(shell pwd)/lib

ifeq ($(PROF_GEN),1)
    DEBUG = 1
endif

ifeq ($(DEBUG),1)
    FFLAGS = $(FFLAGS_DEBUG) -DDEBUG
else
    FFLAGS = $(FFLAGS_RELEASE)
endif

# set up precision
ifeq (,$(PRECISION_REAL))
    PRECISION_REAL = 8
else ifeq (,$(filter 4 8 16,$(PRECISION_REAL)))
    $(error Valid PRECISION_REAL values are '4', '8' or '16')
endif
FFLAGS += -DPRECISION_REAL=$(PRECISION_REAL)
ifeq ($(SUPPORT_QUAD)$(PRECISION_REAL),16)
    $(error Quadruple precision is not supported by $(FC) version $(FCVERSION))
endif

# apply a function to each element of a list
map = $(foreach a,$(2),$(call $(1),$(a)))

# recursively find dependencies
finddeps = $(1) $(if $($(1)),$(call map,finddeps,$($(1))))

# directories to be processed
DIRS = core madcap pacs

# lists of files which are updated in sub-directory Makefile.mk files
MODULES :=
EXECS :=
INCLUDES := -Iinclude -Iinclude/wcslib-4.4.4-Fortran90

PYTHONSOURCES = $(wildcard */src/*.py)
PYTHONTESTS = $(wildcard */test/test_*.py)

.PHONY: all test test-fortran test-python clean distclean python loc
all: lib/libtamasis.so tamasis/tamasisfortran.so python

include $(patsubst %,%/Makefile.mk,$(DIRS))

# libtamasis.so library includes modules created in the sub-directories
lib/libtamasis.so: $(patsubst %,lib/libtamasis%.so,$(DIRS))
	$(FC) -shared -Wl,-soname,$@ -o $@ $^

# rule to create the .mod files without unnecessary recompilations
$(MODULES:.f=.mod):%.mod: %.f %.o
	@true

.SECONDEXPANSION:
$(MODULES:.f=.o) $(EXECS:.f=.o):%.o: %.f $$(addsuffix .mod,$$(sort $$(call map,finddeps,$$($$*))))
	$(FC) $(FFLAGS) $(FFLAGS_PROF) $(FFLAGS_MOD) $(@D) $(INCLUDES) -c -o $@ $<

# tamasisfortran.so should depend on libtamasis.so, but an address mapping problem occurs with the threadprivate variables in module_projection.f
tamasis/tamasisfortran.so: $(join $(DIRS),$(patsubst %,/src/tamasisfortran_%.f90,$(DIRS))) $(MODULES:.f=.o)
	@if ! test -e tamasis; then mkdir tamasis; fi;\
	if test ${PRECISION_REAL} = 4; then echo {\'real\':{\'p\':\'float\'}} > .f2py_f2cmap; fi;\
	if test ${PRECISION_REAL} = 8; then echo {\'real\':{\'p\':\'double\'}} > .f2py_f2cmap; fi;\
	if test ${PRECISION_REAL} = 16; then echo {\'real\':{\'p\':\'long_double\'}} > .f2py_f2cmap; fi;\
	unset LDFLAGS
	f2py --fcompiler=${FCOMPILER} --f90exec=$(FC) --f90flags="$(FFLAGS)" -DF2PY_REPORT_ON_ARRAY_COPY=1 -c $^ -m tamasisfortran $(INCLUDES) $(LDFLAGS); \
	mv -f tamasisfortran.so tamasis/

python:
	@if ! test -e tamasis; then mkdir tamasis; fi;\
	for file in $(PYTHONSOURCES); do \
	ln -sf ../$$file tamasis;\
	done

# test targets
test: test-fortran test-python

test-fortran: $(patsubst %,test-%,$(DIRS))

test-python: tamasis/tamasisfortran.so
	@for test in $(PYTHONTESTS); do \
	echo;\
	echo "Running Python test: "$$test"...";\
	echo "====================";\
	python $$test; \
	done

# clean targets
clean: $(patsubst %,clean-%,$(DIRS))
	@rm -f .f2py_f2cmap

distclean: $(patsubst %,distclean-%,$(DIRS))
	@rm -rf .f2py_f2cmap lib tamasis

# return number of lines of code
loc:
	@find . \( -name .git -prune -or -name include -prune -or -name "Makefile*" -or -name "*f" -or -name "*f90" -or -name "*py" \) -and -not -type l -and -not -type d | sort;\
	find . \( -name .git -prune -or -name include -prune -or -name "Makefile*" -or -name "*f" -or -name "*f90" -or -name "*py" \) -and -not -type l -and -not -type d | xargs cat | wc -l

loc-fortran:
	@find . \( -name .git -prune -or -name include -prune -or -name "*f" -or -name "*f90" \) -and -not -type l -and -not -type d | sort;\
	find . \( -name .git -prune -or -name include -prune -or -name "*f" -or -name "*f90" \) -and -not -type l -and -not -type d | xargs cat | wc -l

loc-python:
	@find . \( -name .git -prune -or -name include -prune -or -name "*py" \) -and -not -type l -and -not -type d | sort;\
	find . \( -name .git -prune -or -name include -prune -or -name "*py" \) -and -not -type l -and -not -type d | xargs cat | wc -l
