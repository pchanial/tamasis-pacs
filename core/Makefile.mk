DIR := core
SDIR := $(DIR)/src
TDIR := $(DIR)/test

MODULESOURCES := $(wildcard $(SDIR)/module_*.f)
TESTSOURCES   := $(filter-out $(TDIR)/test_wcslib%.f,$(wildcard $(TDIR)/test_*.f))

MODULES += $(MODULESOURCES)
EXECS += $(TESTSOURCES)
INCLUDES += -I$(SDIR)

# define module dependencies
$(SDIR)/module_cfitsio := $(SDIR)/module_stdio
$(SDIR)/module_compression := $(SDIR)/module_tamasis
$(SDIR)/module_deglitching := $(SDIR)/module_math $(SDIR)/module_pointingmatrix $(SDIR)/module_string $(SDIR)/module_tamasis
$(SDIR)/module_filtering := $(SDIR)/module_tamasis
$(SDIR)/module_fitstools := $(SDIR)/module_cfitsio $(SDIR)/module_precision $(SDIR)/module_string
$(SDIR)/module_instrument := $(SDIR)/module_string $(SDIR)/module_tamasis
$(SDIR)/module_math := $(SDIR)/module_tamasis
$(SDIR)/module_observation := $(SDIR)/module_fitstools $(SDIR)/module_math $(SDIR)/module_string $(SDIR)/module_tamasis
$(SDIR)/module_optionparser := $(SDIR)/module_string $(SDIR)/module_tamasis
$(SDIR)/module_pointingmatrix := $(SDIR)/module_math $(SDIR)/module_precision $(SDIR)/module_projection $(SDIR)/module_tamasis
$(SDIR)/module_preprocessor := $(SDIR)/module_math $(SDIR)/module_sort $(SDIR)/module_tamasis
$(SDIR)/module_projection := $(SDIR)/module_sort $(SDIR)/module_stack $(SDIR)/module_tamasis
$(SDIR)/module_sort := $(SDIR)/module_math $(SDIR)/module_tamasis
$(SDIR)/module_string := $(SDIR)/module_tamasis
$(SDIR)/module_tamasis := $(SDIR)/module_precision
$(SDIR)/module_wcs := $(SDIR)/module_fitstools $(SDIR)/module_math $(SDIR)/module_string $(SDIR)/module_tamasis $(SDIR)/module_wcslib

.PHONY: core test-core clean-core distclean-core
core: ranlib lib/libtamasiscore.so

lib/libtamasiscore.so: $(MODULESOURCES:.f=.o)
	@if ! test -e lib; then mkdir lib; fi
	$(FC) -shared -Wl,-soname,$@ -o $@ $^

$(TESTSOURCES:.f=):%: lib/libtamasiscore.so %.o
	$(FC) -o $@ $@.o $(LDFLAGS) -ltamasiscore

test-core: lib/libtamasis.so
test-core: $(TESTSOURCES:.f=)
	@for test in $(filter-out %.so,$^); do \
	echo;\
	echo "Running CORE test: "$$test"...";\
	echo "==================";\
	./$$test; \
	done

clean-core:
	@find core \( -perm /u=x -and -type f -and -not -name "*py" \) -exec rm {} ';';\
	find core \( -name '*.o' -or -name "*.mod" -or -name "*~" -or -name "*pyc" \) -exec rm {} ';'

distclean-core: clean-core
	@rm -f lib/libtamasiscore.so