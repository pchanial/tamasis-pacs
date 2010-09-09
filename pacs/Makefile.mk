DIR  := pacs
SDIR := $(DIR)/src
TDIR := $(DIR)/test

MODULESOURCES := $(wildcard $(SDIR)/module_*.f)
EXECSOURCES   := $(SDIR)/pacs_photproject.f
TESTSOURCES   := $(wildcard $(TDIR)/test_*.f)

MODULES += $(MODULESOURCES)
EXECS += $(EXECSOURCES) $(TESTSOURCES)
INCLUDES += -I$(SDIR)

$(SDIR)/module_pacsinstrument := $(SDIR)/module_pacsobservation

.PHONY: pacs test-pacs clean-pacs
pacs: core $(MODULESOURCES:.f=.mod) $(EXECSOURCES:.f=)

lib/libtamasispacs.so: lib/libtamasiscore.so $(MODULESOURCES:.f=.o)
	@if ! test -e lib; then mkdir lib; fi
	$(FC) -shared -Wl,-soname,$@ -o $@ $^

$(EXECSOURCES:.f=) $(TESTSOURCES:.f=):%: lib/libtamasiscore.so lib/libtamasispacs.so %.o
	$(FC) -o $@ $@.o $(LDFLAGS) -ltamasiscore -ltamasispacs

test-pacs: $(TESTSOURCES:.f=)
	@for test in $^; do \
        echo;\
        echo "Running PACS test: "$$test"...";\
        echo "==================";\
        ./$$test; \
        done

clean-pacs:
	@find pacs \( -perm /u=x -and -type f -and -not -name "*py" \) -exec rm {} ';';\
	find pacs \( -name '*.o' -or -name "*.mod" -or -name "*~" -or -name "*pyc" \) -exec rm {} ';'

dist-clean-pacs: clean-pacs
	@rm -f lib/libtamasispacs.so
