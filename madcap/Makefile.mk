DIR  := madcap
SDIR := $(DIR)/src
TDIR := $(DIR)/test

MODULESOURCES := $(wildcard $(SDIR)/module_*.f)
TESTSOURCES   := $(wildcard $(TDIR)/test_*.f)

MODULES += $(MODULESOURCES)
EXECS += $(TESTSOURCES)
INCLUDES += -I$(SDIR)

.PHONY: madcap test-madcap clean-madcap
madcap: core $(MODULESOURCES:.f=.mod) lib/libtamasismadcap.so

lib/libtamasismadcap.so: lib/libtamasiscore.so $(MODULESOURCES:.f=.o)
	@if ! test -e lib; then mkdir lib; fi
	$(FC) -shared -Wl,-soname,$@ -o $@ $^

$(TESTSOURCES:.f=):%: lib/libtamasiscore.so lib/libtamasismadcap.so %.o
	$(FC) -o $@ $@.o $(LDFLAGS) -ltamasiscore -ltamasismadcap

test-madcap: $(TESTSOURCES:.f=)
	@for test in $^; do \
        echo;\
        echo "Running MADCAP test: "$$test"...";\
        echo "==================";\
        ./$$test; \
        done

clean-madcap:
	@find madcap \( -perm /u=x -and -type f -and -not -name "*py" \) -exec rm {} ';';\
	find madcap \( -name '*.o' -or -name "*.mod" -or -name "*~" -or -name "*pyc" \) -exec rm {} ';'

dist-clean-madcap: clean-madcap
	@rm -f lib/libtamasismadcap.so
