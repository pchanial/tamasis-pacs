RANLIBSOURCES = $(wildcard ranlib/src/*.f) $(wildcard ranlib/linpack/*.f)
RANLIBTESTS = $(wildcard ranlib/test/*.f)
RANLIBFFLAGS = -O3 -fPIC -std=legacy

.PHONY: ranlib
ranlib: ranlib/libranlib.so

ranlib/libranlib.so: $(RANLIBSOURCES:.f=.o)
	$(FC) -shared -Wl,-soname,$@ -o $@ $^

ranlib/src/%.o: ranlib/src/%.f
	$(FC) $(RANLIBFFLAGS) -c -o $@ $<

ranlib/linpack/%.o: ranlib/linpack/%.f
	$(FC) $(RANLIBFFLAGS) -c -o $@ $<

ranlib/test/%.o: ranlib/test/%.f
	$(FC) $(RANLIBFFLAGS) -c -o $@ $<

ranlib/test/%: ranlib/test/%.o ranlib/libranlib.so
	$(FC) -o $@ $@.o -Lranlib -lranlib

clean-ranlib:
	@rm -f $(RANLIBTESTS:.f=);\
	find ranlib \( -name '*.o' -or -name "*.mod" -or -name "*.so" -or -name "*~" \) -exec rm {} ';'

tests-ranlib: $(RANLIBTESTS:.f=)
	ranlib/test/tstbot