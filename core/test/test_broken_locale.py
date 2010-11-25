import tamasis

if not tamasis.tmf.test_broken_locale():
    raise Exception("The environment variable LC_NUMERIC should be set to 'POSIX'")
exit() # required for ipython <= 0.10
