platform = ARGUMENTS.get('OS', Platform())

include = "#export/include"
lib = "#export/lib"
bin = "#export/bin"

variant_dir = 'build/' + str(platform)

env = Environment(PLATFORM = platform,
                  BINDIR = bin,
                  INCDIR = include,
                  LIBDIR = lib,
                  PYTHONDIR = '#export/tamasis',
                  FORTRANPPFILESUFFIXES = ['.f', '.f90'],
                 )

Export('env', 'variant_dir')
env.SConscript('SConscript', variant_dir=variant_dir, duplicate=1)
Clean('.', 'build')
