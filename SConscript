import os
import shutil
import string
import sys

Import('env', 'variant_dir')
subdirs = ['core', 'madcap', 'pacs']
module_dir = map(lambda d: os.path.join(d, 'src'), subdirs)
module_dir.append('include')

fortran = None
required_fortran = {'gfortran':'4.5.0', 'ifort':'11.1'}
python_version = sys.version[0:3]
python = 'python' + python_version
required_python = '2.6'
handle_precision_quad = True


#
# PATH variables
#
path = os.environ['PATH']
if os.environ.has_key('LD_LIBRARY_PATH'):
    ld_library_path = filter(lambda x:len(x)>0,os.environ['LD_LIBRARY_PATH'].split(':'))
else:
    ld_library_path = []

#
# Environment
#
fortran = 'gfortran'

env['F2PY'] = 'f2py' + python_version
env['F2PYFCOMPILER'] ='gnu95'
env['FORTRANFLAGS'] = '-fbacktrace -O3 -ffree-form -fopenmp -Wall -fPIC -cpp'
env['CPPDEFINES'] = {'GFORTRAN':None}
env['FORTRAN'] = '/usr/bin/gfortran-4.5'
env['FORTRANPATH'] = module_dir
env['FORTRANMODDIRPREFIX'] = '-J '
#env['FORTRANMODDIR'] = '${TARGET.dir}'
env['FORTRANMODDIR'] = 'include'
env['LIBS'] = ['gomp']
env['LIBPATH'] = ld_library_path
env['SHLINKFLAGS'] = '-shared'
env['PATH'] = path

fortran_version = os.popen(env['FORTRAN'] + ' -dumpversion').read().replace('\n', '')

debug = ARGUMENTS.get('debug', 0)
if int(debug):
    env.Append(FORTRANFLAGS = '-g -fcheck=all')

if fortran_version < '4.6':
    handle_precision_quad = False


#
# Execution environment
#
pythonpath = '.'
if os.environ.has_key('PYTHONPATH'):
    pythonpath = pythonpath + ':' + os.environ['PYTHONPATH']
env['ENV']['PYTHONPATH'] = pythonpath
ld_library_path_env = 'lib' + \
    (':' + string.join(ld_library_path,':') if len(ld_library_path) > 0 else '')
env['ENV']['LD_LIBRARY_PATH'] = ld_library_path_env
env['ENV']['PKG_CONFIG_PATH'] = '/home/pchanial/software/lib/pkgconfig'
env['ENV']['PKG_CONFIG_ALLOW_SYSTEM_CFLAGS'] = 1
env['ENV']['PKG_CONFIG_ALLOW_SYSTEM_LIBS']   = 1

#
# options
#
precision_real = int(ARGUMENTS.get('precision_real', 8))
if precision_real not in (4, 8, 16):
    print('Real precision must be 4 (single), 8 (double) or 16 (quad).')
    Exit(1)
env.Append(CPPDEFINES={'PRECISION_REAL':precision_real})

# the variable CPPDEFINES is not used by the fortran builder...
#env.Append(FORTRANFLAGS=' ' + env.subst('$_CPPDEFFLAGS'))


#
# Configuration
#
def CheckFortran(context):
    context.Message('Checking for Fortran compiler...')
    result = fortran in required_fortran.keys()
    context.Result(result)
    return result

def CheckFortranVersion(context):
    context.Message('Checking for Fortran compiler version...')
    result = fortran_version >= required_fortran[fortran]
    context.Result(result)
    return result

def CheckPythonVersion(context):
    context.Message('Checking for Python version...')
    result = python_version >= required_python and python_version < '3'
    context.Result(result)
    return result

def CheckQuadPrecision(context):
    print 'ici'
    context.Message('Checking for quadruple precision support...')
    context.Result(handle_precision_quad)
    return handle_precision_quad

required_packages = ['numpy',
                     #'scipy',
                     'matplotlib',
                     'fftw3',
                     'pyfits',
                     'kapteyn',
                     'mpi4py']

def CheckPackage(context, package):
    context.Message('Checking for Python package ' + package + '...')
    try:
        __import__(package)
        result = True
    except ImportError:
        result = False
    context.Result(result)
    return result

if not env.GetOption('clean'):
    tests = {}
    tests['CheckFortran'] = CheckFortran
    tests['CheckFortranVersion'] = CheckFortranVersion
    tests['CheckPythonVersion'] = CheckPythonVersion
    tests['CheckQuadPrecision'] = CheckQuadPrecision
    for package in required_packages:
        tests['CheckPackage_' + package] = lambda c: CheckPackage(c, package)
    conf = Configure(env, custom_tests=tests)

    if not conf.CheckFortran():
        if fortran is None:
            print('No Fortran compiler has been found.')
            Exit(1)
        
        print('The compiler ' + fortran + ' is not supported.')
        Exit(1)

    if not conf.CheckFortranVersion():
        print('The version of the compiler ' + env['FORTRAN'] + ' (' +         \
              fortran_version + ') should be ' + 'newer or equal to ' +        \
              required_fortran[fortran] + '.')
        Exit(1)

    if not conf.CheckPythonVersion():
        if python_version >= '3':
            print('Python 3 is not supported.')
        else:
            print('The python version (' + python_version + ') should be ' +   \
                  'newer or equal to ' + required_python + '.')
        Exit(1)

    if precision_real == 16 and not conf.CheckQuadPrecision():
        print('The compiler ' + env['FORTRAN'] + ' (version ' +                \
              fortran_version + ') does not handle quad precision.')
        Exit(1)

    if not conf.CheckLibWithHeader('m', 'math.h', 'c'):
        print('Did not find math library, exiting.')
        Exit(1)

    if not conf.CheckLib('cfitsio'):
        print('Did not find cfitsio library, exiting.')
        Exit(1)

    if not conf.CheckLib('fftw3'):
        print('Did not find fftw3 library, exiting.')
        Exit(1)

    for package in required_packages:
        if not getattr(conf, 'CheckPackage_' + package)():
            print('Python package ' + package + ' is not found.')
            Exit(1)

    env = conf.Finish()

#
# Libraries
#
env['STATIC_AND_SHARED_OBJECTS_ARE_THE_SAME'] = 1
env.ParseConfig("pkg-config cfitsio --cflags")
env.ParseConfig("pkg-config fftw3   --cflags")
env.ParseConfig("pkg-config wcslib  --cflags")

def get_include_path(package):
    includes = env.Split(env.backtick('pkg-config --cflags-only-I ' + package) \
                  .replace('\n', ''))
    if len(includes) == 0:
        raise ValueError('The include path for package ' + package +           \
                         ' cannot be found')
    if len(includes) > 1:
        raise ValueError('The include path for package ' + package +           \
                         ' is ambiguous: ' + str(includes))
    return includes[0][2:]

# FFTW3
env.Append(FORTRANPATH=get_include_path('fftw3'))

# WCSLIB
def f77_to_f90(target, source, env):
    f90name = str(target[0])
    f77name = str(source[0])
    with open(f77name, 'r') as f77:
        f77content = f77.readlines()
        f90content = []
        for i, l in enumerate(f77content):
            if l[0:1] != ' ':
                continue
            if i+1 < len(f77content) and f77content[i+1][5:6] not in ('',' '):
                l = l[0:-1] + ' &\n'
            l = '      ' + l[6:]
            f90content.append(l)
        with open(f90name, 'w') as f90:
            f90.writelines(f90content)
files = Glob(get_include_path('wcslib') + os.sep + '*inc')
for f77File in files:
    f90File = File(os.path.join('include', os.path.basename(str(f77File))))
    env.Command(f90File, f77File, f77_to_f90)


#
# Build part
#
objs = []
libs = []
objsall = []
modsall = []
for subdir in subdirs:
    objs_mods = env.Object(Glob(os.path.join(subdir, 'src', '*.f')))
    # Remove any mod files. These should not be passed to the linker.
    objs = filter(lambda o: str(o)[-4:] != '.mod', objs_mods)
    libs.append(env.SharedLibrary('lib/tamasis'+subdir, objs)[0])
    objsall += objs
    modsall += filter(lambda o: str(o)[-4:] == '.mod', objs_mods)
objall = Flatten(objsall)
libtamasis = env.Library('lib/tamasis', libs)


#
# Python/Fortran interface
#
def f2py_builder_fn(target, source, env):
    ctype = {4:'float', 8:'double', 16:'long_double'}
    os.system("echo {\\'real\\':{\\'p\\':\\'" + ctype[precision_real] + "\\'}} > .f2py_f2cmap")
    for t in target:
        t_ = os.path.splitext(os.path.basename(str(t)))[0]
        cmd = '$F2PY --fcompiler=$F2PYFCOMPILER --f90exec=$FORTRAN'
        cmd += ' --f90flags="$FORTRANFLAGS $CPPFLAGS $_CPPDEFFLAGS"'
        cmd += ' -DF2PY_REPORT_ON_ARRAY_COPY=1'
        cmd += ' -m ' + t_ + ' -c ' + string.join([str(s) for s in source], ' ')
        cmd += ' ' + string.join(map(str, objall))
        cmd += ' -I' + variant_dir + os.sep + 'include'
        cmd += string.join([' -l'+f for f in env['LIBS']], '')
        cmd += string.join([' -L'+f for f in env['LIBPATH']], '')
        result = env.Execute(cmd)
        if result:
            return result
        shutil.move(t_+'.so', str(t))
    return None

f2py_builder = Builder(action=f2py_builder_fn, suffix = '.so')
env.Append(BUILDERS={'F2py' : f2py_builder})
libtamasisfortran = env.F2py(target='lib'+os.sep+'tamasisfortran',
                             source=[os.path.join(s,'src','tamasisfortran_' + \
                                     s + '.f90') for s in subdirs])
Depends(libtamasisfortran, objall)


#
# Installation
#
Install(env['INCDIR'], modsall)
Install(env['LIBDIR'], [libtamasis, libs])
Install(env['PYTHONDIR'], [libtamasisfortran])
for subdir in subdirs:
    files = Glob(os.path.join(subdir, 'src', '*.py'))
    Install(env['PYTHONDIR'], files)
    Clean('.', [str(file)+'c' for file in files])


#
# test, test-fortran, test-python
#
test_fortran_actions = []
env_test = env.Clone()
env_test.Append(LIBS=['wcs','tamasis']) #XXX wcs should not be added
env_test.Append(LIBPATH='lib')
env_test.AppendENVPath('LD_LIBRARY_PATH', variant_dir+os.sep+'lib')
for subdir in subdirs:
    srcs = Glob(os.path.join(subdir, 'test', 'test_*.f'))
    tests = []
    for src in srcs:
        if os.path.basename(str(src)) in ('test_wcslib1.f','test_wcslib2.f'):
            continue
        prog = env_test.Program(src)[0]
        tests.append(prog)
        Depends(prog, libtamasis)
    Depends('test-fortran', tests)
    title = 'Running Fortran ' + subdir.upper() + ' test:'
    test_fortran_actions += ['@echo;echo ' + title + ' '+ str(t) + '...;echo '+\
                             len(title)*'=' + '; ' + variant_dir + os.sep + str(t) for t in tests]

test_python_actions = []
for subdir in subdirs:
    tests = Glob(os.path.join(subdir, 'test', 'test_*py'))
    title = 'Running Python ' + subdir.upper() + ' test:'
    test_python_actions += ['@echo;echo ' + title + ' '+ str(t) + '...;'       \
                            'echo ' + len(title)*'=' + ';  ' + python + ' ' +  \
                            str(t) for t in tests]

test_fortran = env_test.AlwaysBuild(env.Alias('test-fortran', [], test_fortran_actions))
test_python = env.AlwaysBuild(env.Alias('test-python', [], test_python_actions))
env.AlwaysBuild(env.Alias('test', [test_fortran, test_python]))

#
# loc, loc-fortran, loc-python
#
env.AlwaysBuild(env.Alias('loc', [], '@find . \( -name .git -prune -or -name include -prune -or -name "Makefile*" -or -name "*f" -or -name "*f90" -or -name "*py" \) -and -not -type l -and -not -type d | xargs cat | wc -l'))
env.AlwaysBuild(env.Alias('loc-fortran', [], '@find . \( -name .git -prune -or -name include -prune -or -name "*f" -or -name "*f90" \) -and -not -type l -and -not -type d | xargs cat | wc -l'))
env.AlwaysBuild(env.Alias('loc-python', [], '@find . \( -name .git -prune -or -name include -prune -or -name "*py" \) -and -not -type l -and -not -type d | xargs cat | wc -l'))
