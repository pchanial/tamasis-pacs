#! /usr/bin/env python

import glob
import os
import string
import sys

from waflib import Options
from waflib.Build import BuildContext
from waflib.Configure import conf
from waflib.Context import Context

top = '.'
out = 'build'


# List of packages to be built
subdirs = ['core', 'madcap', 'pacs']

# Required libraries
libraries = ['CFITSIO', 'FFTW3', 'WCSLIB', 'OPENMP']

# Required Python packages
required_modules = ['numpy',
                    'scipy',
                    'matplotlib',
                    'fftw3',
                    'pyfits',
                    'kapteyn',
                    'mpi4py']


#
# Read options
#

def options(opt):

    opt.load('compiler_c')
    opt.load('compiler_fc')
    opt.load('python')
    opt.add_option('--precision-real',
                   action='store',
                   default='8',
                   choices=['4', '8', '16'],
                   help='precision of reals: single (4), double (8)or quadruple precision (16)')
    opt.add_option('--debug',
                   action='store_true',
                   default=False,
                   help='use compiler flags for debugging (export symbols, add checks...)')

#
# Configuration
#

def configure(conf):

    conf.load('compiler_c')
    conf.load('compiler_fc')
    conf.load('python')
    if conf.env.FC_NAME == 'GFORTRAN':
        conf.env.FCFLAGS = ['-ffree-form', '-fopenmp', '-Wall', '-fPIC', '-cpp']
        conf.env.FCFLAGS += ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if conf.options.debug else ['-O3']
        conf.env.LIB_OPENMP = ['gomp']
        conf.env.F2PYFCOMPILER = 'gnu95'
    elif conf.env.FC_NAME == 'IFORT':
        conf.env.FCFLAGS = ['-warn all', '-fpp', '-fPIC', '-free', '-openmp', '-ftz', '-fp-model precise', '-ftrapuv']
        conf.env.FCFLAGS = ['-debug', '-check all', '-traceback'] if conf.options.debug else ['-fast']
        conf.env.LIB_OPENMP = ['iomp5']
        conf.env.F2PYFCOMPILER = 'intelem'

    # these two environment variables are required by pkg-config
    os.putenv('PKG_CONFIG_ALLOW_SYSTEM_CFLAGS', '1')
    os.putenv('PKG_CONFIG_ALLOW_SYSTEM_LIBS',   '1')

    conf.check_python_version((2,6))
    for module in required_modules:
        conf.check_python_module(module)
    conf.check_fortran()
    conf.check_fortran_verbose_flag()
    conf.check_fortran_clib()
    conf.check_fortran_dummy_main()
    conf.check_fortran_mangling()
    if conf.options.precision_real == '16':
        conf.check_cc(
            fragment         = """
program test
   integer, parameter :: qp = selected_real_kind(2*precision(1.0d0))
   real(kind=qp)      :: quad
end program test""",
            compile_filename = 'test.f',
            features         = 'fc fcprogram',
            msg              = 'Checking for quadruple precision')
    conf.find_program('f2py' + sys.version[0:3], var='F2PY')
    conf.check_cfg(package='cfitsio', args=['--libs', '--cflags'])
    conf.check_cfg(package='fftw3',   args=['--libs', '--cflags'])
    conf.check_cfg(package='wcslib',  args=['--libs', '--cflags'])

    conf.env.SHAREDIR = os.path.abspath(conf.env.PYTHONDIR + '/../../../share')
    conf.define(conf.env.FC_NAME, 1)
    conf.define('TAMASISDIR', os.path.join(conf.env.SHAREDIR, 'tamasis'))
    conf.define('PRECISION_REAL', int(conf.options.precision_real))

    # Since WCSLIB doesn't provide F90 include files, we generate them
    # from the F77 ones
    f77dir = conf.root.find_node(conf.env['INCLUDES_WCSLIB'][0])
    f90dir = conf.bldnode.make_node('include').mkdir()
    f77files = f77dir.ant_glob('*inc')
    for f77file in f77dir.ant_glob('*inc'):
        f90file = conf.bldnode.make_node('include/'+f77file.name)
        f77_to_f90(f77file, f90file)

    # Write the file .f2py_f2cmap which maps the real type parameter 'p' to a C type
    ctype = {'4':'float', '8':'double', '16':'long_double'}
    os.system("echo {\\'real\\':{\\'p\\':\\'" + ctype[conf.options.precision_real] +              \
              "\\'}} > %s/.f2py_f2cmap" % out)

    print conf.env



#
# Build & Installation
#

def build(bld):

    # Build static libraries libtamasis*.a
    for subdir in subdirs:
        files = bld.srcnode.ant_glob(subdir+'/src/module_*.f')
        bld(features = 'fc fcstlib',
            source   = files,
            target   = 'tamasis' + subdir,
            includes = 'include',
            use      = libraries)
    bld.add_group()
    
    # Python extension tamasisfortran.so
    target = 'tamasisfortran.so'
    source = [bld.srcnode.find_node('%s/src/tamasisfortran_%s.f90' % (s,s)) for s in subdirs]

    #XXX this should be a Task...
    cmd = '${F2PY} --fcompiler=${F2PYFCOMPILER} --f90exec=${FC}'
    cmd += ' --f90flags="${FCFLAGS}"'
    cmd += ' -DF2PY_REPORT_ON_ARRAY_COPY=1'
    cmd += ' -m tamasisfortran -c ' + string.join(map(str, [s.abspath() for s in source]))
    cmd += ' ' + string.join([bld.env.fcstlib_PATTERN%('tamasis'+s) for s in reversed(subdirs)])
    cmd += ' ${FCINCPATH_ST:include} ' + string.join(['${FCINCPATH_ST:INCLUDES_%s}' % u for u in libraries])
    cmd += ' ' + string.join(['${FCLIBPATH_ST:LIBPATH_%s}' % u for u in libraries])
    cmd += ' ' + string.join(['${FCLIB_ST:LIB_%s}' % u for u in libraries])
    cmd += ' ${DEFINES_ST:DEFINES}'

    bld(rule=cmd,
        source=source + [bld.env.fcstlib_PATTERN%('tamasis'+s) for s in subdirs],
        target=target,
        use=libraries)

    # Installation
    bld.install_files('${PYTHONDIR}/tamasis', bld.srcnode.ant_glob('*/src/*py') + ['tamasisfortran.so'])
    for subdir in subdirs:
        datanode = bld.srcnode.find_node(subdir+'/data')
        if datanode is not None:
            bld.install_files('${SHAREDIR}/tamasis/'+subdir, datanode.ant_glob('*'))

        
#
# test, test-fortran, test-python
#

class test(BuildContext):
    """run Fortran and Python test suites"""
    cmd = 'test'
    fun = 'test_fun'

def test_fun(ctx):
    Options.commands = ['build', 'install', '_test-fortran', '_test-python'] + Options.commands

class test_fortran(BuildContext):
    """run Fortran test suite"""
    cmd = 'test-fortran'
    fun = 'test_fortran_fun'

def test_fortran_fun(ctx):
    Options.commands = ['build', 'install', '_test-fortran'] + Options.commands

class _test_fortran(BuildContext):
    """(internal)"""
    cmd = '_test-fortran'
    fun = '_test_fortran_fun'

def _test_fortran_fun(bld):
    for subdir in subdirs:
        files = bld.srcnode.ant_glob(subdir+'/test/test_*.f')
        for file in files:
            if str(file) in ('test_wcslib1.f', 'test_wcslib2.f'): continue
            bld(features = 'fc fcprogram',
                source   = file,
                target   = os.path.splitext(str(file))[0],
                includes = '. include',
                stlib    = ['tamasis' + s for s in reversed(subdirs)],
                libpath  = '.',
                use      = libraries)
    bld.add_group()
    for subdir in subdirs:
        files = bld.path.ant_glob(subdir+'/test/test_*.f')
        for file in files:
            if str(file).startswith('test_wcs'):
                continue
            bld(rule   = os.path.join(out,os.path.splitext(str(file))[0]),
                cwd    = bld.path.abspath(),
                always = True)
            bld.add_group()

class test_python(BuildContext):
    """run Python test suite"""
    cmd = 'test-python'
    fun = 'test_python_fun'

def test_python_fun(ctx):
    Options.commands = ['build', 'install', '_test-python'] + Options.commands

class _test_python(BuildContext):
    """(internal)"""
    cmd = '_test-python'
    fun = '_test_python_fun'

def _test_python_fun(bld):
    for subdir in subdirs:
        files = bld.path.ant_glob(subdir+'/test/test_*.py')
        for file in files:
            bld(rule   = 'python ' + file.abspath(),
                always = True)
            bld.add_group()


#
# loc, loc-fortran, loc-python
#

def loc(ctx):
    """total number of lines of code"""
    os.system('find . \( -name ".waf*" -prune -or -name build -prune -or -name .git -prune -or -name include -prune -or -name "wscript" -or -name "*f" -or -name "*f90" -or -name "*py" \) -and -not -type l -and -not -type d | xargs cat | wc -l')

class loc_fortran(Context):
    """number of lines of Fortran code"""
    cmd = 'loc-fortran'
    fun = 'loc_fortran_fun'

def loc_fortran_fun(ctx):
    os.system('find . \( -name ".waf*" -prune -or -name build -prune -name .git -prune -or -name include -prune -or -name "*f" -or -name "*f90" \) -and -not -type l -and -not -type d | xargs cat | wc -l')

class loc_python(Context):
    """number of lines of Python code"""
    cmd = 'loc-python'
    fun = 'loc_python_fun'

def loc_python_fun(ctx):
    os.system('find . \( -name ".waf*" -prune -or -name build -prune -name .git -prune -or -name include -prune -or -name "*py" \) -and -not -type l -and -not -type d | xargs cat | wc -l')


#
# Miscellaneous helpers
#

def f77_to_f90(source, target):
    f77name = source.abspath()
    f90name = target.abspath()
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