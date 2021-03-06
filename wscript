#! /usr/bin/env python

import glob
import multiprocessing
import numpy as np
import os
import string
import subprocess
import sys

from waflib import Logs, Options, Tools
from waflib.Build import BuildContext
from waflib.Configure import conf
from waflib.Context import Context

APPNAME = 'tamasis'
VERSION = '5.0'

top = '.'
out = 'build'


# List of packages to be built
subdirs = ['core', 'madcap', 'pacs']

# Required libraries
libraries = ['CFITSIO', 'FFTW3', 'OPENMP']

# Required Python packages
required_modules = ['numpy>=1.6',
                    'scipy>=0.9',
                    'pyoperators>=0.7',
                    'pysimulators>=0.2',
                    'matplotlib',
                    'fftw3',
                    'pyfits>=3.0',
                    'kapteyn',
                    'numexpr>=2.0']



#
# Read options
#

def options(opt):

    opt.load('compiler_c')
    opt.load('compiler_fc')
    opt.load('python')
    opt.add_option('--enable-mpi',
                   action='store_true',
                   default=False,
                   help='MPI-enabled installation')
    opt.add_option('--enable-wcslib',
                   action='store_true',
                   default=False,
                   help='enables linking WCSLIB')
    opt.add_option('--precision-real',
                   action='store',
                   default='8',
                   choices=['4', '8', '16'],
                   help='precision of reals: single (4), double (8) or quadrupl'
                        'e precision (16)')
    opt.add_option('--debug',
                   action='store_true',
                   default=False,
                   help='use compiler flags for debugging (export symbols, add '
                        'checks...)')
    opt.add_option('--np',
                   action='store',
                   help='number of MPI copies')



#
# Configuration
#

def configure(conf):

    global required_modules, libraries

    if conf.options.enable_mpi:
        required_modules += ['mpi4py']
        libraries += ['MPI']

    if conf.options.enable_wcslib:
        libraries += ['WCSLIB']

    for platform in 'linux', 'darwin', 'default':
        Tools.compiler_fc.fc_compiler[platform] = ['ifort', 'gfortran']

    conf.load('compiler_c')
    conf.load('compiler_fc')
    conf.load('python')

    conf.check_python_version((2,6))
    for module in required_modules:
        conf.check_python_module(module)
    conf.check_fortran()
    conf.check_fortran_verbose_flag()
    conf.check_fortran_clib()
    conf.check_fortran_dummy_main()
    conf.check_fortran_mangling()

    if conf.env.FC_NAME == 'GFORTRAN':
        conf.env.FCFLAGS = ['-ffree-form', '-Wall', '-fPIC', '-cpp']
        if conf.options.debug:
            conf.env.FCFLAGS += ['-g', '-fcheck=all', '-fbacktrace']
        else:
            conf.env.FCFLAGS += ['-O3']
        if 'OPENMP' in libraries:
            conf.env.FCFLAGS_OPENMP = ['-fopenmp']
            conf.env.LIB_OPENMP = ['gomp']
        conf.env.F2PYFCOMPILER = 'gnu95'
        conf.env.LIB_LAPACK = ['lapack', 'blas']
    elif conf.env.FC_NAME == 'IFORT':
        conf.env.FCFLAGS = ['-fpp', '-fPIC', '-free', '-ftz', '-fp-model',
                            'precise', '-ftrapuv', '-fast', '-warn', 'all']
        if conf.options.debug:
            conf.env.FCFLAGS += ['-debug', '-check', 'all', '-traceback']
        if 'OPENMP' in libraries:
            conf.env.FCFLAGS_OPENMP = ['-openmp']
            conf.env.LIB_OPENMP = ['iomp5']
        conf.env.F2PYFCOMPILER = 'intelem'
        # should add -Wl,--start-group and -Wl,--end-group
        conf.env.LIB_LAPACK = ['mkl_intel_lp64','mkl_intel_thread','mkl_core', 'iomp5']

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

    conf.find_program(['ipython'  + sys.version[0:3], 
                       'ipython-' + sys.version[0:3],
                       'ipython'], var='IPYTHON')
    conf.env.IPYTHON_VERSION = _get_version(conf.env.IPYTHON, flag='-Version')

    conf.find_program(['f2py'  + sys.version[0:3], 
                       'f2py-' + sys.version[0:3],
                       'f2py'], var='F2PY')
        
    conf.find_program(['nosetests'  + sys.version[0:3], 
                       'nosetests-' + sys.version[0:3],
                       'nosetests'], var='NOSETESTS')

    # these two environment variables are required by pkg-config
    os.putenv('PKG_CONFIG_ALLOW_SYSTEM_CFLAGS', '1')
    os.putenv('PKG_CONFIG_ALLOW_SYSTEM_LIBS',   '1')
    conf.check_cfg(package='cfitsio', args=['--libs', '--cflags'])
    conf.check_cfg(package='fftw3',   args=['--libs', '--cflags'])
    conf.check_cc(
        fragment         = """
program test
    include 'fftw3.f'
    integer, parameter :: dp = kind(1.d0)
    integer*8          :: plan
    real(dp)           :: in(10), out(10)
    call dfftw_plan_r2r_1d(plan, 10, in, out, FFTW_R2HC, FFTW_ESTIMATE)
    call dfftw_destroy_plan(plan)
end program test""",
        compile_filename = 'test.f',
        features         = 'fc fcprogram',
        msg              = "Checking for 'fftw3' double precision",
        use              = ['FFTW3'])

    if 'LAPACK' in libraries:
        fragment = """
program test
    integer, parameter :: p = kind(1.d0)
    real(p)            :: u(2), v(2)
    u = 1._p
    v = 0._p
    call daxpy(2, 3._p, u, 1, v, 1)
end program test
"""
        if conf.options.precision_real == '4':
            fragment = fragment.replace('kind(1.d0)', 'kind(1.)') \
                               .replace('daxpy','saxpy')
        conf.check_cc(
            fragment         = fragment,
            compile_filename = 'test.f',
            features         = 'fc fcprogram',
            msg              = "Checking for 'blas'",
            use              = ['LAPACK'])

        fragment = """
program test
    integer, parameter :: p = kind(1.d0)
    integer            :: status
    real(p)            :: cd(2,2), vr(2), vim(2), junk, work(6)
    cd = 0
    call dgeev('N', 'N', 2, cd, 2, vr, vim, junk, 1, junk, 1, work, size(work), status)
end program test
"""
        if conf.options.precision_real == '4':
            fragment = fragment.replace('kind(1.d0)', 'kind(1.)') \
                               .replace('dgeev','sgeev')
        conf.check_cc(
            fragment         = fragment,
            compile_filename = 'test.f',
            features         = 'fc fcprogram',
            msg              = "Checking for 'lapack'",
            use              = ['LAPACK'])

    if conf.options.enable_mpi:
        conf.env.OMPI_FC = conf.env.FC
        conf.find_program('mpif90', var='MPIFC')
        conf.env.FC = conf.env.MPIFC
        conf.check_cc(fragment    = 'program test\nuse mpi\nend',
                      compile_filename = 'test.f',
                      features    = 'fc fcprogram',
                      msg         = 'Checking for the MPI fortran module',
                      define_name = 'HAVE_MPI_MODULE',
                      mandatory   = False)
        if conf.env['HAVE_MPI_MODULE'] == 0:
            conf.check_cc(fragment= "program test\ninclude 'mpif.h'\nend",
                          compile_filename = 'test.f',
                          features    = 'fc fcprogram',
                          msg         = 'Checking for the MPI fortran header',
                          define_name = 'HAVE_MPI_HEADER')
        conf.define('HAVE_MPI', 1)


    if 'WCSLIB' in libraries:
        conf.check_cfg(package='wcslib',  args=['--libs', '--cflags'])
        conf.check_cfg(modversion='wcslib')
        check_wcslib_external(conf.env)
        # Since WCSLIB doesn't provide F90 include files, we generate them
        # from the F77 ones
        f77dir = conf.root.find_node(conf.env['INCLUDES_WCSLIB'][0])
        f90dir = conf.bldnode.make_node('include').mkdir()
        f77files = f77dir.ant_glob('*inc')
        for f77file in f77dir.ant_glob('*inc'):
            f90file = conf.bldnode.make_node('include/'+f77file.name)
            f77_to_f90(f77file, f90file)
        conf.define('HAVE_WCSLIB', 1)

    conf.env.SHAREDIR = os.path.abspath(conf.env.PYTHONDIR + '/../../../share')
    conf.define(conf.env.FC_NAME, 1)
    conf.define('TAMASIS_DIR', os.path.join(conf.env.SHAREDIR, 'tamasis'))
    conf.define('PRECISION_REAL', int(conf.options.precision_real))

    # Write the file .f2py_f2cmap which maps the real or complex type parameter
    # 'p' to a C type
    p = conf.options.precision_real
    ctype = {'4' : 'float',
             '8' : 'double',
             '16': 'long_double'}[p]
    cctype = {'4' : 'complex_float',
              '8' : 'complex_double',
              '16': 'complex_long_double'}[p]
    os.system("echo {\\'real\\':{\\'p\\':\\'" + ctype + \
        "\\'}, \\'complex\\':{\\'p\\':\\'" + cctype + "\\'}} > %s/.f2py_f2cmap"
        % out)



#
# Build & Installation
#

def build(bld):
    global libraries

    # Build static libraries libtamasis*.a
    for subdir in subdirs:
        files = bld.srcnode.ant_glob(subdir+'/src/module_*.f')
        if subdir == 'core' and 'WCSLIB' not in libraries:
            files = [f for f in files if 'wcslib' not in f.abspath()]
        bld(features = 'fc fcstlib',
            source   = files,
            target   = 'tamasis' + subdir,
            includes = 'include',
            use      = libraries)
    bld.add_group()
    
    # Python extension tamasisfortran.so
    target = 'tamasisfortran.so'
    source = [bld.srcnode.find_node('%s/src/tamasisfortran_%s.f90' % (s,s))
              for s in subdirs]
    additional_interfaces = ['math', 'operators', 'pointing', 'processing',
                             'wcsutils']
    if bld.env.HAVE_MPI:
        additional_interfaces += [ 'mpi' ]
    for a in additional_interfaces:
        source += bld.srcnode.ant_glob('core/src/tamasisfortran_core_{0}.f90'.
                                       format(a))

    #XXX this should be a Task...
    cmd = ''
    if bld.env.HAVE_MPI:
        cmd += 'OMPI_FC=${OMPI_FC} '
    cmd = '${F2PY} --fcompiler=${F2PYFCOMPILER} --f90exec=${FC}'
    cmd += ' --f90flags="${FCFLAGS}'
    cmd += ' ${FCFLAGS_OPENMP}"' if 'OPENMP' in libraries else '"'
    cmd += ' --quiet' if bld.options.verbose == 0 else ''
    cmd += ' -DF2PY_REPORT_ON_ARRAY_COPY=1'
    cmd += ' -m tamasisfortran -c ' + string.join(map(str, [s.abspath() for s in source]))
    cmd += ' ' + string.join([bld.env.fcstlib_PATTERN%('tamasis'+s) for s in reversed(subdirs)])
    cmd += ' ${FCINCPATH_ST:include} ' + string.join(['${FCINCPATH_ST:INCLUDES_%s}' % u for u in libraries])
    cmd += ' ' + string.join(['${FCLIBPATH_ST:LIBPATH_%s}' % u for u in libraries])
    cmd += ' ' + string.join(['${FCLIB_ST:LIB_%s}' % u for u in libraries])
    cmd += ' ${DEFINES_ST:DEFINES}'

    bld(rule=cmd,
        color='YELLOW',
        source=source + [bld.env.fcstlib_PATTERN%('tamasis' + s)
                         for s in subdirs],
        target=target,
        use=libraries)

    # include configuration information in Python files
    bld(features= 'subst', # feature 'subst' overrides source/target processing
        source='core/src/var.py.in',
        target='core/src/var.py',
        VERSION=VERSION, # variable to use in the substitution
        HAVE_MPI=str(bld.env.HAVE_MPI == 1),
        HAVE_WCSLIB=str(bld.env.HAVE_WCSLIB == 1))

    # Installation
    files = bld.path.ant_glob('*/src/*py') + \
            [bld.bldnode.find_or_declare('core/src/var.py')]
    pyfiles = [f for f in files if not os.access(f.abspath(), os.X_OK)]
    execfiles = [f for f in files if os.access(f.abspath(), os.X_OK)]
    bld.install_files('${PYTHONDIR}/tamasis', pyfiles + ['tamasisfortran.so'])
    bld.install_files('${PYTHONDIR}/tamasis', execfiles, chmod=0o755)

    for f in execfiles:
        bf = os.path.basename(f.abspath())
        bld.symlink_as('${BINDIR}/' + bf.replace('.py', ''),
                       bld.env['PYTHONDIR']+'/tamasis/' + bf)

    for subdir in subdirs:
        node = bld.srcnode.find_node(subdir+'/src')
        if node is not None:
            bld.install_files('${SHAREDIR}/tamasis/'+subdir,
                              node.ant_glob('*.cfg'))
        node = bld.srcnode.find_node(subdir+'/data')
        if node is not None:
            bld.install_files('${SHAREDIR}/tamasis/'+subdir, node.ant_glob('*'))
        node = bld.srcnode.find_node(subdir+'/example')
        if node is not None:
            bld.install_files('${SHAREDIR}/tamasis/'+subdir, node.ant_glob('*'))

    bld.install_files('${LIBDIR}/jython/tamasishcss', bld.srcnode.ant_glob(
                      'pacs/hcss/tamasishcss/*py'))



#
# test, test-fortran, test-python
#

class test(BuildContext):
    """run Fortran and Python test suites"""
    cmd = 'test'
    fun = 'test_fun'

def test_fun(ctx):
    tests = ['test-fortran', 'test-python']
    if 'MPI' in libraries:
        tests += ['test-mpi']
    Options.commands = tests + Options.commands

class test_fortran(BuildContext):
    """run Fortran test suite"""
    cmd = 'test-fortran'
    fun = 'test_fortran_fun'

def test_fortran_fun(bld):
    build(bld)
    tamasislib = ' '.join(['tamasis' + s for s in reversed(subdirs)])
    for subdir in subdirs:
        files = bld.srcnode.ant_glob(subdir+'/test/test_*.f')
        for file in files:
            if 'WCSLIB' not in libraries and 'wcslib' in str(file):
                continue
            bld(features = 'fc fcprogram',
                source   = file,
                target   = os.path.splitext(str(file))[0],
                includes = '. include',
                libpath  = '.',
                use      = ' '.join(libraries) + ' ' + tamasislib)
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

def test_python_fun(bld):
    if bld.env.IPYTHON_VERSION >= '0.12':
        options = '--no-confirm-exit --i '
    else:
        options = '-noconfirm_exit '
#    bld(rule='${IPYTHON} ' + options + bld.path.find_node(
#        'core/test/test_broken_locale.py').abspath(), always=True)

    for subdir in subdirs:
        if subdir == 'pacs':
            continue
        files = bld.path.ant_glob(subdir + '/test/test_*.py')
        for file in files:
            file = file.abspath()
            if 'test_display' in file: continue
            if 'test_mpi' in file: continue
            bld(rule='${NOSETESTS} ' + file + (' --nocapture'
                if bld.options.verbose else ''), always=True)
            bld.add_group()
    if 'pacs' in subdirs:
        test_pacs_fun(bld)

class test_display(BuildContext):
    """run Python display suite"""
    cmd = 'test-display'
    fun = 'test_display_fun'

def test_display_fun(bld):
    for subdir in subdirs:
        files = bld.path.ant_glob(subdir + '/test/test_display*.py')
        for file in files:
            file = file.abspath()
            bld(rule='${PYTHON} ' + file + (' > /dev/null'
                if not bld.options.verbose else ''), always=True)
            bld.add_group()

class test_pacs(BuildContext):
    """run PACS test suite"""
    cmd = 'test-pacs'
    fun = 'test_pacs_fun'

def test_pacs_fun(bld):
    subdir = 'pacs'
    files = bld.path.ant_glob(subdir+'/test/test_*.py')
    for file in files:
        file = file.abspath()
        if 'test_display' in file: continue
        if bld.env.HAVE_MPI and 'test_mpi' in file: continue
        if '_nl_' in file: continue
        bld(rule='${NOSETESTS} ' + file + (' --nocapture'
            if bld.options.verbose else ''), always=True)
        bld.add_group()

class test_mpi(BuildContext):
    """run MPI test suite"""
    cmd = 'test-mpi'
    fun = 'test_mpi_fun'

def test_mpi_fun(bld):

    if not bld.env.HAVE_MPI:
        Logs.pprint('RED', 'The package is not configured to use MPI. Use --ena'
                    'ble-mpi option.')
        return

    if bld.options.np is None:
        ns = 2**np.arange(int(np.floor(np.log2(multiprocessing.cpu_count())))+1)
    else:
        ns = [ int(bld.options.np) ]
    for subdir in subdirs:
        files = bld.path.ant_glob(subdir+'/test/test_mpi*.py')
        for file in files:
            for n in ns:
                bld(rule='export OMP_NUM_THREADS=1; mpirun -n ' + str(n) +
                    ' ${NOSETESTS} ' + file.abspath() + (' > /dev/null' if bld.
                    options.verbose == 0 else '') + '; unset OMP_NUM_THREADS',
                    always=True)
                bld.add_group()



#
# loc, loc-fortran, loc-python
#

def loc(ctx):
    """
    Python/Fortran/all lines of code.

    """
    print('Python:')
    os.system('find . \( -name ".waf*" -prune -or -name ".bento*" -prune -or -name build -prune -or -name .git -prune -or -name include -prune -or -name "wscript" -or -name "*py" \) -and -not -type l -and -not -type d | xargs cat | wc -l')
    print('Fortran:')
    os.system('find . \( -name ".waf*" -prune -or -name ".bento*" -prune -or -name build -prune -or -name .git -prune -or -name include -prune -or -name "*.f" -or -name "*.f90" \) -and -not -type l -and -not -type d | xargs cat | wc -l')
    print('Total:')
    os.system('find . \( -name ".waf*" -prune -or -name ".bento*" -prune -or -name build -prune -or -name .git -prune -or -name include -prune -or -name "wscript" -or -name "*.f" -or -name "*.f90" -or -name "*py" \) -and -not -type l -and -not -type d | xargs cat | wc -l')


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

def check_wcslib_external(env):
    wme = 'WCSLIB_MISSING_EXTERNAL'
    i = (i for i in range(len(env.DEFINES)) if env.DEFINES[i].startswith('WCSLIB_VERSION')).next()
    version = env.DEFINES[i].split('=')[1].replace('"', '')
    env.DEFINES[i] = wme + '=' + str(int(version < '4.5'))
    env.define_key[i] = wme

def _get_version(progname, flag=None):
    return '0.12'
    p = subprocess.Popen([progname, '--version'], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if stderr != '' and flag is not None:
        p = subprocess.Popen([progname, '-Version'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
    if stderr != '':
        raise RuntimeError(stderr)
    return stdout[:-1]

def _check_git_version(ctx):
    global VERSION
    cmd = 'git log HEAD^..'.split()
    # Python 2.7
    # line = subprocess.check_output(cmd).split('\n')[0]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    (stdout, stderr) = p.communicate()
    commit = stdout.split('\n')[0].split()[1]
    tagdir = os.path.join(ctx.path.abspath(), '.git', 'refs', 'tags')
    tags = os.listdir(os.path.join(ctx.path.abspath(), '.git', 'refs', 'tags'))
    for tag in os.listdir(tagdir):
        if commit == open(os.path.join(tagdir,tag)).read()[0:-1]:
            VERSION = tag
            if VERSION[0] == 'v': VERSION = VERSION[1:]
            break
    else:
        dev = '-dev' in VERSION
        VERSION = VERSION.replace('-dev', '') + '.' + commit[0:8] + \
                  ('-dev' if dev else '')

def dist(ctx):
    _check_git_version(ctx)
