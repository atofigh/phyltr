import os, sys

# We want to seperate build directory from source directory
BuildDir('build', 'src')

# Set up the command line options of the SConstruct file
opts = Options()
opts.AddOptions(
    BoolOption('profile', 'compile with -pg -O2 and without -g', 'False'),
    BoolOption('debug', 'compile with -g without optimizations', 'False'),
    BoolOption('release', 'compile with -O2 without -g', 'False')
    )
env = Environment(options = opts)
Help(opts.GenerateHelpText(env))

# Ensure that at most one of profile, debug, and release has been set.
if env['debug'] + env['profile'] + env['release'] > 1:
    print >> sys.stderr, "At most one of 'debug', 'profile', and 'release' may be set"

# Set the compiler flags
if env['debug']:
    env.Append(CCFLAGS = ' -g')
elif env['profile']:
    env.Append(CCFLAGS = ' -pg -O2', LINKFLAGS = ' -pg')
elif env['release']:
    env.Append(CCFLAGS = ' -O2')
else:
    env.Append(CCFLAGS = ' -g -O2')
env.Append(CCFLAGS = '-Wall -pedantic -ansi -Wno-long-long')

# Make sure you are using SYSDIR if it is set
env.Replace(SYSDIR = os.environ['SYSDIR'])

if len(env['SYSDIR']) != 0:
    env.Append(CPPPATH = '$SYSDIR/include',
               LIBPATH = '$SYSDIR/lib')
    env['ENV']['LD_LIBRARY_PATH'] = env.subst('$SYSDIR/lib')


# Declare the programs and their source files
common_objs = env.Object('build/common.cc')

phyltr_dp = env.Program(target = 'bin/phyltr-dp',
                        source = ['build/phyltr-dp.cc'] + common_objs,
                        LIBS = ['NHparser', 'boost_program_options'])

phyltr_fpt = env.Program(target = 'bin/phyltr-fpt',
                         source = ['build/phyltr-fpt.cc'] + common_objs,
                         LIBS = ['NHparser', 'boost_program_options'])

phyltr_gen_stree = env.Program(target = 'bin/phyltr-gen-stree',
                               source = ['build/phyltr-gen-stree.cc'] + common_objs,
                               LIBS = ['boost_program_options'])

phyltr_gen_gtree = env.Program(target = 'bin/phyltr-gen-gtree',
                               source = ['build/phyltr-gen-gtree.cc'] + common_objs,
                               LIBS = ['NHparser', 'boost_program_options'])

# Declare a test target that runs the test script(s)
test_target = env.Command('test', None, '@cd tests; python test.py')
env.Depends(test_target, [phyltr_dp, phyltr_fpt])

# Set the default targets
env.Default([phyltr_fpt, phyltr_dp, phyltr_gen_stree, phyltr_gen_gtree])
