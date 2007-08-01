import os

# We want to seperate build directory from source directory
BuildDir('build', 'src')

# Set up the environment
env = Environment()
env.Append(CCFLAGS = '-g -Wall -pedantic -ansi')

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

env.Default([phyltr_fpt, phyltr_dp])

test_target = env.Command('test', None, '@cd tests; python test.py')
env.Clean([phyltr_dp, phyltr_fpt],
          Split('tests/s tests/g tests/sigma tests/res-fpt tests/res-dp'))
env.Depends(test_target, [phyltr_dp, phyltr_fpt])
