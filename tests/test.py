#! /usr/bin/env python

import sys, os, tempfile

sys.path.append('../scripts')

from phyltr import *


# We do some testing of small trees

print " Test: small random trees ".center(80, '=')

s_file = tempfile.NamedTemporaryFile()
g_file = tempfile.NamedTemporaryFile()
sigma_file = tempfile.NamedTemporaryFile()
res_dp_file = tempfile.NamedTemporaryFile()
res_fpt_file = tempfile.NamedTemporaryFile()

file_names = ' '.join([s_file.name, g_file.name, sigma_file.name])

for i in range(10):
    for s_size in range(2, 10):
        for g_size in range(s_size, 10):
            (s, g, sigma) = create_random_input(s_size, g_size)
            s_file.truncate(0)
            s_file.seek(0)
            g_file.truncate(0)
            g_file.seek(0)
            sigma_file.truncate(0)
            sigma_file.seek(0)
            res_dp_file.truncate(0)
            res_dp_file.seek(0)
            res_fpt_file.truncate(0)
            res_fpt_file.seek(0)

            print >> s_file, s + ";"
            print >> g_file, g + ";"
            for m in sigma:
                print >> sigma_file, m[0], m[1]

            s_file.flush()
            g_file.flush()
            sigma_file.flush()

            infile = os.popen('../bin/phyltr-dp -c ' + file_names, 'r')
            opt_cost = int(infile.read())
            infile.close()

            infile = os.popen("../bin/phyltr-dp " + file_names)
            res_dp_file.writelines(infile.readlines())
            infile.close()
            res_dp_file.flush()

            infile = os.popen("../bin/phyltr-fpt " + file_names + " " +
                              str(opt_cost) + " " + str(opt_cost))
            res_fpt_file.writelines(infile.readlines())
            infile.close()
            res_fpt_file.flush()

            diff = os.system("diff " +
                             res_dp_file.name + " " +
                             res_fpt_file.name)
            if diff != 0:
                raise RuntimeError, "Test Failed"
    sys.stderr.write('.')
print


print " PASSED ".center(80, '=')

