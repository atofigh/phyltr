#! /usr/bin/env python

import sys
import os

sys.path.append('../scripts')

from phyltr import *


# We do some testing of small trees

print " Test: small random trees ".center(80, '=')

for s_size in range(2, 10):
    for g_size in range(s_size, 10):
        (s, g, sigma) = create_random_input(s_size, g_size)
        s_file = open("s", "w")
        g_file = open("g", "w")
        sigma_file = open("sigma", "w")

        print >> s_file, s + ";"
        print >> g_file, g + ";"
        for m in sigma:
            print >> sigma_file, m[0], m[1]

        s_file.close()
        g_file.close()
        sigma_file.close()
            
        infile = os.popen('../bin/phyltr-dp -c s g sigma', 'r')

        opt_cost = int(infile.read())
        infile.close()

        os.system("../bin/phyltr-dp s g sigma > res-dp")
        os.system("../bin/phyltr-fpt s g sigma "
                  + str(opt_cost) + " " + str(opt_cost) + " "
                  + "> res-fpt")

        diff = os.system("diff res-dp res-fpt")
        if diff != 0:
            raise RuntimeError, "Test Failed"


print " PASSED ".center(80, '=')

os.remove("s")
os.remove("g")
os.remove("sigma")
os.remove("res-dp")
os.remove("res-fpt")
