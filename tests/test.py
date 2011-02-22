#! /usr/bin/env python

# Copyright (C) 2010, 2011 Ali Tofigh
#
# This file is part of PhylTr, a package for phylogenetic analysis
# using duplications and transfers.
#
# PhylTr is released under the terms of the license contained in the
# file LICENSE.

import sys, os, os.path, tempfile

if len(sys.argv) != 3:
    print >> sys.stderr, "Usage: python test.py <scripts dir> <bin dir>"
    sys.exit(1)

scriptsdir = sys.argv[1]
bindir = sys.argv[2]
sys.path.append(scriptsdir)

from phyltr import *

def create_and_write_random_input(s_size, g_size, s_file, g_file, sigma_file):
    (s, g, sigma) = create_random_input(s_size, g_size)
    s_file.truncate(0)
    s_file.seek(0)
    g_file.truncate(0)
    g_file.seek(0)
    sigma_file.truncate(0)
    sigma_file.seek(0)
    
    print >> s_file, s + ";"
    print >> g_file, g + ";"
    for m in sigma:
        print >> sigma_file, m[0], m[1]
        
    s_file.flush()
    g_file.flush()
    sigma_file.flush()

    return (s, g, sigma)
    

s_file = tempfile.NamedTemporaryFile()
g_file = tempfile.NamedTemporaryFile()
sigma_file = tempfile.NamedTemporaryFile()
res_dp_file = tempfile.NamedTemporaryFile()
res_fpt_file = tempfile.NamedTemporaryFile()

file_names = ' '.join([s_file.name, g_file.name, sigma_file.name])

print " Test: phyltr-event-combinations with small random tree ".center(80, '=')
for s_size in range(2, 15):
    for g_size in range(s_size, 20):
        for i in range(10):
            (s, g, sigma) = create_and_write_random_input(s_size, g_size,
                                                          s_file, g_file, sigma_file)

            infile = os.popen(os.path.join(bindir, 'phyltr-event-combinations ') + file_names, 'r')
            combo_lines = infile.readlines()
            infile.close();

            combo_lines = [l.split() for l in combo_lines]
            plus_lines = [l for l in combo_lines if l[0] == '+']

            # For each event pair that is optimal for more than one
            # value, test a cost in the middle of the range to ensure
            # that that is the only event pair that is optimal there.
            for l in plus_lines:
                dups = int(l[1])
                transfers = int(l[2])
                lower_numer = int(l[3].split('/')[0])
                upper_numer = int(l[4].split('/')[0])
                denominator = int(l[3].split('/')[1])
                # set the cost of duplication and transfer as integers
                cost_dup = lower_numer + upper_numer
                cost_transfer = denominator - lower_numer + denominator - upper_numer

                infile = os.popen(os.path.join(bindir, 'phyltr-dp') + ' -d ' + str(cost_dup)
                                  + ' -t ' + str(cost_transfer) + ' '
                                  + file_names, 'r')
                dp_lines = [l.split() for l in infile if len(l.strip()) != 0]
                infile.close()

                transfer_lines = [l for l in dp_lines if l[0] == 'Transfer']
                dup_lines = [l for l in dp_lines if l[0] == 'Duplications:']
                for i in range(len(dup_lines)):
                    if len(transfer_lines[i]) - 2 != transfers or len(dup_lines[i]) - 1 != dups:
                        print s
                        print g
                        for a in sigma:
                            print a[0], a[1]
                        print 'dup cost:\t', cost_dup
                        print 'transfer cost:\t', cost_transfer
                        print
                        raise RuntimeError, "Test Failed"

            # For each lower and upper limit of a cost range, check
            # that exactly the right event pairs are optimal

            # events[r], where r is a rational number, is a set that
            # contains all event pairs optimal when D = r and T = 1-r
            # according to phyltr-event-combinations
            events = {}
            for l in combo_lines:
                if l[0] == '-':
                    continue
                (numerator, denominator) = l[3].split('/')
                if numerator != denominator and numerator != 0:
                    lower = Rational(int(numerator), int(denominator))
                    a_set = events.get(lower, set())
                    a_set.add( (int(l[1]), int(l[2])) )
                    events[lower] = a_set

                (numerator, denominator) = l[4].split('/')
                if numerator != denominator and numerator != 0:
                    upper = Rational(int(numerator), int(denominator))
                    a_set = events.get(upper, set())
                    a_set.add( (int(l[1]), int(l[2])) )
                    events[upper] = a_set

            # For each cost, test if phyltr-dp agrees with
            # phyltr-event-combinations
            for (key, events_set) in events.iteritems():
                cost_dup = key.numerator
                cost_transfer = key.denominator - key.numerator
                infile = os.popen(os.path.join(bindir, 'phyltr-dp') + ' -d ' + str(cost_dup)
                                  + ' -t ' + str(cost_transfer) + ' '
                                  + file_names, 'r')
                
                dp_lines = [l.split() for l in infile if len(l.strip()) != 0]
                infile.close()

                transfer_lines = [l for l in dp_lines if l[0] == 'Transfer']
                dup_lines = [l for l in dp_lines if l[0] == 'Duplications:']
                dp_set = set()
                for i in range(len(dup_lines)):
                    dp_set.add( (len(dup_lines[i]) - 1, len(transfer_lines[i]) - 2) )
                    
                if events_set != dp_set:
                    print s
                    print g
                    for a in sigma:
                        print a[0], a[1]
                    print 'dup cost:\t', cost_dup
                    print 'transfer cost:\t', cost_transfer
                    print 'events:\t', events_set
                    print 'dp:\t', dp_set
                    print
                    raise RuntimeError, "Test Failed"
                    
                
    sys.stderr.write('.')
print
print " PASSED ".center(80, '=')


print " Test: small random trees with phyltr-dp and phyltr-fpt ".center(80, '=')
for s_size in range(2, 15):
    for g_size in range(s_size, 15):
        for i in range(10):
            (s, g, sigma) = create_and_write_random_input(s_size, g_size,
                                                          s_file, g_file, sigma_file)

            infile = os.popen(os.path.join(bindir, 'phyltr-dp ') + '-c ' + file_names, 'r')
            opt_cost = int(infile.read())
            infile.close()

            infile = os.popen(os.path.join(bindir, "phyltr-dp ") + file_names)
            res_dp_file.writelines(infile.readlines())
            infile.close()
            res_dp_file.flush()

            infile = os.popen(os.path.join(bindir, "phyltr-fpt ") + file_names + " " +
                              str(opt_cost) + " " + str(opt_cost))
            res_fpt_file.writelines(infile.readlines())
            infile.close()
            res_fpt_file.flush()

            diff = os.system("diff " +
                             res_dp_file.name + " " +
                             res_fpt_file.name)
            if diff != 0:
                print s
                print g
                for a in sigma:
                    print a[0], a[1],
                print
                raise RuntimeError, "Test Failed"
    sys.stderr.write('.')
print
print " PASSED ".center(80, '=')
