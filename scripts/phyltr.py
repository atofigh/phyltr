"""Utility functions useful when working with the PhylTr-package

Description will follow later...
"""

import random as _random
import types as _types

_LEAF_MARK = 'X'
_alphabet = 'abcdefghijklmnopqrstuvwxyz'
_rng = _random.Random()
_Z_table = [0, 1, 1]

def _isqrt(n):
    xn = 1
    xn1 = (xn + n/xn)/2
    while abs(xn1 - xn) > 1:
        xn = xn1
        xn1 = (xn + n/xn)/2
    while xn1*xn1 > n:
        xn1 -= 1
    return xn1

def _W(x, i, j):
    assert 0 <= i <= j <= x-1, (i, j, x)

    return x*(x+1)/2 - (x-i)*(x-i + 1)/2 + j - i


def _inv_W(w, x):
    assert type(x) == _types.IntType or type(x) == _types.LongType
    assert x >= 1
    
    import math

    r = (2*x - 1)**2 - 8*(w-x+1)
    isqrt_r = _isqrt(r)
    s = (2*x + 1)**2 - 8*w
    isqrt_s = _isqrt(s)
    
    i1 = 2*x - 1 - isqrt_r
    if isqrt_r*isqrt_r != r:
        i1 -= 1
    i1 /= 2

    i2 = 2*x + 1 - isqrt_s
    if isqrt_s*isqrt_s != s:
        i2 -= 1
    i2 /= 2

    i = max(i1, i2)
    j = w - x*(x+1)/2 + (x-i)*(x-i + 1)/2 + i

    return i, j

def _compute_Z_table(n):
    assert type(n) == _types.IntType or type(n) == _types.LongType

    for k in range(len(_Z_table), n+1):
        Z_tmp = 0
        i = 1
        while 2*i < k:
            Z_tmp += _Z_table[i] * _Z_table[k-i]
            i += 1

        if i*2 == k:
            Z_tmp += _Z_table[i]*(_Z_table[i] + 1) / 2

        _Z_table.append(Z_tmp)
        

def _Z(n):
    assert type(n) == _types.IntType or type(n) == _types.LongType

    if n >= len(_Z_table):
        _compute_Z_table(n)

    return _Z_table[n]

def _inv_rank(rank, n):
    assert rank <= _Z(n), (n, rank, _Z(n))

    if n == 1:
        return _LEAF_MARK
    
    a = 0
    p = 0
    while p < rank:
        a += 1
        p += _Z(a)*_Z(n-a)
    p -= _Z(a)*_Z(n-a)

    V = rank - p - 1
    if 2*a == n:
        left, right = _inv_W(V, _Z(a))
    else:
        left = V / _Z(n-a)
        right = V % _Z(n-a)

#    print "rank:", rank, "\tn:", n, "\ta:", a, "\tp:", p, "\tV:", V, "\tleft:", left, "\tright:", right

    return "(" + _inv_rank(left + 1, a) + ", " + _inv_rank(right + 1, n-a) + ")"

    
def _get_random_tree(n):
    """Returns a rooted full binary tree of size n"""

    assert (type(n) == _types.IntType or type(n) == _types.LongType) and n >= 2

    rank = _rng.randint(1, _Z(n))
    return _inv_rank(rank, n)


def create_random_input(species, genes):
    def label(n):
        assert type(n) == _types.IntType or type(n) == _types.LongType

        size = len(_alphabet)
        l = ""
        while n >= 0:
            nr = n % size
            l += _alphabet[nr]
            n /= size
            n -= 1

        return l[-1::-1]

    if type(species) != _types.IntType and type(species) != _types.LongType:
        raise TypeError, "Expected an integer, got a " + str(type(species))
    if type(genes) != _types.IntType and type(genes) != _types.LongType:
        raise TypeError, "Expected an integer, got a " + str(type(species))
    if species > genes:
        raise ValueError, "The number of species exceeds the number of genes"
    if species < 2 or genes < 2:
        raise ValueError, "genes and species must be at least 2"

    
    gene_tree = _get_random_tree(genes)
    species_tree = _get_random_tree(species)

    species_parts = species_tree.split(_LEAF_MARK)
    gene_parts = gene_tree.split(_LEAF_MARK)

    labels = []
    for i in range(species):
        labels.append(label(i))

    new_species_parts = [species_parts[0]]
    for i in range(species):
        new_species_parts.append(labels[i])
        new_species_parts.append(species_parts[i+1])

    excess_labels = []
    for i in range(genes-species):
        excess_labels.append(_random.choice(labels))

    labels += excess_labels
    labels.sort()

    gene_labels = []
    sigma = []
    i = 0
    while i < len(labels):
        j = 1
        gene_labels.append(labels[i] + str(j))
        sigma.append((labels[i] + str(j), labels[i]))
        while i+1 < len(labels) and labels[i] == labels[i+1]:
            i += 1
            j += 1
            gene_labels.append(labels[i] + str(j))
            sigma.append((labels[i] + str(j), labels[i]))
        i += 1

    _random.shuffle(gene_labels)

    new_gene_parts = [gene_parts[0]]
    for i in range(genes):
        new_gene_parts.append(gene_labels[i])
        new_gene_parts.append(gene_parts[i+1])
        

    return "".join(new_species_parts), "".join(new_gene_parts), sigma

