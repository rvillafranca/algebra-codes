#!/usr/bin/env sage
# Encoding: utf-8

load('weightdistr.sage')

# Example 1 --------------------------------------------------------------------

# Generator Matrix
G = [[1, 0, 0, 0, 2, 0, 2, 0, 0, 0, 1, 0],
     [0, 1, 0, 1, 0, 0, 0, 2, 0, 2, 0, 0],
     [2, 0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0],
     [0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 2],
     [1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
     [2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0],
     [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]

p = 3       # Field characteristic
n = 12      # Code length
K = GF(p)

B = [[K(x) for x in v] for v in G]      # Convert G to a K-basis
wd = weightdist(p, n, B)

print('Distribuição de pesos:\n\n'
      '\tWeight  # words\n'
      '\t------  -------')
for w in range(n + 1):
    if (wd[0][w] != 0):
        print('\t{:5}   {:6}'.format(w, wd[0][w]))
print('\n\n{} words computed\n'.format(wd[2]))

# Example 2 --------------------------------------------------------------------

# Group algebra of C14 over GF(2^2)
q = 4
d = 2       # Degree of field extension
n = 14
K.<z> = GF(q)
G.<ag> = AbelianGroup([n])
KG = GroupAlgebra(G, K)

a = KG(ag^7)    # Generator of subgroup of order 2
g = KG(ag^8)    # Generator of subgroup of order 7

# Idempotents of KG
e1 = 1 + g + g^2 + g^4
e2 = 1 + g^3 + g^5 + g^6
e3 = 1 + g + g^2 + g^3 + g^4 + g^5 + g^6

# Cyclic code of length 14 over GF(2^2)
basis = [e1, g^3*e1, g^5*e1, (a - 1)*e1, (a - 1)* g^3*e1, (a - 1)*g^5*e1,
         (a - 1)*e2, (a - 1)* g^2*e2, (a - 1)*g^6*e2]

# Compute the weight distribution
B = [tuple(b.coefficient(ag^k) for k in range(n)) for b in basis]
wd = weightdist(q, n, B)

# wd[0] has the weight distribution of the code seen as over K
print('\n\nWeight distribution over GF(4):\n\n'
      '\tWeight  # words\n'
      '\t------  -------')
for w in range(n + 1):
    if (wd[0][w] != 0):
        print('\t{:5}   {:6}'.format(w, wd[0][w]))

# If K is not a prime field, wd[1] has the weight distribution of the code seen
# as over the prime field of K
print('\n\nWeight distribution over GF(2):\n\n'
      '\tWeight  # words\n'
      '\t------  -------')
for w in range(n * d + 1):
    if (wd[1][w] != 0):
        print('\t{:5}   {:6}'.format(w, wd[1][w]))

# wd[2] has the total number of computed words
print('\n\n{} words computed\n'.format(wd[2]))
