import ctypes as ct

C = ct.CDLL("./exhaustive.so")

def weightdist(q, n, basis, minstop=0):
    p = prime_divisors(q)
    if (len(p) > 1):
        return None

    basis_p = list()
    e = log(q, p[0])
    for j in range(len(basis)):
        for i in range(e):
            b = [0] * n * e
            for k in range(n):
                if (e > 1):
                    x = basis[j][k].integer_representation()
                else:
                    x = int(basis[j][k])
                y = [0] * e
                for m in range(e):
                    x, y[m] = divmod(x, p[0])
                for m in range(e):
                    b[k*e + m] = y[(m + i) % e]
            basis_p.append(b)

    k = len(basis_p)
    B = tuple(int(basis_p[i][j]) for i in range(k) for j in range(n*e))

    B_t = ct.c_int * len(B)
    C.para_wd.argtypes = [ct.c_int, ct.c_int, ct.c_int, ct.c_int,
                          B_t, ct.c_int]
    C.para_wd.restype = ct.POINTER(ct.c_long)
    d = C.para_wd(p[0], e, n * e, k, B_t(*B), minstop)
    if (e == 1):
        ret_tuple = d[0:(n+1)], None, d[n+1]
    else:
        ret_tuple = d[0:(n+1)], d[(n+1):(n*(e+1))+2], d[n*(e+1)+2]

    C.free_memory()

    return ret_tuple
