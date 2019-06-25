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
    c0 = (0,) * n*e

    B_t = ct.c_int * len(B)
    c0_t = ct.c_int * len(c0)
    if (e == 1):
        C.para_wdp.argtypes = [ct.c_int, ct.c_int, ct.c_int, B_t, ct.c_int]
        C.para_wdp.restype = ct.POINTER(ct.c_long)
        d = C.para_wdp(p[0], n, k, B_t(*B), minstop)
        ret_tuple = d[0:(n+1)], None, d[n+1]
    else:
        C.wde.argtypes = [ct.c_int, ct.c_int, ct.c_int, ct.c_int,
                          B_t, c0_t, ct.c_int]
        C.wde.restype = ct.POINTER(ct.c_long)
        d = C.wde(p[0], e, n*e, k, B_t(*B), c0_t(*c0), minstop)
        ret_tuple = d[0:(n+1)], d[(n+1):(n*(e+1))+2], d[n*(e+1)+2]

    C.free_memory()

    return ret_tuple
