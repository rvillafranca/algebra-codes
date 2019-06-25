#include <stdlib.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>

long *distr = NULL;     // Weight distributions, used as return from the
                        // functions "wdp" and "wde". Should be freed after
                        // used.

/******************************************************************************

    Weight distribution over prime finite fields

    Parameters:
      p     field characteristic;
      n     code length;
      k     code dimension;
      B     base of the code / generator matrix:
            base vectors (matrix rows) must be concatenated as a single
            vector of length n*k;
      c0    starting vector (to divide de code in cosets);
      min   stop criterium:
            if this parameter is positive the execution is interrupted
            if a word of weight less than its value is found.

    Return:
      a vector "distr" of length n+2:
      - for 0 <= i <= n, distr[i] contains the number of codewords of
        weight i found;
      - distr[n+1] contains the total number of codewords computed.


******************************************************************************/
long *wdp(int p, int n, int k, int *B, int *c0, int min)
{
    int *c, *a, *s;
    int i, j, w;

    // Weight distribution and number of words computed
    distr = (long*)calloc(n + 2, sizeof(long));

    // Codewords computed
    c = (int*)malloc(n * sizeof(int));
    for (i = 0; i < n; i++)
        c[i] = c0[i];

    // Coefficient and upward/downward vectors
    a = (int*)calloc(k, sizeof(int));
    s = (int*)calloc(k, sizeof(int));

    i = 0;
    while (i < k) {
        w = 0;
        for (j = 0; j < n; j++)
            if ((c[j] % p) != 0)
                w++;
        distr[w]++;
        distr[n + 1]++;

        if (min && w > 0 && w < min)
            break;

        while (i < k) {
            if (s[i] && a[i] > 0) {
                a[i]--;
                for (j = 0; j < n; j++)
                    c[j] = (c[j] - B[i * n + j]);
                i = 0;
                break;
            }else if (!s[i] && a[i] < (p - 1)) {
                a[i]++;
                for (j = 0; j < n; j++)
                    c[j] = (c[j] + B[i * n + j]);
                i = 0;
                break;
            }else {
                s[i] = !s[i];
                i++;
            }
        }
    }

    free(s);
    free(a);
    free(c);

    return distr;
}


/******************************************************************************

    Weight distribution over non-prime finite fields, i.e., GF(p^e)

    Parameters:
      p     field characteristic;
      e     extension degree;
      n     code length (as seen over the prime field GF(p));
      k     code dimension (as seen over the prime field GF(p));
      B     base of the code / generator matrix:
            - base vectors (matrix rows) must be concatenated as a single
              vector of length n*k;
            - each coefficient in GF(p^e) should be split into e consecutive
              coefficients in GF(p);
      c0    starting vector (to divide de code in cosets);
      min   stop criterium:
            if this parameter is positive the execution is interrupted
            if a word of weight less than its value is found.

    Return:
      a vector "distr" of length n/e + n + 3:
      - for 0 <= i <= n/e, distr[i] contains the number of codewords of
        weight i found, as seen over GF(p^e);
      - for 0 <= i <= n, distr[n/e + 1 + i] contains the number of codewords of
        weight i found, as seen over GF(p);
      - distr[n/e + n + 2] contains the total number of codewords computed.


******************************************************************************/
long *wde(int p, int e, int n, int k, int *B, int *c0, int min)
{
    int *c, *a, *s;
    int i, j, wp, we, z, ne;

    if (n % e != 0)
        return NULL;
    ne = n / e;

    // Weight distribution and number of words computed
    distr = (long*)calloc((ne + n + 3), sizeof(long));

    // Codewords computed
    c = (int*)malloc(n * sizeof(int));
    for (i = 0; i < n; i++)
        c[i] = c0[i];

    // Coefficient and upward/downward vectors
    a = (int*)calloc(k, sizeof(int));
    s = (int*)calloc(k, sizeof(int));

    i = 0;
    while (i < k) {
        wp = we = z = 0;
        for (j = 0; j < n; j++) {
            if ((c[j] % p) != 0) {
                wp++;
                z = 1;
            }
            if ((j + 1) % e == 0) {
                we += z;
                z = 0;
            }
        }
        distr[we]++;
        distr[ne + 1 + wp]++;
        distr[ne + n + 2]++;

        if (min && we > 0 && we < min)
            break;

        while (i < k) {
            if (s[i] && a[i] > 0) {
                a[i]--;
                for (j = 0; j < n; j++)
                    c[j] = (c[j] - B[i * n + j]);
                i = 0;
                break;
            } else if (!s[i] && a[i] < (p - 1)) {
                a[i]++;
                for (j = 0; j < n; j++)
                    c[j] = (c[j] + B[i * n + j]);
                i = 0;
                break;
            } else {
                s[i] = !s[i];
                i++;
            }
        }
    }

    free(s);
    free(a);
    free(c);

    return distr;
}


/******************************************************************************

    Free memory allocated for results by wdp, wde or para_wd.

******************************************************************************/
void free_memory()
{
    if (distr != NULL) {
        free(distr);
        distr = NULL;
    }
}


/******************************************************************************

    Weight distribution over finite fields (parallel processing)

    Parameters:
      p     field characteristic;
      e     extension degree;
      n     code length (as seen over the prime field GF(p));
      k     code dimension (as seen over the prime field GF(p));
      B     base of the code / generator matrix:
            - base vectors (matrix rows) must be concatenated as a single
              vector of length n*k;
            - each coefficient in GF(p^e) should be split into e consecutive
              coefficients in GF(p);
      min   stop criterium:
            if this parameter is positive the execution is interrupted
            if a word of weight less than its value is found.

    Return:
    If e == 1:
      a vector "distr" of length n+2:
      - for 0 <= i <= n, distr[i] contains the number of codewords of
        weight i found;
      - distr[n+1] contains the total number of codewords computed.
    If e > 1:
      a vector "distr" of length n/e + n + 3:
      - for 0 <= i <= n/e, distr[i] contains the number of codewords of
        weight i found, as seen over GF(p^e);
      - for 0 <= i <= n, distr[n/e + 1 + i] contains the number of codewords of
        weight i found, as seen over GF(p);
      - distr[n/e + n + 2] contains the total number of codewords computed.

******************************************************************************/
long *para_wd(int p, int e, int n, int k, int *B, int min)
{
    pid_t pid;
    int a, *c0, N;
    long *partial_wd, *aux;

    N = (e == 1) ? (n + 1) : (n / e + n + 2);

    // Auxiliary memory shared by processes
    aux = (long*)mmap(NULL, (N + 1) * sizeof(long), PROT_READ | PROT_WRITE,
                      MAP_SHARED | MAP_ANONYMOUS, -1, 0);

    c0 = (int*)malloc(n * sizeof(int));
    for (a = 0; a < p; a++) {
        pid = fork();
        if (pid == 0) {
            for (int i = 0; i < n; i++)
                c0[i] = B[i] * a;
            if (e == 1)
                partial_wd = wdp(p, n, k - 1, &B[n], c0, min);
            else
                partial_wd = wde(p, e, n, k - 1, &B[n], c0, min);
            for (int w = 0; w <= N; w++)
                aux[w] += partial_wd[w];
            free(partial_wd);
            exit(EXIT_SUCCESS);
        }
    }

    // Wait for children to finish their computation
    while (wait(NULL) > 0);

    // Copy shared memory to new allocated vector
    distr = (long*)calloc(N + 1, sizeof(long));
    for (int w = 0; w <= N; w++)
        distr[w] = aux[w];

    // Release shared memory
    munmap(aux, (N + 1) * sizeof(long));

    return distr;
}


/* TODO:
    interprocess communication for 'min';
    Galois rings routines;
*/
