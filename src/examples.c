#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "exhaustive.h"

int main()
{
    long *wd;               // Pointer to results
    clock_t start, end;     // Elapsed time measuring

    // EX 1: Code over prime finite field ==> use "wdp"
    // Cyclic code of length 12 and dimension 7 over GF(3)
    // ===================================================

    // Generator matrix (as vector)
    int code1[] = { 1, 0, 0, 0, 2, 0, 2, 0, 0, 0, 1, 0,
                    0, 1, 0, 1, 0, 0, 0, 2, 0, 2, 0, 0,
                    2, 0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0,
                    0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 2,
                    1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
                    2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    // Starting vector (unless dividing code for parallel processing,
    // this should be the null vector of appropriate size)
    // int nullvec12[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    // Parameters
    int p = 3;      // Field characteristic
    int n = 12;     // Code length
    int k = 7;      // Code dimension

    // Call weight distribution exhaustive calculator
    start = clock();
    // wd = wdp(p, n, k, code1, nullvec12, 0);
    wd = para_wdp(p, n, k, code1, 0);
    end = clock();
    if (wd == NULL)
        exit(EXIT_FAILURE);

    // Show results
    printf("The linear code generated by the matrix\n");
    for (int i = 0; i < k; i++) {
        printf("\n\t");
        for (int j = 0; j < n; j++)
            printf("%d  ", code1[i * n + j]);
    }
    printf("\n\nover GF(3) has the following weight distribution:\n\n");
    printf("\tweight   # of words\n\t------   ----------\n");
    for (int w = 0; w <= n; w++)
        if (wd[w] != 0)
            printf("\t  %2d        %4ld\n", w, wd[w]);
    printf("\t-------------------\n");
    printf("(%d words computed in %.3lfs)\n\n", wd[n + 1],
           (double)(end - start) / CLOCKS_PER_SEC);

    // Free memory allocated for results
    free(wd);


    // EX 2: Code over finite field extension ==> use "wde"
    // Cyclic code of length 14 and dimension 9 over GF(2^2)
    // =====================================================

    // Generator matrix (as vector)
    int code2[] = { 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                    1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                    1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                    1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
                    1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
                    1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0,
                    0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0,
                    0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0,
                    0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                    0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                    0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
                    0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
                    0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1,
                    0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1,
                    0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1,
                    0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0 };

    // Starting vector (unless dividing code for parallel processing,
    // this should be the null vector of appropriate size)
    int nullvec28[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    // Parameters
    p = 2;      // Field characteristic
    int e = 2;  // Extension degree
    n = 28;     // Code length over GF(2)
    k = 18;     // Code dimension over GF(2)

    // Call weight distribution exhaustive calculator
    start = clock();
    wd = wde(p, e, n, k, code2, nullvec28, 0);
    end = clock();
    if (wd == NULL)
        exit(EXIT_FAILURE);

    // Show results
    printf("The linear code generated by the matrix\n");
    for (int i = 0; i < k / e; i++) {
        printf("\n    ");
        for (int j = 0; j < n; j += 2)
            printf("%d+%dz ", code2[i * n + j], code2[i * n + j + 1]);
    }
    printf("\n\nover GF(4)=GF(2)[z] has the following weight distribution:\n\n");
    printf("\tweight   # of words\n\t------   ----------\n");
    for (int w = 0; w <= n / e; w++)
        if (wd[w] != 0)
            printf("\t  %2d        %5ld\n", w, wd[w]);
    printf("\t-------------------\n\n");

    printf("Seen as a code over GF(2), it has the following weight distribution:\n\n");
    printf("\tweight   # of words\n\t------   ----------\n");
    for (int w = 0; w <= n; w++)
        if (wd[n / e + 1 + w] != 0)
            printf("\t  %2d        %5ld\n", w, wd[n / e + 1 + w]);
    printf("\t-------------------\n");


    printf("(%d words computed in %.3lfs)\n\n", wd[n / e + n + 2],
           (double)(end - start) / CLOCKS_PER_SEC);

    // Free memory allocated for results
    free(wd);
}
