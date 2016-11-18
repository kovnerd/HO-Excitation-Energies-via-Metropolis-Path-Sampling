/*
 * Allocate (and free) two-dimensional array.
 * 
 * See http://c-faq.com/aryptr/dynmuldimary.html for a discussion.
 */

#include <stdlib.h>
#include "matrixmem.h"

int **matrix_allocate_int (int lx, int ly)
{
    int **s = (int **) malloc((unsigned long) lx * sizeof(int *));
    if (s == NULL)
    {
        return s;
    }

    s[0] = (int *) malloc((unsigned long) (lx * ly) * sizeof(int));
    if (s[0] == NULL)
    {
        return NULL;
    }

    for (int i = 1; i < lx; i++)
    {
        s[i] = s[0] + i*ly;
    }
    return s;
}

double **matrix_allocate_double (int lx, int ly)
{
    double **s = (double **) malloc((unsigned long) lx * sizeof(double *));
    if (s == NULL)
    {
        return s;
    }

    s[0] = (double *) malloc((unsigned long) (lx * ly) * sizeof(double));
    if (s[0] == NULL)
    {
        return NULL;
    }

    for (int i = 1; i < lx; i++)
    {
        s[i] = s[0] + i*ly;
    }
    return s;
}

void matrix_free_int (int **s)
{
    free (s[0]);
    free (s);
}

void matrix_free_double (double **s)
{
    free (s[0]);
    free (s);
}
