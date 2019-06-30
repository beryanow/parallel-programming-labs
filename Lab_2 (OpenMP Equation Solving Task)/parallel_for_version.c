#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>

#define N 12000

double *multiplying_matrixes(double **A, double *b) {
    double *final = (double *) malloc(N * sizeof(double));
    
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        double sum = 0;
#pragma omp parallel for
        for (int j = 0; j < N; j++) {
            sum += A[i][j] * b[i];
        }
        final[i] = sum;
    }
    return final;
}

double *multuplying_by_scalar(double *A, double t) {
    double *final = (double *) malloc(N * sizeof(double));
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        final[i] = A[i] * t;
    }
    return final;
}

double *matrixes_subtraction(double *A, double *b) {
    double *final = (double *) malloc(N * sizeof(double));
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        final[i] = A[i] - b[i];
    }
    return final;
}

double module(double *a) {
    double sum = 0;
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        sum += a[i] * a[i];
    }
    return sqrt(sum);
}

int main() {
    double **A = (double **) malloc(N * sizeof(double *));
    
    omp_set_num_threads(1);
    for (int i = 0; i < N; i++) {
        A[i] = (double *) malloc(N * sizeof(double));
    }
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j)
                A[i][j] = 1.0;
            else
                A[i][j] = 2.0;
        }
    }
    
    double *b = (double *) malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        b[i] = N + 1;
    }
    
    double *x = (double *) malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        x[i] = 0;
    }
    
    double e = 0.00001;
    double t = 0.0001;
    
    struct timeval tv1, tv2, dtv;
    struct timezone tz;
    gettimeofday(&tv1, &tz);
    
    while (module(matrixes_subtraction(multiplying_matrixes(A, x), b)) / module(b) > e) {
        x = matrixes_subtraction(x, multuplying_by_scalar(matrixes_subtraction(multiplying_matrixes(A, x), b), t));
    }
    
    gettimeofday(&tv2, &tz);
    dtv.tv_sec = tv2.tv_sec - tv1.tv_sec;
    dtv.tv_usec = tv2.tv_usec - tv1.tv_usec;
    if (dtv.tv_usec < 0) {
        dtv.tv_sec--;
        dtv.tv_usec += 1000000;
    }
    printf("%f seconds\n", ((double) dtv.tv_sec * 1000 + dtv.tv_usec / 1000) / 1000);
    
    //    for (int i = 0; i < N; i++) {
    //        printf("%f ", x[i]);
    //    }
    return 0;
}

