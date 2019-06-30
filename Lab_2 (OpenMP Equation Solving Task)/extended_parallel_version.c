#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>

#define N 12000

double *multiplying_matrixes(double **A, double *b, int k, int start) {
    double *final = (double *) malloc(k * sizeof(double));
    for (int i = 0; i < k; i++) {
        final[i] = 0;
    }
    int y = start;
    for (int i = 0; i < k; i++) {
        double sum = 0;
        for (int j = 0; j < N; j++) {
            sum += A[i][j] * b[y];
        }
        y++;
        final[i] = sum;
    }
    return final;
}

double *multuplying_by_scalar(double *A, double t, int k) {
    double *final = (double *) malloc(k * sizeof(double));
    for (int i = 0; i < k; i++) {
        final[i] = A[i] * t;
    }
    return final;
}

double *matrixes_subtraction(double *A, double *b, int k, int start) {
    double *final = (double *) malloc(k * sizeof(double));
    for (int i = 0; i < k; i++) {
        final[i] = 0;
    }
    
    int y = start;
    for (int i = 0; i < k; i++) {
        final[i] = A[i] - b[y];
        y++;
    }
    return final;
}

int main(int argc, char *argv[]) {
    double **A = (double **) malloc(N * sizeof(double *));
    
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
    
    double b_sqrt = sqrt((double) (N + 1) * (N + 1) * N);
    
    double *x = (double *) malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        x[i] = 0;
    }
    
    double e = 0.00001;
    double t = 0.0001;
    double k;
    int y1 = -1;
    
    int l = 1;
    int kint, kint1;
    double sum = 0;
    double s;
    
    struct timeval tv1, tv2, dtv;
    struct timezone tz;
    gettimeofday(&tv1, &tz);
    omp_set_num_threads(4);
#pragma omp parallel
    {
        int size = omp_get_num_threads();
        int rank = omp_get_thread_num();
        
        if (N < size) {
            printf("Reduce the number of processes\n");
            exit(0);
        } else {
            k = (double) N / size;
            
            kint = (int) k;
            kint1 = N - kint * (size - 1);
            
            if (k - (int) k > 0) {
                y1 = 0;
            }
        }
        
        double **A1 = (double **) malloc(kint * sizeof(double *));
        double **A2 = (double **) malloc(kint1 * sizeof(double *));
        
        for (int i = 0; i < kint; i++) {
            A1[i] = (double *) malloc(N * sizeof(double));
        }
        
        for (int i = 0; i < kint1; i++) {
            A2[i] = (double *) malloc(N * sizeof(double));
        }
        
        if (y1 != -1) {
            if (rank != size - 1) {
                for (int i = 0; i < kint; i++) {
                    for (int j = 0; j < N; j++) {
                        A1[i][j] = A[i + rank * kint][j];
                    }
                }
            } else {
                for (int i = 0; i < kint1; i++) {
                    for (int j = 0; j < N; j++) {
                        A2[i][j] = A[N - kint1 + i][j];
                    }
                }
            }
        } else {
            for (int i = 0; i < kint; i++) {
                for (int j = 0; j < N; j++) {
                    A1[i][j] = A[i + rank * kint][j];
                }
            }
        }
        
        double *A_new;
        if (y1 != -1) {
            if (rank != size - 1) {
                A_new = (double *) malloc(kint * sizeof(double));
            } else {
                A_new = (double *) malloc(kint1 * sizeof(double));
            }
        } else {
            A_new = (double *) malloc(kint * sizeof(double));
        }
        
        
        while (l == 1) {
            s = 0;
            if (y1 != -1) {
                if (rank != size - 1) {
                    int g1 = rank * kint;
                    A_new = multiplying_matrixes(A1, x, kint, g1);
                    A_new = matrixes_subtraction(A_new, b, kint, g1);
                    
                } else {
                    int g1 = N - kint1;
                    A_new = multiplying_matrixes(A2, x, kint1, g1);
                    A_new = matrixes_subtraction(A_new, b, kint1, g1);
                }
            } else {
                int g1 = rank * kint;
                A_new = multiplying_matrixes(A1, x, kint, g1);
                A_new = matrixes_subtraction(A_new, b, kint, g1);
            }
            
            if ((y1 != -1) && (rank == size - 1)) {
                for (int i = 0; i < kint1; i++) {
                    s += A_new[i] * A_new[i];
                }
            } else {
                for (int i = 0; i < kint; i++) {
                    s += A_new[i] * A_new[i];
                }
            }
            
            sum += s;
            
            if (y1 != -1) {
                if (rank != size - 1) {
                    A_new = multuplying_by_scalar(A_new, t, kint);
                } else {
                    A_new = multuplying_by_scalar(A_new, t, kint1);
                }
            } else {
                A_new = multuplying_by_scalar(A_new, t, kint);
            }
            
            double *everybody;
            if (y1 != -1) {
                if (rank != size - 1) {
                    everybody = (double *) malloc(kint * sizeof(double));
                } else {
                    everybody = (double *) malloc(kint1 * sizeof(double));
                }
            } else {
                everybody = (double *) malloc(kint * sizeof(double));
            }
            
            for (int i = 0; i < size; i++) {
                if (rank == i) {
                    if (y1 != -1) {
                        if (rank != size - 1) {
                            int g1 = 0;
                            everybody = matrixes_subtraction(x, A_new, kint, g1);
                        } else {
                            int g1 = 0;
                            everybody = matrixes_subtraction(x, A_new, kint1, g1);
                        }
                    } else {
                        int g1 = 0;
                        everybody = matrixes_subtraction(x, A_new, kint, g1);
                    }
                }
            }
            
            if (l == 0) {
                break;
            }
            
#pragma omp barrier
            if (y1 != -1) {
                if (rank != size - 1) {
                    for (int i = 0; i < kint; i++) {
                        x[i + rank * kint] = everybody[i];
                    }
                } else {
                    for (int i = 0; i < kint1; i++) {
                        x[N - kint1 + i] = everybody[i];
                    }
                }
            } else {
                for (int i = 0; i < kint; i++) {
                    x[i + rank * kint] = everybody[i];
                }
            }
            
            if (rank == 0) {
                if ((sqrt(sum) / b_sqrt) < e) {
                    l = 0;
                } else sum = 0;
            }
        }
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

