#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>

#define N 3

double *multiplying_matrixes(double **A, double *b, int k, int start, int k_size) {
    double *final = (double *) malloc(k_size * sizeof(double));
    for (int i = 0; i < k_size; i++) {
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

double *matrixes_subtraction(double *A, double *b, int k, int start, int k_size) {
    double *final = (double *) malloc(k_size * sizeof(double));
    for (int i = 0; i < k_size; i++) {
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
    int size;
    int rank;
    
    struct timeval tv1, tv2, dtv;
    struct timezone tz;
    if (rank == 0)
        gettimeofday(&tv1, &tz);
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
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
    int *l = (int *) malloc(sizeof(int));
    
    l[0] = 1;
    int kint, kint1;
    
    if (N < size) {
        MPI_Finalize();
        if (rank == 0)
            printf("Reduce the number of processes\n");
        return 0;
        
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
    
    double *A_new = (double *) malloc(kint1 * sizeof(double));
    while (l[0] == 1) {
        MPI_Barrier(MPI_COMM_WORLD);
        
        for (int i = 0; i < kint1; i++) {
            A_new[i] = 0;
        }
        
        if (y1 != -1) {
            if (rank != size - 1) {
                int g1 = rank * kint;
                A_new = multiplying_matrixes(A1, x, kint, g1, kint1);
                A_new = matrixes_subtraction(A_new, b, kint, g1, kint1);
                
            } else {
                int g1 = N - kint1;
                A_new = multiplying_matrixes(A2, x, kint1, g1, kint1);
                A_new = matrixes_subtraction(A_new, b, kint1, g1, kint1);
                
            }
        } else {
            int g1 = rank * kint;
            A_new = multiplying_matrixes(A1, x, kint, g1, kint1);
            A_new = matrixes_subtraction(A_new, b, kint, g1, kint1);
        }
        
        double s = 0, sum;
        for (int i = 0; i < kint1; i++) {
            s += A_new[i] * A_new[i];
        }
        
        MPI_Reduce(&s, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        A_new = multuplying_by_scalar(A_new, t, kint1);
        
        double *everybody = (double *) malloc(kint1 * sizeof(double));
        for (int i = 0; i < kint1; i++) {
            everybody[i] = 0;
        }
        
        for (int i = 0; i < size; i++) {
            if (rank == i) {
                if (y1 != -1) {
                    if (rank != size - 1) {
                        int g1 = 0;
                        everybody = matrixes_subtraction(x, A_new, kint, g1, kint1);
                    } else {
                        int g1 = 0;
                        everybody = matrixes_subtraction(x, A_new, kint1, g1, kint1);
                    }
                } else {
                    int g1 = 0;
                    everybody = matrixes_subtraction(x, A_new, kint, g1, kint1);
                }
            }
        }
        
        double *new_x = (double *) malloc(size * kint1 * sizeof(double));
        
        MPI_Gather(everybody, kint1, MPI_DOUBLE, new_x, kint1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if (rank == 0) {
            if (y1 != -1) {
                int h = 0;
                int y = 0;
                while (h != N - kint1) {
                    for (int k = 0; k < kint; k++) {
                        x[h] = new_x[k + y];
                        h++;
                    }
                    y += kint1;
                }
                for (int u = y; u < kint1 * size; u++) {
                    x[h] = new_x[u];
                    h++;
                }
            } else {
                int h = 0;
                int y = 0;
                while (h != N) {
                    for (int k = 0; k < kint; k++) {
                        x[h] = new_x[k + y];
                        h++;
                    }
                    y += kint1;
                }
            }
        }
        
        MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if (rank == 0) {
            if ((sqrt(sum) / b_sqrt) < e) {
                l[0] = 0;
            }
        }
        MPI_Bcast(l, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    
    //    if (rank == 0) {
    //        gettimeofday(&tv2, &tz);
    //        dtv.tv_sec = tv2.tv_sec - tv1.tv_sec;
    //        dtv.tv_usec = tv2.tv_usec - tv1.tv_usec;
    //        if (dtv.tv_usec < 0) {
    //            dtv.tv_sec--;
    //            dtv.tv_usec += 1000000;
    //        }
    //        printf("%f seconds\n", ((double) dtv.tv_sec * 1000 + dtv.tv_usec / 1000)/ 1000);
    //    }
    
    if (rank == 0)
        for (int i = 0; i < N; i++) {
            printf("%f ", x[i]);
        }
    
    MPI_Finalize();
    return 0;
}

