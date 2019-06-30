#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>

#define N1 4
#define N2 4
#define N3 4

#define P1 4
#define P2 4

double **multiplying_matrixes(double **A, double **B, int overall_frame_y, int overall_frame_x) {
    double **C_sub = (double **) malloc(overall_frame_y * sizeof(double));
    int i;
    for (i = 0; i < overall_frame_y; i++) {
        C_sub[i] = (double *) malloc(overall_frame_x * sizeof(double));
    }

    for (i = 0; i < overall_frame_y; i++) {
        int k;
        for (k = 0; k < overall_frame_x; k++) {
            double sum = 0;
            int j;
            for (j = 0; j < N2; j++) {
                sum += A[i][j] * B[j][k];
            }
            C_sub[i][k] = sum;
        }
    }

    return C_sub;
}

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    int dims[2] = {P1, P2}, periods[2] = {0, 0}, coords[2], reorder = 1;
    int size, rank;

    MPI_Comm comm2d;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Dims_create(size, 2, dims);

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm2d);
    MPI_Comm_rank(comm2d, &rank);
    MPI_Cart_get(comm2d, 2, dims, periods, coords);

    double **A = (double **) malloc(N1 * sizeof(double *));
    int i;
    int j;
    for (i = 0; i < N1; i++) {
        A[i] = (double *) malloc(N2 * sizeof(double));
    }

    double **B = (double **) malloc(N2 * sizeof(double *));
    for (i = 0; i < N2; i++) {
        B[i] = (double *) malloc(N3 * sizeof(double));
    }

    if (rank == 0) {
        for (i = 0; i < N1; i++) {
            for (j = 0; j < N2; j++) {
                if (i != j)
                    A[i][j] = 0;
                else
                    A[i][j] = 1.0;
            }
        }

        for (i = 0; i < N2; i++) {
            for (j = 0; j < N3; j++) {
                if (i == 0)
                    B[i][j] = 1;
                else
                    B[i][j] = 0;
            }
        }
    }

    struct timeval tv1, tv2, dtv;
    struct timezone tz;
    if (rank == 0)
        gettimeofday(&tv1, &tz);


    int overall_frame_y = N1 / P1;
    int overall_frame_x = N3 / P2;

    double **A_overall = (double **) malloc(overall_frame_y * sizeof(double *));
    for (i = 0; i < overall_frame_y; i++) {
        A_overall[i] = (double *) malloc(N2 * sizeof(double));
    }

    double **B_overall = (double **) malloc(N2 * sizeof(double *));
    for (i = 0; i < N2; i++) {
        B_overall[i] = (double *) malloc(overall_frame_x * sizeof(double));
    }

    if (rank == 0) {
        int y;
        for (y = 1; y < P1 * P2; y++) {
            int divided_i = y / P2;
            int divided_j = y % P2;

            for (i = 0; i < overall_frame_y; i++) {
                for (j = 0; j < N2; j++) {
                    A_overall[i][j] = A[divided_i * overall_frame_y + i][j];
                }
            }

            for (i = 0; i < N2; i++) {
                for (j = 0; j < overall_frame_x; j++) {
                    B_overall[i][j] = B[i][overall_frame_x * divided_j + j];
                }
            }

            for (i = 0; i < overall_frame_y; i++) {
                MPI_Send(A_overall[i], N2, MPI_DOUBLE, y, i, comm2d);
            }

            for (i = 0; i < N2; i++) {
                MPI_Send(B_overall[i], overall_frame_x, MPI_DOUBLE, y, i, comm2d);
            }
        }

        y = 0;
        int divided_i = y / P2;
        int divided_j = y % P2;

        for (i = 0; i < overall_frame_y; i++) {
            for (j = 0; j < N2; j++) {
                A_overall[i][j] = A[divided_i * overall_frame_y + i][j];
            }
        }
        for (i = 0; i < N2; i++) {
            for (j = 0; j < overall_frame_x; j++) {
                B_overall[i][j] = B[i][overall_frame_x * divided_j + j];
            }
        }

    }

    int y;
    for (y = 1; y < P1 * P2; y++) {
        if (rank == y) {
            for (i = 0; i < overall_frame_y; i++) {
                MPI_Recv(A_overall[i], N2, MPI_DOUBLE, 0, i, comm2d, MPI_STATUS_IGNORE);
            }
            for (i = 0; i < N2; i++) {
                MPI_Recv(B_overall[i], overall_frame_x, MPI_DOUBLE, 0, i, comm2d, MPI_STATUS_IGNORE);
            }
        }
    }

    int amount = overall_frame_y * overall_frame_x;
    double *C_all = (double *) malloc(amount * P1 * P2 * sizeof(double));

    double **C_overall = (double **) malloc(N1 * sizeof(double *));
    for (i = 0; i < N1; i++) {
        C_overall[i] = (double *) malloc(N3 * sizeof(double));
    }
    double **C_sub = multiplying_matrixes(A_overall, B_overall, overall_frame_y, overall_frame_x);

    int k = 0;
    double *C = (double *) malloc(amount * sizeof(double));
    for (i = 0; i < overall_frame_y; i++) {
        for (j = 0; j < overall_frame_x; j++) {
            C[k] = C_sub[i][j];
            k++;
        }
    }

    MPI_Gather(C, amount, MPI_DOUBLE, C_all, amount, MPI_DOUBLE, 0, comm2d);
    int p = 0;
    int r = 0;
    if (rank == 0) {
        while (r != P1 * P2) {
            double *part = (double *) malloc(amount * sizeof(double));
            for (i = 0; i < amount; i++) {
                part[i] = C_all[p];
                p++;
            }
            int divided_i = r / P2;
            int divided_j = r % P2;

            int u = 0;
            for (i = 0; i < overall_frame_y; i++) {
                for (j = 0; j < overall_frame_x; j++) {
                    C_overall[i + overall_frame_y * divided_i][j + overall_frame_x * divided_j] = part[u];
                    u++;
                }
            }
            r++;
        }
    }

    if (rank == 0) {
        gettimeofday(&tv2, &tz);
        dtv.tv_sec = tv2.tv_sec - tv1.tv_sec;
        dtv.tv_usec = tv2.tv_usec - tv1.tv_usec;
        if (dtv.tv_usec < 0) {
            dtv.tv_sec--;
            dtv.tv_usec += 1000000;
        }
        printf("%f seconds\n", ((double) dtv.tv_sec * 1000 + dtv.tv_usec / 1000) / 1000);
    }

    if (rank == 0)
        for (i = 0; i < N1; i++) {
            for (j = 0; j < N3; j++) {
                printf("%f ", C_overall[i][j]);
            }
            printf("\n");
        }

    MPI_Finalize();
    return 0;
}

