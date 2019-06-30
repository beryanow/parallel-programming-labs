#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>

#define a 100000.0
#define e 0.00000001

double f(double x, double y, double z) {
    return x * x + y * y + z * z;
}

double p(double x, double y, double z) {
    return 6 - a * f(x, y, z);
}

double x_new(double x0, double i, double hx) {
    return x0 + i * hx;
}

double y_new(double y0, double j, double hy) {
    return y0 + j * hy;
}

double z_new(double z0, double k, double hz) {
    return z0 + k * hz;
}

void count_step(int size, int rank, double Ny, double Nz, double ***temp_addition_second, double ***array_divided,
                double ***addition_second, double ***temp_addition_first, double ***addition_first, int kint, int kint1,
                double ***array_new, double hx, double hy, double hz, double x0, double y0, double z0, int y1) {
    int k, j, i;
    if (size != 1) {
        for (i = 0; i < size; i++) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == i) {
                if (i == 0) {
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            temp_addition_second[0][j][k] = array_divided[kint - 1][j][k];
                        }
                        MPI_Send(temp_addition_second[0][j], (int) Nz, MPI_DOUBLE, i + 1, 123, MPI_COMM_WORLD);
                    }
                } else if (i != size - 1) {
                    for (j = 0; j < Ny; j++) {
                        MPI_Recv(temp_addition_second[0][j], (int) Nz, MPI_DOUBLE, i - 1, 123, MPI_COMM_WORLD,
                                 MPI_STATUS_IGNORE);
                    }
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            addition_second[0][j][k] = temp_addition_second[0][j][k];
                        }
                    }
                    
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            temp_addition_second[0][j][k] = array_divided[kint - 1][j][k];
                        }
                        MPI_Send(temp_addition_second[0][j], (int) Nz, MPI_DOUBLE, i + 1, 123, MPI_COMM_WORLD);
                    }
                } else {
                    for (j = 0; j < Ny; j++) {
                        MPI_Recv(temp_addition_second[0][j], (int) Nz, MPI_DOUBLE, i - 1, 123, MPI_COMM_WORLD,
                                 MPI_STATUS_IGNORE);
                    }
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            addition_second[0][j][k] = temp_addition_second[0][j][k];
                        }
                    }
                }
            }
        }
    }
    
    if (size != 1) {
        for (i = size - 1; i >= 0; i--) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == i) {
                if (i == size - 1) {
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            temp_addition_first[0][j][k] = array_divided[0][j][k];
                        }
                        MPI_Send(temp_addition_first[0][j], (int) Nz, MPI_DOUBLE, i - 1, 123, MPI_COMM_WORLD);
                    }
                } else if (i != 0) {
                    for (j = 0; j < Ny; j++) {
                        MPI_Recv(temp_addition_first[0][j], (int) Nz, MPI_DOUBLE, i + 1, 123, MPI_COMM_WORLD,
                                 MPI_STATUS_IGNORE);
                    }
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            addition_first[0][j][k] = temp_addition_first[0][j][k];
                        }
                    }
                    
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            temp_addition_first[0][j][k] = array_divided[0][j][k];
                        }
                        MPI_Send(temp_addition_first[0][j], (int) Nz, MPI_DOUBLE, i - 1, 123, MPI_COMM_WORLD);
                    }
                } else {
                    for (j = 0; j < Ny; j++) {
                        MPI_Recv(temp_addition_first[0][j], (int) Nz, MPI_DOUBLE, i + 1, 123, MPI_COMM_WORLD,
                                 MPI_STATUS_IGNORE);
                    }
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            addition_first[0][j][k] = temp_addition_first[0][j][k];
                        }
                    }
                }
            }
        }
    }
    
    if (size == 1) {
        for (i = 1; i < kint - 1; i++) {
            for (j = 1; j < Ny - 1; j++) {
                for (k = 1; k < Nz - 1; k++) {
                    array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                    (((array_divided[i + 1][j][k] + array_divided[i - 1][j][k]) /
                      (hx * hx)) +
                     ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                      (hy * hy)) +
                     ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                      (hz * hz)) -
                     p(x_new(x0, i, hx), y_new(y0, j, hy), z_new(z0, k, hz)));
                }
            }
        }
    } else {
        if (kint == 1) {
            if (size != 2) {
                if ((rank != 0) && (rank != size - 1)) {
                    for (j = 1; j < Ny - 1; j++) {
                        for (k = 1; k < Nz - 1; k++) {
                            array_new[kint - 1][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                            (((addition_first[0][j][k] + addition_second[0][j][k]) /
                              (hx * hx)) +
                             ((array_divided[kint - 1][j + 1][k] +
                               array_divided[kint - 1][j - 1][k]) /
                              (hy * hy)) +
                             ((array_divided[kint - 1][j][k + 1] +
                               array_divided[kint - 1][j][k - 1]) /
                              (hz * hz)) -
                             p(x_new(x0, rank, hx), y_new(y0, j, hy),
                               z_new(z0, k, hz)));
                        }
                    }
                }
                if ((rank == size - 1) && (kint1 != 1)) {
                    if (y1 != -1) {
                        for (i = 1; i < kint1 - 1; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                                    (((array_divided[i + 1][j][k] +
                                       array_divided[i - 1][j][k]) /
                                      (hx * hx)) +
                                     ((array_divided[i][j + 1][k] +
                                       array_divided[i][j - 1][k]) /
                                      (hy * hy)) +
                                     ((array_divided[i][j][k + 1] +
                                       array_divided[i][j][k - 1]) /
                                      (hz * hz)) -
                                     p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                                       z_new(z0, k, hz)));
                                }
                            }
                        }
                        int i = 0;
                        for (j = 1; j < Ny - 1; j++) {
                            for (k = 1; k < Nz - 1; k++) {
                                array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                                (((array_divided[i + 1][j][k] + addition_second[0][j][k]) /
                                  (hx * hx)) +
                                 ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                                  (hy * hy)) +
                                 ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                                  (hz * hz)) -
                                 p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                                   z_new(z0, k, hz)));
                            }
                        }
                    } else {
                        for (i = 1; i < kint - 1; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                                    (((array_divided[i + 1][j][k] +
                                       array_divided[i - 1][j][k]) /
                                      (hx * hx)) +
                                     ((array_divided[i][j + 1][k] +
                                       array_divided[i][j - 1][k]) /
                                      (hy * hy)) +
                                     ((array_divided[i][j][k + 1] +
                                       array_divided[i][j][k - 1]) /
                                      (hz * hz)) -
                                     p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                                       z_new(z0, k, hz)));
                                }
                            }
                        }
                        int i = 0;
                        for (j = 1; j < Ny - 1; j++) {
                            for (k = 1; k < Nz - 1; k++) {
                                array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                                (((array_divided[i + 1][j][k] + addition_second[0][j][k]) /
                                  (hx * hx)) +
                                 ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                                  (hy * hy)) +
                                 ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                                  (hz * hz)) -
                                 p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                                   z_new(z0, k, hz)));
                            }
                        }
                    }
                }
            }
        } else {
            if (y1 != -1) {
                if (rank == 0) {
                    for (i = 1; i < kint - 1; i++) {
                        for (j = 1; j < Ny - 1; j++) {
                            for (k = 1; k < Nz - 1; k++) {
                                array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                                (((array_divided[i + 1][j][k] + array_divided[i - 1][j][k]) /
                                  (hx * hx)) +
                                 ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                                  (hy * hy)) +
                                 ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                                  (hz * hz)) -
                                 p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                                   z_new(z0, k, hz)));
                            }
                        }
                    }
                    int i = kint - 1;
                    for (j = 1; j < Ny - 1; j++) {
                        for (k = 1; k < Nz - 1; k++) {
                            array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                            (((addition_first[0][j][k] + array_divided[i - 1][j][k]) /
                              (hx * hx)) +
                             ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                              (hy * hy)) +
                             ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                              (hz * hz)) -
                             p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                               z_new(z0, k, hz)));
                        }
                    }
                } else if (rank == size - 1) {
                    for (i = 1; i < kint1 - 1; i++) {
                        for (j = 1; j < Ny - 1; j++) {
                            for (k = 1; k < Nz - 1; k++) {
                                array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                                (((array_divided[i + 1][j][k] + array_divided[i - 1][j][k]) /
                                  (hx * hx)) +
                                 ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                                  (hy * hy)) +
                                 ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                                  (hz * hz)) -
                                 p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                                   z_new(z0, k, hz)));
                            }
                        }
                    }
                    int i = 0;
                    for (j = 1; j < Ny - 1; j++) {
                        for (k = 1; k < Nz - 1; k++) {
                            array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                            (((array_divided[i + 1][j][k] + addition_second[0][j][k]) /
                              (hx * hx)) +
                             ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                              (hy * hy)) +
                             ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                              (hz * hz)) -
                             p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                               z_new(z0, k, hz)));
                        }
                    }
                } else {
                    for (i = 1; i < kint - 1; i++) {
                        for (j = 1; j < Ny - 1; j++) {
                            for (k = 1; k < Nz - 1; k++) {
                                array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                                (((array_divided[i + 1][j][k] + array_divided[i - 1][j][k]) /
                                  (hx * hx)) +
                                 ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                                  (hy * hy)) +
                                 ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                                  (hz * hz)) -
                                 p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                                   z_new(z0, k, hz)));
                            }
                        }
                    }
                    int i = 0;
                    for (j = 1; j < Ny - 1; j++) {
                        for (k = 1; k < Nz - 1; k++) {
                            array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                            (((array_divided[i + 1][j][k] + addition_second[0][j][k]) /
                              (hx * hx)) +
                             ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                              (hy * hy)) +
                             ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                              (hz * hz)) -
                             p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                               z_new(z0, k, hz)));
                        }
                    }
                    i = kint - 1;
                    for (j = 1; j < Ny - 1; j++) {
                        for (k = 1; k < Nz - 1; k++) {
                            array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                            (((addition_first[0][j][k] + array_divided[i - 1][j][k]) /
                              (hx * hx)) +
                             ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                              (hy * hy)) +
                             ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                              (hz * hz)) -
                             p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                               z_new(z0, k, hz)));
                        }
                    }
                }
            } else {
                if (rank == 0) {
                    for (i = 1; i < kint - 1; i++) {
                        for (j = 1; j < Ny - 1; j++) {
                            for (k = 1; k < Nz - 1; k++) {
                                array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                                (((array_divided[i + 1][j][k] + array_divided[i - 1][j][k]) /
                                  (hx * hx)) +
                                 ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                                  (hy * hy)) +
                                 ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                                  (hz * hz)) -
                                 p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                                   z_new(z0, k, hz)));
                            }
                        }
                    }
                    int i = kint - 1;
                    for (j = 1; j < Ny - 1; j++) {
                        for (k = 1; k < Nz - 1; k++) {
                            array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                            (((addition_first[0][j][k] + array_divided[i - 1][j][k]) /
                              (hx * hx)) +
                             ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                              (hy * hy)) +
                             ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                              (hz * hz)) -
                             p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                               z_new(z0, k, hz)));
                        }
                    }
                } else if (rank == size - 1) {
                    for (i = 1; i < kint - 1; i++) {
                        for (j = 1; j < Ny - 1; j++) {
                            for (k = 1; k < Nz - 1; k++) {
                                array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                                (((array_divided[i + 1][j][k] + array_divided[i - 1][j][k]) /
                                  (hx * hx)) +
                                 ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                                  (hy * hy)) +
                                 ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                                  (hz * hz)) -
                                 p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                                   z_new(z0, k, hz)));
                            }
                        }
                    }
                    int i = 0;
                    for (j = 1; j < Ny - 1; j++) {
                        for (k = 1; k < Nz - 1; k++) {
                            array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                            (((array_divided[i + 1][j][k] + addition_second[0][j][k]) /
                              (hx * hx)) +
                             ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                              (hy * hy)) +
                             ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                              (hz * hz)) -
                             p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                               z_new(z0, k, hz)));
                        }
                    }
                } else {
                    for (i = 1; i < kint - 1; i++) {
                        for (j = 1; j < Ny - 1; j++) {
                            for (k = 1; k < Nz - 1; k++) {
                                array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                                (((array_divided[i + 1][j][k] + array_divided[i - 1][j][k]) /
                                  (hx * hx)) +
                                 ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                                  (hy * hy)) +
                                 ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                                  (hz * hz)) -
                                 p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                                   z_new(z0, k, hz)));
                            }
                        }
                    }
                    int i = 0;
                    for (j = 1; j < Ny - 1; j++) {
                        for (k = 1; k < Nz - 1; k++) {
                            array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                            (((array_divided[i + 1][j][k] + addition_second[0][j][k]) /
                              (hx * hx)) +
                             ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                              (hy * hy)) +
                             ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                              (hz * hz)) -
                             p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                               z_new(z0, k, hz)));
                        }
                    }
                    i = kint - 1;
                    for (j = 1; j < Ny - 1; j++) {
                        for (k = 1; k < Nz - 1; k++) {
                            array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                            (((addition_first[0][j][k] + array_divided[i - 1][j][k]) /
                              (hx * hx)) +
                             ((array_divided[i][j + 1][k] + array_divided[i][j - 1][k]) /
                              (hy * hy)) +
                             ((array_divided[i][j][k + 1] + array_divided[i][j][k - 1]) /
                              (hz * hz)) -
                             p(x_new(x0, i + rank * kint, hx), y_new(y0, j, hy),
                               z_new(z0, k, hz)));
                        }
                    }
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    
    int j, k, i;
    double x0 = -1, y0 = -1, z0 = -1, Dx = 2, Dy = 2, Dz = 2;
    double Nx = 5, Ny = 3, Nz = 5;
    
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    double ***array;
    if (rank == 0) {
        array = (double ***) malloc((int) Nx * sizeof(double **));
        for (i = 0; i < Nx; i++) {
            array[i] = (double **) malloc((int) Ny * sizeof(double *));
            for (j = 0; j < Ny; j++) {
                array[i][j] = (double *) malloc((int) Nz * sizeof(double));
            }
        }
    }
    
    struct timeval tv1, tv2, dtv;
    struct timezone tz;
    if (rank == 0)
    gettimeofday(&tv1, &tz);
    
    
    double t;
    int kint, kint1;
    int y1 = -1;
    
    if ((int) Nx < size) {
        MPI_Finalize();
        if (rank == 0)
        printf("Reduce the number of processes\n");
        return 0;
        
    } else {
        t = Nx / size;
        
        kint = (int) t;
        kint1 = (int) Nx - kint * (size - 1);
        
        if (t - (int) t > 0) {
            y1 = 0;
        }
    }
    
    double ***array_divided;
    
    // создание раздельных кусков
    
    if ((y1 != -1) && (rank == size - 1)) {
        array_divided = (double ***) malloc(kint1 * sizeof(double **));
        for (i = 0; i < kint1; i++) {
            array_divided[i] = (double **) malloc((int) Ny * sizeof(double *));
            for (j = 0; j < Ny; j++) {
                array_divided[i][j] = (double *) malloc((int) Nz * sizeof(double));
            }
        }
    } else {
        array_divided = (double ***) malloc(kint * sizeof(double **));
        for (i = 0; i < kint; i++) {
            array_divided[i] = (double **) malloc((int) Ny * sizeof(double *));
            for (j = 0; j < Ny; j++) {
                array_divided[i][j] = (double *) malloc((int) Nz * sizeof(double));
            }
        }
    }
    
    if ((y1 != -1) && (rank == size - 1)) {
        for (i = 0; i < kint1; i++) {
            for (j = 0; j < Ny; j++) {
                for (k = 0; k < Nz; k++) {
                    array_divided[i][j][k] = 0;
                }
            }
        }
    } else {
        for (i = 0; i < kint; i++) {
            for (j = 0; j < Ny; j++) {
                for (k = 0; k < Nz; k++) {
                    array_divided[i][j][k] = 0;
                }
            }
        }
    }
    
    // конец создания раздельных кусков
    
    double hx = Dx / (Nx - 1);
    double hy = Dy / (Nx - 1);
    double hz = Dz / (Nx - 1);
    
    // заполнение верхних/нижних границ
    
    if (y1 != -1) {
        if (rank == 0) {
            for (j = 0; j < Ny; j++) {
                for (k = 0; k < Nz; k++) {
                    array_divided[0][j][k] = f(x_new(x0, 0, hx), y_new(y0, j, hy), z_new(z0, k, hz));
                }
            }
        } else if (rank == size - 1) {
            for (j = 0; j < Ny; j++) {
                for (k = 0; k < Nz; k++) {
                    array_divided[kint1 - 1][j][k] = f(x_new(x0, Nx - 1, hx), y_new(y0, j, hy), z_new(z0, k, hz));
                }
            }
        }
    } else {
        if (size == 1) {
            for (j = 0; j < Ny; j++) {
                for (k = 0; k < Nz; k++) {
                    array_divided[0][j][k] = f(x_new(x0, 0, hx), y_new(y0, j, hy), z_new(z0, k, hz));
                    array_divided[kint - 1][j][k] = f(x_new(x0, Nx - 1, hx), y_new(y0, j, hy), z_new(z0, k, hz));
                }
            }
        } else {
            if (rank == 0) {
                for (j = 0; j < Ny; j++) {
                    for (k = 0; k < Nz; k++) {
                        array_divided[0][j][k] = f(x_new(x0, 0, hx), y_new(y0, j, hy), z_new(z0, k, hz));
                    }
                }
            } else if (rank == size - 1) {
                for (j = 0; j < Ny; j++) {
                    for (k = 0; k < Nz; k++) {
                        array_divided[kint - 1][j][k] = f(x_new(x0, Nx - 1, hx), y_new(y0, j, hy), z_new(z0, k, hz));
                    }
                }
            }
        }
    }
    
    // конец заполнения верхних/нижних границ
    
    // заполнение границ со сторон
    
    if (y1 != -1) {
        if (rank == 0) {
            for (i = 1; i < kint; i++) {
                for (j = 0; j < Ny; j++) {
                    array_divided[i][j][(int) Nz - 1] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy),
                                                          z_new(z0, Nz - 1, hz));
                    array_divided[i][j][0] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy), z_new(z0, 0, hz));
                }
                
                for (k = 0; k < Nz; k++) {
                    array_divided[i][(int) Ny - 1][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, Ny - 1, hy),
                                                          z_new(z0, k, hz));
                    array_divided[i][0][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, 0, hy), z_new(z0, k, hz));
                }
            }
        } else if (rank == size - 1) {
            for (i = 0; i < kint1 - 1; i++) {
                for (j = 0; j < Ny; j++) {
                    array_divided[i][j][(int) Nz - 1] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy),
                                                          z_new(z0, Nz - 1, hz));
                    array_divided[i][j][0] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy), z_new(z0, 0, hz));
                }
                
                for (k = 0; k < Nz; k++) {
                    array_divided[i][(int) Ny - 1][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, Ny - 1, hy),
                                                          z_new(z0, k, hz));
                    array_divided[i][0][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, 0, hy), z_new(z0, k, hz));
                }
            }
        } else {
            for (i = 0; i < kint; i++) {
                for (j = 0; j < Ny; j++) {
                    array_divided[i][j][(int) Nz - 1] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy),
                                                          z_new(z0, Nz - 1, hz));
                    array_divided[i][j][0] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy), z_new(z0, 0, hz));
                }
                
                for (k = 0; k < Nz; k++) {
                    array_divided[i][(int) Ny - 1][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, Ny - 1, hy),
                                                          z_new(z0, k, hz));
                    array_divided[i][0][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, 0, hy), z_new(z0, k, hz));
                }
            }
        }
    } else {
        if (size == 1) {
            for (i = 1; i < kint - 1; i++) {
                for (j = 0; j < Ny; j++) {
                    array_divided[i][j][(int) Nz - 1] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy),
                                                          z_new(z0, Nz - 1, hz));
                    array_divided[i][j][0] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy), z_new(z0, 0, hz));
                }
                
                for (k = 0; k < Nz; k++) {
                    array_divided[i][(int) Ny - 1][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, Ny - 1, hy),
                                                          z_new(z0, k, hz));
                    array_divided[i][0][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, 0, hy), z_new(z0, k, hz));
                }
            }
        } else {
            if (rank == 0) {
                for (i = 1; i < kint; i++) {
                    for (j = 0; j < Ny; j++) {
                        array_divided[i][j][(int) Nz - 1] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy),
                                                              z_new(z0, Nz - 1, hz));
                        array_divided[i][j][0] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy),
                                                   z_new(z0, 0, hz));
                    }
                    
                    for (k = 0; k < Nz; k++) {
                        array_divided[i][(int) Ny - 1][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, Ny - 1, hy),
                                                              z_new(z0, k, hz));
                        array_divided[i][0][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, 0, hy),
                                                   z_new(z0, k, hz));
                    }
                }
            } else if (rank == size - 1) {
                for (i = 0; i < kint - 1; i++) {
                    for (j = 0; j < Ny; j++) {
                        array_divided[i][j][(int) Nz - 1] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy),
                                                              z_new(z0, Nz - 1, hz));
                        array_divided[i][j][0] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy),
                                                   z_new(z0, 0, hz));
                    }
                    
                    for (k = 0; k < Nz; k++) {
                        array_divided[i][(int) Ny - 1][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, Ny - 1, hy),
                                                              z_new(z0, k, hz));
                        array_divided[i][0][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, 0, hy),
                                                   z_new(z0, k, hz));
                    }
                }
            } else {
                for (i = 0; i < kint; i++) {
                    for (j = 0; j < Ny; j++) {
                        array_divided[i][j][(int) Nz - 1] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy),
                                                              z_new(z0, Nz - 1, hz));
                        array_divided[i][j][0] = f(x_new(x0, i + kint * rank, hx), y_new(y0, j, hy), z_new(z0, 0, hz));
                    }
                    
                    for (k = 0; k < Nz; k++) {
                        array_divided[i][(int) Ny - 1][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, Ny - 1, hy),
                                                              z_new(z0, k, hz));
                        array_divided[i][0][k] = f(x_new(x0, i + kint * rank, hx), y_new(y0, 0, hy), z_new(z0, k, hz));
                    }
                }
            }
        }
    }
    
    // конец заполнения границ со сторон
    
    double ***array_new;
    
    if ((y1 != -1) && (rank == size - 1)) {
        array_new = (double ***) malloc(kint1 * sizeof(double **));
        for (i = 0; i < kint1; i++) {
            array_new[i] = (double **) malloc((int) Ny * sizeof(double *));
            for (j = 0; j < Ny; j++) {
                array_new[i][j] = (double *) malloc((int) Nz * sizeof(double));
            }
        }
    } else {
        array_new = (double ***) malloc(kint * sizeof(double **));
        for (i = 0; i < kint; i++) {
            array_new[i] = (double **) malloc((int) Ny * sizeof(double *));
            for (j = 0; j < Ny; j++) {
                array_new[i][j] = (double *) malloc((int) Nz * sizeof(double));
            }
        }
    }
    
    // создание граничных плоскостей, которые будут передаваться для подсчёта нового шага
    
    double ***temp_addition_first = (double ***) malloc(sizeof(double **));
    temp_addition_first[0] = (double **) malloc((int) Ny * sizeof(double *));
    for (j = 0; j < Ny; j++) {
        temp_addition_first[0][j] = (double *) malloc((int) Nz * sizeof(double));
    }
    
    double ***temp_addition_second = (double ***) malloc(sizeof(double **));
    temp_addition_second[0] = (double **) malloc((int) Ny * sizeof(double *));
    for (j = 0; j < Ny; j++) {
        temp_addition_second[0][j] = (double *) malloc((int) Nz * sizeof(double));
    }
    
    double ***addition_first = (double ***) malloc(sizeof(double **));
    addition_first[0] = (double **) malloc((int) Ny * sizeof(double *));
    for (j = 0; j < Ny; j++) {
        addition_first[0][j] = (double *) malloc((int) Nz * sizeof(double));
    }
    
    double ***addition_second = (double ***) malloc(sizeof(double **));
    addition_second[0] = (double **) malloc((int) Ny * sizeof(double *));
    for (j = 0; j < Ny; j++) {
        addition_second[0][j] = (double *) malloc((int) Nz * sizeof(double));
    }
    
    for (j = 0; j < Ny; j++) {
        for (k = 0; k < Nz; k++) {
            addition_first[0][j][k] = 0;
            addition_second[0][j][k] = 0;
        }
    }
    
    for (j = 0; j < Ny; j++) {
        for (k = 0; k < Nz; k++) {
            temp_addition_first[0][j][k] = 0;
            temp_addition_second[0][j][k] = 0;
        }
    }
    
    // конец создания граничных плоскостей, которые будут передаваться для подсчёта нового шага
    
    // пересылка границ
    
    if (size != 1) {
        for (i = 0; i < size; i++) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == i) {
                if (i == 0) {
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            temp_addition_second[0][j][k] = array_divided[kint - 1][j][k];
                        }
                        MPI_Send(temp_addition_second[0][j], (int) Nz, MPI_DOUBLE, i + 1, 123, MPI_COMM_WORLD);
                    }
                } else if (i != size - 1) {
                    for (j = 0; j < Ny; j++) {
                        MPI_Recv(temp_addition_second[0][j], (int) Nz, MPI_DOUBLE, i - 1, 123, MPI_COMM_WORLD,
                                 MPI_STATUS_IGNORE);
                    }
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            addition_second[0][j][k] = temp_addition_second[0][j][k];
                        }
                    }
                    
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            temp_addition_second[0][j][k] = array_divided[kint - 1][j][k];
                        }
                        MPI_Send(temp_addition_second[0][j], (int) Nz, MPI_DOUBLE, i + 1, 123, MPI_COMM_WORLD);
                    }
                } else {
                    for (j = 0; j < Ny; j++) {
                        MPI_Recv(temp_addition_second[0][j], (int) Nz, MPI_DOUBLE, i - 1, 123, MPI_COMM_WORLD,
                                 MPI_STATUS_IGNORE);
                    }
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            addition_second[0][j][k] = temp_addition_second[0][j][k];
                        }
                    }
                }
            }
        }
    }
    
    if (size != 1) {
        for (i = size - 1; i >= 0; i--) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == i) {
                if (i == size - 1) {
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            temp_addition_first[0][j][k] = array_divided[0][j][k];
                        }
                        MPI_Send(temp_addition_first[0][j], (int) Nz, MPI_DOUBLE, i - 1, 123, MPI_COMM_WORLD);
                    }
                } else if (i != 0) {
                    for (j = 0; j < Ny; j++) {
                        MPI_Recv(temp_addition_first[0][j], (int) Nz, MPI_DOUBLE, i + 1, 123, MPI_COMM_WORLD,
                                 MPI_STATUS_IGNORE);
                    }
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            addition_first[0][j][k] = temp_addition_first[0][j][k];
                        }
                    }
                    
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            temp_addition_first[0][j][k] = array_divided[0][j][k];
                        }
                        MPI_Send(temp_addition_first[0][j], (int) Nz, MPI_DOUBLE, i - 1, 123, MPI_COMM_WORLD);
                    }
                } else {
                    for (j = 0; j < Ny; j++) {
                        MPI_Recv(temp_addition_first[0][j], (int) Nz, MPI_DOUBLE, i + 1, 123, MPI_COMM_WORLD,
                                 MPI_STATUS_IGNORE);
                    }
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            addition_first[0][j][k] = temp_addition_first[0][j][k];
                        }
                    }
                }
            }
        }
    }
    
    // конец пересылки границ
    
    count_step(size, rank, Ny, Nz, temp_addition_second, array_divided, addition_second, temp_addition_first,
               addition_first, kint, kint1, array_new, hx, hy, hz, x0, y0, z0, y1);
    
    double *max_array = (double *) malloc(size * sizeof(double));
    
    while (1) {
        double max = -1;
        
        // подсчёт максимума для каждого процесса
        
        if (size == 1) {
            for (i = 1; i < kint - 1; i++) {
                for (j = 1; j < Ny - 1; j++) {
                    for (k = 1; k < Nz - 1; k++) {
                        double value;
                        if ((value = array_new[i][j][k] - array_divided[i][j][k]) < 0) {
                            value *= -1;
                        }
                        
                        if (value > max) {
                            max = value;
                        }
                    }
                }
            }
        } else {
            if (kint == 1) {
                if (size != 2) {
                    if ((rank != 0) && (rank != size - 1)) {
                        for (j = 1; j < Ny - 1; j++) {
                            for (k = 1; k < Nz - 1; k++) {
                                double value;
                                if ((value = array_new[kint - 1][j][k] - array_divided[kint - 1][j][k]) < 0) {
                                    value *= -1;
                                }
                                
                                if (value > max) {
                                    max = value;
                                }
                            }
                        }
                    }
                }
            } else {
                if (y1 != -1) {
                    if (rank == 0) {
                        for (i = 1; i < kint; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    double value;
                                    if ((value = array_new[i][j][k] - array_divided[i][j][k]) < 0) {
                                        value *= -1;
                                    }
                                    
                                    if (value > max) {
                                        max = value;
                                    }
                                }
                            }
                        }
                    } else if (rank == size - 1) {
                        for (i = 0; i < kint1 - 1; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    double value;
                                    if ((value = array_new[i][j][k] - array_divided[i][j][k]) < 0) {
                                        value *= -1;
                                    }
                                    
                                    if (value > max) {
                                        max = value;
                                    }
                                }
                            }
                        }
                    } else {
                        for (i = 0; i < kint; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    double value;
                                    if ((value = array_new[i][j][k] - array_divided[i][j][k]) < 0) {
                                        value *= -1;
                                    }
                                    
                                    if (value > max) {
                                        max = value;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (rank == 0) {
                        for (i = 1; i < kint; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    double value;
                                    if ((value = array_new[i][j][k] - array_divided[i][j][k]) < 0) {
                                        value *= -1;
                                    }
                                    
                                    if (value > max) {
                                        max = value;
                                    }
                                }
                            }
                        }
                    } else if (rank == size - 1) {
                        for (i = 0; i < kint - 1; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    double value;
                                    if ((value = array_new[i][j][k] - array_divided[i][j][k]) < 0) {
                                        value *= -1;
                                    }
                                    
                                    if (value > max) {
                                        max = value;
                                    }
                                }
                            }
                        }
                    } else {
                        for (i = 0; i < kint; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    double value;
                                    if ((value = array_new[i][j][k] - array_divided[i][j][k]) < 0) {
                                        value *= -1;
                                    }
                                    
                                    if (value > max) {
                                        max = value;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // конец подсчёта максимума для каждого процесса
        
        // копирование полученных значений
        
        if (size == 1) {
            for (i = 1; i < kint - 1; i++) {
                for (j = 1; j < Ny - 1; j++) {
                    for (k = 1; k < Nz - 1; k++) {
                        array_divided[i][j][k] = array_new[i][j][k];
                    }
                }
            }
        } else {
            if (kint == 1) {
                if (size != 2) {
                    if ((rank != 0) && (rank != size - 1)) {
                        for (j = 1; j < Ny - 1; j++) {
                            for (k = 1; k < Nz - 1; k++) {
                                array_divided[kint - 1][j][k] = array_new[kint - 1][j][k];
                            }
                        }
                    }
                    if ((rank == size - 1) && (kint1 != 1)) {
                        if (y1 != -1) {
                            for (i = 0; i < kint1 - 1; i++) {
                                for (j = 1; j < Ny - 1; j++) {
                                    for (k = 1; k < Nz - 1; k++) {
                                        array_divided[i][j][k] = array_new[i][j][k];
                                    }
                                }
                            }
                        } else {
                            for (i = 0; i < kint - 1; i++) {
                                for (j = 1; j < Ny - 1; j++) {
                                    for (k = 1; k < Nz - 1; k++) {
                                        array_divided[i][j][k] = array_new[i][j][k];
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (y1 != -1) {
                    if (rank == 0) {
                        for (i = 1; i < kint; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    array_divided[i][j][k] = array_new[i][j][k];
                                }
                            }
                        }
                    } else if (rank == size - 1) {
                        for (i = 0; i < kint1 - 1; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    array_divided[i][j][k] = array_new[i][j][k];
                                }
                            }
                        }
                    } else {
                        for (i = 0; i < kint; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    array_divided[i][j][k] = array_new[i][j][k];
                                }
                            }
                        }
                    }
                } else {
                    if (rank == 0) {
                        for (i = 1; i < kint; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    array_divided[i][j][k] = array_new[i][j][k];
                                }
                            }
                        }
                    } else if (rank == size - 1) {
                        for (i = 0; i < kint - 1; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    array_divided[i][j][k] = array_new[i][j][k];
                                }
                            }
                        }
                    } else {
                        for (i = 0; i < kint; i++) {
                            for (j = 1; j < Ny - 1; j++) {
                                for (k = 1; k < Nz - 1; k++) {
                                    array_divided[i][j][k] = array_new[i][j][k];
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // конец копирования полученных значений
        
        MPI_Gather(&max, 1, MPI_DOUBLE, max_array, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if (rank == 0) {
            max = -1;
            for (i = 0; i < size; i++) {
                if (max_array[i] > max) {
                    max = max_array[i];
                }
            }
        }
        
        MPI_Bcast(&max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if (max < e) {
            break;
        }
        
        count_step(size, rank, Ny, Nz, temp_addition_second, array_divided, addition_second, temp_addition_first,
                   addition_first, kint, kint1, array_new, hx, hy, hz, x0, y0, z0, y1);
    }
    
    // сборка разделённых частей из всех процессов в изначальный куб в нулевом процессе
    
    double ***array_divided_new = (double ***) malloc(kint1 * sizeof(double **));
    for (i = 0; i < kint1; i++) {
        array_divided_new[i] = (double **) malloc((int) Ny * sizeof(double *));
        for (j = 0; j < Ny; j++) {
            array_divided_new[i][j] = (double *) malloc((int) Nz * sizeof(double));
        }
    }
    
    if ((y1 != -1) && (rank == size - 1)) {
        for (i = 0; i < kint1; i++) {
            for (j = 0; j < Ny; j++) {
                for (k = 0; k < Nz; k++) {
                    array_divided_new[i][j][k] = array_divided[i][j][k];
                }
            }
        }
        
        for (i = 0; i < kint1; i++) {
            for (j = 0; j < Ny; j++) {
                MPI_Send(array_divided_new[i][j], (int) Nz, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD);
            }
        }
    } else {
        if (rank != 0) {
            for (i = 0; i < kint; i++) {
                for (j = 0; j < Ny; j++) {
                    MPI_Send(array_divided[i][j], (int) Nz, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD);
                }
            }
        }
    }
    
    if (rank == 0) {
        if (y1 != -1) {
            int y = 0;
            for (i = 0; i < kint; i++) {
                for (j = 0; j < Ny; j++) {
                    for (k = 0; k < Nz; k++) {
                        array[i + y * kint][j][k] = array_divided[i][j][k];
                    }
                }
            }
            
            for (y = 1; y < size - 1; y++) {
                for (i = 0; i < kint; i++) {
                    for (j = 0; j < Ny; j++) {
                        MPI_Recv(array_divided[i][j], (int) Nz, MPI_DOUBLE, y, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
                for (i = 0; i < kint; i++) {
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            array[i + y * kint][j][k] = array_divided[i][j][k];
                        }
                    }
                }
            }
            y = size - 1;
            
            for (i = 0; i < kint1; i++) {
                for (j = 0; j < Ny; j++) {
                    MPI_Recv(array_divided_new[i][j], (int) Nz, MPI_DOUBLE, y, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            for (i = 0; i < kint1; i++) {
                for (j = 0; j < Ny; j++) {
                    for (k = 0; k < Nz; k++) {
                        array[(size - 1) * kint + i][j][k] = array_divided_new[i][j][k];
                    }
                }
            }
            
        } else {
            int y = 0;
            for (i = 0; i < kint; i++) {
                for (j = 0; j < Ny; j++) {
                    for (k = 0; k < Nz; k++) {
                        array[i + y * kint][j][k] = array_divided[i][j][k];
                    }
                }
            }
            for (y = 1; y < size; y++) {
                for (i = 0; i < kint; i++) {
                    for (j = 0; j < Ny; j++) {
                        MPI_Recv(array_divided[i][j], (int) Nz, MPI_DOUBLE, y, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
                for (i = 0; i < kint; i++) {
                    for (j = 0; j < Ny; j++) {
                        for (k = 0; k < Nz; k++) {
                            array[i + y * kint][j][k] = array_divided[i][j][k];
                        }
                    }
                }
            }
        }
    }
    
    // конец сборки разделённых частей из всех процессов в изначальный куб в нулевом процессе
    
//        if (rank == 0) {
//            gettimeofday(&tv2, &tz);
//            dtv.tv_sec = tv2.tv_sec - tv1.tv_sec;
//            dtv.tv_usec = tv2.tv_usec - tv1.tv_usec;
//            if (dtv.tv_usec < 0) {
//                dtv.tv_sec--;
//                dtv.tv_usec += 1000000;
//            }
//            printf("%f seconds\n", ((double) dtv.tv_sec * 1000 + dtv.tv_usec / 1000) / 1000);
//        }
    
    if (rank == 0) {
        for (i = 0; i < Nx; i++) {
            for (j = 0; j < Ny; j++) {
                for (k = 0; k < Nz; k++) {
                    printf(" %f ", array[i][j][k]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }
    
    MPI_Finalize();
    
    return 0;
}

