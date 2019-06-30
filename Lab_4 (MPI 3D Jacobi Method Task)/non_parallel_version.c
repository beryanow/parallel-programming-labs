#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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

int main() {
    double x0 = -1, y0 = -1, z0 = -1, Dx = 2, Dy = 2, Dz = 2;
    double Nx = 5, Ny = 3, Nz = 5;

    double ***array = (double ***) malloc((int) Nx * sizeof(double **));
    for (int i = 0; i < Nx; i++) {
        array[i] = (double **) malloc((int) Ny * sizeof(double *));
        for (int j = 0; j < Ny; j++) {
            array[i][j] = (double *) malloc((int) Nz * sizeof(double));
        }
    }

    double ***array_new = (double ***) malloc((int) Nx * sizeof(double **));
    for (int i = 0; i < Nx; i++) {
        array_new[i] = (double **) malloc((int) Ny * sizeof(double *));
        for (int j = 0; j < Ny; j++) {
            array_new[i][j] = (double *) malloc((int) Nz * sizeof(double));
        }
    }

    double hx = Dx / (Nx - 1);
    double hy = Dy / (Nx - 1);
    double hz = Dz / (Nx - 1);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                array[i][j][k] = 0;
            }
        }
    }

    for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nz; k++) {
            array[(int) Nx - 1][j][k] = f(x_new(x0, Nx - 1, hx), y_new(y0, j, hy), z_new(z0, k, hz));
            array[0][j][k] = f(x_new(x0, 0, hx), y_new(y0, j, hy), z_new(z0, k, hz));
        }
    }

    for (int i = 0; i < Nx; i++) {
        for (int k = 0; k < Nz; k++) {
            array[i][(int) Ny - 1][k] = f(x_new(x0, i, hx), y_new(y0, Ny - 1, hy), z_new(z0, k, hz));
            array[i][0][k] = f(x_new(x0, i, hx), y_new(y0, 0, hy), z_new(z0, k, hz));
        }
    }

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            array[i][j][(int) Nz - 1] = f(x_new(x0, i, hx), y_new(y0, j, hy), z_new(z0, Nz - 1, hz));
            array[i][j][0] = f(x_new(x0, i, hx), y_new(y0, j, hy), z_new(z0, 0, hz));
        }
    }

    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                                     (((array[i + 1][j][k] + array[i - 1][j][k]) / (hx * hx)) +
                                      ((array[i][j + 1][k] + array[i][j - 1][k]) / (hy * hy)) +
                                      ((array[i][j][k + 1] + array[i][j][k - 1]) / (hz * hz)) -
                                      p(x_new(x0, i, hx), y_new(y0, j, hy), z_new(z0, k, hz)));
            }
        }
    }

    int l = 0;
    while (l == 0) {
        double max = -1;

        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int k = 1; k < Nz - 1; k++) {
                    double value;
                    if ((value = array_new[i][j][k] - array[i][j][k]) < 0) {
                        value *= -1;
                    }

                    if (value > max) {
                        max = value;
                    }
                }
            }
        }

        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int k = 1; k < Nz - 1; k++) {
                    array[i][j][k] = array_new[i][j][k];
                }
            }
        }

        if (max < e) {
            l = 1;
            break;
        }

        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int k = 1; k < Nz - 1; k++) {
                    array_new[i][j][k] = (1 / ((2 / hx * hx) + (2 / hy * hy) + (2 / hz * hz) + a)) *
                                         (((array[i + 1][j][k] + array[i - 1][j][k]) / (hx * hx)) +
                                          ((array[i][j + 1][k] + array[i][j - 1][k]) / (hy * hy)) +
                                          ((array[i][j][k + 1] + array[i][j][k - 1]) / (hz * hz)) -
                                          p(x_new(x0, i, hx), y_new(y0, j, hy), z_new(z0, k, hz)));
                }
            }
        }
    }

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                printf(" %f ", array[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
    return 0;
}
