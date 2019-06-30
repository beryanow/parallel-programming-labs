#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>
#include <mpi.h>

#define pi 3.14
#define N 1500000
#define LX 100
#define NX 10
#define ITERATIONS 50

double field_value(double xk) {
    return sin(2 * pi * xk / NX);
}

int main(int argc, char *argv[]) {
    int size;
    int rank;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double *overall_coordinates;
    
    double *coordinates;
    coordinates = (double *) malloc(sizeof(double) * N);
    double *speeds;
    speeds = (double *) malloc(sizeof(double) * N);
    double *speeds_last;
    speeds_last = (double *) malloc(sizeof(double) * N);
    
    double *charges_nodes;
    charges_nodes = (double *) malloc(sizeof(double) * NX);
    double *field_values_nodes;
    field_values_nodes = (double *) malloc(sizeof(double) * NX);
    double *charges;
    charges = (double *) malloc(sizeof(double) * N);
    double *weights;
    weights = (double *) malloc(sizeof(double) * N);
    double *nodes_coordinates;
    nodes_coordinates = (double *) malloc(sizeof(double) * NX);
    double *field_coordinates;
    field_coordinates = (double *) malloc(sizeof(double) * N);
    
    double step_for_coordinates = (double) LX / N;
    double step_for_grid = (double) LX / NX;
    
    int divided = N / size;
    
    for (int i = 0; i < divided; i++) {
        coordinates[i] = (i + rank * divided) * step_for_coordinates;
        speeds[i] = 0;
        charges[i] = 0;
    }
    
    //    int j = 0;
    //    while (j != divided / 2) {
    //        int y = rand() % divided;
    //        if (charges[y] == 0) {
    //            charges[y] = -1;
    //            weights[y] = 1;
    //            j++;
    //        }
    //    }
    
    for (int i = 0; i < divided; i++) {
        if (charges[i] == 0) {
            charges[i] = 1;
            weights[i] = 100;
        }
    }
    
    for (int i = 0; i < NX; i++) {
        nodes_coordinates[i] = i * step_for_grid;
        field_values_nodes[i] = field_value(i * step_for_grid);
    }
    
    struct timeval tv1, tv2, dtv;
    struct timezone tz;
    if (rank == 0) {
        gettimeofday(&tv1, &tz);
    }
    
    int f = 0;
    int g = 0;
    while (g != ITERATIONS) {
        
        for (int i = 0; i < NX; i++) {
            double y = i * step_for_grid;
            double y1 = (i - 1) * step_for_grid;
            double y2 = (i + 1) * step_for_grid;
            
            double y_start = y - ((y - y1) / 2);
            double y_end = y + ((y2 - y) / 2);
            if (y1 < 0) {
                y_start = 0;
            }
            if (y2 > LX) {
                y_end = LX;
            }
            
            double charge_sum = 0;
            for (int v = 0; v < divided; v++) {
                if ((coordinates[v] >= y_start) && (coordinates[v] <= y_end)) {
                    charge_sum += charges[v];
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            
            double sum = 0;
            
            MPI_Reduce(&charge_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            
            if (rank == 0) {
                charges_nodes[i] = sum / N;
            }
        }
        
        for (int i = 0; i < divided; i++) {
            int d;
            for (d = 0; d < NX - 1; d++) {
                if ((nodes_coordinates[d] <= coordinates[i]) && (nodes_coordinates[d + 1] >= coordinates[i])) {
                    break;
                }
            }
            field_coordinates[i] =
            ((coordinates[i] - nodes_coordinates[d]) * field_values_nodes[d + 1] /
             step_for_grid) +
            (coordinates[i] - nodes_coordinates[d + 1]) / step_for_grid *
            field_values_nodes[d];
            speeds[i] = charges[i] / weights[i] * field_coordinates[i];
        }
        
        if (f == 1) {
            for (int i = 0; i < divided; i++) {
                coordinates[i] = speeds_last[i] + coordinates[i];
            }
        }
        
        for (int i = 0; i < divided; i++) {
            speeds_last[i] = speeds[i];
        }
        
        //        if (rank == 0) {
        //            overall_coordinates = (double *) malloc(sizeof(double) * N);
        //        }
        //        MPI_Gather(coordinates, divided, MPI_DOUBLE, overall_coordinates, divided, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //
        //        if (rank == 0) {
        //            char name[11];
        //            for (int i = 0; i < 11; i++) {
        //                name[i] = '\0';
        //            }
        //            sprintf(name, "%d", g);
        //            int h;
        //            for (h = 0; h < 11; h++) {
        //                if (name[h] == '\0') {
        //                    break;
        //                }
        //            }
        //            name[h] = '.';
        //            name[h + 1] = 't';
        //            name[h + 2] = 'x';
        //            name[h + 3] = 't';
        //
        //            FILE *file = fopen(name, "wb");
        //
        //            for (int i = 0; i < N; i++) {
        //                fprintf(file, "Number: %d Coordinate: %f\n", i, overall_coordinates[i]);
        //            }
        //
        //            fprintf(file, "\n");
        //
        //            for (int i = 0; i < NX; i++) {
        //                fprintf(file, "Node number: %d Density: %f\n", i, charges_nodes[i]);
        //            }
        //
        //        }
        f = 1;
        g++;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        gettimeofday(&tv2, &tz);
        dtv.tv_sec = tv2.tv_sec - tv1.tv_sec;
        dtv.tv_usec = tv2.tv_usec - tv1.tv_usec;
        if (dtv.tv_usec < 0) {
            dtv.tv_sec--;
            dtv.tv_usec += 1000000;
        }
        
        printf("Time: %f seconds\n", ((double) dtv.tv_sec * 1000 + dtv.tv_usec / 1000) / 1000);
    }
    MPI_Finalize();
    
    return 0;
}

