#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define TASK_AMOUNT 100
#define ITERATIONS 1
#define MIN_BALANCE 20

double every_iteration_local_time = 0;
double global_res = 0;

int tasks[TASK_AMOUNT];
int rank = 0;
int size = 0;
int iteration = 0;

int tasks_implemented_amount = 0;
int total_tasks_implemented_amount = 0;
int task_to_do = 0;
int current_task = 0;
const int request = -1;
const int finish = -2;

pthread_mutex_t mutex;

void *receiving(void *me) {
    int tmp = 0;
    int tasks_to_send = 0;
    MPI_Status status;
    while (1) {
        MPI_Recv(&tmp, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        if (tmp == finish)
            pthread_exit(0);

        pthread_mutex_lock(&mutex);
        if (task_to_do >= MIN_BALANCE) {
            tasks_to_send = task_to_do / 2;
            MPI_Send(&tasks_to_send, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
            MPI_Send(tasks + current_task + 1, tasks_to_send, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
            task_to_do -= tasks_to_send;
            current_task += tasks_to_send;
        } else {
            tasks_to_send = 0;
            MPI_Send(&tasks_to_send, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
        }
        pthread_mutex_unlock(&mutex);
    }
}

void implementing(int tasks_number) {
    pthread_mutex_lock(&mutex);
    current_task = 0;
    task_to_do = tasks_number;
    while (task_to_do != 0) {
        task_to_do--;
        int local_current_task = current_task;
        pthread_mutex_unlock(&mutex);
        tasks_implemented_amount++;
        for (int j = 0; j < tasks[local_current_task]; j++) {
            global_res += sqrt(j);
        }
        pthread_mutex_lock(&mutex);
        current_task++;
    }
    pthread_mutex_unlock(&mutex);
}

int main(int argc, char **argv) {
    int level;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &level);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    pthread_t receiver;
    pthread_attr_t attrs;

    pthread_mutex_init(&mutex, NULL);
    pthread_attr_init(&attrs);
    pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);
    pthread_create(&receiver, &attrs, receiving, NULL);
    pthread_attr_destroy(&attrs);

    int remained_tasks;
    double time_start, time_finish;

    while (iteration != ITERATIONS) {
        time_start = MPI_Wtime();
        tasks_implemented_amount = 0;

        for (int i = 0; i < TASK_AMOUNT; ++i) {
            tasks[i] = abs(rank - (iteration % size)) * 10000;
        }

        implementing(TASK_AMOUNT);
        remained_tasks = 0;
        bool tasks_available = true;
        while (tasks_available) {
            tasks_available = false;
            for (int i = (rank + 1) % size; i != rank; i = (i + 1) % size) {
                MPI_Send(&request, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Recv(&remained_tasks, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (remained_tasks != 0) {
                    MPI_Recv(tasks, remained_tasks, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    implementing(remained_tasks);
                    tasks_available = true;
                }
            }
        }
        time_finish = MPI_Wtime();
        every_iteration_local_time += time_finish - time_start;
        total_tasks_implemented_amount += tasks_implemented_amount;
        MPI_Barrier(MPI_COMM_WORLD);
        iteration++;
    }
    MPI_Send(&finish, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);

    pthread_join(receiver, NULL);
    pthread_mutex_destroy(&mutex);

    double tmp;
    MPI_Reduce(&global_res, &tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    printf("Total time by process %d: %lf Total tasks: %d\n", rank, every_iteration_local_time,
           total_tasks_implemented_amount);
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (rank == 0)
//        printf("Final result: %f\n", tmp);
    MPI_Finalize();
    return 0;
}
