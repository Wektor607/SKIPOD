#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>

void get_walltime(double *wcTime);
int data_set(double *a, double *b, double *c, int NUM_THREADS, long long int N);
int serial_mult(double *a, double *b, double *c, int NUM_THREADS, long long int N);
int mp_mult(double *a, double *b, double *c, int NUM_THREADS, long long int N);

int mp_mult(double *a, double *b, double *c, int NUM_THREADS, long long int N)
{
    int i, j, k;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel shared(a, b, c) private(i, j, k)
    {
    #pragma omp for schedule(dynamic)
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                for (k = 0; k < N; k++)
                {
                    c[i * N + j] += a[i * N + k] * b[k * N + j];
                }
            }
        }
    }
    return 0;
}

int serial_mult(double *a, double *b, double *c, int NUM_THREADS, long long int N)
{
    int i, j, k;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            for (k = 0; k < N; k++)
            {
                c[i * N + j] += a[i * N + k] * b[k * N + j];
            }
        }
    }
    return 0;
}

int data_set(double *a, double *b, double *c, int NUM_THREADS, long long int N)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            a[i * N + j] = 2.0;
            b[i * N + j] = 3.0;
            c[i * N + j] = 0.0;
        }
    }
    return 0;
}

void get_walltime(double *wcTime)
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    *wcTime = (double)(tp.tv_sec + tp.tv_usec / 1000000.0);
}

int main(int argv, char **argc)
{
    if (argv != 3)
    {
        printf("ERROR: \nPlease enter NUM_THREADS and NODE_NUMBERA: \n./run NUM_THREADS NODE_NUMBERS\n");
        return -1;
    }

    int NUM_THREADS = atoi(argc[1]);
    long long int N = atoi(argc[2]);

    double *a = (double *)malloc(N * N * sizeof(double));
    double *b = (double *)malloc(N * N * sizeof(double));
    double *c = (double *)malloc(N * N * sizeof(double));

    double d_s, d_e;
    double d_s_mp, d_e_mp;
    double d_s_mpi, d_e_mpi;
    int i, j;

    data_set(a, b, c, NUM_THREADS, N);
    get_walltime(&d_s);
    serial_mult(a, b, c, NUM_THREADS, N);
    get_walltime(&d_e);
    data_set(a, b, c, NUM_THREADS, N);
    printf("openMP:\n");
    get_walltime(&d_s_mp);
    mp_mult(a, b, c, NUM_THREADS, N);
    get_walltime(&d_e_mp);

    printf("-----------------\nN=%lld, NUM_THREADS=%d\n", N, NUM_THREADS);
    printf("serial time: %lf\n", d_e - d_s);
    printf("openMP time: %lf\n-----------------\n\n\n", d_e_mp - d_s_mp);
    free(a);
    free(b);
    free(c);

    return 0;
}





