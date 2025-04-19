#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>


void writeToFile(char *filename,double* Xvec,long converging_iteration,long n){

    printf("#Writing to file: %s\n\n", filename);

    FILE *_file;
    char line[100];
    double *vec = NULL;
    long size = 0;
    _file = fopen(filename, "r");
    _file = fopen(filename, "w");
    fprintf(_file, "Converging iteration : %d\n", converging_iteration);
    for (long i = 0; i < n; i++)
    {
        fprintf(_file, "x%d : %f\n", i + 1, Xvec[i]);
    }
    fclose(_file);
}

double* readFile_To_vec(char *filename,long *n){

    printf("#Reading from file: %s\n", filename);

    FILE *_file;
    char line[100];
    double *vec = NULL;
    long size = 0;
    _file = fopen(filename, "r");
    if (_file == NULL)
    {
        printf(" Failed to open the Kmat_file.\n");
    }
    while (fgets(line, sizeof(line), _file))
    {
        double value = atof(line);
        size++;
        vec = realloc(vec, size * sizeof(double));
        vec[size - 1] = value;
    }
    *n = size;
    fclose(_file);
    return vec;
}

void multiply_Ax(long n, double *A, double *x, double *result){
    for(long i=0;i<n;i++){
        result[i]=0;
        for(long j=0;j<n;j++){
            result[i] += A[n*i+j]*x[j];
        }
    }
}

double multiply_xy(long n, double *x, double *y){
    double sum=0;
    for(long i=0;i<n;i++){
        sum += x[i]*y[i];
    }
    return sum;
}

void subtract_xy(long n, double *x, double *y, double *result){
    for(long i=0;i<n;i++){
        result[i] = x[i]-y[i];
    }
}

void add_xy(long n, double *x, double *y, double *result){
    for(long i=0;i<n;i++){
        result[i] = x[i]+y[i];
    }
}

void MinimumResidual(long n, double *A, double *b, double *x, long max_iterations, double tolerance, long *converging_iteration, double *avg_convergence_rate) {
    double *r = (double *)malloc(n * sizeof(double));
    double *Ax = (double *)malloc(n * sizeof(double));
    double *p = (double *)malloc(n * sizeof(double));
    double *Ap = (double *)malloc(n * sizeof(double));
    double *convergence_rate = (double *)calloc(max_iterations, sizeof(double)); 

    // if (r == NULL || Ax == NULL || p == NULL || Ap == NULL || convergence_rate == NULL) {
    //     fprintf(stderr, "Memory allocation error\n");
    //     return;
    // }

    double alpha = 0;
    double rTr_old = 0;
    double rTr_new = 0; 
    double residual_norm = 0;
    double prev_residual_norm = 0;

    // Initialize x to 0
    for(long i = 0; i < n; i++) {
        x[i] = 0;
    }

    multiply_Ax(n, A, x, Ax);    // Compute Ax
    subtract_xy(n, b, Ax, r);    // r = b - Ax
    rTr_old = multiply_xy(n, r, r);
    prev_residual_norm = sqrt(rTr_old);

    for(long k = 0; k < max_iterations; k++) {
        multiply_Ax(n, A, r, p);    // p = A * r
        double p_dot_p = multiply_xy(n, p, p);
        
        if (p_dot_p == 0) {
            *converging_iteration = k + 1;
            break;
        }

        alpha = multiply_xy(n, p, r) / p_dot_p;

        for(long i = 0; i < n; i++) {
            x[i] += alpha * r[i];   
            r[i] -= alpha * p[i];   
        }

        rTr_new = multiply_xy(n, r, r);   
        residual_norm = sqrt(rTr_new);
        convergence_rate[k] = residual_norm / prev_residual_norm;

        if (residual_norm < tolerance) {
            *converging_iteration = k + 1;
            break;
        }

        rTr_old = rTr_new;
        prev_residual_norm = residual_norm;
    }
    // Free allocated memory
    free(p);
    free(Ap);
    free(convergence_rate);
    free(r);
    free(Ax);
}


int main(){

    long max_iterations = 100000;
    double tolerance = 1e-6;
    long n;
    long converging_iteration;
    double Spec_rad_value;
    double avg_convergence_rate;
    FILE *Fvec_file;
    FILE *Kmat_file;
    FILE *Xvec_file;
    FILE *Spec_rad;

    char line[100];
    double *Fvec = NULL;
    double *Kmat = NULL;
    long size = 0;

    Fvec = readFile_To_vec("/home/amaan/Documents/matrix-solvers-2/MinimumResidual/input_mat/Fvec.txt", &n);
    Kmat = readFile_To_vec("/home/amaan/Documents/matrix-solvers-2/MinimumResidual/input_mat/Kmat.txt", &size);

    double* newKmat = (double *)malloc(size * sizeof(double));
    newKmat = Kmat;
    double *Xvec = (double *)malloc(n * sizeof(double));
    MinimumResidual(n,Kmat,Fvec,Xvec,max_iterations,tolerance,&converging_iteration,&avg_convergence_rate);
    writeToFile("/home/amaan/Documents/matrix-solvers-2/MinimumResidual/MinimumResidual-Sol.txt",Xvec,converging_iteration,n);

    free(Fvec);
    free(Xvec);
    free(newKmat);

    return 0;
}
