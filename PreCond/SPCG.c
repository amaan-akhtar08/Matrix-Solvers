#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>

// Helper function declarations
void writeToFile(char *filename, double* Xvec, long converging_iteration, long n);
void execute_spec_rad_py(bool* flag);
double* readFile_To_vec(char *filename, long *n);
void multiply_Ax(long n, double *A, double *x, double *result);
double multiply_xy(long n, double *x, double *y);
void subtract_xy(long n, double *x, double *y, double *result);
void add_xy(long n, double *x, double *y, double *result);

// Function to compute the optimal omega value
double optimum_omega(long n, double Spec_rad_value) {
    return 2 / (1 + sqrt(1 - pow(Spec_rad_value, 2)));
}

// Function to generate matrix (using external C script)
void generate_mat(bool* flag, long dx) {
    char commandline[300];
    sprintf(commandline, "/home/amaan/Documents/matrix-solvers-2/FDMmat_generator.out %ld 0 1 0 0", dx-1);
    system("g++ /home/amaan/Documents/matrix-solvers-2/FDMmat_generator.cpp -o /home/amaan/Documents/matrix-solvers-2/FDMmat_generator.out");
    int result = system(commandline);
    if (result == 0) {
        printf("\n#Matrix generator script executed successfully.\n");
        *flag = true;
    } else {
        printf("#Matrix generator script execution failed.\n");
    }
}

// PreconditionSolve for SPCG (split diagonal approximation)
void preconditionSolve(long n, double *D, double *L, double *U, double *r, double *z) {
    // Perform forward and backward substitution using split preconditioner
    // Diagonal (D) preconditioner
    for (long i = 0; i < n; i++) {
        z[i] = r[i] / D[i];  // Diagonal preconditioner approximation
    }
}

// Function to perform the Split Preconditioned Conjugate Gradient algorithm
void splitPrecondConjugateGrad(long n, double *A, double *b, double *x, long max_iterations, double tolerance, long *converging_iteration) {
    // Split A into D (diagonal), L (lower off-diagonal), and U (upper off-diagonal)
    double *D = (double *)malloc(n * sizeof(double));  // Diagonal part of A
    double *L = (double *)malloc(n * n * sizeof(double));  // Lower part of A
    double *U = (double *)malloc(n * n * sizeof(double));  // Upper part of A
    for(long i = 0; i < n; i++) {
        printf("diag - %lf\n", A[i * n + i]);
        if (A[i * n + i] == 0) {
            
            printf("Zero diagonal element at row %ld in matrix A. Exiting...\n", i);
            exit(EXIT_FAILURE);
        }
    }

    // Extract D, L, U from A
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            if (i == j) {
                D[i] = A[i * n + j];
            } else if (i > j) {
                L[i * n + j] = A[i * n + j];  // Lower off-diagonal
                U[i * n + j] = 0;
            } else {
                U[i * n + j] = A[i * n + j];  // Upper off-diagonal
                L[i * n + j] = 0;
            }
        }
    }

    // Allocate memory for vectors
    double *r = (double *)malloc(n * sizeof(double));  // Residual vector
    double *z = (double *)malloc(n * sizeof(double));  // Preconditioned residual
    double *p = (double *)malloc(n * sizeof(double));  // Search direction
    double *q = (double *)malloc(n * sizeof(double));  // Temporary vector
    double *Ap = (double *)malloc(n * sizeof(double)); // Matrix-vector product A*p

    // Initial residual r = b - A * x
    multiply_Ax(n, A, x, Ap);
    for (long i = 0; i < n; i++) {
        r[i] = b[i] - Ap[i];
    }

    // Preconditioned residual z = M^(-1) * r
    preconditionSolve(n, D, L, U, r, z);

    writeToFile("/home/amaan/Documents/matrix-solvers-2/PreCond/z.txt", z, *converging_iteration, n);



    // Initial search direction p = z
    for (long i = 0; i < n; i++) {
        p[i] = z[i];
    }

    double rz_old = multiply_xy(n, r, z);  // r^T * z
    double rz_new, alpha, beta;

    for (long k = 0; k < max_iterations; k++) {
        // Compute A * p
        multiply_Ax(n, A, p, Ap);

        // Compute step size alpha = (r^T * z) / (p^T * A * p)
        double pAp = multiply_xy(n, p, Ap);
        alpha = rz_old / pAp;

        // Update solution x = x + alpha * p
        for (long i = 0; i < n; i++) {
            x[i] += alpha * p[i];
        }

        // Update residual r = r - alpha * A * p
        for (long i = 0; i < n; i++) {
            r[i] -= alpha * Ap[i];
        }

        // Check for convergence ||r|| < tolerance
        double norm_r = multiply_xy(n, r, r);
        if (norm_r < tolerance) {
            *converging_iteration = k + 1;
            break;
        }

        // Solve for the new preconditioned residual z = M^(-1) * r
        preconditionSolve(n, D, L, U, r, z);

        // Compute new rz_new = r^T * z
        rz_new = multiply_xy(n, r, z);

        // Compute direction update factor beta = (rz_new / rz_old)
        beta = rz_new / rz_old;

        // Update search direction p = z + beta * p
        for (long i = 0; i < n; i++) {
            p[i] = z[i] + beta * p[i];
        }

        // Update rz_old = rz_new for the next iteration
        rz_old = rz_new;
    }

    // Clean up memory
    free(D);
    free(L);
    free(U);
    free(r);
    free(z);
    free(p);
    free(q);
    free(Ap);
}

// Main function for executing SPCG with matrix generation
int main() {
    struct timeval start, end_program, start_solver, end_mat_gen;
    long seconds, useconds;
    double elapsed;

    bool flag = false;
    long dx;
    printf("Enter the value of dx: ");
    scanf("%ld", &dx);

    gettimeofday(&start, NULL);

    // Generate matrix
    generate_mat(&flag, dx);
    if (!flag) return 1;
    gettimeofday(&end_mat_gen, NULL);
    seconds = end_mat_gen.tv_sec - start.tv_sec;
    useconds = end_mat_gen.tv_usec - start.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    printf("Matrix Generation Time: %f seconds\n\n", elapsed);

    // Read the matrix and vectors from file
    long n;
    long size;
    double *Fvec = readFile_To_vec("/home/amaan/Documents/matrix-solvers-2/PreCond/input_mat/Fvec.txt", &n);
    double *Kmat = readFile_To_vec("/home/amaan/Documents/matrix-solvers-2/PreCond/input_mat/Kmat.txt", &size);

    // Solution vector
    double *Xvec = (double *)calloc(n, sizeof(double));  // Solution initialized to 0

    // Parameters
    long max_iterations = 1000;
    double tolerance = 1e-6;
    long converging_iteration = max_iterations;

    // Solve the system using SPCG
    gettimeofday(&start_solver, NULL);
    splitPrecondConjugateGrad(n, Kmat, Fvec, Xvec, max_iterations, tolerance, &converging_iteration);
    gettimeofday(&end_program, NULL);

    // Calculate time taken
    seconds = end_program.tv_sec - start_solver.tv_sec;
    useconds = end_program.tv_usec - start_solver.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    printf("SPCG Solver Time: %f seconds\n", elapsed);

    // Write results to file
    writeToFile("/home/amaan/Documents/matrix-solvers-2/output/SPCG-Sol.txt", Xvec, converging_iteration, n);

    //printing parameters:
    printf("#NxN = %ldx%ld\n", n, n);
    printf(" converging iterations=%d\n", converging_iteration);
    printf(" tolerance=%e\n\n", tolerance);

    // Clean up
    free(Fvec);
    free(Kmat);
    free(Xvec);

    return 0;
}

// Other utility functions are the same as in your original code.

void execute_spec_rad_py(bool* flag){
    // Calling Python interpreter to run a Python script

    int result = system("python3 spec-rad.py");

    if (result == 0){
        printf("Python script executed successfully, ");
        *flag = true;
    }
    else printf("Python script execution failed.\n");
    return;
}

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

void convertToSymmetric(long n, double* A, double* K){
    long initial_sum=n+1;
    for(long i=0;i<n;i++){
        long sum = initial_sum;
        for(long j=i*(n+1)+1;j<n*(i+1);j++){
            K[j]=(A[j]+A[sum-j])/2;
            K[sum-j]=K[j];
            // printf("%lf : %ld : %ld\n",K[j],j,i);
            sum+=n+1;
        }
        initial_sum+=2*(n+1);
    }
    for(long i=0;i<n;i++){    //n*i+i
        K[n*i+i] = A[n*i+i];
    }
    return;
}

bool checkSymmetry(long n,double *A){
    long initial_sum=n+1;
    for(long i=0;i<n;i++){
        long sum = initial_sum;
        for(long j=i*(n+1)+1;j<n*(i+1);j++){
            if(A[j]!=A[sum-j]){
                return 0;
            }
            sum+=n+1;
        }
        initial_sum+=2*(n+1);
    }
    return 1;
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

