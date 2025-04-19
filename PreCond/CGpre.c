#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>

void writeToFile(char *filename,double* Xvec,long converging_iteration,long n);
void execute_spec_rad_py(bool* flag);
double* readFile_To_vec(char *filename,long *n);
void convertToSymmetric(long n, double* A, double* K);
bool checkSymmetry(long n,double *A);
void multiply_Ax(long n, double *A, double *x, double *result);
double multiply_xy(long n, double *x, double *y);
void subtract_xy(long n, double *x, double *y, double *result);
void add_xy(long n, double *x, double *y, double *result);
double optimum_omega(long n, double Spec_rad_value)
{
    double pi = 3.141592;
    double w = 2 / (1 + sqrt(1 - pow(Spec_rad_value, 2)));
    return (w);
}

void generate_mat(bool* flag,long dx){
    char commandline[300];
    sprintf(commandline, "/home/amaan/Documents/matrix-solvers-2/FDMmat_generator.out %ld 0 1 0 0", dx-1);    //R, Tb, Tt, Tl, Tr
    system("g++ /home/amaan/Documents/matrix-solvers-2/FDMmat_generator.cpp -o /home/amaan/Documents/matrix-solvers-2/FDMmat_generator.out");
    int result = system(commandline);  

    if (result == 0){
        printf("\n#C script for matrix generator executed successfully, ");
        *flag = true;
    }
    else printf("#C script execution failed.\n");
    return;
}
void preconditionSolve(long n, double **M, double *r, double *z);

void SORPreconditioner(long n, double *A, double **M, double omega) {
    // Assume M is an identity matrix initially
    for (long i = 0; i < n; i++) {
        M[i] = (double *)calloc(n, sizeof(double));
        M[i][i] = 1.0;
    }
    
    // Modify M based on the SOR method
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            if (i != j) {
                M[i][j] = -1*omega * A[i * n + j];
            }
        }
        M[i][i] = 1.0 / (omega * A[i * n + i]); // Relaxation factor
    }
    
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            if(M[i][j]==0){
                M[i][j]+=0.000001;
            }
        }
    }
    double *Marray = (double *)malloc(n * n * sizeof(double));
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            Marray[i * n + j] = M[i][j];
        }
    }
    writeToFile("/home/amaan/Documents/matrix-solvers-2/PreCond/M.txt", Marray, 0, n);
}
void conjugateGradPreconditioned(double n, double *newKmat, double *Fvec, double *Xvec, double max_iterations, double tolerance,
                                long *converging_iteration, double *avg_convergence_rate, double **preconditioner)
{
    // Allocate memory for vectors
    double *r = (double *)malloc(n * sizeof(double));  // Residual vector
    double *z = (double *)malloc(n * sizeof(double));  // Preconditioned residual
    double *p = (double *)malloc(n * sizeof(double));  // Search direction
    double *q = (double *)malloc(n * sizeof(double));  // Temporary vector for A * p
    double *Ap = (double *)malloc(n * sizeof(double)); // Matrix-vector product A*p

    // Initial residual r = Fvec - A * Xvec
    multiply_Ax(n, newKmat, Xvec, Ap);
    for (long i = 0; i < n; i++) {
        r[i] = Fvec[i] - Ap[i];
    }

    // Solve for the preconditioned residual z = M^(-1) * r
    preconditionSolve(n, preconditioner, r, z);  // M^(-1) * r

    writeToFile("/home/amaan/Documents/matrix-solvers-2/PreCond/z.txt", z, *converging_iteration, n);

    // Initial search direction p = z
    for (long i = 0; i < n; i++) {
        p[i] = z[i];
    }

    double rz_old = multiply_xy(n, r, z);  // Initial residual inner product
    double rz_new, alpha, beta;

    for (long k = 0; k < max_iterations; k++) {
        // Compute A * p
        multiply_Ax(n, newKmat, p, Ap);

        // Compute step size alpha = (r^T * z) / (p^T * A * p)
        double pAp = multiply_xy(n, p, Ap);
        // if (k==0)printf("pAp = %lf\n", pAp);
        alpha = rz_old / pAp;

        // Update solution Xvec = Xvec + alpha * p
        for (long i = 0; i < n; i++) {
            Xvec[i] += alpha * p[i];
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
        preconditionSolve(n, preconditioner, r, z);

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
    free(r);
    free(z);
    free(p);
    free(q);
    free(Ap);
}
void preconditionSolve(long n, double **M, double *r, double *z) {
    // Perform forward and backward substitution to solve M * z = r
    // This would depend on the structure of M (SGS or another preconditioner)

    // Forward substitution for lower triangular system (D + L)
    for (long i = 0; i < n; i++) {
        z[i] = r[i];
        for (long j = 0; j < i; j++) {
            z[i] -= M[i][j] * z[j];
        }
        z[i] /= M[i][i];
    }

    // Backward substitution for upper triangular system (D + U)
    for (long i = n - 1; i >= 0; i--) {
        for (long j = i + 1; j < n; j++) {
            z[i] -= M[i][j] * z[j];
        }
        z[i] /= M[i][i];
    }
    printf("z[0] = %lf\n", z[0]);
    
}


void conjugateGrad(long n, double *A, double*b, double *x, long max_iterations, double tolerance, long *converging_iteration,double omega) {
    double **precond = (double **)malloc(n*sizeof(double *));
    double avg_convergence_rate;
    SORPreconditioner(n, A, precond, omega);
    conjugateGradPreconditioned(n, A, b, x, max_iterations, tolerance, converging_iteration, &avg_convergence_rate, precond);
    for(long i=0;i<n;i++){
        free(precond[i]);
    }
    free(precond);
}
int main()
{
    struct timeval start, end_program,start_solver, start_python, end_mat_gen, end_python, end_solver;
    long seconds, useconds;
    double elapsed;

    bool flag=false;
    long dx;
    printf("Enter the value of dx: ");
    scanf("%ld",&dx);
    // Record the start time
    gettimeofday(&start, NULL);

    generate_mat(&flag,dx);
    if(flag==false) return 1;
    gettimeofday(&end_mat_gen, NULL);
    seconds = end_mat_gen.tv_sec - start.tv_sec;
    useconds = end_mat_gen.tv_usec - start.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    printf("Time taken: %f seconds\n\n", elapsed);


    gettimeofday(&start_python, NULL);
    flag=false;
    execute_spec_rad_py(&flag);  //to get the spectral radius
    if(flag==false) return 1;
    gettimeofday(&end_python, NULL);
    seconds = end_python.tv_sec - start_python.tv_sec;
    useconds = end_python.tv_usec - start_python.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    printf("Time taken: %f seconds\n\n", elapsed);

    // parameters :
    long max_iterations = 30;
    double tolerance = 1e-6;
    //auto parameters :
    long n;
    double omega;       //omega automatically calculated in SOR optimised
    long converging_iteration = max_iterations;
    double Spec_rad_value;
    // end of parameters//

    FILE *Spec_rad;
    char line[100];
    double *Fvec = NULL;
    double *Kmat = NULL;
    long size = 0;

    Fvec = readFile_To_vec("/home/amaan/Documents/matrix-solvers-2/PreCond/input_mat/Fvec.txt", &n);
    Kmat = readFile_To_vec("/home/amaan/Documents/matrix-solvers-2/PreCond/input_mat/Kmat.txt", &size);
    
    Spec_rad = fopen("spectral_radius_TJ.txt", "r");
    if (Spec_rad == NULL)
    {
        printf("Failed to open the spectral_radius_TJ.\n");
        return 1;
    }
    while (fgets(line, sizeof(line), Spec_rad))
    {
        double value = atof(line);
        Spec_rad_value = value;
    }
    fclose(Spec_rad);

    double *Xvec = (double *)malloc(n * sizeof(double));  //solution vector
    omega = optimum_omega(n,Spec_rad_value);


    gettimeofday(&start_solver, NULL);
    conjugateGrad(n,Kmat,Fvec,Xvec,max_iterations,tolerance,&converging_iteration,omega);
    gettimeofday(&end_solver, NULL);
    seconds = end_solver.tv_sec - start_solver.tv_sec;
    useconds = end_solver.tv_usec - start_solver.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    printf("SOR Solver Time taken: %f seconds\n", elapsed);

    writeToFile("/home/amaan/Documents/matrix-solvers-2/output/PreCondCG-Sol.txt",Xvec,converging_iteration,n);
    //printing parameters:
    printf("NxN = %ldx%ld\n", n,n);
    printf("converging iterations=%ld\n", converging_iteration);
    printf("omega=%lf\n", omega);
    printf("tolerance=%e\n\n", tolerance);
    //end of printing parameters

    free(Fvec);
    free(Kmat);
    free(Xvec);

    gettimeofday(&end_program, NULL);
    seconds = end_program.tv_sec - start.tv_sec;
    useconds = end_program.tv_usec - start.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    printf("Total Time taken: %f seconds\n", elapsed);

    return 0;
}

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

