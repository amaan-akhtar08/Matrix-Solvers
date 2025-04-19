#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>

void writeToFile(char *filename,double* Xvec,long converging_iteration,long n);
void execute_py(bool* flag);
double* readFile_To_vec(char *filename,long *n);
void convertToSymmetric(long n, double* A, double* K);
bool checkSymmetry(long n,double *A);
void multiply_Ax(long n, double *A, double *x, double *result);
double multiply_xy(long n, double *x, double *y);
void subtract_xy(long n, double *x, double *y, double *result);
void add_xy(long n, double *x, double *y, double *result);

void generate_mat(bool* flag,long dx){
    char commandline[100];
    sprintf(commandline, "D:/matrix-solvers-2/FDMmat_generator.exe %ld 0 1 0 0", dx-1);    //R, Tb, Tt, Tl, Tr
    system("g++ D:/matrix-solvers-2/FDMmat_generator.cpp -o D:/matrix-solvers-2/FDMmat_generator.exe");
    int result = system(commandline);  

    if (result == 0){
        printf("\n#C script for matrix generator executed successfully, ");
        *flag = true;
    }
    else printf("#C script execution failed.\n");
    return;
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
            printf("p_dot_p = 0\n");
            break;
        }

        alpha = multiply_xy(n, p, r) / p_dot_p;

        for(long i = 0; i < n; i++) {
            x[i] += alpha * r[i];   
            r[i] -= alpha * p[i];   
        }
        printf("x%ld = %lf\n",1,x[0]);

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

    // Computing average convergence rate
    double sum = 0;
    for(int i = 0; i < *converging_iteration; i++) {
        sum += convergence_rate[i];
    }
    
    *avg_convergence_rate = sum / (*converging_iteration);

    // Free allocated memory
    free(p);
    free(Ap);
    free(convergence_rate);
    free(r);
    free(Ax);
}


int main(){
    struct timeval start, end_program,start_solver, start_python, end_mat_gen, end_python, end_solver;
    long seconds, useconds;
    double elapsed;
    long dx;
    printf("#Enter the value of dx :");
    scanf("%ld",&dx);

    gettimeofday(&start, NULL);
    bool flag=false;
    generate_mat(&flag,dx);
    if(flag==false) return 1;
    gettimeofday(&end_mat_gen, NULL);
    seconds = end_mat_gen.tv_sec - start.tv_sec;
    useconds = end_mat_gen.tv_usec - start.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    printf("Time taken: %f seconds\n\n", elapsed);

    // parameters :
    long max_iterations = 10000;
    double tolerance = 1e-6;
    //auto parameters :
    long n;
    long converging_iteration = max_iterations;
    double Spec_rad_value;
    double avg_convergence_rate;
    // end of parameters//

    FILE *Fvec_file;
    FILE *Kmat_file;
    FILE *Xvec_file;
    FILE *Spec_rad;

    char line[100];
    double *Fvec = NULL;
    double *Kmat = NULL;
    // double *newKmat = NULL;
    long size = 0;

    Fvec = readFile_To_vec("D:/matrix-solvers-2/MinimumResidual/input_mat/Fvec.txt", &n);
    Kmat = readFile_To_vec("D:/matrix-solvers-2/MinimumResidual/input_mat/Kmat.txt", &size);

    double* newKmat = (double *)malloc(size * sizeof(double));
    // for(long i=0;i<size;i++) printf("%lf\n",Kmat[i]);
    
    flag = checkSymmetry(n,Kmat);
    if(flag==1){
        printf("\n#The matrix is symmetric!!!\n\n");
        newKmat = Kmat;
    }else{
        printf("#The matrix is not symmetric!!!\n");
        printf(" Converting it into (A+A^T)/2 -->\n\n");
        convertToSymmetric(n,Kmat,newKmat);
        free(Kmat);
    }
    // for(long i=0;i<size;i++) printf("%lf\n",newKmat[i]);
    double *Xvec = (double *)malloc(n * sizeof(double));  //solution vector


    gettimeofday(&start_solver, NULL);
    MinimumResidual(n,Kmat,Fvec,Xvec,max_iterations,tolerance,&converging_iteration,&avg_convergence_rate);
    gettimeofday(&end_solver, NULL);
    seconds = end_solver.tv_sec - start_solver.tv_sec;
    useconds = end_solver.tv_usec - start_solver.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    printf("#Minimum Residual Solver Time taken: %f seconds\n\n", elapsed);
    
    
    writeToFile("D:/matrix-solvers-2/output/MinimumResidual-Sol.txt",Xvec,converging_iteration,n);

    //printing parameters:
    printf("#NxN = %ldx%ld\n", n, n);
    printf(" converging iterations=%d\n", converging_iteration);
    printf(" tolerance=%e\n", tolerance);
    printf(" Average convergence rate = %lf\n\n", avg_convergence_rate);

    // gettimeofday(&start_python, NULL);
    // flag=false;
    // execute_py(&flag);  //to get the convergence rate
    // if(flag==false) return 1;
    // gettimeofday(&end_python, NULL);
    // seconds = end_python.tv_sec - start_python.tv_sec;
    // useconds = end_python.tv_usec - start_python.tv_usec;
    // elapsed = seconds + useconds / 1000000.0;
    // printf("Time taken: %f seconds\n\n", elapsed);
    // //end of printing parameters

    free(Fvec);
    free(Xvec);
    free(newKmat);
    gettimeofday(&end_program, NULL);
    seconds = end_program.tv_sec - start.tv_sec;
    useconds = end_program.tv_usec - start.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    printf("#Total Time taken: %f seconds\n", elapsed);

    return 0;
}

void execute_py(bool* flag){

    int result = system("python3 D:/matrix-solvers-2/MinimumResidual/ConvRate.py");

    if (result == 0){
        printf("#Python script executed successfully, ");
        *flag = true;
    }
    else printf("#Python script execution failed.\n");
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