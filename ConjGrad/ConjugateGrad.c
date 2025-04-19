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
void conjugateGrad(long n, double *A, double*b, double *x, long max_iterations, double tolerance, long *converging_iteration, double *avg_convergence_rate) {
    double *r = (double *)malloc(n * sizeof(double));
    double *Ax = (double *)malloc(n * sizeof(double));
    double *Ar = (double *)malloc(n * sizeof(double));
    double *p = (double *)malloc(n * sizeof(double));
    double *Ap = (double *)malloc(n * sizeof(double));
    double *convergence_rate = (double *)calloc(n, sizeof(double));
    // Initialize residual r = b - Ax
    multiply_Ax(n, A, x, Ap); 
    subtract_xy(n, b, Ap, r);
    
    // Initialize search direction p = r
    for (int i = 0; i < n; ++i) {
        p[i] = r[i];
    }

    double rTr = 0.0;
    double rTr_new = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    
    // Compute initial rTr
    rTr = multiply_xy(n, r, r);
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        multiply_Ax(n, A, p, Ap); // Ap = A*p
        
        // Compute alpha
        double pAp = 0;
        pAp = multiply_xy(n, p, Ap);
        alpha = rTr / pAp;
        
        // Update x and r
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        
        rTr_new = multiply_xy(n, r, r);

        convergence_rate[iter]= rTr_new/rTr;
        if (sqrt(rTr_new) < tolerance) {
            *converging_iteration = iter + 1;
            break;
        }
        beta = rTr_new / rTr;
        for (int i = 0; i < n; ++i) {
            p[i] = r[i] + beta * p[i];
        }
        rTr = rTr_new;
    }
    double sum=0;
    for(int i=0;i<*converging_iteration;i++){
        sum+=convergence_rate[i];
    }
    *avg_convergence_rate = sum/(*converging_iteration);

    free(r);
    free(Ax);
    free(Ar);
    free(p);
    free(Ap);
    free(convergence_rate);
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
    long max_iterations = 1000;
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

    Fvec = readFile_To_vec("D:/matrix-solvers-2/ConjGrad/input_mat/Fvec.txt", &n);
    Kmat = readFile_To_vec("D:/matrix-solvers-2/ConjGrad/input_mat/Kmat.txt", &size);

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
    double *Xvec = (double *)calloc(n , sizeof(double));  //solution vector


    gettimeofday(&start_solver, NULL);
    conjugateGrad(n,newKmat,Fvec,Xvec,max_iterations,tolerance,&converging_iteration,&avg_convergence_rate);
    gettimeofday(&end_solver, NULL);
    seconds = end_solver.tv_sec - start_solver.tv_sec;
    useconds = end_solver.tv_usec - start_solver.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    printf("#Conjugate Gradient Solver Time taken: %f seconds\n\n", elapsed);
    
    writeToFile("D:/matrix-solvers-2/output/ConjGrad-Sol.txt",Xvec,converging_iteration,n);

    //printing parameters:
    printf("#NxN = %ldx%ld\n", n, n);
    printf(" converging iterations=%d\n", converging_iteration);
    printf(" tolerance=%e\n", tolerance);
    printf(" Average convergence rate = %lf\n\n", avg_convergence_rate);

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

    int result = system("python3 D:/matrix-solvers-2/ConjGrad/ConvRate.py");

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