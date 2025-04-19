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

void BiCGSTAB(long n, double *A, double *b, double *x, long max_iterations, double tolerance, long *converging_iteration, double *avg_convergence_rate) {
    double *r = (double *)malloc(n * sizeof(double));
    double *r_old = (double *)malloc(n * sizeof(double));
    double *r0 = (double *)malloc(n * sizeof(double));
    double *p = (double *)malloc(n * sizeof(double));
    double *Ap = (double *)malloc(n * sizeof(double));
    double *Ax = (double *)malloc(n * sizeof(double));
    double *s = (double *)malloc(n * sizeof(double));
    double *As = (double *)malloc(n * sizeof(double));

    double alpha = 0, omega =0,beta =0;

    multiply_Ax(n, A, x, Ax);    // Compute Ax
    subtract_xy(n, b, Ax, r);    // r = b - Ax
    for(long i=0;i<n;i++){
        p[i]=r[i];  //p0 = r0
    }
    for(long i=0;i<n;i++){
        r0[i]=r[i];  //r0*(arbitrary) = r
    }
    for(long i=0;i<n;i++){
        r_old[i]=r[i];  //r_old = r
    }
    for(long k=0;k<max_iterations;k++){
        printf("Iteration : %ld\n",k);
        multiply_Ax(n,A,p,Ap);   //Ap = A*p
        alpha = multiply_xy(n, r, r0) / multiply_xy(n, Ap, r0);
        for(long i=0;i<n;i++){
            s[i] = r[i] - alpha*Ap[i];
        }
        multiply_Ax(n,A,s,As);   //As = A*s
        omega = multiply_xy(n, As, s) / multiply_xy(n, As, As);
        for(long i=0;i<n;i++){
            x[i] = x[i] + alpha*p[i] + omega*s[i];
            r[i] = s[i] - omega*As[i];
        }
        if(multiply_xy(n,r,r)<tolerance){
            *converging_iteration = k+1;
            break;
        }
        beta = multiply_xy(n,r,r0) / multiply_xy(n,r_old,r0);
        for(long i=0;i<n;i++){
            p[i] = r[i] + beta*(p[i] - omega*Ap[i]);
        }
        for(long i=0;i<n;i++){
            r_old[i] = r[i];
        }
    }
    free(r);
    free(r_old);
    free(r0);
    free(p);
    free(Ap);
    free(Ax);
    free(s);
    free(As);
}


int main(){
    struct timeval start, end_program,start_solver, start_python, end_mat_gen, end_python, end_solver;
    long seconds, useconds;
    double elapsed;
    long dx;
    // printf("#Enter the value of dx :");
    // scanf("%ld",&dx);

    // gettimeofday(&start, NULL);
    bool flag=false;
    // generate_mat(&flag,dx);
    // if(flag==false) return 1;
    // gettimeofday(&end_mat_gen, NULL);
    // seconds = end_mat_gen.tv_sec - start.tv_sec;
    // useconds = end_mat_gen.tv_usec - start.tv_usec;
    // elapsed = seconds + useconds / 1000000.0;
    // printf("Time taken: %f seconds\n\n", elapsed);

    // parameters :
    long max_iterations = 10000;
    double tolerance = 1e-9;
    //auto parameters :
    long n;
    long converging_iteration;
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

    Fvec = readFile_To_vec("/home/amaan/Documents/matrix-solvers-2/BiCGSTAB/input_mat/Fvec.txt", &n);
    Kmat = readFile_To_vec("/home/amaan/Documents/matrix-solvers-2/BiCGSTAB/input_mat/Kmat.txt", &size);

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
    BiCGSTAB(n,Kmat,Fvec,Xvec,max_iterations,tolerance,&converging_iteration,&avg_convergence_rate);
    gettimeofday(&end_solver, NULL);
    seconds = end_solver.tv_sec - start_solver.tv_sec;
    useconds = end_solver.tv_usec - start_solver.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    printf("#Bi-Conjugate Gradient Stabalised Solver Time taken: %f seconds\n\n", elapsed);
    
    writeToFile("/home/amaan/Documents/matrix-solvers-2/BiCGSTAB/BiCGSTAB-Sol.txt",Xvec,converging_iteration,n);

    //printing parameters:
    printf("#NxN = %ldx%ld\n", n, n);
    printf(" converging iterations=%d\n", converging_iteration);
    printf(" tolerance=%e\n", tolerance);
    //printf(" Average convergence rate = %lf\n\n", avg_convergence_rate);

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

    int result = system("python3 /home/amaan/Documents/matrix-solvers-2/BiCGSTAB/ConvRate.py");

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