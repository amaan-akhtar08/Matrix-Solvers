#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>

void execute_spec_rad_py(bool* flag);
void generate_mat(bool* flag,long dx);
double optimum_omega(long n, double Spec_rad_value);
void SOR(double *Kmat, double *Fvec, double *Xvec, long n, long iterations, double *omega, double tolerance, long *convereging_iteration);
void writeToFile(char *filename,double* Xvec,long converging_iteration,long n);
double* readFile_To_vec(char *filename,long *n);

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
    long max_iterations = 10000;
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

    Fvec = readFile_To_vec("D:/matrix-solvers-2/SORopti/input_mat/Fvec.txt",&n);
    Kmat = readFile_To_vec("D:/matrix-solvers-2/SORopti/input_mat/Kmat.txt",&size);
    
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
    SOR(Kmat, Fvec, Xvec, n, max_iterations, &omega, tolerance, &converging_iteration);
    gettimeofday(&end_solver, NULL);
    seconds = end_solver.tv_sec - start_solver.tv_sec;
    useconds = end_solver.tv_usec - start_solver.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    printf("SOR Solver Time taken: %f seconds\n", elapsed);

    writeToFile("D:/matrix-solvers-2/output/SORopti-Sol.txt",Xvec,converging_iteration,n);
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
double optimum_omega(long n, double Spec_rad_value)
{
    double pi = 3.141592;
    double w = 2 / (1 + sqrt(1 - pow(Spec_rad_value, 2)));
    return (w);
}
void SOR(double *Kmat, double *Fvec, double *Xvec, long n, long iterations, double *omega, double tolerance, long *convereging_iteration)
{
    double **Amat = (double **)malloc(n * sizeof(double *));
    double *Xvec_old = (double *)malloc(n * sizeof(double));
    for (long i = 0; i < n; i++)
    {
        Amat[i] = (double *)malloc(n * sizeof(double));
    }
    for (long i = 0; i < n; i++)
    {
        for (long j = 0; j < n; j++)
        {
            Amat[i][j] = Kmat[i * n + j];
        }
    }

    // Initial guess xi = 0
    for (long i = 0; i < n; i++)
    {
        Xvec_old[i] = 0;
        Xvec[i] = 0;
    }

    for (long iter = 0; iter < iterations; iter++)
    {
        for (long i = 0; i < n; i++)
        {
            double sum = 0;
            for (long j = 0; j < n; j++)
            {
                if (i != j)
                {
                    sum += Amat[i][j] * Xvec[j];
                }
            }

            double newValue = (Fvec[i] - sum) / Amat[i][i];
            Xvec[i] = (1 - (*omega)) * Xvec[i] + (*omega) * newValue;
        }

        // Check convergence
        double norm_diff = 0;
        for (long i = 0; i < n; i++)
        {
            norm_diff += (Xvec[i] - Xvec_old[i]) * (Xvec[i] - Xvec_old[i]);
        }
        norm_diff = sqrt(norm_diff);
        if (norm_diff < tolerance)
        {
            *convereging_iteration = iter + 1;
            break;
        }

        for (long i = 0; i < n; i++)
        {
            Xvec_old[i] = Xvec[i];
        }
    }
    for (long i = 0; i < n; i++)
    {
        free(Amat[i]);
    }
    free(Xvec_old);
    free(Amat);
    
}

