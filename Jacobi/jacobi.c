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
void jacobi(double* Kmat, double* Fvec, double* Xvec,int n,int iterations,double tolerance, int* convereging_iteration){
    double **Amat = (double **)malloc(n * sizeof(double *));
    double * Xvec_old = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        Amat[i] = (double *)malloc(n * sizeof(double));
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            Amat[i][j] = Kmat[i*n+j];
        }
    }
    //Initial guess xi = 0
    for(int i=0;i<n;i++){
        Xvec_old[i] = 0;
        Xvec[i] = 0;
    }

    for(int iter=0;iter<iterations;iter++){
        for(int i=0;i<n;i++){
            double sum = 0;
            for(int j=0;j<n;j++){
                if(i!=j){
                    sum += Amat[i][j]*Xvec_old[j];
                }
            }
            Xvec[i] = (Fvec[i] - sum)/Amat[i][i];
        }
        // Check convergence
        double norm_diff = 0;
        for (int i = 0; i < n; i++) {
            norm_diff += (Xvec[i] - Xvec_old[i]) * (Xvec[i] - Xvec_old[i]);
        }
        norm_diff = sqrt(norm_diff);
        
        if (norm_diff < tolerance) {
            *convereging_iteration = iter + 1;
            break;
        }

        for(int i=0;i<n;i++){
            Xvec_old[i] = Xvec[i];
        }
    }
    for (int i = 0; i < n; i++) {
        free(Amat[i]);
    }
    free(Amat);
    free(Xvec_old);
}
int main() {

    //parameters :
    int n;
    int iterations = 1000;
    int convereging_iteration=1000;
    double tolerance = 1e-6;
    //end of parameters//
    long dx;
    printf("# Enter the value of dx :");
    scanf("%ld",&dx);

    n = (dx-1)*(dx-1);

    bool flag =false;
    generate_mat(&flag,dx);
    if(flag==false) return 1;

    FILE *Fvec_file;
    FILE *Kmat_file;
    FILE *Xvec_file;

    char line[100];
    double *Fvec = NULL;
    double *Kmat = NULL;
    double * Xvec = (double *)malloc(n * sizeof(double));
    int size = 0;

    Fvec_file = fopen("input_mat/Fvec.txt", "r");
    if (Fvec_file == NULL) {
        printf("Failed to open the Fvec_file.\n");
        return 1;
    }

    while (fgets(line, sizeof(line), Fvec_file)) {
        double value = atof(line);
        size++;
        Fvec = realloc(Fvec, size * sizeof(double));
        Fvec[size - 1] = value;
    }
    fclose(Fvec_file);

    // for (int i = 0; i < size; i++) {
    //     printf("%d : %.2f\n", i+1, Fvec[i]);
    // }
    
//
    
    size = 0;

    Kmat_file = fopen("input_mat/Kmat.txt", "r");
    if (Kmat_file == NULL) {
        printf("Failed to open the Kmat_file.\n");
        return 1;
    }

    while (fgets(line, sizeof(line), Kmat_file)) {
        double value = atof(line);
        size++;
        Kmat = realloc(Kmat, size * sizeof(double));
        Kmat[size - 1] = value;
    }
    fclose(Kmat_file);
//
    // for (int i = 0; i < size; i++) {
    //     printf("%d : %.2f\n", i+1, Kmat[i]);
    // }

    jacobi(Kmat, Fvec, Xvec, n, iterations, tolerance, &convereging_iteration);

    Xvec_file = fopen("D:/matrix-solvers-2/output/Jacobi-Sol.txt", "w");
    if (Xvec_file == NULL) {
        printf("Failed to open the Xvec_output file.\n");
        return 1;
    }
    fprintf(Xvec_file, "Converging iteration : %d\n", convereging_iteration);
    for (int i = 0; i < n; i++) {
        fprintf(Xvec_file, "x%d : %f\n", i + 1, Xvec[i]);
    }
    fclose(Xvec_file);

    for (int i = 0; i < n; i++) {
        printf("x%d : %f\n", i+1, Xvec[i]);
    }

    free(Fvec);
    free(Kmat);
    free(Xvec);

    return 0;
}