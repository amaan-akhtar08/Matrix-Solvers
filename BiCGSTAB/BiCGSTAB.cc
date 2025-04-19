#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>

using namespace std;

void writeToFile(const string &filename, const vector<double>& Xvec, long converging_iteration);
void execute_py(bool& flag);
vector<double> readFile_To_vec(const string &filename);
void convertToSymmetric(long n, const vector<double>& A, vector<double>& K);
bool checkSymmetry(long n, const vector<double>& A);
void multiply_Ax(long n, const vector<double>& A, const vector<double>& x, vector<double>& result);
double multiply_xy(long n, const vector<double>& x, const vector<double>& y);
void subtract_xy(long n, const vector<double>& x, const vector<double>& y, vector<double>& result);
void add_xy(long n, const vector<double>& x, const vector<double>& y, vector<double>& result);
void generate_mat(bool* flag,long dx){
    char commandline[100];
    sprintf(commandline, "D:/matrix-solvers-2/FDMmat_generator.exe %ld 0 1 0 0", dx-1);    //R, Tb, Tt, Tl, Tr
    system("g++ D:/matrix-solvers-2/FDMmat_generator.cpp -o D:/matrix-solvers-2/FDMmat_generator.exe");
    int result = system(commandline);  

    if (result == 0){
        printf("\n#C script for matrix generator executed successfully, \n");
        *flag = true;
    }
    else printf("#C script execution failed.\n");
    return;
}
void BiCGSTAB(long n, const vector<double>& A, const vector<double>& b, vector<double>& x, long max_iterations, double tolerance, long &converging_iteration) {
    vector<double> r(n), r_old(n), r0(n), p(n), Ap(n), Ax(n), s(n), As(n);

    double alpha = 0, omega = 0, beta = 0;

    multiply_Ax(n, A, x, Ax);    // Compute Ax
    subtract_xy(n, b, Ax, r);    // r = b - Ax
    for(long i = 0; i < n; i++) {
        p[i] = r[i];  // p0 = r0
    }
    for(long i = 0; i < n; i++) {
        r0[i] = r[i];  // r0*(arbitrary) = r
    }                      
    for(long i = 0; i < n; i++) {
        r_old[i] = r[i];  // r_old = r
    }                

    for (long k = 0; k < max_iterations; ++k) {
        //cout << "Iteration : " << k+1 << endl;
        multiply_Ax(n, A, p, Ap);   // Ap = A*p
        alpha = multiply_xy(n, r, r0) / multiply_xy(n, Ap, r0);
        for (long i = 0; i < n; ++i) {
            s[i] = r[i] - alpha * Ap[i];
        }
        multiply_Ax(n, A, s, As);   // As = A*s
        omega = multiply_xy(n, As, s) / multiply_xy(n, As, As);
        for (long i = 0; i < n; ++i) {
            x[i] = x[i] + alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * As[i];
        }
        if (multiply_xy(n, r, r) < tolerance) {
            converging_iteration = k + 1;
            break;
        }
        beta = multiply_xy(n, r, r0) / multiply_xy(n, r_old, r0);
        for (long i = 0; i < n; ++i) {
            p[i] = r[i] + beta * (p[i] - omega * Ap[i]);
        }
        for (long i = 0; i < n; ++i) {
            r_old[i] = r[i];
        }
    }
}

int main() {
    struct timeval start, end_program, start_solver, end_solver;
    long seconds, useconds;
    double elapsed;
    long dx;bool flag=false;

    printf("#Enter the value of dx :");
    scanf("%ld",&dx);
    generate_mat(&flag,dx);
    if(flag==false) return 1;

    // Parameters:
    long max_iterations = 1000;
    double tolerance = 1e-6;
    long n;
    long converging_iteration = max_iterations;

    vector<double> Fvec = readFile_To_vec("D:/matrix-solvers-2/BiCGSTAB/input_mat/Fvec.txt");
    vector<double> Kmat = readFile_To_vec("D:/matrix-solvers-2/BiCGSTAB/input_mat/Kmat.txt");

    n = sqrt(Kmat.size());

    vector<double> newKmat(n * n);
    newKmat = Kmat;

    vector<double> Xvec(n, 0);  // Solution vector

    gettimeofday(&start_solver, NULL);
    BiCGSTAB(n, newKmat, Fvec, Xvec, max_iterations, tolerance, converging_iteration);
    gettimeofday(&end_solver, NULL);
    seconds = end_solver.tv_sec - start_solver.tv_sec;
    useconds = end_solver.tv_usec - start_solver.tv_usec;
    elapsed = seconds + useconds / 1000000.0;
    cout << "\n# Converging Iteration :"<< converging_iteration <<endl;
    cout << "#Bi-Conjugate Gradient Stabilized Solver Time taken: " << elapsed << " seconds\n\n";

    writeToFile("D:/matrix-solvers-2/output/BiCGSTAB-Sol.txt", Xvec, converging_iteration);

    return 0;
}

void writeToFile(const string &filename, const vector<double>& Xvec, long converging_iteration) {
    cout << "#Writing to file: " << filename << endl;

    ofstream outfile(filename);
    outfile << "Converging iteration : " << converging_iteration << endl;
    for (size_t i = 0; i < Xvec.size(); ++i) {
        double value = Xvec[i];
        value = round(value * 1000000.0) / 1000000.0;
        outfile << "x" << i + 1 << " : " << value << endl;
    }
    outfile.close();
}

vector<double> readFile_To_vec(const string &filename) {
    cout << "#Reading from file: " << filename << endl;

    ifstream infile(filename);
    vector<double> vec;
    string line;
    while (getline(infile, line)) {
        vec.push_back(stod(line));
    }
    infile.close();
    return vec;
}

void convertToSymmetric(long n, const vector<double>& A, vector<double>& K) {
    for (long i = 0; i < n; ++i) {
        for (long j = i; j < n; ++j) {
            K[i * n + j] = (A[i * n + j] + A[j * n + i]) / 2;
            K[j * n + i] = K[i * n + j];
        }
    }
}

bool checkSymmetry(long n, const vector<double>& A) {
    for (long i = 0; i < n; ++i) {
        for (long j = i + 1; j < n; ++j) {
            if (A[i * n + j] != A[j * n + i]) {
                return false;
            }
        }
    }
    return true;
}

void multiply_Ax(long n, const vector<double>& A, const vector<double>& x, vector<double>& result) {
    for (long i = 0; i < n; ++i) {
        result[i] = 0;
        for (long j = 0; j < n; ++j) {
            result[i] += A[i * n + j] * x[j];
        }
    }
}

double multiply_xy(long n, const vector<double>& x, const vector<double>& y) {
    double sum = 0;
    for (long i = 0; i < n; ++i) {
        sum += x[i] * y[i];
    }
    return sum;
}

void subtract_xy(long n, const vector<double>& x, const vector<double>& y, vector<double>& result) {
    for (long i = 0; i < n; ++i) {
        result[i] = x[i] - y[i];
    }
}

void add_xy(long n, const vector<double>& x, const vector<double>& y, vector<double>& result) {
    for (long i = 0; i < n; ++i) {
        result[i] = x[i] + y[i];
    }
}
