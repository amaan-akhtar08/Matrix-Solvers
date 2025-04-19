#include<stdio.h>
#include<stdlib.h>
#include <bits/stdc++.h>
using namespace std;

int main(int argc, char *argv[]){
    if (argc != 6) {
        printf("Usage: %s R Tb Tt Tl Tr\n", argv[0]);  //
        return 1; 
    }

    long R = atof(argv[1]);
    double Tb = atof(argv[2]);
    double Tt = atof(argv[3]);
    double Tl = atof(argv[4]);
    double Tr = atof(argv[5]);

    long MESH[R][R];
    long ii,jj;
    for(jj=0;jj<R;jj++){
        for(ii=0;ii<R;ii++){
            MESH[ii][jj]=jj*R+ii;
        }
    }
    // printf("*\n");
    long n;
    n=R*R;
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);
    long i,j,k;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            A[i][j]=0;
            // printf("i=%ld j=%ld\n",i,j);
        }
        b[i]=0;
    }
    // printf("*\n");
    for(jj=1;jj<(R-1);jj++){
        for(ii=1;ii<(R-1);ii++){
            A[MESH[ii][jj]][MESH[ii+1][jj]]=1;
            A[MESH[ii][jj]][MESH[ii-1][jj]]=1;
            A[MESH[ii][jj]][MESH[ii][jj+1]]=1;
            A[MESH[ii][jj]][MESH[ii][jj-1]]=1;
            A[MESH[ii][jj]][MESH[ii][jj]]=-4;
        }
    }
    // printf("*\n");

    for(jj=1;jj<R-1;jj++){
        A[MESH[0][jj]][MESH[1][jj]]=1;
        A[MESH[0][jj]][MESH[0][jj+1]]=1;
        A[MESH[0][jj]][MESH[0][jj-1]]=1;
        A[MESH[0][jj]][MESH[0][jj]]=-4;
        b[MESH[0][jj]]=-Tl;
    }
    for(jj=1;jj<R-1;jj++){
        A[MESH[R-1][jj]][MESH[R-2][jj]]=1;
        A[MESH[R-1][jj]][MESH[R-1][jj+1]]=1;
        A[MESH[R-1][jj]][MESH[R-1][jj-1]]=1;
        A[MESH[R-1][jj]][MESH[R-1][jj]]=-4;
        b[MESH[R-1][jj]]=-Tr;
    }
    for(ii=1;ii<R-1;ii++){
        A[MESH[ii][0]][MESH[ii+1][0]]=1;
        A[MESH[ii][0]][MESH[ii-1][0]]=1;
        A[MESH[ii][0]][MESH[ii][1]]=1;
        A[MESH[ii][0]][MESH[ii][0]]=-4;
        b[MESH[ii][0]]=-Tb;
    }
    for(ii=1;ii<R-1;ii++){
        A[MESH[ii][R-1]][MESH[ii+1][R-1]]=1;
        A[MESH[ii][R-1]][MESH[ii-1][R-1]]=1;
        A[MESH[ii][R-1]][MESH[ii][R-2]]=1;
        A[MESH[ii][R-1]][MESH[ii][R-1]]=-4;
        b[MESH[ii][R-1]]=-Tt;
    }
    A[MESH[0][0]][MESH[0][1]]=1;
    A[MESH[0][0]][MESH[1][0]]=1;
    A[MESH[0][0]][MESH[0][0]]=-4;
    b[MESH[0][0]]=-Tl-Tb;
    A[MESH[R-1][0]][MESH[R-2][0]]=1;
    A[MESH[R-1][0]][MESH[R-1][1]]=1;
    A[MESH[R-1][0]][MESH[R-1][0]]=-4;
    b[MESH[R-1][0]]=-Tr-Tb;
    A[MESH[R-1][R-1]][MESH[R-1][R-2]]=1;
    A[MESH[R-1][R-1]][MESH[R-2][R-1]]=1;
    A[MESH[R-1][R-1]][MESH[R-1][R-1]]=-4;
    b[MESH[R-1][R-1]]=-Tr-Tt;
    A[MESH[0][R-1]][MESH[0][R-2]]=1;
    A[MESH[0][R-1]][MESH[1][R-1]]=1;
    A[MESH[0][R-1]][MESH[0][R-1]]=-4;
    b[MESH[0][R-1]]=-Tl-Tt;

    double x0[n],x1[n],r[n],p[n];

    FILE *FDM_file;
    FDM_file = fopen("input_mat/Kmat.txt", "w");
    for (long i = 0; i < n; i++){
        for(j=0;j<n;j++){
            fprintf(FDM_file, "%lf\n",A[i][j]); //
        }
    }
    fclose(FDM_file);
    FILE *FDMb_file;
    FDMb_file = fopen("input_mat/Fvec.txt", "w");
    for (long i = 0; i < n; i++){
        fprintf(FDM_file, "%lf\n",b[i]); //
    }
    fclose(FDMb_file);
    return 0;
}