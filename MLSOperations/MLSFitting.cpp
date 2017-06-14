/*
 Files Required:
 * 
 * theta_input.bin
 * Hsin_output.bin
 * Hsout_output.bin
 * index_Hsin.bin
 * index_Hsout.bin
 * 
 */


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include "eigen-eigen-3.0.3/Eigen/Dense"

using namespace std;
using namespace Eigen;

int MLSFitting(double val1, double val2, double val3, double val4, double val5, char* inputdir, char* outputdir, char* outfile1, char* outfile2 ) {
    
    string inputDir = inputdir;
    string outputDir = outputdir;
    
    VectorXf theta_new (5);
    theta_new << val1, val2, val3, val4, val5;
    
    // load high-fidelity theta matrix
    long nNum;
    string infile = "theta_input.bin";
    FILE *fin = fopen(string( inputDir + infile).c_str(), "rb");
    if (fin == NULL) {
        printf("Theta file was not opened");
        return 0;
    } else {
        fseek ( fin , 0 , SEEK_END );
        nNum = ftell (fin) / sizeof(float);
        fseek ( fin , 0 , SEEK_SET );
    }
    MatrixXf theta_Matrix (5, nNum/5);
    float* theta = new float[nNum];
    fread(theta, sizeof(*theta), nNum * 3, fin);
    for (int i = 0; i < nNum; i++) {
            theta_Matrix(i % 5, (int)(i / 5)) = theta[i];
    }
    fclose(fin);
    delete theta;

    // calculate mean and std for each theta row
    VectorXf v_std (5);
    VectorXf v_mean (5);
    for (int i = 0; i < 5; i++) {
        v_mean(i) = theta_Matrix.row(i).mean();
        VectorXf temp_std_vector = theta_Matrix.row(i);
        for (int j = 0; j < nNum/5; j++) {
            temp_std_vector (j) = pow((float)(temp_std_vector (j) - v_mean(i)) , 2);
        }
        v_std (i) = sqrt(temp_std_vector.sum()/(nNum/5 - 1));
    }

    // nomalization of the high fidelity theta matrix
    MatrixXf yb (5, nNum/5);
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < nNum/5; j++) {
            yb (i, j) = (theta_Matrix(i, j) - v_mean(i))/v_std(i);
        }
    }

    int n = yb.rows();  // number of rows of the high fidelity theta matrix
    int ns = yb.cols();  // number of columns of the high fidelity theta matrix
    // Initialize B matrix.
    MatrixXf B (ns, 15); //1+5+9 
    for (int i = 0; i < ns; i++) {
        B (i, 0) = 1;   
        B (i, 1) = yb (0, i);
        B (i, 2) = yb (1, i);
        B (i, 3) = yb (2, i);
        B (i, 4) = yb (3, i);
        B (i, 5) = yb (4, i);
        B (i, 6) = yb (0, i) * yb (0, i);
        B (i, 7) = yb (0, i) * yb (1, i);
        B (i, 8) = yb (0, i) * yb (2, i);
        B (i, 9) = yb (0, i) * yb (3, i);
        B (i, 10) = yb (1, i) * yb (1, i);
        B (i, 11) = yb (1, i) * yb (2, i);
        B (i, 12) = yb (1, i) * yb (3, i);
        B (i, 13) = yb (2, i) * yb (2, i);
        B (i, 14) = yb (3, i) * yb (3, i);
    }

    float c = 0.4;
    int k = 2;
    int n_d = 200;  // number of high fidelity data used for fitting. Distance will be Set based on this number
    MatrixXf v (5 ,5);
    v.setIdentity();
    v(0,0) = 2;
    v(1,1) = 2;

    VectorXf y (5);
    y = (theta_new - v_mean);
    for (int i = 0; i < 5; i++) {
        y(i) = y(i)/v_std(i);
    }

    MatrixXf rep_matrix (5, nNum/5);
    VectorXf temp_row (nNum/5);
    temp_row.setOnes();
    for (int i = 0; i < 5; i++) {
        rep_matrix.row(i) = temp_row * y(i);
    }
    
    VectorXf ro (nNum/5);
    for (int i = 0; i < nNum/5; i++) {
        ro(i) = (v*(rep_matrix.col(i) - yb.col(i))).transpose()*(v*(rep_matrix.col(i) - yb.col(i)));
        ro(i) = sqrt(ro(i));
    }

    VectorXf ro_index_number (nNum/5);
    for (int i = 0; i < nNum/5 - 1; i++) {
        ro_index_number (i) = i;
    }
    
    for (int i = 0; i < nNum/5 - 1; i++) {
        for (int j = 0; j < nNum/5 - i - 1; j++) {
            if (ro(j) > ro(j+1)) {
                float temp = ro(j);
                ro(j) = ro(j+1);
                ro(j+1) = temp;

                temp = ro_index_number(j);
                ro_index_number(j) = ro_index_number(j+1);
                ro_index_number(j+1) = temp;
            }
        } 
    }

    float D = ro(n_d - 1);
    D = 3.1519;
    MatrixXf W (n_d, n_d);
    W.setIdentity();
    for (int i = 0; i < n_d; i++) {
        W(i, i) = (exp(-1*pow(1*ro(i)/D/c, k)) - exp(-1*pow(1/c, k)))/(1 - exp(-1*pow(1/c,k)));
    }

    MatrixXf Ba (n_d, 15);
    for (int i = 0; i < n_d; i++) {
        Ba.row(i) = B.row(ro_index_number(i));
    }

    VectorXf b (15);
    b(0) = 1;
    b(1) = y(0);
    b(2) = y(1);
    b(3) = y(2);
    b(4) = y(3);
    b(5) = y(4);
    b(6) = y(0) * y(0);
    b(7) = y(0) * y(1);
    b(8) = y(0) * y(2);
    b(9) = y(0) * y(3);
    b(10) = y(1) * y(1);
    b(11) = y(1) * y(2);
    b(12) = y(1) * y(3);
    b(13) = y(2) * y(2);
    b(14) = y(3) * y(3);

    MatrixXf M (15, 15);
    M = Ba.transpose()*W*Ba;

    MatrixXf L (15, n_d);
    L = Ba.transpose()*W;

    MatrixXf aux;
    aux = b.transpose()*M.inverse()*L;

    long nNum_Hsin;
    infile = "Hsin_output.bin";
    FILE *fin_Hsin = fopen(string( inputDir + infile).c_str(), "rb");
    if (fin_Hsin == NULL) {
        printf("Hsin file was not opened");
        return 0;
    } else {
        fseek ( fin_Hsin , 0 , SEEK_END );
        nNum_Hsin = ftell (fin_Hsin) / sizeof(float);
        fseek ( fin_Hsin , 0 , SEEK_SET );
    }
    MatrixXf Hsin_matrix (nNum/5, nNum_Hsin/(nNum/5));
    float* Hsin = new float[nNum_Hsin];
    fread(Hsin, sizeof(*Hsin), nNum_Hsin * 3, fin_Hsin);
    for (int i = 0; i < nNum_Hsin; i++) {
            Hsin_matrix(i % (nNum/5), (int)(i / (nNum/5))) = Hsin[i];
    }
    fclose(fin_Hsin);
    delete Hsin;

    long nNum_Hsout;
    infile = "Hsout_output.bin"; 
    FILE *fin_Hsout = fopen(string( inputDir + infile).c_str(), "rb");
    if (fin_Hsout == NULL) {
        printf("Hsout file was not opened");
        return 0;
    } else {
        fseek ( fin_Hsout , 0 , SEEK_END );
        nNum_Hsout = ftell (fin_Hsout) / sizeof(float);
        fseek ( fin_Hsout , 0 , SEEK_SET );
    }
    MatrixXf Hsout_matrix (nNum/5, nNum_Hsout/(nNum/5));
    float* Hsout = new float[nNum_Hsout];
    fread(Hsout, sizeof(*Hsout), nNum_Hsout * 3, fin_Hsout);
    for (int i = 0; i < nNum_Hsout; i++) {
            Hsout_matrix(i % (nNum/5), (int)(i / (nNum/5))) = Hsout[i];
    }
    fclose(fin_Hsout);
    delete Hsout;

    MatrixXf Hsin_matrix_nd (n_d, nNum_Hsin/(nNum/5));
    for (int i = 0; i < n_d; i++) {
        Hsin_matrix_nd.row(i) = Hsin_matrix.row(ro_index_number(i));
    }
    VectorXf Fin (nNum_Hsin/(nNum/5));
    for (int i = 0; i < nNum_Hsin/(nNum/5); i++){
        float temp_sum = 0;
        for (int j = 0; j < 200; j++) {
            temp_sum = temp_sum + aux(j)*log(Hsin_matrix_nd(j,i));
        }

        Fin(i) = exp(temp_sum);
    }

    MatrixXf Hsout_matrix_nd (n_d, nNum_Hsout/(nNum/5));
    for (int i = 0; i < n_d; i++) {
        Hsout_matrix_nd.row(i) = Hsout_matrix.row(ro_index_number(i));
    }
    VectorXf Fout (nNum_Hsout/(nNum/5));
    for (int i = 0; i < nNum_Hsout/(nNum/5); i++){
        float temp_sum = 0;
        for (int j = 0; j < 200; j++) {
            temp_sum = temp_sum + aux(j)*log(Hsout_matrix_nd(j,i));
        }

        Fout(i) = exp(temp_sum);
    }
//    Fout = exp(aux * log(Hsout_matrix_nd));

    int nx_Hsin = 100;
    int ny_Hsin = 150;
    int nx_Hsout = 80;
    int ny_Hsout = 120;
    
    infile = "index_Hsin.bin";
    FILE *fin_index_Hsin = fopen(string( inputDir + infile).c_str(), "rb");
    if (fin_index_Hsin == NULL) {
        printf("index_Hsin file was not opened");
        return 0;
    }
    int* index_Hsin = new int [nNum_Hsin/(nNum/5)];
    fread(index_Hsin, sizeof(*index_Hsin), nNum_Hsin/(nNum/5) * 3, fin_index_Hsin);
    
    infile = "index_Hsout.bin";
    FILE *fin_index_Hsout = fopen(string( inputDir + infile).c_str(), "rb");
    if (fin_index_Hsout == NULL) {
        printf("index_Hsin file was not opened");
        return 0;
    }
    int* index_Hsout = new int [nNum_Hsout/(nNum/5)];
    fread(index_Hsout, sizeof(*index_Hsout), nNum_Hsout/(nNum/5) * 3, fin_index_Hsout);


    MatrixXf results_Hsin (nx_Hsin, ny_Hsin);
    results_Hsin.setZero();
    for (int k = 0; k < nNum_Hsin/(nNum/5); k++) {
        int x = (index_Hsin[k] - 1) % nx_Hsin;
        int y = (int)((index_Hsin[k] - 1)/nx_Hsin);
        results_Hsin (x, y) = Fin(k);
    }
    
//    string outfile = "results_Hsin.csv";
    FILE *fin_results_Hsin = fopen( string( outputDir + outfile1).c_str(), "w");
    if (fin_results_Hsin == NULL) {
        printf("results_Hsin file was not opened");
        return 0;
    }
    for (int i = 0; i < nx_Hsin; i ++) {
        for (int j = 0; j < ny_Hsin; j++) {
           fprintf (fin_results_Hsin, "%7.4f   ,", results_Hsin(i,j));
        }
        fprintf (fin_results_Hsin, "\n");
    }
    fclose(fin_results_Hsin);

    MatrixXf results_Hsout (nx_Hsout, ny_Hsout);
    results_Hsout.setZero();
    for (int k = 0; k < nNum_Hsout/(nNum/5); k++) {
        int x = (index_Hsout[k] - 1) % nx_Hsout;
        int y = (int)((index_Hsout[k] - 1)/nx_Hsout);
        results_Hsout (x, y) = Fout(k);
    }
    
//    outfile = "results_Hsout.csv";
    FILE *fin_results_Hsout = fopen(string( outputDir + outfile2).c_str(), "w");
    if (fin_results_Hsout == NULL) {
        printf("results_Hsout file was not opened");
        return 0;
    }
    for (int i = 0; i < nx_Hsout; i ++) {
        for (int j = 0; j < ny_Hsout; j++) {
           fprintf (fin_results_Hsout, "%7.4f   ,", results_Hsout(i,j));
        }
        fprintf (fin_results_Hsout, "\n");
    }
    fclose(fin_results_Hsout);

    delete index_Hsin;
    delete index_Hsout;
    return 0;
    
}