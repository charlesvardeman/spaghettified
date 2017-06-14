/* 
 * File:   main.cpp
 * Author: Cheng
 *
 * Created on July 17, 2013, 4:11 PM
 */

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream> 
#include <string.h>
#include <eigen3/Eigen/Dense>
#include <math.h>

#include "./json/json.h"
#include "cmdline.h"

using namespace std;
using namespace Eigen;

#define PI 3.1415926535897932384626433832795
#define EPS 2.2204e-16

char* inputDir;
char* outputDir;
bool DEBUG_MODE, DETERMINISTIC_ONLY;

//char dir[] = "/Users/Cheng/Documents/CRC projects/Projects/Cyber-Eye/newMatlabCodes2.0/Hawaii_data/newData1.0/data/";
//char dir[] = "/vagrant/HAKOUv3/Hawaii_data/newData1.0/data/";

/*
 * read binary data into matrix, float version
 */
int readDataIntoMatrix_float (char* rootDir, char* fileName, MatrixXf& data_Matrix) {
    long nNum;
    int rows, cols;
    char filePath[200];
    // get the full path of the data
    strcpy(filePath, rootDir);
    strcat(filePath, "/");
    strcat(filePath, fileName);
    
    FILE *fin = fopen(filePath, "rb");
    if (fin == NULL) {
        cout << "Cannot open file: " << filePath << endl;
        return 0;
    } 
    
    fread(&rows, sizeof(int), 1 , fin); // get the rows of the matrix
    fread(&cols, sizeof(int), 1 , fin); // get the columns of the matrix
    nNum = rows * cols; // get the total number of data

    data_Matrix.resize(rows, cols);
  
    float* data = new float[nNum];
    fread(data, sizeof(*data), nNum, fin);
    for (int i = 0; i < nNum; i++) {
            data_Matrix(i % rows, (int)(i / rows)) = data[i];   // matrix initialization
    }
    delete data;

    fclose(fin);
    
    return 1;
}

/*
 * read binary data into matrix, double version
 */
int readDataIntoMatrix_double (char* rootDir, char* fileName, MatrixXd& data_Matrix) {
    long nNum;
    int rows, cols;
    char filePath[200];
    // get the full path of the data
    strcpy(filePath, rootDir);
    strcat(filePath, "/");
    strcat(filePath, fileName);
    
    FILE *fin = fopen(filePath, "rb");
    if (fin == NULL) {
        cout << "Cannot open file: " << filePath << endl;
        return 0;
    } 
    
    fread(&rows, sizeof(int), 1 , fin); // get the rows of the matrix
    fread(&cols, sizeof(int), 1 , fin); // get the columns of the matrix
    nNum = rows * cols; // get the total number of data

    data_Matrix.resize(rows, cols);
  
    double* data = new double[nNum];
    fread(data, sizeof(*data), nNum, fin);
    for (int i = 0; i < nNum; i++) {
            data_Matrix(i % rows, (int)(i / rows)) = data[i];   // matrix initialization
    }
    delete data;

    fclose(fin);
    
    return 1;
}
    
/*
 * Read binary data into vector, float version
 */
int readDataIntoVector_float (char* rootDir, char* fileName, VectorXf& data_Vector) {
    long nNum;
    int rows, cols;
    char filePath[200];
    // get the full path of the data
    strcpy(filePath, rootDir);
    strcat(filePath, "/");
    strcat(filePath, fileName);
    
    FILE *fin = fopen(filePath, "rb");
    if (fin == NULL) {
        cout << "Cannot open file: " << filePath << endl;
        return 0;
    } 
    
    fread(&rows, sizeof(int), 1 , fin); // get the rows of the matrix
    fread(&cols, sizeof(int), 1 , fin); // get the columns of the matrix
    nNum = rows * cols; // get the total number of data in this vector
    
    data_Vector.resize(nNum);
   
    float* data = new float[nNum];
    fread(data, sizeof(*data), nNum, fin);
    for (int i = 0; i < nNum; i++) {
            data_Vector(i) = data[i];   // vector initialization
    }
    delete data;

    fclose(fin);
    
    return 1;
}
    
/*
 * Read binary data into vector, double version
 */
int readDataIntoVector_double (char* rootDir, char* fileName, VectorXd& data_Vector) {
    long nNum;
    int rows, cols;
    char filePath[200];
    // get the full path of the data
    strcpy(filePath, rootDir);
    strcat(filePath, "/");
    strcat(filePath, fileName);
    
    FILE *fin = fopen(filePath, "rb");
    if (fin == NULL) {
        cout << "Cannot open file: " << filePath << endl;
        return 0;
    } 
    
    fread(&rows, sizeof(int), 1 , fin); // get the rows of the matrix
    fread(&cols, sizeof(int), 1 , fin); // get the columns of the matrix
    nNum = rows * cols; // get the total number of data in this vector
    
    data_Vector.resize(nNum);
   
    double* data = new double[nNum];
    fread(data, sizeof(*data), nNum, fin);
    for (int i = 0; i < nNum; i++) {
            data_Vector(i) = data[i];   // vector initialization
    }
    delete data;

    fclose(fin);
    
    return 1;
}


/*
 * convert number to string
 */
string convertNumberToString (int inputNumber){
    std::ostringstream ostr;
    ostr << inputNumber;
    std::string returnString = ostr.str();
    return returnString;
}

// write vectors into outfile file, row by row
void writeVectorToFile (FILE *fin, VectorXd& input) {
    for (int i = 0; i < input.size(); i++) {
        fprintf (fin, "%7.4f   ,", input(i));
    }
    fprintf (fin, "\n");
}

// write json string into json file
void jsonWriter (std::ofstream& fin_results, vector<VectorXd> output_vector, vector<string> output_names) {
    
    Json::Value event;  
    Json::FastWriter fast_writer; 

    for (int i = 0; i < output_names.size(); i ++) {
        Json::Value temp(Json::arrayValue);
        
        for (int j = 0; j < output_vector.at(i).size(); j ++) {
           temp.append(Json::Value(output_vector.at(i)(j))); 
        }
        
        event[output_names.at(i)] = temp;
    }
    
    fin_results << fast_writer.write(event);
}


void calculate_core(VectorXd& b, VectorXd& r, VectorXd& y, MatrixXd& v, MatrixXd& yb, MatrixXf& B, VectorXd& par_s, VectorXd& index, MatrixXd& rand_v) {
    int ns = yb.cols();
    int N = rand_v.cols();
    // calculate ro
    MatrixXd ro(yb.rows(), ns);
    for (int i = 0; i < ns; i++) {
        ro.col(i) = y;
    }   
    ro = ro - yb;
    
    // end of calculating ro

    // calculate correlation td
    MatrixXd td(yb.rows(), ns);
    
    for (int i = 0; i < td.rows(); i++) {
        for (int j = 0; j < ns; j++) {
            td(i, j) = pow(abs(ro(i, j)), par_s(0)) * -1 * v(i, 1);
        }
    }
    
    // end of calculating td
    
    // calculate r
    r.resize(ns);
    for (int i = 0; i < ns; i++) {
        r(i) = exp(td.col(i).sum());
    }
    // end of calculating r
    
    // force release some memories
    td.resize(0,0);

    // evaluation of  basis functions at new point
    b(0) = 1;
    for (int i = 0; i < y.size(); i++)
        b(i + 1) = y (i);

    //        b(0) = 1;
    //        b(1) = y(0);
    //        b(2) = y(1);
    //        b(3) = y(2);
    //        b(4) = y(3);
    //        b(5) = y(4);
    //        b(6) = y(0) * y(0);
    //        b(7) = y(0) * y(1);
    //        b(8) = y(0) * y(2);
    //        b(9) = y(0) * y(3);
    //        b(10) = y(0) * y(4);
    //        b(11) = y(1) * y(1);
    //        b(12) = y(1) * y(2);
    //        b(13) = y(1) * y(3);
    //        b(14) = y(1) * y(4);
    //        b(15) = y(2) * y(2);
    //        b(16) = y(2) * y(3);
    //        b(17) = y(2) * y(4);
    //        b(18) = y(3) * y(3);
    //        b(19) = y(3) * y(4);
    //        b(20) = y(4) * y(4);
    int index_tracker = 0;
    for (int i = 0; i < index.size(); i++) {
        for (int j = i; j < index.size(); j++) {
            b(index.size() + 1 + index_tracker) = y(i) * y(j);
            index_tracker ++;
        }
    }
    
}

void calculate_output(MatrixXf& f, MatrixXf& fp, VectorXd& analysis, VectorXd& theta, MatrixXf& F, MatrixXd& v, MatrixXd& yb, MatrixXf& B, VectorXd& par_s, VectorXd& index, MatrixXd& rand_v, double time) {
    int ns = yb.cols();
    int N = rand_v.cols();
    
    
    if (analysis(0) == 1) {
        cout << "deterministic calculation" << endl;
        f.resize(1, F.cols());
        // calculate y. Normalize input and avoid extrapolations
        MatrixXd v_diag(v.rows(), v.cols());
        v_diag.setZero();
        for (int i = 0; i < v.rows(); i++) {
            v_diag(i, i) = v(i, 0);
        }
        VectorXd min_vector(v.rows());
        for (int i = 0; i < v.rows(); i++) {
            double temp;
            temp = v(i, 2) > theta(i) ? v(i, 2) : theta(i);
            min_vector(i) = temp > v(i, 3) ? v(i, 3) : temp;
        }
        
        VectorXd y;
        y = v_diag * min_vector;
        
        // end of y calculation

        // initialize vector b, matrix M and L
        VectorXd b(B.rows());

        // define sr and sri;
        VectorXd r;
        
        // perform the core-MLE, get b M and L
        calculate_core( b, r, y, v, yb, B, par_s, index, rand_v);
              
        f.row(0) = b.transpose().cast<float>()*B + r.transpose().cast<float>()*F;

    }
    
    if (analysis(1) == 1) {
        cout << "probabilistic calculation" << endl;
  
        fp.resize(N, F.cols());
        double taux = 1 + (time - 12)/12;
        fp.resize(N, F.cols());
        for (int n =0; n < N; n++) {
//            cout << n << endl;
            // calculate sample for Monte Carlo, normalize input and avoid extrapolations
            VectorXd min_vector(v.rows());
            for (int i = 0; i < v.rows(); i++) {
                double temp;
                temp = v(i, 2) > theta(i) + taux * v(i, 4) * rand_v(i, n) ? v(i, 2) : theta(i) + taux * v(i, 4) * rand_v(i, n);
                min_vector(i) = temp > v(i, 3) ? v(i, 3) : temp;
            }
            
            VectorXd y(v.rows());
            for (int i = 0; i < v.rows(); i++)
                y(i) = v(i, 0) * min_vector(i);
            // end of y calculation
            
            // initialize vector b, matrix M and L
            VectorXd b(B.rows());
            // define r;
            VectorXd r;

            // perform the core-MLE, get b M and L
            calculate_core( b, r, y, v, yb, B, par_s, index, rand_v);

            fp.row(n) = b.transpose().cast<float>()*B + r.transpose().cast<float>()*F;
        }
    }

}

double diffclock(clock_t clock1,clock_t clock2)
{
    double diffticks=clock1-clock2;
    double diffms=(diffticks)/CLOCKS_PER_SEC;
    return diffms;
} 

/*
 * Main function to calculate MLE
 */
int main(int argc, char** argv) {
    float flongitude, angle, cp, vel, Rmax, time;

    try {
        // Parse arguments
        gengetopt_args_info args_info;

        if (cmdline_parser (argc, argv, &args_info) != 0)
            exit(1);

        // flags
        DEBUG_MODE = args_info.debug_flag;
        DETERMINISTIC_ONLY = args_info.detonly_flag;

        if( args_info.inputs_num == 8 ) {
            int idx = 0;
            inputDir = args_info.inputs[idx++];
            outputDir = args_info.inputs[idx++];
            flongitude = atof(args_info.inputs[idx++]);
            angle = atof(args_info.inputs[idx++]);
            cp = atof(args_info.inputs[idx++]);
            vel = atof(args_info.inputs[idx++]);
            Rmax = atof(args_info.inputs[idx++]);
            time = atof(args_info.inputs[idx++]);


        } else {
            cout << "Requires 8 args." << endl;
            exit(1);
        }


    } catch (...) {
        cout << "There was an error parsing the command line arguments and initializing variables.";
        exit(1);
    }
    
//    VectorXf test(5);
//    test << 2, 1, 4, -6, 9;
//    sort(test.data(), test.data()+test.size()); 
//    cout << test;
//    exit(1);
    
//    VectorXd test(5);
//    VectorXd test2(5);
//    VectorXf final;
//    test << 2, 1, 4, -6, 9;
//    test2 << -3,-2, 9, -4, 2;
//    final = test.cast<float>() + test2.cast<float>();
//    cout << final;
    
//    exit(1);
 
//    float angle = 210;
//    float time = 30;
//    float vel = 16;
//    float cp = 950;
//    float Rmax = 40;
//    int cs = 4;
//    int sc = 1;
    // initialize theta
    VectorXd theta(5); 
    theta << flongitude, angle, cp, vel, Rmax;
    
    // choose output style. Must be integers. First is for deterministic (1 is run, 0 is not), second is for probabilistic
    VectorXd analysis(2); 
    analysis << 1, int(!DETERMINISTIC_ONLY);
    
    // initialize par_s
    VectorXd par_s; 
    readDataIntoVector_double(inputDir, (char *)"par_s.bin", par_s);
     
    // load high-fidelity data matrix
    MatrixXf F;
    readDataIntoMatrix_float(inputDir, (char *)"F.bin", F);

    // load weight matrix
    MatrixXd v;
    readDataIntoMatrix_double(inputDir, (char *)"v.bin", v);
    
    // load rand_v matrix
    MatrixXd rand_v;
    readDataIntoMatrix_double(inputDir, (char *)"rand_v.bin", rand_v);
    
    // load support points yb matrix
    MatrixXd yb;
    readDataIntoMatrix_double(inputDir, (char *)"yb.bin", yb);
    
    // load support points B matrix
    MatrixXf B;
    readDataIntoMatrix_float(inputDir, (char *)"B.bin", B);
    
    // load index vector
    VectorXd index;
    readDataIntoVector_double(inputDir, (char *)"index.bin", index);
    
    // calculate deterministic and probabilistic output
    MatrixXf f;
    MatrixXf fp;
    calculate_output(f, fp, analysis, theta, F, v, yb, B, par_s, index, rand_v, time);
    
    // force release some memories
    F.resize(0,0);
    v.resize(0,0);
//    rand_v.resize(0,0);       // don't know why this cannot be force released. When I did it, got an error
    yb.resize(0,0);
    B.resize(0,0);
    
    // PCA matrix loading for both cases
    MatrixXf aC;
    readDataIntoMatrix_float(inputDir, (char *)"aC.bin", aC);
    
    MatrixXd am;
    readDataIntoMatrix_double(inputDir, (char *)"am.bin", am);
    
    // transform the deterministic output back to original format
    MatrixXf Res;
    Res = f * aC + am.cast<float>();
    
    // transform the probabilistic output back to original format
    MatrixXf Resp;
    // if probabilistic calculation needed
    if (analysis[1] == 1) {
        
        Resp = fp * aC;
        
        for (int i = 0; i < fp.rows(); i++) {
            Resp.row(i) = Resp.row(i) + am.cast<float>();
        }
        
        // sort probabilistic results before saving it
        for (int i = 0; i < Resp.cols(); i++) {
            sort(Resp.col(i).data(), Resp.col(i).data()+Resp.col(i).size());
        }
    }
    
    // force release some memories
    f.resize(0,0);
    fp.resize(0,0);
    am.resize(0,0);
    aC.resize(0,0);
    
    // output results
    VectorXd longitude;
    readDataIntoVector_double(inputDir, (char *)"long_grid.bin", longitude);

    VectorXd latitude;
    readDataIntoVector_double(inputDir, (char *)"lat_grid.bin", latitude);
    
    // save all the output vectors within one big data frame
    vector<VectorXd> output_vectors;
    vector<string> output_names;
    output_vectors.push_back(theta);
    output_names.push_back("theta");

    output_vectors.push_back(par_s);
    output_names.push_back("par_s");

    output_vectors.push_back(longitude);
    output_names.push_back("long");

    output_vectors.push_back(latitude);
    output_names.push_back("lat");
    
    VectorXd latlongSize(2); 
    latlongSize << latitude.size(), longitude.size();
    output_vectors.push_back(latlongSize);
    output_names.push_back("lat_long_size");
    
    // write json string to output file
    std::ofstream fin_results(string(string(outputDir)+"results.json").c_str());
    jsonWriter (fin_results, output_vectors, output_names);
    fin_results.close();
    
    // output deterministic results
    int res_row = 1;
    int res_col = Res.size();
    std::ofstream fout(string(string(outputDir)+"result_deter.bin").c_str(), ios::out|ios::binary);
    if(!fout) {
        cout << "Cannot open result_deter.bin.";
        exit(1);
    }

    fout.write((char *) (&res_row), sizeof(res_row));
    fout.write((char *) (&res_col), sizeof(res_col));
    fout.write((char *) Res.data(), Res.rows() * Res.cols() * sizeof(float));
    fout.close();

    if (analysis[1] == 1) {
        int resp_row = rand_v.cols();
        int resp_col = (int)Resp.size()/resp_row;
        std::ofstream fout(string(string(outputDir)+"result_prob.bin").c_str(), ios::out|ios::binary);
        if(!fout) {
            cout << "Cannot open result_prob.bin.";
            exit(1);
        }
        fout.write((char *) (&resp_row), sizeof(resp_row));
        fout.write((char *) (&resp_col), sizeof(resp_col));
        fout.write((char *) Resp.data(), Resp.rows() * Resp.cols() * sizeof(float));
        fout.close();
    }
  
    return 0;
}



