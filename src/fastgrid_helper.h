#ifndef fastgrid_H
#define fastgrid_H


#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

#include <Rdefines.h>
#include <Rinternals.h>

#include <float.h> //DBL_EPSILON
#include <R_ext/Lapack.h>
#include <Rmath.h>

#define RUNIF runif
#define PRINTF Rprintf
#define MAX(A,B)    ((A) > (B) ? (A) : (B))
#define MIN(A,B)    ((A) < (B) ? (A) : (B))


using namespace scythe;
using namespace std;

//struct Eigen {
//   Matrix<> values;
//   Matrix<> vectors;
//};
//Eigen eigen (const Matrix<>& A, bool vectors=true);

void SampleReplace(int k, int n, int *y);
void SampleNoReplace(int k, int n, int *y, int *x);
void make_symmetric(double* matrix, int rows);
Matrix<> crossprod1(const Matrix<>& A);
Matrix<> tcrossprod1(const Matrix<>& A);
Matrix<> myqr_getQ (const Matrix<>& A);
Matrix<> qr_solve (const Matrix<>& A, const Matrix<>& b);

//struct Eigen {
//   Matrix<> values;
//   Matrix<> vectors;
//};
//
//Eigen eigen (const Matrix<>& A, bool vectors=true);

void _preprocess(Matrix<double, Row>& Z, Matrix<double, Row>& Y);

double Mstep_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x,vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds,
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, 
    double* logliks);

// upper hinge         
double M10_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds, 
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    double* logliks);

double M20_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds,
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    vector<double>& x2cusum, vector<double>& x3cusum, vector<double>& xrcusum, Matrix<double,Row>& xBcusum, 
    double* logliks);

double M20c_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds,
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    vector<double>& x2cusum, vector<double>& x3cusum, vector<double>& xrcusum, Matrix<double,Row>& xBcusum, 
    double* logliks);

double M22_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds, 
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    vector<double>& x2cusum, vector<double>& x3cusum, vector<double>& xrcusum, Matrix<double,Row>& xBcusum, 
    Matrix<double,Row>& BcusumR, vector<double>& rcusumR, vector<double>& WcusumR, vector<double>& xcusumR, 
    vector<double>& x2cusumR, vector<double>& x3cusumR, vector<double>& xrcusumR, Matrix<double,Row>& xBcusumR, 
    double* logliks);

double M22c_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, 
    int nThresholds, 
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    vector<double>& x2cusum, vector<double>& x3cusum, vector<double>& xrcusum, Matrix<double,Row>& xBcusum, 
    Matrix<double,Row>& BcusumR, vector<double>& rcusumR, vector<double>& WcusumR, vector<double>& xcusumR, 
    vector<double>& x2cusumR, vector<double>& x3cusumR, vector<double>& xrcusumR, Matrix<double,Row>& xBcusumR, 
    double* logliks);

double M30_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds,
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    vector<double>& x2cusum, vector<double>& x3cusum, vector<double>& x4cusum, vector<double>& x5cusum, 
    vector<double>& xrcusum, vector<double>& x2rcusum, Matrix<double,Row>& xBcusum, Matrix<double,Row>& x2Bcusum, 
    double* logliks);

double M33c_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds,
    // forward cumulative sum
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    vector<double>& x2cusum, vector<double>& x3cusum, vector<double>& x4cusum, vector<double>& x5cusum, 
    vector<double>& xrcusum, vector<double>& x2rcusum, Matrix<double,Row>& xBcusum, Matrix<double,Row>& x2Bcusum, 
    // reverse cumulative sum
    Matrix<double,Row>& BcusumR, vector<double>& rcusumR, vector<double>& WcusumR, vector<double>& xcusumR, 
    vector<double>& x2cusumR, vector<double>& x3cusumR, vector<double>& x4cusumR, vector<double>& x5cusumR, 
    vector<double>& xrcusumR, vector<double>& x2rcusumR, Matrix<double,Row>& xBcusumR, Matrix<double,Row>& x2BcusumR, 
    double* logliks);

double M111_search(
    Matrix<double,Row>& B, Matrix<double,Row>& r, vector<double>& x, vector<double>& w, 
    int * thresholdIdx, vector<double>& thresholds, int nThresholds,
    // forward cumulative sum
    Matrix<double,Row>& Bcusum, vector<double>& rcusum, vector<double>& Wcusum, vector<double>& xcusum, 
    vector<double>& x2cusum,// vector<double>& x3cusum, vector<double>& x4cusum, vector<double>& x5cusum, 
    double* ehatvec, int verbose);


#endif /* fastgrid_H */
