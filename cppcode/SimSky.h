#pragma once
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
using namespace Eigen;


typedef SparseMatrix<double> SpM;
typedef SparseVector<double> SpV;

VectorXd SimSky(const SpM& mat, const MatrixXd& W, int m1, int m2, int k, double c, const SpV& query_vec, bool reorthog);

