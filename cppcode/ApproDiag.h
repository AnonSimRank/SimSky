#pragma once
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
using namespace Eigen;

typedef SparseMatrix<double> SpM;

MatrixXd ApproDiag(const SpM& mat, double c, int k, const MatrixXd& W);