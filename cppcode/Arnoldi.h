#pragma once
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
using namespace Eigen;

typedef SparseMatrix<double> SpM;
typedef SparseVector<double> SpV;

struct OrthoHess
{
	MatrixXd U, T;
	int m;
};

OrthoHess Arnoldi(const SpM& mat,  int m, const SpV& vec, bool reorthog);
