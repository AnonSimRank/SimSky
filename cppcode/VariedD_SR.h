#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
using namespace Eigen;

typedef SparseMatrix<double> SpM;


class VariedD
{
public:
	VariedD(SpM P_, int k_, double c_) : P(P_), k(k_), c(c_) {}

	MatrixXd Compute_Diagonal_Matrix();
	VectorXd Optimized_Single_Source(int query);

	inline int Matrix_Rows()
	{
		return P.rows();
	}

private:
	SpM P;
	int k;
	double c;
};


