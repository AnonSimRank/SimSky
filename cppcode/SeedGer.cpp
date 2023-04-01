#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include "SeedGer.h"

using namespace std;
using namespace Eigen;

VectorXd SeedGer(SpM W, SpV v, int kmax, double c)
{
	int mat_row_num = W.rows();
	MatrixXd U(mat_row_num, kmax + 1);
	U.setZero();

	MatrixXd V(mat_row_num, kmax + 1);
	V.setZero();

	U.col(0) = v;
	for (int L = 1; L <= kmax; L++)
	{
		U.col(L) = W * U.col(L - 1);
	}

	V.col(0) = U.col(kmax);
	for (int L = 1; L <= kmax; L++)
	{
		V.col(L) = c * W.transpose() * V.col(L - 1) + U.col(kmax - L);
	}

	VectorXd ss = (1 - c) * V.col(kmax);
	return ss;
}
