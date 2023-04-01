#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;

#include "VariedD_SR.h"
 
MatrixXd VariedD::Compute_Diagonal_Matrix()
{
	int row_num = Matrix_Rows();

	VectorXd D_0(row_num);
	D_0.setOnes();

	MatrixXd D(row_num, k + 1);
	D.setZero();
	D.col(0) = D_0;

	for (int j = 1; j <= k; j++)
	{
		for (int i = 0; i < row_num; i++)
		{
			SparseVector<double> e_i(row_num);
			e_i.coeffRef(i) = 1.0;
			double t = 0.0;
			
			MatrixXd temp_mat(row_num, j + 1);
			temp_mat.setZero();
			temp_mat.col(0) = e_i;

			for (int L = 1; L <= j; L++)
			{
				temp_mat.col(L) = sqrt(c) * P * temp_mat.col(L-1);
				t = t + temp_mat.col(L).cwiseProduct(temp_mat.col(L)).dot(D.col(j - L));
			}
			D.coeffRef(i, j) = 1 - t;
		}
	}
	return D;

}


VectorXd VariedD::Optimized_Single_Source(int query)
{
	MatrixXd D = Compute_Diagonal_Matrix();
	int row_num = Matrix_Rows();

	MatrixXd mat_x(row_num, k + 1);
	mat_x.setZero();
	SparseVector<double> query_vec(row_num);
	query_vec.coeffRef(query) = 1;
	mat_x.col(0) = query_vec;

	for (int L = 1; L <= k; L++)
	{
		mat_x.col(L) = P * mat_x.col(L - 1);
	}

	MatrixXd mat_y(row_num, k + 1);
	mat_y.setZero();
	VectorXd y_0 = D.col(0).cwiseProduct(mat_x.col(k));
	mat_y.col(0) = y_0;

	for (int L = 1; L <= k; L++)
	{
		mat_y.col(L) = c * P.transpose() * mat_y.col(L - 1) + D.col(L).cwiseProduct(mat_x.col(k - L));
	}

	return mat_y.col(k);

}
