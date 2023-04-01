#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
using namespace Eigen;

#include "ApproDiag.h"
#include <iostream>
#include <string>
using namespace std;

MatrixXd ApproDiag(SpM mat, double c, int k, MatrixXd W)
{
	const int col_num = mat.cols();
	const int n = mat.rows();
	const int W_col_num = W.rows();
	const int s = W.cols();

	VectorXd one_vec_s = VectorXd::Ones(s);
	MatrixXd W_W = W.cwiseProduct(W);
	VectorXd denom = W_W * one_vec_s;

	MatrixXd D = MatrixXd::Zero(n,k + 1);
	VectorXd one_vec_n = VectorXd::Ones(n);
	D.col(0) = one_vec_n;


	/*for (int j = 0; j < k; j++)
	{
		VectorXd nume = VectorXd::Zero(n);
		for (int i = 0; i < s; i++)
		{
			VectorXd w_i = W.col(i);
			MatrixXd x_vec_cell = MatrixXd::Zero(W_col_num, j + 1);
			x_vec_cell.col(0) = w_i;
			for (int a = 1; a <= j; a++)
			{
				x_vec_cell.col(a) = mat * x_vec_cell.col(a-1);
			}
			MatrixXd y_vec_cell = MatrixXd::Zero(W_col_num, j + 1);
			y_vec_cell.col(0) = D.col(0).cwiseProduct(x_vec_cell.col(j));
			if (j == 0) {
				y_vec_cell.col(j+1) = c * mat.transpose() * y_vec_cell.col(j);
			}
			else
			{
				for (int b = 1; b <= j; b++)
				{
					y_vec_cell.col(b) = D.col(b).cwiseProduct(x_vec_cell.col(j+1 - b)) + c * mat.transpose() * y_vec_cell.col(b - 1);
				}
				y_vec_cell.col(j) = c * mat.transpose() * y_vec_cell.col(j-1);
			}
			nume = nume + w_i.cwiseProduct(y_vec_cell.col(j));
		}
		D.col(j + 1) = one_vec_n - nume.cwiseQuotient(denom);
	}
	return D;*/

	for (int j = 1; j <= k; j++)
	{
		VectorXd nume = VectorXd::Zero(n);
		for (int i = 0; i < s; i++) 
		{
			VectorXd w_i = W.col(i);
			MatrixXd x_vec_cell = MatrixXd::Zero(W_col_num, j + 1);
			x_vec_cell.col(0) = w_i;
			for (int a = 1; a <= j; a++)
			{
				x_vec_cell.col(a) = mat * x_vec_cell.col(a - 1);
			}
			MatrixXd y_vec_cell = MatrixXd::Zero(W_col_num, j + 1);
			y_vec_cell.col(0) = D.col(0).cwiseProduct(x_vec_cell.col(j));
			if (j == 1)
			{
				y_vec_cell.col(j) = c * mat.transpose() * y_vec_cell.col(j - 1);
			}
			else
			{
				for (int b = 1; b <= j - 1; b++)
				{
					y_vec_cell.col(b) = c * mat.transpose() * y_vec_cell.col(b - 1) + D.col(b).cwiseProduct(x_vec_cell.col(j - b));
				}
				y_vec_cell.col(j) = c * mat.transpose() * y_vec_cell.col(j - 1);

			}
			nume = nume + w_i.cwiseProduct(y_vec_cell.col(j));
		}
		
		D.col(j) = one_vec_n - nume.cwiseQuotient(denom);
	}
	return D;
}

