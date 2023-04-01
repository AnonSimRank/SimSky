#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
using namespace Eigen;

#include "ApproDiag.h"
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <chrono>
using namespace std;

MatrixXd ApproDiag(const SpM& mat, double c, int k, const MatrixXd& W)
{
	auto tp1 = chrono::system_clock::now();


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

	auto cmt = c * mat.transpose();


	auto tp2 = chrono::system_clock::now();

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

	unordered_map<int, MatrixXd> x_vec_cell_map(s * 2);

	for (int i = 0; i < s; i++) 
	{
		x_vec_cell_map.emplace(i, MatrixXd::Zero(W_col_num, k + 1));
	}

#pragma omp parallel for
	for (int i = 0; i < s; i++)
	{
		const VectorXd& w_i = W.col(i);
		MatrixXd& x_vec_cell = x_vec_cell_map[i];

		x_vec_cell.col(0) = w_i;

		for (int a = 1; a <= k; a++)
		{
			x_vec_cell.col(a) = mat * x_vec_cell.col(a - 1);
		}
	}


	auto tp3 = chrono::system_clock::now();

#pragma omp parallel for
	for (int j = 1; j <= k; j++)
	{
		std::vector<VectorXd> numes(s);

		VectorXd nume = VectorXd::Zero(n);

		for (int i = 0; i < s; i++) 
		{
			const VectorXd& w_i = W.col(i);
			const MatrixXd& x_vec_cell = x_vec_cell_map[i];

			MatrixXd y_vec_cell = MatrixXd::Zero(W_col_num, j + 1);
			y_vec_cell.col(0) = D.col(0).cwiseProduct(x_vec_cell.col(j));
			if (j == 1)
			{
				y_vec_cell.col(j) = cmt * y_vec_cell.col(j - 1);
			}
			else
			{
				for (int b = 1; b < j; b++)
				{
					y_vec_cell.col(b) = cmt * y_vec_cell.col(b - 1) + D.col(b).cwiseProduct(x_vec_cell.col(j - b));
				}

				y_vec_cell.col(j) = cmt * y_vec_cell.col(j - 1);
			}

			nume += w_i.cwiseProduct(y_vec_cell.col(j));
		}

		D.col(j) = one_vec_n - nume.cwiseQuotient(denom);
	}

	auto tp4 = chrono::system_clock::now();

	double dur1 = chrono::duration_cast<chrono::microseconds>(tp2 - tp1).count();
	double dur2 = chrono::duration_cast<chrono::microseconds>(tp3 - tp2).count();
	double dur3 = chrono::duration_cast<chrono::microseconds>(tp4 - tp3).count();

	//double aaa = chrono::microseconds::period::num / chrono::microseconds::period::den;

	//cout << dur1 * chrono::microseconds::period::num / chrono::microseconds::period::den << "\t\t" 
	//	<< dur2 * chrono::microseconds::period::num / chrono::microseconds::period::den << "\t\t" 
	//	<< dur3 * chrono::microseconds::period::num / chrono::microseconds::period::den << endl;


	return D;
}

