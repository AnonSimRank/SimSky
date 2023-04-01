#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/MatrixFunctions>
using namespace Eigen;

#include <iostream>
#include <chrono>
using namespace std;


#include "SimSky.h"
#include "Arnoldi.h"
#include "ApproDiag.h"

VectorXd SimSky(const SpM& mat, const MatrixXd& W, int m1, int m2, int k, double c, const SpV& query_vec, bool reorthog)
{
	const int n = mat.rows();
	MatrixXd V = MatrixXd::Zero((k + 1) * n, m2 + 1);
	MatrixXd H = MatrixXd::Zero(m2 + 1, m2);

	MatrixXd Q = MatrixXd::Zero(n, m2 + 1);
	MatrixXd P = MatrixXd::Zero(n, m2 + 1);

	// 	OrthoHess orth_hess = Arnoldi(mat, m1, query_vec, reorthog);


	auto start = chrono::system_clock::now();
	OrthoHess orth_hess = Arnoldi(mat, m1, query_vec, reorthog);
	auto end = chrono::system_clock::now();
	double dur = chrono::duration_cast<chrono::microseconds>(end - start).count();
	//cout << dur * chrono::microseconds::period::num / chrono::microseconds::period::den << endl;

	MatrixXd U_pre = orth_hess.U;
	MatrixXd T_pre = orth_hess.T;
	const int m1_pre = orth_hess.m;

	MatrixXd U_m1 = U_pre.block(0, 0, n, m1);
	MatrixXd T_m1 = T_pre.block(0, 0, m1, m1);



	// MatrixXd D = ApproDiag(mat, c, k, W);
	auto start1 = chrono::system_clock::now();
	MatrixXd D = ApproDiag(mat, c, k, W);
	auto end1 = chrono::system_clock::now();

	double dur1 = chrono::duration_cast<chrono::microseconds>(end1 - start1).count();

	//cout << dur1 * chrono::microseconds::period::num / chrono::microseconds::period::den << endl;

	SparseVector<double> e1(m1_pre);
	e1.setZero();
	e1.coeffRef(0) = 1;

	MatrixXd Tm_e1 = MatrixXd::Zero(m1, k);
	Tm_e1.col(0) = e1;

	for (int ii = 1; ii < k; ii++)
	{
		Tm_e1.col(ii) = T_m1 * Tm_e1.col(ii - 1);
	}

	MatrixXd U_Mul_T = U_pre * T_pre;
	VectorXd beta_vec = T_pre * Tm_e1.col(k - 1);
	double beta = beta_vec.norm();

	V.block(0, 0, n, 1) = 1 / beta * D.col(0).cwiseProduct(U_Mul_T * Tm_e1.col(k - 1));
	V.block(n, 0, n, 1) = 1 / beta * D.col(k).cwiseProduct(U_m1 * e1);

	MatrixXd U_T_Tm_e1 = MatrixXd::Zero(n, k - 1);

#pragma omp parallel for
	for (int jj = 0; jj < k - 1; jj++)
	{
		U_T_Tm_e1.col(jj) = U_Mul_T * Tm_e1.col(jj);

		V.block((jj + 2) * n, 0, n, 1) = 1 / beta * D.col(k - 1 - jj).cwiseProduct(U_T_Tm_e1.col(jj));
	}

	double eps = 1.0e-16;
	double tol = n * eps;

	for (int i = 0; i < m2; i++)
	{
		VectorXd V_i = V.col(i);
		VectorXd r = c * mat.transpose() * V_i.head(n) + V_i.tail(n);
		VectorXd t = V_i.head(k * n);
		VectorXd s((k + 1) * n);
		s.setZero();
		s << r, t;
		double ow = r.norm();
		//#pragma omp parallel for
		for (int j = 0; j <= i; j++)
		{
			VectorXd V_j = V.col(j);
			double temp = r.transpose() * V_j.head(n);
			s = s - temp * V_j;
			H(j, i) = temp;
		}
		if (reorthog == true) {
			for (int j = 0; j <= i; j++)
			{
				VectorXd V_j = V.col(j);
				double temp = s.head(n).transpose() * V_j.head(n);
				H(j, i) = H(j, i) + temp;
				s = s - temp * V_j;
			}
		}
		H(i + 1, i) = s.head(n).norm();
		if (H(i + 1, i) <= tol * ow)
		{
			m2 = i;
			H = H.block(0, 0, m2 + 1, m2);
			V = V.block(0, 0, (k + 1) * n, m2 + 1);
			Q = V.block(0, 0, n, m2 + 1);
			P = V.block(k * n, 0, n, m2 + 1);
		}
		V.col(i + 1) = s / H(i + 1, i);
	}
	Q = V.block(0, 0, n, m2 + 1);
	P = V.block(k * n, 0, n, m2 + 1);
	MatrixXd Hm = H.block(0, 0, m2, m2);

	SparseVector<double> e2(m2);
	e2.setZero();
	e2.coeffRef(0) = 1;

	MatrixPower<MatrixXd> Apow(Hm);
	MatrixXd Hm_power = Apow(k - 1);
	VectorXd Soar_ss = beta * Q * H * Hm_power * e2;
	return Soar_ss;
}