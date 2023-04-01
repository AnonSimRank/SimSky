#include <iostream>
using namespace std;

#include "arnoldi.h"


OrthoHess Arnoldi(const SpM& mat,  int m, const SpV& vec, bool reorthog)
{
	size_t col_num = mat.cols();
	size_t row_num = mat.rows();

	OrthoHess s;
	MatrixXd Ortho(row_num, m + 1), Hess(m + 1, m);
	Ortho.setZero();
	Hess.setZero();
	VectorXd w(col_num);
	w.setZero();

	Ortho.col(0) = vec / vec.norm();
	double eps = 1.0e-16;
	for (int k = 0; k < m; k++) 
	{
		w = mat * Ortho.col(k);
		double ow = w.norm();
		for (int i = 0; i <= k; i++)
		{
			Hess(i, k) = w.transpose() * Ortho.col(i);
			w = w - Hess(i, k) * Ortho.col(i);
		}
		if (reorthog == true)
		{
			for (int i = 0; i <= k; i++)
			{
				double tmp = w.transpose() * Ortho.col(i);
				w = w - tmp * Ortho.col(i);
				Hess(i, k) = Hess(i, k) + tmp;
			}
		}

		Hess(k + 1, k) = w.norm();
		double tol = row_num * eps;
		if (Hess(k + 1, k) <= tol * ow)
		{
			m = k;
			Ortho = Ortho.block(0, 0, row_num, m+1);
			Hess = Hess.block(0, 0, m+1, m);
			s.m = m;
			s.U = Ortho;
			s.T = Hess;
			return s;
		}
		Ortho.col(k+1) = w / Hess(k+1, k);
	}

	s.m = m;
	s.U = Ortho;
	s.T = Hess;
	return s;
}