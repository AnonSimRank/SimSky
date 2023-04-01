#include <iostream>
#include <string>
#include <chrono>
using namespace std;

#include <iomanip>
#include "ReadData.h"
#include "Arnoldi.h"
#include "ApproDiag.h"
#include "SimSky.h"
#include "Txt2Vec.h"


#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/MatrixFunctions>
using namespace Eigen;

typedef SparseMatrix<double> SpM;
typedef SparseVector<double> SpV;


int main() {

	string file_txt = "./test.txt";
	SpM sparse_mat = read_txt(file_txt, true);
	const int row_num = sparse_mat.rows();
	SpV vec(row_num);
	vec.setZero();
	vec.coeffRef(1) = 1;
	int m1 = 5;
	int m2 = 5;
	int k = 5;
	double c = 0.8;
	int s = row_num;
	MatrixXd W(row_num, s);
	W.setZero();
	W.block(0, 0, s, s) = MatrixXd::Identity(s, s);
	VectorXd t(row_num - s);
	t.setZero();
	for (int i = 0; i < row_num - s; i++)
	{
		if (i % 2 == 0) {
			t.coeffRef(i) = -1;
		}
		else {
			t.coeffRef(i) = 1;
		}
	}
	W.block(s, s - 1, row_num - s, 1) = t;

	
	VectorXd soar_ss = SimSky(sparse_mat, W, m1, m2, k, c, vec, false);
	cout << soar_ss << endl;
	


	

	system("pause");
}




