#include <Eigen/Dense>
using namespace Eigen;

#include <stdio.h>
#include "DiagDom.h"


// Determine matrix is strictly diagonally dominant whether or not
void DiagDom(MatrixXd P)
{
	int rows_num = P.rows();
	int col_num = P.cols();

	P = P.cwiseAbs();

	for (int i = 0; i < rows_num; i++)
	{
		if (2*P.coeffRef(i,i) >= P.row(i).sum())
		{
			continue;
		}
		else
		{
			printf("matrix P is not diagonally dominant\n");
		}
	}

}