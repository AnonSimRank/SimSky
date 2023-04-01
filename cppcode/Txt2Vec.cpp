#include "Txt2Vec.h"

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
using namespace Eigen;

#include <iostream>

VectorXd Txt2Vec(string file_name)
{
	int row_num = 0;
	int MAX_DIM = 100000;
	VectorXd ss(MAX_DIM);
	ss.setZero();
	string simrank_score;
	ifstream infile;
	infile.open(file_name);
	assert(infile.is_open());
	while (getline(infile, simrank_score)) 
	{	
		ss.coeffRef(row_num) = std::stod(simrank_score);
		row_num++;
	}
	return ss.head(row_num);
}