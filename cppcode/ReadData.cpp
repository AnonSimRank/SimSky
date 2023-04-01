#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include "ReadData.h"

using namespace std;
using namespace Eigen;

//.mtx includes node i, node j, weight w
SpM read_mtx_lib(string fname, SpM W)
{
	loadMarket(W, fname);
	for (int i = 0; i < W.outerSize(); ++i)
	{
		double sum = W.col(i).sum();
		for (SpM::InnerIterator it(W, i); it; ++it)
		{
			it.valueRef() /= sum;
		}
	}
	return W;
}

//.mtx includes node i,node j,node start from 1
SpM read_mtx(string fname)
{
	ifstream fin(fname);
	int M, N, L;
	while (fin.peek() == '%')
	{
		fin.ignore(2048, '\n');
	}

	fin >> M >> N >> L;
	SpM mat1(M, N);
	mat1.reserve(L);

	vector<Eigen::Triplet<double>> triple;
	for (int i = 0; i < L; ++i)
	{
		int m, n;
		fin >> m >> n;
		triple.push_back(Triplet<double>(m-1, n-1, 1));
	}
	fin.close();
	mat1.setFromTriplets(triple.begin(), triple.end());

	for (int i = 0; i < mat1.outerSize(); ++i)
	{
		double sum = mat1.col(i).sum();
		for (SpM::InnerIterator it(mat1, i); it; ++it)
		{
			it.valueRef() /= sum;
		}
	}
	return mat1;
}


//txt file contains FromNode i, EndNode j
// argument first_ determine node count start from 0 or 1 
//if first_ is true, node count start from 0
//if first_ is false, node count start from 1
// add node_number node_number edge_number manually after "# FromNodeId	ToNodeId"
SpM read_txt(string fname, bool first_)
{
	ifstream fin(fname);
	int M, N, L;
	while (fin.peek() == '#')
	{
		fin.ignore(2048, '\n');
	}

	fin >> M >> N >> L;
	SpM mat(M, N);
	mat.reserve(L);

	vector<Eigen::Triplet<double>> triple;
	for (int i = 0; i < L; ++i)
	{
		int m, n;
		fin >> m >> n;

		if (first_ == true)
		{
			triple.push_back(Triplet<double>(m, n, 1));
		}
		else
		{
			triple.push_back(Triplet<double>(m-1, n-1, 1));
		}
	}
	fin.close();
	mat.setFromTriplets(triple.begin(), triple.end());
	for (int i = 0; i < mat.outerSize(); ++i)
	{
		double sum = mat.col(i).sum();
		for (SpM::InnerIterator it(mat, i); it; ++it)
		{
			it.valueRef() /= sum;
		}
	}
	return mat;
}


