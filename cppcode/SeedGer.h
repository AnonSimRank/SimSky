#pragma once
#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;

// implement Theorem 1 of "Efficient Partial-Pairs SimRank Search on Large Graphs"

typedef SparseMatrix<double> SpM;
typedef SparseVector<double> SpV;

VectorXd SeedGer(SpM W, SpV v, int kmax, double c);


