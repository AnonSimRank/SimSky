#pragma once
#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SpM;

SpM read_mtx_lib(string fname, SpM W);
SpM read_mtx(string fname);
SpM read_txt(string fname, bool first_);

