#pragma once
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
using namespace Eigen;

#include <iostream>
using namespace std;

VectorXd Txt2Vec(string file_name);
