#ifndef QUADRICOPTIMIZER_H_
#define QUADRICOPTIMIZER_H_

#include <iostream>
#include <Eigen/Sparse>

class QuadricOptimizer
{
public:
	QuadricOptimizer();
	~QuadricOptimizer();

	//solve 1/2ftQf subject to 
    void solve(Eigen::SparseMatrix<double> &Q, Eigen::VectorXd &c, Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b, Eigen::VectorXd &x);

	/* data */
	int NUMCON;
	int NUMVAR;
	int NUMANZ;
	int NUMQNZ;
};

#endif
