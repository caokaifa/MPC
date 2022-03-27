#include <iostream>
// osqp-eigen
#include "OsqpEigen/OsqpEigen.h"
// eigen
#include <Eigen/Dense>
 
using namespace std;
 
int initMat(OsqpEigen::Solver& solver)
{
    unsigned int numOfVar = 3;
    unsigned int numOfCons = 4;
    solver.data()->setNumberOfVariables(numOfVar);
    solver.data()->setNumberOfConstraints(numOfCons);
 
    Eigen::SparseMatrix<double> hessian;
    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double> linearMatrix;
    Eigen::VectorXd lowerBound;
    Eigen::VectorXd upperBound;
 
    hessian.resize(numOfVar, numOfVar);
    gradient.resize(numOfVar);
    linearMatrix.resize(numOfCons, numOfVar);
    lowerBound.resize(numOfCons);
    upperBound.resize(numOfCons);
 
 
    hessian.insert(0, 0) = 1;
    hessian.insert(0, 1) = -1;
    hessian.insert(0, 2) = 1;
    hessian.insert(1, 0) = -1;
    hessian.insert(1, 1) = 2;
    hessian.insert(1, 2) = -2;
    hessian.insert(2, 0) = 1;
    hessian.insert(2, 1) = -2;
    hessian.insert(2, 2) = 4;
   //cout << "hessian" << hessian << endl;
   /* hessian << 1, -1, 1,
              -1, 2, -2,
               1, -2, 4;*/
 
    gradient << 2, -3, 1;
    
    linearMatrix.insert(0, 0) = 1;
    linearMatrix.insert(1, 1) = 1;
    linearMatrix.insert(2, 2) = 1;
    linearMatrix.insert(3, 0) = 1;
    linearMatrix.insert(3, 1) = 1;
    linearMatrix.insert(3, 2) = 1;
    /*linearMatrix << 1, 0, 0,
                    0, 1, 0,
                    0, 0, 1,
                    1, 1, 1;*/
 
    lowerBound << 0, 0, 0, 0.5;
    upperBound << 1, 1, 1, 0.5;
 
    if (!solver.data()->setHessianMatrix(hessian)) return false;
    if (!solver.data()->setGradient(gradient)) return false;
    if (!solver.data()->setLinearConstraintsMatrix(linearMatrix)) return false;
    if (!solver.data()->setLowerBound(lowerBound)) return false;
    if (!solver.data()->setUpperBound(upperBound)) return false;
 
    return true;
}
 
int main() {
    OsqpEigen::Solver solver;
 
    // set the solver
    solver.settings()->setWarmStart(true);
 
    // instantiate the solver
    if (initMat(solver)) {
        if (!solver.initSolver()) return 1;
    }
    else {
        cout << "initilize QP solver failed" << endl;
        return 1;
    }
 
    // solve
    solver.solve();
 
    Eigen::VectorXd QPSolution;
    QPSolution = solver.getSolution();
    
    cout << "x1 = " << QPSolution[0] << endl
        << "x2 = " << QPSolution[1] << endl
        << "x3 = " << QPSolution[2] << endl;
}
