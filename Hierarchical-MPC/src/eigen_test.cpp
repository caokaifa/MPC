#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;

int main()

{

  MatrixXd m(2,2);

  m(0,0) = 3;

  m(1,0) = 2.5;

  m(0,1) = -1;

  m(1,1) = m(1,0) + m(0,1);

  std::cout << "Here is the matrix m:\n" << m << std::endl;

  VectorXd v(2);

  v(0) = 4;

  v(1) = v(0) - 1;

  std::cout << "Here is the vector v:\n" << v << std::endl;


  // v.segment(0, 1)<<9;
  
  // std::cout << "---Here is the vector v:\n" << v << std::endl;
  v.head(1)<<-1;
  std::cout<<"----88888-----"<< v <<std::endl;

  // Matrix<double,3,3> A;
  // Matrix<double,3,2> B;
  // Matrix<double,3,3> A, M12;
  // Matrix<double,3,3> B;


  // std::cout<<"----666---"<<Matrix3d::Identity()<<std::endl;
  // A <<   0.0, 0.0, 1,
  //           0.0, 0.0, 2,
  //           0.0, 0.0, 0.0;

  // B <<   4, 0.0,
  //           5, 0.0,
  //           6, 7;

  // Matrix<double,6,6> aux, M;
  // aux.setZero();
  // aux.block<3,3>(0,0) << A;//矩阵快操作，

  // // cout<<"999999999999999"<<nx<<endl;
  // aux.block<3,3>(0, 3) << Matrix3d::Identity();

  // std::cout<<"----8888---"<<aux<<std::endl;

  // M = (aux*0.1).exp();
  // M12 = M.block<3,3>(0,3);


std::cout<<Vector2d(1, -1).normalized()<<std::endl;

}
