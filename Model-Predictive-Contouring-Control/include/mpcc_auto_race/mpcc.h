//
// Created by yuwei on 11/19/19.
//
// ros
#include <ros/ros.h>
#include <ros/package.h>
#include <sensor_msgs/LaserScan.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/PointStamped.h>
#include <geometry_msgs/Pose.h>
#include <geometry_msgs/Point.h>
#include <ackermann_msgs/AckermannDriveStamped.h>
#include <nav_msgs/OccupancyGrid.h>
#include <tf/transform_listener.h>

// standard
#include <math.h>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <random>
#include <mpcc_auto_race/spline.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

#include <nav_msgs/Odometry.h>

#include <mpcc_auto_race/track.h>
#include <Eigen/Sparse>
#include "OsqpEigen/OsqpEigen.h"

using namespace std;
using namespace Eigen;

const int nx = 4;
const int nu = 2;
const double CAR_LENGTH = 0.35;


enum rviz_id{
    CENTERLINE,
    CENTERLINE_SPLINE,
    THETA_EST,
    PREDICTION,
    BORDERLINES
};


typedef struct Stage{
    Matrix<double, nx, nx> Q;
    Matrix<double, nu, nu> R;
    Matrix<double, nx, 1> fx;
    Matrix<double, nu, 1> fu;
    Matrix<double, nx, nx> Ad;
    Matrix<double, nx, nu> Bd;
    Matrix<double, nx, 1> dv_cost_dx;
}Stage;


class MPCC{

public:
    MPCC(ros::NodeHandle& nh);

private:
    ros::NodeHandle nh_;
    ros::Publisher track_viz_pub_;
    ros::Publisher mpcc_viz_pub_;
    ros::Publisher drive_pub_;
    ros::Subscriber odom_sub_;

    Track track_;

    Eigen::VectorXd QPSolution_;
    geometry_msgs::Pose car_pose_;

   // OsqpEigen::Solver solver_;
    bool reset_flag;   // reset linearization point to x0
    double theta_prev_;

    /* Paramaters */
    string pose_topic;

    double q_c;
    double q_l;
    double q_v;
    double r_v;
    double r_steer;

    double last_steer_cmd;
    int N;  // horizon
    double Ts;

    double SPEED_MAX;
    double SPEED_MIN;
    double STEER_MAX;
    string drive_topic;

    double speed_m_;
    vector<geometry_msgs::Point> border_lines;


    void getParams(ros::NodeHandle& nh);

    void odom_callback(const nav_msgs::Odometry::ConstPtr &odom_msg);
    void visualize_centerline();

    void getMatrices(Matrix<double,nx,nx>& Q, Matrix<double,nx, 1>& fx, Matrix<double, nx, 1>& dv_cost_dx,  Matrix<double,nu, 1>& fu, Matrix<double,nx,nx>& Ad,
            Matrix<double,nx, nu>& Bd, Matrix<double,nx,1>& x_prev, Matrix<double,nu,1>& u_prev);

//    void constructQPHessian(const Matrix<double,nx, nx>& Q, const DiagonalMatrix<double,4>& R, int horizon, SparseMatrix<double> &hessianMatrix);
//
//    void constructConstraintMatrix(const Matrix<double,nx, nx>& Q, int horizon, Matrix<double,nx,nx>& Ad,
//            Matrix<double,nx,nu>& Bd, SparseMatrix<double>& constraintMatrix);
//
//    void constructConstraintVectors(const Matrix<double,nx,1>& xMax, const Matrix<double,nx,1>& xMin, const Matrix<double,nx,1>& x0,
//            int horizon, VectorXd& lowerBound, VectorXd& upperBound);

    void initialize_solution(Eigen::VectorXd& x0);

    void solveQP(Eigen::VectorXd& QPSolution_shifted, Eigen::VectorXd& QP_x0, bool reset);

    Eigen::VectorXd simulate_dynamics(Eigen::VectorXd& state, Eigen::VectorXd& input, double dt);

    void applyControl();
    void visualize_mpcc_solution();
};