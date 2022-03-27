//
// Created by yuwei on 12/16/19.
//
//
// Created by yuwei on 11/19/19.
//
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
#include <tf2_ros/transform_listener.h>
#include <cmath> 
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
#include <hmpc_auto_race/spline.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

#include <nav_msgs/Odometry.h>

#include <hmpc_auto_race/track.h>
#include <Eigen/Sparse>
#include "OsqpEigen/OsqpEigen.h"
#include <hmpc_auto_race/occupancy_grid.h>
#include <mutex>

using namespace std;
using namespace Eigen;

const int nx = 3;
const int nu = 2;
const double CAR_LENGTH = 0.35;


enum rviz_id{
    CENTERLINE,
    CENTERLINE_POINTS,
    THETA_EST,
    PREDICTION,
    BORDERLINES,
    TRAJECTORY_REF,
    TRAJECTORIES,
    MAX_THETA,
    DEBUG
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


class HMPC{

public:
    HMPC(ros::NodeHandle& nh);

private:
    ros::NodeHandle nh_;
    ros::Publisher track_viz_pub_;
    ros::Publisher trajectories_viz_pub_;
    ros::Publisher hmpc_viz_pub_;
    ros::Publisher drive_pub_;
    ros::Publisher debugger_pub_;
    ros::Publisher obstacle_pub_;

    ros::Subscriber odom_sub_;
    ros::Subscriber rrt_sub_;
    ros::Subscriber map_sub_;
    ros::Subscriber scan_sub_;
    /*Paramaters*/
    string pose_topic;
    string drive_topic;
    string scan_topic;
    double Ts;
    int speed_num;
    int steer_num;
    int N;
    double SPEED_MAX;
    double STEER_MAX;
    double ACCELERATION_MAX;
    double DECELERATION_MAX;
    double MAP_MARGIN;
    double SPEED_THRESHOLD;
    // MPC params
    double q_x;
    double q_y;
    double q_yaw;
    double r_v;
    double r_steer;
    double q_s;
    Matrix<double, nx, nx> Q;
    Matrix<double, nu, nu> R;


    mutex m;
    Track track_;
    tf::Transform tf_;
    tf::Vector3 car_pos_;
    double yaw_;
    double car_theta_;
    double speed_m_;
    vector<vector<vector<Vector3d>>> trajectory_table_;
    vector<Vector3d> trajectory_ref_;
    Vector2d input_ref_;
    double last_speed_cmd_;

    vector<geometry_msgs::Point> rrt_path_;
    nav_msgs::OccupancyGrid map_;
    nav_msgs::OccupancyGrid map_updated_;

    vector<geometry_msgs::Point> border_lines;

    geometry_msgs::TransformStamped tf_laser_to_map_;
    geometry_msgs::TransformStamped tf_baselink_to_map_;

    tf2_ros::Buffer tf_buffer_;
    tf2_ros::TransformListener tf2_listener_;

    double map_x_global;
    double map_y_global;

    std::vector<double> obstacle_x;
    std::vector<double> obstacle_y;

    std::vector<double> obstacle_x_inter;
    std::vector<double> obstacle_y_inter;
 
    void getParameters(ros::NodeHandle& nh);
    void init_occupancy_grid();
    void visualize_centerline();

    void obstacle();
    void compute_trajectory_table();
    void odom_callback(const nav_msgs::Odometry::ConstPtr &odom_msg);
    void select_trajectory();
    void simulate_dynamics(Vector3d& state, Vector2d& input, double dt, Vector3d& new_state);
    void execute_MPC();
    void get_linearized_dynamics(Matrix<double,nx,nx>& Ad, Matrix<double,nx, nu>& Bd, Matrix<double,nx,1>& hd, Matrix<double,nx,1>& x_op, Matrix<double,nu,1>& u_op);
    void applyControl(VectorXd& QPSolution);
    void visualize_trajectories(int low, int high);
    void visualize_mpc_solution(VectorXd& QPSolution);

    void rrt_path_callback(const visualization_msgs::Marker::ConstPtr &path_msg);
    void map_callback(const nav_msgs::OccupancyGrid::Ptr &map_msg);
    void scan_callback(const sensor_msgs::LaserScan::ConstPtr& scan_msg);




};