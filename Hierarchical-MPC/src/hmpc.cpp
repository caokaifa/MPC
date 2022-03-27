//
// Created by yuwei on 12/16/19.

#include <ros/ros.h>
#include <Eigen/Dense>
#include <hmpc_auto_race/hmpc.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <csignal>


//const string file_name = "/home/yuwei/rcws/logs/final1030.csv";
const string file_name = "/home/robert/catkin_ws/src/Hierarchical-MPC/src/yuwei_wp.csv";
const double RRT_INTERVAL = 0.1;

HMPC::HMPC(ros::NodeHandle &nh): nh_(nh), track_(Track(file_name, 0.1)),tf2_listener_(tf_buffer_){
   

    getParameters(nh_);
    init_occupancy_grid();

    odom_sub_ = nh_.subscribe(pose_topic, 10, &HMPC::odom_callback, this);
    drive_pub_ = nh_.advertise<ackermann_msgs::AckermannDriveStamped>(drive_topic, 1);

    scan_sub_ = nh_.subscribe(scan_topic, 1, &HMPC::scan_callback, this);
    // rrt_sub_ = nh_.subscribe("path_found", 1, &HMPC::rrt_path_callback, this);
    // map_sub_ = nh_.subscribe("map_updated", 1, &HMPC::map_callback, this);
    
    obstacle_pub_ =nh_.advertise<visualization_msgs::MarkerArray>("obstacle", 10);
    track_viz_pub_ = nh_.advertise<visualization_msgs::MarkerArray>("track_centerline", 10);
    trajectories_viz_pub_ = nh_.advertise<visualization_msgs::MarkerArray>("trajectories", 10);

    hmpc_viz_pub_ = nh_.advertise<visualization_msgs::MarkerArray>("HMPC", 10);
    debugger_pub_ = nh_.advertise<visualization_msgs::Marker>("Debugger", 10);

    compute_trajectory_table();
}

void HMPC::getParameters(ros::NodeHandle &nh) {
    nh.getParam("pose_topic", pose_topic);
    nh.getParam("drive_topic", drive_topic);
    nh.getParam("scan_topic", scan_topic);
    nh.getParam("speed_num", speed_num);
    nh.getParam("steer_num", steer_num);
    nh.getParam("N",N);
    nh.getParam("Ts",Ts);
    nh.getParam("ACCELERATION_MAX", ACCELERATION_MAX);
    nh.getParam("DECELERATION_MAX", DECELERATION_MAX);
    nh.getParam("SPEED_MAX", SPEED_MAX);
    nh.getParam("STEER_MAX", STEER_MAX);

    nh.getParam("q_x",q_x);
    nh.getParam("q_y",q_y);
    nh.getParam("q_yaw",q_yaw);
    nh.getParam("r_v",r_v);
    nh.getParam("r_steer",r_steer);
    nh.getParam("q_s",q_s);

    Q.setZero(); R.setZero();
    Q.diagonal() << q_x, q_y, q_yaw;
    R.diagonal() << r_v, r_steer;

    nh.getParam("MAP_MARGIN",MAP_MARGIN);
}

void HMPC::init_occupancy_grid(){
    boost::shared_ptr<nav_msgs::OccupancyGrid const> map_ptr;
    map_ptr = ros::topic::waitForMessage<nav_msgs::OccupancyGrid>("map", ros::Duration(5.0));
    if (map_ptr == nullptr){ROS_INFO("No map received");}
    else{
        map_ = *map_ptr;
        map_updated_ = map_;
        ROS_INFO("Map received");
    }
    ROS_INFO("Initializing occupancy grid for map ...");
    occupancy_grid::inflate_map(map_, MAP_MARGIN);
}

void HMPC::obstacle(){
    // plot waypoints
    visualization_msgs::Marker dots;
    dots.header.frame_id = "map";
    dots.id = rviz_id::CENTERLINE_POINTS;
    dots.ns = "obstacle";
    dots.type = visualization_msgs::Marker::POINTS;
    dots.scale.x = dots.scale.y = 0.08;
    dots.scale.z = 0.04;
    dots.action = visualization_msgs::Marker::ADD;
    dots.pose.orientation.w = 1.0;
    dots.color.g = 1.0;
    dots.color.a = 1.0;
    //dots.lifetime = ros::Duration();
    
    for (int i = 0; i < obstacle_x.size(); i++) {
        geometry_msgs::Point p;
        p.x = obstacle_x[i];
        p.y = obstacle_y[i];
        dots.points.push_back(p);
    }


    visualization_msgs::MarkerArray markers;
    markers.markers.push_back(dots);

    obstacle_pub_.publish(markers);
}


void HMPC::visualize_centerline(){
    // plot waypoints
    visualization_msgs::Marker dots;
    dots.header.frame_id = "map";
    dots.id = rviz_id::CENTERLINE_POINTS;
    dots.ns = "centerline_points";
    dots.type = visualization_msgs::Marker::POINTS;
    dots.scale.x = dots.scale.y = 0.08;
    dots.scale.z = 0.04;
    dots.action = visualization_msgs::Marker::ADD;
    dots.pose.orientation.w = 1.0;
    dots.color.r = 1.0;
    dots.color.a = 1.0;
    //dots.lifetime = ros::Duration();
    
    for (int i = 0; i < track_.centerline.size(); i++) {
        geometry_msgs::Point p;
        p.x = track_.centerline.at(i).x;
        p.y = track_.centerline.at(i).y;
        dots.points.push_back(p);
    }

    visualization_msgs::Marker line_dots;
    line_dots.header.stamp = ros::Time::now();
    line_dots.header.frame_id = "map";
    line_dots.id = rviz_id::CENTERLINE;
    line_dots.ns = "centerline";
    line_dots.type = visualization_msgs::Marker::LINE_STRIP;
    line_dots.scale.x = line_dots.scale.y = 0.02;
    line_dots.scale.z = 0.02;
    line_dots.action = visualization_msgs::Marker::ADD;
    line_dots.pose.orientation.w = 1.0;
    line_dots.color.b = 1.0;
    line_dots.color.a = 1.0;
    line_dots.points = dots.points;

    visualization_msgs::MarkerArray markers;
    markers.markers.push_back(dots);
    markers.markers.push_back(line_dots);

    track_viz_pub_.publish(markers);
}





void HMPC::odom_callback(const nav_msgs::Odometry::ConstPtr &odom_msg){
    /* process pose info */
    visualize_centerline();

    speed_m_ = odom_msg->twist.twist.linear.x;

    //speed_m_ = 3.2;

    tf::Quaternion q_tf;
    tf::quaternionMsgToTF(odom_msg->pose.pose.orientation, q_tf); //转化为tf形式

   // std::cout<<"----q_tf----"<<odom_msg->pose.pose.orientation<<" "<<q_tf<<std::endl;
    tf::Matrix3x3 rot(q_tf);
    double roll, pitch;
    rot.getRPY(roll, pitch, yaw_);  //转化为欧拉角

    
    car_pos_ = tf::Vector3(odom_msg->pose.pose.position.x, odom_msg->pose.pose.position.y, 0.0);

    // std::cout<<"------1.car_pos_"<<car_pos_[0]<<car_pos_[1]<<std::endl;
    // std:;cout<<"-----1.yaw_"<<yaw_<<std::endl;
    // tf_.setOrigin(car_pos_);
    // tf_.setRotation(q_tf);

    // std::cout<<"------2.car_pos_"<<car_pos_[0]<<car_pos_[1]<<std::endl;
    //  std::cout<<"-----2.yaw_"<<yaw_<<std::endl;
    /*select reference trajectory*/
    double now = ros::Time::now().toSec();
    select_trajectory();
    std::cout<<"-------time diff--------="<<ros::Time::now().toSec()-now<<std::endl;
    /*low level MPC to track the reference*/
    execute_MPC();
    std::cout<<"-------time diff-MPC-------="<<ros::Time::now().toSec()-now<<std::endl;
}


void HMPC::scan_callback(const sensor_msgs::LaserScan::ConstPtr &scan_msg)
{
    
   m.lock();
   std::cout<<"receiving scan date..................."<<std::endl;
   obstacle_x={};
   obstacle_y={};
    try
    {
        tf_laser_to_map_ = tf_buffer_.lookupTransform("map", "laser", ros::Time(0));
        tf_baselink_to_map_ = tf_buffer_.lookupTransform("map", "base_link", ros::Time(0));

    }
    catch (tf::TransformException& ex)
    {
        ROS_ERROR("%s",ex.what());
        ros::Duration(0.1).sleep();

        std::cout<<"----error------"<<std::endl;
    }
    const auto translation = tf_laser_to_map_.transform.translation;

    // std::cout<<"translation.x="<<translation.x<<" "<<"translation.y"<<translation.y<<std::endl;

    // std::cout<<"base_link_to_map"<<tf_baselink_to_map_.transform.translation.x<<std::endl;

    const double yaw = tf::getYaw(tf_laser_to_map_.transform.rotation);

    map_x_global=translation.x;
    map_y_global=translation.y;

    const auto start = static_cast<int>(2*scan_msg->ranges.size()/6);

    
    const auto end = static_cast<int>(4*scan_msg->ranges.size()/6);
    const double angle_increment = scan_msg->angle_increment;
    double theta = scan_msg->angle_min + angle_increment*(start-1);

    for(int i=start; i<end; ++i)
    {
        theta+=angle_increment;

        const double hit_range = scan_msg->ranges[i];
        if(std::isinf(hit_range) || std::isnan(hit_range)) continue;

     
        const double x_base_link = hit_range*cos(theta);
        const double y_base_link = hit_range*sin(theta);

        if(fabs(x_base_link) > 0.6|| fabs(y_base_link) > 0.6) continue;
        // if(fabs(x_base_link) < 0.3|| fabs(y_base_link) <0.3) continue;

        // laser hit x, y in base_link frame
        const double x_map = x_base_link*cos(yaw) - y_base_link*sin(yaw) + translation.x;
        const double y_map = x_base_link*sin(yaw) + y_base_link*cos(yaw) + translation.y;

        obstacle_x.push_back(x_map);
        obstacle_y.push_back(y_map);
       
        
        // std::vector<int> index_of_expanded_obstacles = get_expanded_row_major_indices(x_map, y_map);

        // for(const auto& index: index_of_expanded_obstacles)
        // {
        //     if(input_map_.data[index] != 100)
        //     {
        //         input_map_.data[index] = 100;
        //         new_obstacles_.emplace_back(index);
        //     }
        // }
    }

    m.unlock();

    obstacle();
    
    // clear_obstacles_count_++;
    // if(clear_obstacles_count_ > 10)
    // {
    //     for(const auto index: new_obstacles_)
    //     {
    //         input_map_.data[index] = 0;
    //     }
    //     new_obstacles_.clear();
    //     clear_obstacles_count_ = 0;
    //     ROS_INFO("Obstacles Cleared");
    // }

    // dynamic_map_pub_.publish(input_map_);
    // ROS_INFO("Map Updated");
   
}





void HMPC::map_callback(const nav_msgs::OccupancyGrid::Ptr &map_msg){
    visualization_msgs::Marker dots;
    dots.header.frame_id = "map";
    dots.id = rviz_id::DEBUG;
    dots.ns = "debug_points";
    dots.type = visualization_msgs::Marker::POINTS;
    dots.scale.x = dots.scale.y = 0.1;
    dots.scale.z = 0.1;
    dots.action = visualization_msgs::Marker::ADD;
    dots.pose.orientation.w = 1.0;
    dots.color.b = 1.0;
    dots.color.g = 0.5;
    dots.color.a = 1.0;

    vector<geometry_msgs::Point> path_processed;

    if(rrt_path_.empty()) return;
    for (int i=0; i< rrt_path_.size()-1; i++){
        path_processed.push_back(rrt_path_[i]);
        double dist = sqrt(pow(rrt_path_[i+1].x-rrt_path_[i].x, 2)
                           +pow(rrt_path_[i+1].y-rrt_path_[i].y, 2));
        if (dist < RRT_INTERVAL) continue;
        int num = static_cast<int>(ceil(dist/RRT_INTERVAL));
        for(int j=1; j< num; j++){
            geometry_msgs::Point p;
            p.x = rrt_path_[i].x + j*((rrt_path_[i+1].x - rrt_path_[i].x)/num);
            p.y = rrt_path_[i].y + j*((rrt_path_[i+1].y - rrt_path_[i].y)/num);
            path_processed.push_back(p);
        }
    }

    for (int i=0; i<path_processed.size(); i++){
        double theta = track_.findTheta(path_processed[i].x, path_processed[i].y, 0, true);
        Vector2d p_path(path_processed[i].x,path_processed[i].y);
        Vector2d p_proj(track_.x_eval(theta), track_.y_eval(theta));
        Vector2d p1, p2;
        int t=0;
        // search one direction until hit obstacle
        while(true){
            double x = (p_path + t*map_msg->info.resolution*(p_path-p_proj).normalized())(0);
            double y = (p_path + t*map_msg->info.resolution*(p_path-p_proj).normalized())(1);
            if(occupancy_grid::is_xy_occupied(*map_msg, x, y)){
                p1(0) = x; p1(1) = y;
                geometry_msgs::Point point;
                point.x = x; point.y = y;
                dots.points.push_back(point);
                break;
            }
            t++;
        }
        t=0;
        //search the other direction until hit obstacle
        while(true){
            double x = (p_path - t*map_msg->info.resolution*(p_path-p_proj).normalized())(0);
            double y = (p_path - t*map_msg->info.resolution*(p_path-p_proj).normalized())(1);
            if(occupancy_grid::is_xy_occupied(*map_msg, x, y)){
                p2(0) = x; p2(1) = y;
                geometry_msgs::Point point;
                point.x = x; point.y = y;
                dots.points.push_back(point);
                break;
            }
            t++;
        }
        double dx_dtheta = track_.x_eval_d(theta);
        double dy_dtheta = track_.y_eval_d(theta);
        double right_width = Vector2d(dy_dtheta, -dx_dtheta).dot(p1-p_proj)>0 ? (p1-p_proj).norm() : -(p1-p_proj).norm();
        double left_width = Vector2d(-dy_dtheta, dx_dtheta).dot(p2-p_proj)>0 ? (p2-p_proj).norm() : -(p2-p_proj).norm();
//        right_width = -0.15;
//        left_width =  0.4;
        cout<<"p1: "<< p1<<endl;
        cout<<"p2: "<< p2<<endl;

        track_.setHalfWidth(theta, left_width, right_width);
    }
    debugger_pub_.publish(dots);
}


void HMPC:: rrt_path_callback(const visualization_msgs::Marker::ConstPtr &path_msg){
    rrt_path_ = path_msg->points;
}

void HMPC::select_trajectory(){
    double low = max(0.3, speed_m_-DECELERATION_MAX*Ts);
    double high = min(SPEED_MAX, speed_m_ + 2.0*ACCELERATION_MAX*Ts);
    double dv = SPEED_MAX/speed_num;
    double ds = 2*STEER_MAX/steer_num;
    int low_ind = static_cast<int>(low/dv);
    int high_ind = static_cast<int>(high/dv);
    double max_theta = 0.0;
    car_theta_ = track_.findTheta(car_pos_.x(), car_pos_.y(), 0, true);
    double pos_in_map_x, pos_in_map_y;
    vector<Vector3d> traj;

    std::cout<<"------select_trajectory------"<<std::endl;
   obstacle_x_inter=obstacle_x;
   obstacle_y_inter=obstacle_y;
   double cost_obstacle_distance=0;
   visualize_trajectories(low_ind, high_ind);

    using namespace occupancy_grid;

    for (int i=0; i< steer_num+1; i++) {
        for (int j = low_ind; j <= high_ind; j++) {
            // first check if the whole trajectory lies inside track
            traj = trajectory_table_[i][j];
            bool break_flag = false;
            for (int k=0; k<traj.size(); k++){
                pos_in_map_x = cos(yaw_)*traj[k](0) - sin(yaw_)*traj[k](1) + car_pos_.x();
                pos_in_map_y = sin(yaw_)*traj[k](0) + cos(yaw_)*traj[k](1) + car_pos_.y();
                // update with poses transformed into map frame
                traj[k](0) = pos_in_map_x;   // x
                traj[k](1) = pos_in_map_y;   // y
                traj[k](2) += yaw_;          //yaw
                if (is_xy_occupied(map_, pos_in_map_x, pos_in_map_y)){
                    break_flag = true;
                    break;
                }
                
                // std::cout<<"obstacles is checking8888888="<<obstacle_x_inter.size()<<std::endl;
                // if(!obstacle_x_inter.empty())
                // {
                //     m.lock();
                //     for(int k=0;k<obstacle_x_inter.size();k++)
                //     {
                        
                //         double diff_x=obstacle_x_inter[k]-pos_in_map_x;
                //         double diff_y=obstacle_y_inter[k]-pos_in_map_y;
                //         double distance_obstacle=sqrt(pow(diff_x,2)+pow(diff_y,2));
                //         cost_obstacle_distance=cost_obstacle_distance+distance_obstacle;
                //         //std::cout<<"distance_obstacle--------="<<distance_obstacle<<std::endl;
                //         if(distance_obstacle<0.35)
                //         {
                //             cost_obstacle_distance=0;
                //             break_flag=true;
                //             break;
                //             std::cout<<"obstacles is runing on the road"<<std::endl;
                //         }

                //     }
                //     m.unlock();
                // }
               
            }

         
            if (break_flag) break;
            // compare to max progress so far
            double theta = track_.findTheta(traj.back()(0), traj.back()(1), 0, true);
            double theta_first=track_.findTheta(traj[1](0), traj[1](1), 0, true);

            //新增

            double dx_dtheta = track_.x_eval_d(theta);
            double dy_dtheta = track_.y_eval_d(theta);
            double d2x_dtheta2 = track_.x_eval_dd(theta);
            double d2y_dtheta2 = track_.y_eval_dd(theta);
            double numer = dx_dtheta*d2y_dtheta2 - dy_dtheta*d2x_dtheta2;
            double denom = pow(dx_dtheta,2) + pow(dy_dtheta, 2);

            double dphi_dtheta = numer/denom;
            double phi = atan2(dy_dtheta, dx_dtheta);//在区间[-pi, +pi]弧度


            //新增

            double dx_dtheta_first = track_.x_eval_d(theta_first);
            double dy_dtheta_first = track_.y_eval_d(theta_first);
            double d2x_dtheta2_first = track_.x_eval_dd(theta_first);
            double d2y_dtheta2_first = track_.y_eval_dd(theta_first);
            double numer_first = dx_dtheta_first*d2y_dtheta2_first - dy_dtheta_first*d2x_dtheta2_first;
            double denom_first = pow(dx_dtheta_first,2) + pow(dy_dtheta_first, 2);

            double dphi_dtheta_first = numer_first/denom_first;
            double phi_first = atan2(dy_dtheta_first, dx_dtheta_first);//在区间[-pi, +pi]弧度

        
            //std::cout<<"-----------phi-----------"<<phi<<"    "<<phi_first<<"    "<<phi-phi_first<<std::endl;

            
            //std::cout<<"----i="<<i<<" "<<"------j="<<j<<" "<<"----theta"<<theta<<std::endl;
            //account for discontinuity at joint of starting and goal point
            // if (theta-car_theta_ < -track_.length/2){ theta += track_.length; }
            // if (theta-car_theta_ > track_.length/2){ theta -= track_.length; }


            // if(phi-phi_first>0.5)
            // {
            //     theta=theta- 0.2*j*dv+0.1*i;

            //     std::cout<<"9999999999999999999-"<<std::endl;
            // }
            if (theta>max_theta){
                //std::cout<<"----finial i="<<i<<" "<<"------finial j="<<j<<" "<<-STEER_MAX + i*ds<<std::endl;
                input_ref_(0) = j*dv;  //speed_cmd
                input_ref_(1) = -STEER_MAX + i*ds;  // steer_cmd
                trajectory_ref_ = traj;
                max_theta = theta;
            }
        }

    }

}

void HMPC::simulate_dynamics(Vector3d& state, Eigen::Vector2d& input, double dt, Eigen::Vector3d& new_state){
    VectorXd dynamics(state.size());
    dynamics(0) = input(0)*cos(state(2));
    dynamics(1) = input(0)*sin(state(2));
    dynamics(2) = tan(input(1))*input(0)/CAR_LENGTH;
    new_state = state + dynamics * dt;
};

void HMPC::compute_trajectory_table(){
    trajectory_table_.clear();
    double dv = SPEED_MAX/speed_num;
    double ds = 2*STEER_MAX/steer_num;
    Vector3d state;
    Vector3d new_state;
    Vector2d input;

    vector<vector<Vector3d>> temp;
    vector<Vector3d> trajectory;

    for (int i=0; i<steer_num+1; i++){
        temp.clear();
        for(int j=0; j<speed_num+1; j++){
            trajectory.clear();

            double steer = -STEER_MAX + i*ds;
            double speed = j*dv;
            state.setZero();   // zero in car's body frame
            new_state.setZero();
            input(0) = speed; input(1) = steer;
            trajectory.push_back(state);   // add initial state
            for(int k=0; k<N+1; k++){
                simulate_dynamics(state, input, Ts, new_state);
                trajectory.push_back(new_state);
                state = new_state;
            }
            temp.push_back(trajectory);
        }
        trajectory_table_.push_back(temp);
    }
}

void HMPC::execute_MPC(){

    SparseMatrix<double> HessianMatrix((N+1)*(nx+nu) + (N+1), (N+1)*(nx+nu) + (N+1));
    SparseMatrix<double> constraintMatrix((N+1)*nx+ 2*(N+1) +(N+1)*nu + (N+1) + N, (N+1)*(nx+nu) + (N+1));

    VectorXd gradient((N+1)*(nx+nu) + (N+1));
    gradient.setZero();

    VectorXd lower((N+1)*nx+ 2*(N+1)+(N+1)*nu + (N+1) + N);
    VectorXd upper((N+1)*nx+ 2*(N+1)+(N+1)*nu + (N+1) + N);

    border_lines.clear();

    Matrix<double,nx,1> x_k_ref;
    Matrix<double,nx,1> hd;

    Matrix<double,nu,1> u_k_ref;
    u_k_ref = input_ref_;
    Matrix<double,nx,nx> Ad;
    Matrix<double,nx,nu> Bd;

    for (int i=0; i<N+1; i++){        //0 to N
        x_k_ref = trajectory_ref_[i];
        double theta_ref = track_.findTheta(x_k_ref(0), x_k_ref(1), 0, true);
        if (theta_ref-car_theta_ < -track_.length/2){ theta_ref += track_.length; }
        if (theta_ref-car_theta_ > track_.length/2){ theta_ref -= track_.length; }
        get_linearized_dynamics(Ad, Bd, hd, x_k_ref, u_k_ref);
        /* form Hessian entries*/
        // cost does not depend on x0, only 1 to N
        if (i>0) {
            // populate Qs
            for (int row = 0; row < nx; row++) {
                HessianMatrix.insert(i*nx + row, i*nx + row) = Q(row, row);
                
            }

            
            // populate Rs. every R is diagonal
            for(int row=0; row< nu; row++){
                HessianMatrix.insert((N+1)*nx+i*nu+row, (N+1)*nx+i*nu+row) = R(row,row);
            }
            HessianMatrix.insert((N+1)*(nx+nu)+i, (N+1)*(nx+nu)+i) = q_s;
            /* form gradient vector*/
            gradient.segment<nx>(i*nx) << -Q*x_k_ref;
            gradient.segment<nu>((N+1)*nx+ i*nu) << -R*u_k_ref;
            // std::cout<<"HessianMatrix="<<HessianMatrix<<std::endl;

            // std::cout<<"gradient="<<gradient<<std::endl;
        }
        /* form constraint matrix */
        if (i<N){
            // Ad
            for (int row=0; row<nx; row++){
                for(int col=0; col<nx; col++){
                    constraintMatrix.insert((i+1)*nx+row, i*nx+col) = Ad(row,col);
                }
            }
            // Bd
            for (int row=0; row<nx; row++){
                for(int col=0; col<nu; col++){
                    //constraintMatrixTripletList.push_back(T(i*nx+row, (N+1)*nx+ i*nu+col, s.Bd(row,col)));
                    constraintMatrix.insert((i+1)*nx+row, (N+1)*nx+ i*nu+col) = Bd(row,col);
                }
            }
            lower.segment<nx>((i+1)*nx) = -hd;
            upper.segment<nx>((i+1)*nx) = -hd;
        }
        /* track boundary constraints */
        double dx_dtheta = track_.x_eval_d(theta_ref);
        double dy_dtheta = track_.y_eval_d(theta_ref);

        // std::cout<<"-----dx_dtheta----"<<dx_dtheta<<"----dy_dtheta"<<dy_dtheta<<std::endl;

        constraintMatrix.insert((N+1)*nx+ 2*i, i*nx) = -dy_dtheta;      // a*x
        constraintMatrix.insert((N+1)*nx+ 2*i, i*nx+1) = dx_dtheta;     // b*y
        constraintMatrix.insert((N+1)*nx+ 2*i, (N+1)*(nx+nu) +i) = 1.0;   // min(C1,C2) <= a*x + b*y + s_k <= inf

        constraintMatrix.insert((N+1)*nx+ 2*i+1, i*nx) = -dy_dtheta;      // a*x
        constraintMatrix.insert((N+1)*nx+ 2*i+1, i*nx+1) = dx_dtheta;     // b*y
        constraintMatrix.insert((N+1)*nx+ 2*i+1, (N+1)*(nx+nu) +i) = -1.0;   // -inf <= a*x + b*y - s_k <= max(C1,C2)

        //get upper line and lower line
        Vector2d left_tangent_p, right_tangent_p, center_p;
        Vector2d right_line_p1, right_line_p2, left_line_p1, left_line_p2;
        geometry_msgs::Point r_p1, r_p2, l_p1, l_p2;

        center_p << track_.x_eval(theta_ref), track_.y_eval(theta_ref);
        right_tangent_p = center_p + track_.getRightHalfWidth(theta_ref) * Vector2d(dy_dtheta, -dx_dtheta).normalized();
        left_tangent_p  = center_p + track_.getLeftHalfWidth(theta_ref) * Vector2d(-dy_dtheta, dx_dtheta).normalized();

        right_line_p1 = right_tangent_p + 0.15*Vector2d(dx_dtheta, dy_dtheta).normalized();
        right_line_p2 = right_tangent_p - 0.15*Vector2d(dx_dtheta, dy_dtheta).normalized();
        left_line_p1 = left_tangent_p + 0.15*Vector2d(dx_dtheta, dy_dtheta).normalized();
        left_line_p2 = left_tangent_p - 0.15*Vector2d(dx_dtheta, dy_dtheta).normalized();

        r_p1.x = right_line_p1(0);  r_p1.y = right_line_p1(1);
        r_p2.x = right_line_p2(0);  r_p2.y = right_line_p2(1);
        l_p1.x = left_line_p1(0);   l_p1.y = left_line_p1(1);
        l_p2.x = left_line_p2(0);   l_p2.y = left_line_p2(1);

        border_lines.push_back(r_p1);  border_lines.push_back(r_p2);
        border_lines.push_back(l_p1); border_lines.push_back(l_p2);

        double C1 =  - dy_dtheta*right_tangent_p(0) + dx_dtheta*right_tangent_p(1);
        double C2 = - dy_dtheta*left_tangent_p(0) + dx_dtheta*left_tangent_p(1);

        //lower((N+1)*nx+i) = -OsqpEigen::INFTY;//min(C1, C2);
        //upper((N+1)*nx+i) = OsqpEigen::INFTY;//max(C1, C2);
        lower((N+1)*nx+ 2*i) = min(C1, C2);
        upper((N+1)*nx+ 2*i) = OsqpEigen::INFTY;

        lower((N+1)*nx+ 2*i+1) = -OsqpEigen::INFTY;
        upper((N+1)*nx+ 2*i+1) = max(C1, C2);
        // -I for each x_k+1
        for (int row=0; row<nx; row++) {
            constraintMatrix.insert(i*nx+row, i*nx+row) = -1.0;
        }
        // u_min < u < u_max
        for (int row=0; row<nu; row++){
            constraintMatrix.insert((N+1)*nx+ 2*(N+1) +i*nu+row, (N+1)*nx+i*nu+row) = 1.0;
        }
        // s_k >= 0
        constraintMatrix.insert((N+1)*nx+ 2*(N+1) +(N+1)*nu+ i, (N+1)*(nx+nu)+ i) = 1.0;
        // acceleration constraint
        if(i<N){
            constraintMatrix.insert((N+1)*nx+ 2*(N+1) +(N+1)*nu + (N+1) +i, (N+1)*nx+i*nu) = -1;
            constraintMatrix.insert((N+1)*nx+ 2*(N+1) +(N+1)*nu + (N+1) +i, (N+1)*nx+(i+1)*nu) = 1;
            lower((N+1)*nx+ 2*(N+1) +(N+1)*nu + (N+1) +i) = -OsqpEigen::INFTY;
            upper((N+1)*nx+ 2*(N+1) +(N+1)*nu + (N+1) +i) = ACCELERATION_MAX *Ts;
        }
        // input bounds: speed and steer
        lower((N+1)*nx+ 2*(N+1) +i*nu) = 0.0;
        upper((N+1)*nx+ 2*(N+1)+i*nu) = SPEED_MAX;
        lower((N+1)*nx+ 2*(N+1)+i*nu+1) = -STEER_MAX;
        upper((N+1)*nx+ 2*(N+1)+i*nu+1) = STEER_MAX;
        // s_k >= 0
        lower((N+1)*nx+ 2*(N+1) +(N+1)*nu+ i) = 0.0;
        upper((N+1)*nx+ 2*(N+1) +(N+1)*nu+ i) = OsqpEigen::INFTY;
    }
    // enforcing inital conditions

    lower.head(nx) = -trajectory_ref_[0];  //x0
    upper.head(nx) = -trajectory_ref_[0];
    lower((N+1)*nx+ 2*(N+1)) = max(0.0, last_speed_cmd_- DECELERATION_MAX *Ts);   //u0
    upper((N+1)*nx+ 2*(N+1)) = min(SPEED_MAX, last_speed_cmd_ + ACCELERATION_MAX/2 *Ts);

    SparseMatrix<double> H_t = HessianMatrix.transpose();
    SparseMatrix<double> sparse_I((N+1)*(nx+nu) + (N+1),(N+1)*(nx+nu) + (N+1));
    sparse_I.setIdentity();
    HessianMatrix = 0.5*(HessianMatrix + H_t) + 0.000001*sparse_I;

    OsqpEigen::Solver solver;
    solver.settings()->setWarmStart(true);
    solver.data()->setNumberOfVariables((N+1)*(nx+nu) + (N+1));
    solver.data()->setNumberOfConstraints((N+1)*nx+ 2*(N+1) +(N+1)*nu + (N+1) + N);

    if (!solver.data()->setHessianMatrix(HessianMatrix)) throw "fail set Hessian";
    if (!solver.data()->setGradient(gradient)){throw "fail to set gradient";}
    if (!solver.data()->setLinearConstraintsMatrix(constraintMatrix)) throw"fail to set constraint matrix";
    if (!solver.data()->setLowerBound(lower)){throw "fail to set lower bound";}
    if (!solver.data()->setUpperBound(upper)){throw "fail to set upper bound";}

    if(!solver.initSolver()){ cout<< "fail to initialize solver"<<endl;}

    if(!solver.solve()) {
        return;
    }
    VectorXd QPSolution = solver.getSolution();

    // cout<<"Solution: "<<endl;
    // cout<<QPSolution<<endl;

    applyControl(QPSolution);
    visualize_mpc_solution(QPSolution);

    solver.clearSolver();

}

void HMPC::get_linearized_dynamics(Matrix<double,nx,nx>& Ad, Matrix<double,nx, nu>& Bd, Matrix<double,nx,1>& hd, Matrix<double,nx,1>& x_op, Matrix<double,nu,1>& u_op){
    double yaw = x_op(2);
    double v = u_op(0);
    double steer = u_op(1);

    Vector3d dynamics, h;
    dynamics(0) = u_op(0)*cos(x_op(2));
    dynamics(1) = u_op(0)*sin(x_op(2));
    dynamics(2) = tan(u_op(1))*u_op(0)/CAR_LENGTH;


    Matrix<double,nx,nx> A, M12;
    Matrix<double,nx,nu> B;

    A <<   0.0, 0.0, -v*sin(yaw),
            0.0, 0.0,  v*cos(yaw),
            0.0, 0.0,      0.0;

    B <<   cos(yaw), 0.0,
            sin(yaw), 0.0,
            tan(steer)/CAR_LENGTH, v/(cos(steer)*cos(steer)*CAR_LENGTH);

    Matrix<double,nx+nx,nx+nx> aux, M;
    aux.setZero();
    aux.block<nx,nx>(0,0) << A;
    aux.block<nx,nx>(0, nx) << Matrix3d::Identity();
    M = (aux*Ts).exp();
    M12 = M.block<nx,nx>(0,nx);
    h = dynamics - (A*x_op + B*u_op);

    //Discretize
    Ad = (A*Ts).exp();
    Bd = M12*B;
    hd = M12*h;
}

void HMPC::visualize_trajectories(int low, int high){
    double dv = SPEED_MAX/speed_num;
    //int v_ind = min(speed_num, max(0, static_cast<int>(speed_cmd/dv)));
    visualization_msgs::MarkerArray traj_list;
    visualization_msgs::Marker traj;
    geometry_msgs::Point p;

    traj.header.frame_id = "base_link";
    traj.id = rviz_id::TRAJECTORIES;
    traj.type = visualization_msgs::Marker::LINE_STRIP;
    traj.scale.x = traj.scale.y = 0.01;
    traj.scale.z = 0.02;
    traj.action = visualization_msgs::Marker::ADD;
    traj.pose.orientation.w = 1.0;
    traj.color.r = 1.0;
    traj.color.a = 1.0;

    for (int i=0; i<steer_num+1; i++) {
        for (int v_ind = low; v_ind <= high; v_ind++) {
            traj.points.clear();
            traj.id += 9;
            traj.color.b += 0.1;

            for (int j = 0; j < trajectory_table_[i][v_ind].size(); j++) {
                p.x = trajectory_table_[i][v_ind][j](0);
                p.y = trajectory_table_[i][v_ind][j](1);
                traj.points.push_back(p);
            }
            traj_list.markers.push_back(traj);
        }
    }
    trajectories_viz_pub_.publish(traj_list);
}

void HMPC::visualize_mpc_solution(VectorXd& QPSolution){
    visualization_msgs::Marker traj_ref;
    geometry_msgs::Point p;

    traj_ref.header.stamp = ros::Time::now();
    traj_ref.header.frame_id = "map";
    traj_ref.id = rviz_id::TRAJECTORY_REF;
    traj_ref.ns = "reference_trajectory";
    traj_ref.type = visualization_msgs::Marker::LINE_STRIP;
    traj_ref.scale.x = traj_ref.scale.y = 0.1;
    traj_ref.scale.z = 0.1;
    traj_ref.action = visualization_msgs::Marker::ADD;
    traj_ref.pose.orientation.w = 1.0;
    traj_ref.color.g = 1.0;
    traj_ref.color.a = 1.0;

    for (int i=0; i<trajectory_ref_.size(); i++){
        p.x = trajectory_ref_[i](0);
        p.y =  trajectory_ref_[i](1);
        traj_ref.points.push_back(p);
    }
    visualization_msgs::Marker  max_theta;
    max_theta.header.frame_id = "map";
    max_theta.id = rviz_id::MAX_THETA;
    max_theta.ns = "max_theta";
    max_theta.type = visualization_msgs::Marker::SPHERE;
    max_theta.scale.x = max_theta.scale.y = max_theta.scale.z = 0.5;
    max_theta.action = visualization_msgs::Marker::ADD;
    max_theta.pose.orientation.w = 1.0;
    // max_theta.color.g = 1.0;
    max_theta.color.b = 1.0;
    max_theta.color.a = 1.0;
    double theta = track_.findTheta(p.x, p.y, 0, true);
    max_theta.pose.position.x = track_.x_eval(theta);  max_theta.pose.position.y = track_.y_eval(theta);

    visualization_msgs::Marker pred_dots;
    pred_dots.header.frame_id = "map";
    pred_dots.id = rviz_id::PREDICTION;
    pred_dots.ns = "predicted_positions";
    pred_dots.type = visualization_msgs::Marker::POINTS;
    pred_dots.scale.x = pred_dots.scale.y = pred_dots.scale.z = 0.04;
    pred_dots.action = visualization_msgs::Marker::ADD;
    pred_dots.pose.orientation.w = 1.0;
    // pred_dots.color.g = 1.0;
    pred_dots.color.r = 1.0;
    pred_dots.color.a = 1.0;
    for (int i=0; i<N+1; i++){
        geometry_msgs::Point p;
        p.x = QPSolution(i*nx);
        p.y = QPSolution(i*nx+1);
        pred_dots.points.push_back(p);
    }

    visualization_msgs::Marker borderlines;
    borderlines.header.frame_id = "map";
    borderlines.id = rviz_id::BORDERLINES;
    borderlines.ns = "borderlines";
    borderlines.type = visualization_msgs::Marker::LINE_LIST;
    borderlines.scale.x = 0.1;
    borderlines.action = visualization_msgs::Marker::ADD;
    borderlines.pose.orientation.w = 1.0;
    borderlines.color.b = 1.0;
    borderlines.color.a = 1.0;

    borderlines.points = border_lines;

    visualization_msgs::MarkerArray mpc_markers;
    mpc_markers.markers.push_back(traj_ref);
    mpc_markers.markers.push_back(pred_dots);
    mpc_markers.markers.push_back(borderlines);
    mpc_markers.markers.push_back(max_theta);

    hmpc_viz_pub_.publish(mpc_markers);
}

void HMPC::applyControl(VectorXd& QPSolution) {

    float speed = QPSolution((N+1)*nx);
    float steer = QPSolution((N+1)*nx+1);
    cout<<"speed_cmd: "<<speed<<endl;
    cout<<"steer_cmd: "<<steer<<endl;
    cout<<"input_ref: "<<input_ref_<<endl;

    steer = min(steer, 0.41f);
    steer = max(steer, -0.41f);

    ackermann_msgs::AckermannDriveStamped ack_msg;
    ack_msg.drive.speed = speed;
    ack_msg.drive.steering_angle = steer;
    ack_msg.drive.steering_angle_velocity = 1.0;
    drive_pub_.publish(ack_msg);
    last_speed_cmd_ = speed;
}

void SignalHandler(int signal){
    if(ros::isInitialized() && ros::isStarted() && ros::ok() && !ros::isShuttingDown()){
        ros::shutdown();
    }
}


int main(int argc, char **argv){
    ros::init(argc, argv, "hmpc_node");
    signal(SIGINT, SignalHandler);
    signal(SIGTERM,SignalHandler);
    ros::NodeHandle nh;
    HMPC HMPC(nh);
    ros::Rate rate(25);
    while(ros::ok()){
        ros::spinOnce();
        rate.sleep();
    }
    return 0;
}

