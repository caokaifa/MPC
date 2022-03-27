//
// Created by yuwei on 11/19/19.
#include <ros/ros.h>
#include <Eigen/Dense>
#include <mpcc_auto_race/mpcc.h>
#include <unsupported/Eigen/MatrixFunctions>

const string file_name = "/home/robert/catkin_ws/src/Model-Predictive-Contouring-Control/src/yuwei_wp.csv";

MPCC::MPCC(ros::NodeHandle &nh): nh_(nh), track_(Track(file_name, 0.8)){

    getParams(nh_);

    odom_sub_ = nh_.subscribe(pose_topic, 10, &MPCC::odom_callback, this);
    drive_pub_ = nh_.advertise<ackermann_msgs::AckermannDriveStamped>(drive_topic, 1);

    track_viz_pub_ = nh_.advertise<visualization_msgs::MarkerArray>("track_centerline", 10);
    mpcc_viz_pub_ = nh_.advertise<visualization_msgs::MarkerArray>("mpcc", 10);

    visualize_centerline();
    QPSolution_.setZero();

    reset_flag = true;
    cout<<"caokaifa"<<endl;
}

void MPCC::getParams(ros::NodeHandle &nh) {
    nh.getParam("q_c",q_c);
    nh.getParam("q_l",q_l);
    nh.getParam("q_v",q_v);
    nh.getParam("r_v",r_v);
    nh.getParam("r_steer",r_steer);
    nh.getParam("N",N);
    nh.getParam("Ts",Ts);
    nh.getParam("SPEED_MAX",SPEED_MAX);
    nh.getParam("STEER_MAX",STEER_MAX);
    nh.getParam("pose_topic", pose_topic);
    nh.getParam("drive_topic", drive_topic);
}

void MPCC::visualize_centerline() {
    // plot waypoints
    visualization_msgs::Marker dots;
    dots.header.stamp = ros::Time::now();
    dots.header.frame_id = "map";
    dots.id = rviz_id::CENTERLINE;
    dots.ns = "centerline";
    dots.type = visualization_msgs::Marker::POINTS;
    dots.scale.x = dots.scale.y = 0.08;
    dots.scale.z = 0.04;
    dots.action = visualization_msgs::Marker::ADD;
    dots.pose.orientation.w = 1.0;
    dots.color.r = 1.0;
    dots.color.a = 1.0;
    //dots.lifetime = ros::Duration();

    vector<Point_ref> centerline = track_.centerline;

    for (int i=0; i<centerline.size(); i++){
            geometry_msgs::Point p;
            p.x = centerline.at(i).x;
            p.y = centerline.at(i).y;
            dots.points.push_back(p);
    }

    visualization_msgs::Marker spline_dots;
    spline_dots.header.stamp = ros::Time::now();
    spline_dots.header.frame_id = "map";
    spline_dots.id = rviz_id::CENTERLINE_SPLINE;
    spline_dots.ns = "centerline";
    spline_dots.type = visualization_msgs::Marker::LINE_STRIP;
    spline_dots.scale.x = spline_dots.scale.y = 0.02;
    spline_dots.scale.z = 0.02;
    spline_dots.action = visualization_msgs::Marker::ADD;
    spline_dots.pose.orientation.w = 1.0;
    spline_dots.color.b = 1.0;
    spline_dots.color.a = 1.0;
   // spline_dots.lifetime = ros::Duration();

    for (float t=0.0; t<track_.length; t+=0.05){
        geometry_msgs::Point p;
        p.x = track_.x_eval(t);
        p.y = track_.y_eval(t);
        spline_dots.points.push_back(p);
    }


    visualization_msgs::MarkerArray markers;
    markers.markers.push_back(dots);
    markers.markers.push_back(spline_dots);
    track_viz_pub_.publish(markers);

}

void MPCC::odom_callback(const nav_msgs::Odometry::ConstPtr &odom_msg){
    visualize_centerline();

    car_pose_ = odom_msg->pose.pose;
    tf::Quaternion quat;
    tf::quaternionMsgToTF(car_pose_.orientation, quat);
    tf::Matrix3x3 rot(quat);
    double roll, pitch, yaw;
    rot.getRPY(roll, pitch, yaw);
    // current states as x0: x, y, yaw, theta
    VectorXd x0(nx);
    x0(0) = car_pose_.position.x+0.1;
    x0(1) = car_pose_.position.y+0.1;
    x0(2) = yaw;
    x0(3) = track_.findTheta(x0(0), x0(1), 0, true);//获取最新的车的位置

    cout<<"--------x0: "<<x0<<endl;
    speed_m_ = odom_msg->twist.twist.linear.x;

    VectorXd QPSolution_shifted((N+1)*(nx+nu));
    VectorXd QP_x0((N+1)*(nx+nu));

    for (int i=0; i<N+1; i++){
        QP_x0.segment<nx>(i*nx) = x0;
        QP_x0.segment<nu>((N+1)*nx+i*nu) << 1.0, 0.0; //QP的初始化值
    }


    if(reset_flag||QPSolution_.isZero()){
        /*first run*/
        // assume a 1 m/s initial speed along the centerline, direction aligned with centerline
        // yaw equal to phi at all time
        // use x0, u=0 as linearization point for the entire horizon
        initialize_solution(x0);
        QPSolution_shifted = QP_x0;

        cout<<"QP_shifted: "<<QPSolution_shifted<<endl;
        solveQP(QPSolution_shifted, QP_x0, true);
        theta_prev_ = x0(3);
        reset_flag = false;
    }
    else{

        QPSolution_shifted.head(N*nx) = QPSolution_.segment(nx, N*nx);
        QPSolution_shifted.head(nx) = x0;
        QPSolution_shifted.segment((N+1)*nx, N*nu) = QPSolution_.segment((N+1)*nx+nu, N*nu);
        VectorXd last_state(nx); VectorXd last_input(nu);
        last_state = QPSolution_.segment<nx>(N*nx);
        last_input=  QPSolution_.tail(nu);

        QPSolution_shifted.segment<nx>(N*nx) = simulate_dynamics(last_state, last_input, Ts);
        QPSolution_shifted.tail(nu) = last_input;
        if(abs(x0(3) - theta_prev_) > track_.length/2){
            reset_flag = true;
        }

        solveQP(QPSolution_shifted, QP_x0, true);
        theta_prev_ = x0(3);
    }

}

void MPCC::initialize_solution(VectorXd& x0){

    VectorXd solution_init = VectorXd::Zero((N+1)*(nx+nu));
    solution_init.head(nx) << x0;
    double vx0 = 0.2;
    double theta = x0(3);
    double phi_prev = x0(2);
    for (int i=1; i<N+1; i++){
        theta += vx0 * Ts;
        // to prevent discontinuity in phi
        double phi = track_.getPhi(theta);
        if (phi - phi_prev < -M_PI){
            phi += 2*M_PI;
        }
        if (phi - phi_prev > M_PI){
            phi -= 2*M_PI;
        }
        double x = track_.x_eval(theta);
        double y = track_.y_eval(theta);
        solution_init.segment<nx>(i*nx) << x, y, phi, theta;
        solution_init.segment<nu>((N+1)*nx+i*nu)<< vx0, 0.0;
    }
    QPSolution_ = solution_init;
}

VectorXd MPCC::simulate_dynamics(VectorXd& state, VectorXd& input, double dt){
    VectorXd dynamics(state.size());
    double phi = track_.getPhi(state(3));

    dynamics(0) = input(0)*cos(state(2));
    dynamics(1) = input(0)*sin(state(2));
    dynamics(2) = tan(input(1))*input(0)/CAR_LENGTH;
    dynamics(3) = input(0)*cos(state(2)-phi);

    VectorXd new_state(state.size());
    new_state = state + dynamics * dt;

    return new_state;
};

typedef Eigen::Triplet<double> T;
void MPCC::solveQP(VectorXd& QPSolution_shifted, VectorXd& QP_x0, bool reset){
    // IMPORTANT: z = [x0, ..., x_N, u0, ..., uN]   length: (horizon+1)*(nx+nu)
    vector<T> HessianTripletList;
    vector<T> constraintMatrixTripletList;
    SparseMatrix<double> HessianMatrix((N+1)*(nx+nu),(N+1)*(nx+nu));
    SparseMatrix<double> constraintMatrix((N+1)*nx+(N+1)+(N+1)*nu + 2*N, (N+1)*(nx+nu));

    VectorXd gradient((N+1)*(nx+nu));

    VectorXd lower((N+1)*nx+(N+1)+(N+1)*nu + 2*N);
    VectorXd upper((N+1)*nx+(N+1)+(N+1)*nu + 2*N);

    lower.head(nx) = -QP_x0.head(nx);
    upper.head(nx) = -QP_x0.head(nx);

    border_lines.clear();
   // lower.setZero();
    //upper.setZero();
//
//    lower.head(nx) = Vector4d::Zero();
//    upper.head(nx) = Vector4d::Zero();
//    simple_lower.head(nx) = Vector4d::Zero();
//    simple_upper.head(nx) = Vector4d::Zero();

    for (int i=0; i<N+1; i++){        //0 to N
        Stage s;
        Matrix<double,nx,1> x_prev;
        Matrix<double,nu,1> u_prev;

        if(reset){
            x_prev = QP_x0.segment<nx>(i*nx);
            u_prev = QP_x0.segment<nu>((N+1)*nx + i*nu);
        }
        else {
             x_prev = QPSolution_shifted.segment<nx>(i*nx);
             u_prev = QPSolution_shifted.segment<nu>((N +1)*nx + i*nu);
        }
//        cout<< "x_prev: "<<x_prev<<endl;
//        cout<<"u_prev: "<<u_prev<<endl;
        getMatrices(s.Q, s.fx, s.dv_cost_dx, s.fu, s.Ad, s.Bd, x_prev, u_prev);
        s.R << r_v, 0.0, 0.0, r_steer;
        /* form Hessian entries*/
        // cost does not depend on x0, only 1 to N

        if (i>0) {
            // populate Qs
            for (int row = 0; row < nx; row++) {
                for (int col = 0; col < nx; col++) {
                    HessianMatrix.insert(i*nx + row, i*nx + col) = s.Q(row, col);
                }
            }
            // populate Rs. every R is diagonal
            for(int row=0; row< nu; row++){
                HessianMatrix.insert((N+1)*nx+i*nu+row, (N+1)*nx+i*nu+row) = s.R(row,row);
            }

            /* form gradient vector*/
            gradient.segment<nx>(i*nx) << s.fx;
            gradient.segment<nu>((N+1)*nx+i*nu) << s.fu;
        }
        /* form constraint matrix */
        if (i<N){
            Matrix<double,nx,1> x_kp1_prev;
            Matrix<double,nu,1> u_kp1_prev;
            if(reset){
                x_kp1_prev = QP_x0.segment<nx>((i+1)*nx);
                u_kp1_prev = QP_x0.segment<nu>((N+1)*nx + (i+1)*nu);
            }
            else {
                x_kp1_prev = QPSolution_shifted.segment<nx>((i+1)*nx);
                u_kp1_prev = QPSolution_shifted.segment<nu>((N+1)*nx + (i+1)*nu);
            }
            // Ad
            for (int row=0; row<nx; row++){
                for(int col=0; col<nx; col++){
                    //constraintMatrixTripletList.push_back(T(i*nx+row, (i-1)*nx+col, s.Ad(row,col)));
                    constraintMatrix.insert((i+1)*nx+row, i*nx+col) = s.Ad(row,col);
                }
            }
            // Bd
            for (int row=0; row<nx; row++){
                for(int col=0; col<nu; col++){
                    //constraintMatrixTripletList.push_back(T(i*nx+row, (N+1)*nx+ i*nu+col, s.Bd(row,col)));
                    constraintMatrix.insert((i+1)*nx+row, (N+1)*nx+ i*nu+col) = s.Bd(row,col);
                }
            }
            lower.segment<nx>((i+1)*nx) = -x_kp1_prev + s.Ad*x_prev + s.Bd*u_prev;
            upper.segment<nx>((i+1)*nx) = -x_kp1_prev + s.Ad*x_prev + s.Bd*u_prev;

        }

        /* track boundary constraints */

        Matrix<double,nx,1> x_prev_b;
        x_prev_b = QPSolution_shifted.segment<nx>(i*nx);
        double dx_dtheta = track_.x_eval_d(x_prev_b(3));
        double dy_dtheta = track_.y_eval_d(x_prev_b(3));

        constraintMatrix.insert((N+1)*nx+i, i*nx) = -dy_dtheta;
        constraintMatrix.insert((N+1)*nx+i, i*nx+1) = dx_dtheta;

         //get upper line and lower line
         Vector2d left_tangent_p, right_tangent_p, center_p;
         Vector2d right_line_p1, right_line_p2, left_line_p1, left_line_p2;
         geometry_msgs::Point r_p1, r_p2, l_p1, l_p2;

         center_p << track_.x_eval(x_prev_b(3)), track_.y_eval(x_prev_b(3));
         right_tangent_p = center_p + 0.5* Vector2d(dy_dtheta, -dx_dtheta).normalized();
         left_tangent_p  = center_p + track_.getHalfWidth(x_prev_b(3)) * Vector2d(-dy_dtheta, dx_dtheta).normalized();

         right_line_p1 = right_tangent_p + 0.8*Vector2d(dx_dtheta, dy_dtheta).normalized();
         right_line_p2 = right_tangent_p - 0.8*Vector2d(dx_dtheta, dy_dtheta).normalized();
         left_line_p1 = left_tangent_p + 0.8*Vector2d(dx_dtheta, dy_dtheta).normalized();
         left_line_p2 = left_tangent_p - 0.8*Vector2d(dx_dtheta, dy_dtheta).normalized();

         r_p1.x = right_line_p1(0);  r_p1.y = right_line_p1(1);
         r_p2.x = right_line_p2(0);  r_p2.y = right_line_p2(1);
         l_p1.x = left_line_p1(0);   l_p1.y = left_line_p1(1);
         l_p2.x = left_line_p2(0);   l_p2.y = left_line_p2(1);

        border_lines.push_back(r_p1);  border_lines.push_back(r_p2);
        border_lines.push_back(l_p1); border_lines.push_back(l_p2);


        double C1 =  - dy_dtheta*right_tangent_p(0) + dx_dtheta*right_tangent_p(1);
        double C2 = - dy_dtheta*left_tangent_p(0) + dx_dtheta*left_tangent_p(1);

        lower((N+1)*nx+i) = min(C1, C2);
        upper((N+1)*nx+i) = max(C1, C2);

        // -I for each x_k+1
        for (int row=0; row<nx; row++) {
            constraintMatrix.insert(i*nx+row, i*nx+row) = -1.0;
        }
        // u_min < u < u_max
        for (int row=0; row<nu; row++){
            constraintMatrix.insert((N+1)*nx+(N+1)+i*nu+row, (N+1)*nx+i*nu+row) = 1.0;
        }
        // input bounds: speed and steer
        lower((N+1)*nx+(N+1)+i*nu) = 0.0;
        upper((N+1)*nx+(N+1)+i*nu) = SPEED_MAX;
        lower((N+1)*nx+(N+1)+i*nu+1) = -STEER_MAX;
        upper((N+1)*nx+(N+1)+i*nu+1) = STEER_MAX;

        if(i<N){
            // theta must not decrease
            constraintMatrix.insert((N+1)*nx+(N+1)+(N+1)*nu+i, i*nx+3) = -1;
            constraintMatrix.insert((N+1)*nx+(N+1)+(N+1)*nu+i, (i+1)*nx+3) = 1;
            lower((N+1)*nx+(N+1)+(N+1)*nu+i) = 0.0;
            upper((N+1)*nx+(N+1)+(N+1)*nu+i) = OsqpEigen::INFTY;

            constraintMatrix.insert((N+1)*nx+(N+1)+(N+1)*nu+N+i, (N+1)*nx+i*nu) = -1;
            constraintMatrix.insert((N+1)*nx+(N+1)+(N+1)*nu+N+i, (N+1)*nx+(i+1)*nu) = 1;
            lower((N+1)*nx+(N+1)+(N+1)*nu+N+i) = -OsqpEigen::INFTY;
            upper((N+1)*nx+(N+1)+(N+1)*nu+N+i) =  0.2*Ts;
        }

    }
    gradient(N*nx+3) = -q_v;

    SparseMatrix<double> H_t = HessianMatrix.transpose();
    SparseMatrix<double> sparse_I((N+1)*(nx+nu),(N+1)*(nx+nu));
    sparse_I.setIdentity();
    HessianMatrix = 0.5*(HessianMatrix + H_t) + 0.000001*sparse_I;

//    cout<< "Hessian: "<< endl;
//     cout<<HessianMatrix << endl;
//    cout<< "constraint: "<<endl;
//    cout<< constraintMatrix << endl;
//    cout<<"lower: "<<endl;
//     cout<<lower<<endl;
//    cout<<"upper: "<<endl;
//    cout<<upper<<endl;
//     cout<<"gradient: "<<endl;
//     cout<<gradient<<endl;
 
    OsqpEigen::Solver solver;
    solver.settings()->setWarmStart(true);
    solver.data()->setNumberOfVariables((N+1)*(nx+nu));
    solver.data()->setNumberOfConstraints((N+1)*nx+(N+1)+(N+1)*nu + 2*N);

    if (!solver.data()->setHessianMatrix(HessianMatrix)) throw "fail set Hessian";
    if (!solver.data()->setGradient(gradient)){throw "fail to set gradient";}
    if (!solver.data()->setLinearConstraintsMatrix(constraintMatrix)) throw"fail to set constraint matrix";
    if (!solver.data()->setLowerBound(lower)){throw "fail to set lower bound";}
    if (!solver.data()->setUpperBound(upper)){throw "fail to set upper bound";}

    if(!solver.initSolver()){ cout<< "fail to initialize solver"<<endl;}

    if(!solver.solve()) {
        return;
    }

    QPSolution_ = solver.getSolution();
    // cout<<"solution:"<<endl;
    // cout<<QPSolution_<<endl;

    applyControl();
    visualize_mpcc_solution();

    solver.clearSolver();
}

void MPCC::getMatrices(Matrix<double,nx,nx>& Q, Matrix<double,nx, 1>& fx, Matrix<double, nx, 1>& dv_cost_dx, Matrix<double,nu, 1>& fu, Matrix<double,nx,nx>& Ad,
                                Matrix<double,nx, nu>& Bd, Matrix<double,nx,1>& x_prev, Matrix<double,nu,1>& u_prev){
    double x = x_prev(0);
    double y = x_prev(1);
    double yaw = x_prev(2);
    double theta = x_prev(3);
    double v = u_prev(0);
    double steer = u_prev(1);


    double dx_dtheta = track_.x_eval_d(theta);
    double dy_dtheta = track_.y_eval_d(theta);
    double d2x_dtheta2 = track_.x_eval_dd(theta);
    double d2y_dtheta2 = track_.y_eval_dd(theta);

    double numer = dx_dtheta*d2y_dtheta2 - dy_dtheta*d2x_dtheta2;
    double denom = pow(dx_dtheta,2) + pow(dy_dtheta, 2);

    double dphi_dtheta = numer/denom;

    double phi = atan2(dy_dtheta, dx_dtheta);

    cout<<"x_prev: "<<x_prev<<endl;
    cout<<"u_prev: "<<u_prev<<endl;
    cout<<"phi: "<<phi<<endl;
    cout<<"yaw: "<<yaw<<endl;

    double x_diff = x - track_.x_eval(theta);
    double y_diff = y - track_.y_eval(theta);

    Vector2d temp1, temp2;
    temp1<< dphi_dtheta, 1.0;
    temp2<< cos(phi), sin(phi);

    Matrix2d MC, ML;
    MC << x_diff, y_diff, dy_dtheta, -dx_dtheta;
    ML << -y_diff, x_diff, dx_dtheta, dy_dtheta;

    double deC_dtheta = temp1.transpose()* MC * temp2;
    double deL_dtheta = temp1.transpose()* ML * temp2;

    Matrix<double, nx, 2>  error_grad;
    Matrix<double, nx, 1>  eC_grad, eL_grad;

    Matrix2d weights;
    weights.setZero();
    weights.diagonal() << q_c, q_l;

    double my_answer = (x-track_.x_eval(theta))*cos(phi)*dphi_dtheta + (y-track_.y_eval(theta))*sin(phi)*dphi_dtheta
                       -sin(phi)*dx_dtheta + cos(phi)*dy_dtheta;


    eC_grad << -sin(phi), cos(phi), 0.0, deC_dtheta;
    //eL_grad << -cos(phi), -sin(phi), 0.0, deL_dtheta;
    eL_grad << 0.0, 0.0, 0.0, 0.0;

    error_grad << eC_grad, eL_grad;

    Q = error_grad * weights * error_grad.transpose();

    Matrix4d Q_t = Q.transpose();
    Q = 0.5*(Q+Q_t) + 0.00001* Matrix4d::Identity();  // force to be symetric


    double eC = -sin(phi)*x_diff + cos(phi)*y_diff;
   // double eL =  -cos(phi)*x_diff - sin(phi)*y_diff;
    double eL =  0.0;


    Vector2d e;
    e<< eC, eL;

    dv_cost_dx << 0.0, 0.0, q_v*v*sin(yaw-phi), -q_v*v*sin(yaw-phi)*dphi_dtheta;

    fx = q_c *2.0 * (eC-eC_grad.dot(x_prev)) * eC_grad  +  q_l*2.0 * (eL-eL_grad.dot(x_prev)) * eL_grad  ;
    //fx = 2*error_grad*weights*e - 2*error_grad*weights*error_grad.transpose()*x_prev;
    fx +=   dv_cost_dx;
    fu << -q_v * cos(yaw-phi),  0.0;

    Matrix<double,nx,nx> A;
    Matrix<double,nx,nu> B;
    /*Linearized Dynamics */
    A <<   0.0, 0.0, -v*sin(yaw), 0.0,
            0.0, 0.0,  v*cos(yaw), 0.0,
            0.0, 0.0,      0.0,      0.0,
            0.0, 0.0, -v*sin(yaw-phi), v*sin(yaw-phi)*dphi_dtheta;

    B <<   cos(yaw), 0.0,
            sin(yaw), 0.0,
            tan(steer)/CAR_LENGTH, v/(cos(steer)*cos(steer)*CAR_LENGTH), //v*(tan(steer)*tan(steer)+1.0)/CAR_LENGTH
            cos(yaw-phi), 0.0;


    //Discretize with Euler approximation
    Ad = Matrix4d::Identity() + A*Ts;
    //Ad = (A*Ts).exp();
    Bd = Ts*B;

    fx.setZero();
    fu.setZero();
}

void MPCC::applyControl() {
    float speed =  QPSolution_((N+1)*nx+nu);    // u_0 speed
    float steer = QPSolution_((N+1)*nx+nu+1);  // u_0 steer

    cout<<"speed_cmd: "<<speed<<endl;
    cout<<"steer_cmd: "<<steer<<endl;

    steer = min(steer, 0.41f);
    steer = max(steer, -0.41f);

    ackermann_msgs::AckermannDriveStamped ack_msg;
    ack_msg.header.stamp = ros::Time::now();
    ack_msg.drive.speed = speed;
    ack_msg.drive.steering_angle = steer;
    ack_msg.drive.steering_angle_velocity = 1.0;
    drive_pub_.publish(ack_msg);
}

void MPCC::visualize_mpcc_solution() {
    visualization_msgs::MarkerArray markers;

    // Plot estimated theta
    visualization_msgs::Marker dots;
    dots.header.stamp = ros::Time::now();
    dots.header.frame_id = "map";
    dots.id = rviz_id::THETA_EST;
    dots.ns = "theta_est";
    dots.type = visualization_msgs::Marker::POINTS;
    dots.scale.x = dots.scale.y = 0.05;
    dots.scale.z = 0.02;
    dots.action = visualization_msgs::Marker::ADD;
    dots.pose.orientation.w = 1.0;
    dots.color.r = 1.0;
    dots.color.a = 1.0;
    //dots.lifetime = ros::Duration();
    for (int i=0; i<N+1; i++){
        geometry_msgs::Point p;
        double theta_est = QPSolution_(i*nx+3);
        p.x = track_.x_eval(theta_est);
        p.y = track_.y_eval(theta_est);
        dots.points.push_back(p);
    }

    markers.markers.push_back(dots);

    visualization_msgs::Marker pred_dots;
    pred_dots.header.stamp = ros::Time::now();
    pred_dots.header.frame_id = "map";
    pred_dots.id = rviz_id::PREDICTION;
    pred_dots.ns = "predicted_positions";
    pred_dots.type = visualization_msgs::Marker::POINTS;
    pred_dots.scale.x = pred_dots.scale.y = pred_dots.scale.z = 0.08;
    pred_dots.action = visualization_msgs::Marker::ADD;
    pred_dots.pose.orientation.w = 1.0;
    pred_dots.color.g = 1.0;
    pred_dots.color.a = 1.0;
    // spline_dots.lifetime = ros::Duration();
    for (int i=0; i<N+1; i++){
        geometry_msgs::Point p;
        p.x = QPSolution_(i*nx);
        p.y = QPSolution_(i*nx+1);
        pred_dots.points.push_back(p);
    }
    markers.markers.push_back(pred_dots);

    visualization_msgs::Marker borderlines;
    borderlines.header.stamp = ros::Time::now();
    borderlines.header.frame_id = "map";
    borderlines.id = rviz_id::BORDERLINES;
    borderlines.ns = "borderlines";
    borderlines.type = visualization_msgs::Marker::LINE_LIST;
    borderlines.scale.x = 0.03;
    borderlines.action = visualization_msgs::Marker::ADD;
    borderlines.pose.orientation.w = 1.0;
    borderlines.color.b = 1.0;
    borderlines.color.a = 1.0;

    borderlines.points = border_lines;
    markers.markers.push_back(borderlines);

    mpcc_viz_pub_.publish(markers);
}

int main(int argc, char **argv){
    ros::init(argc, argv, "mpcc_node");
    ros::NodeHandle nh;
    MPCC mpcc(nh);
    ros::Rate rate(50);
    while(ros::ok()){
        ros::spinOnce();
        rate.sleep();
    }
    return 0;
}

