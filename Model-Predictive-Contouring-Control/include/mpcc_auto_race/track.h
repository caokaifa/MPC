//
// Created by yuwei on 11/20/19.
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
#include <mpcc_auto_race/CSVReader.h>
#include <Eigen/Dense>

const int SEARCH_RANGE = 10;
const double HALF_WIDTH_MAX = 1.1;
using namespace std;
typedef struct Point_ref{
    double x;
    double y;
    double theta;  // theta represents progress along centerline
    double half_width= 0.8;
}Point_ref;


class Track{
public:
    vector<Point_ref> centerline;
    tk::spline X_spline;
    tk::spline Y_spline;
    double length;
    double space;

    Track(string file_name, float space) : space(space){
        centerline.clear();
        vector<geometry_msgs::Point> waypoints;
        CSVReader reader(file_name);
        // Get the data from CSV File
        std::vector<std::vector<std::string> > dataList = reader.getData();
        // Print the content of row by row on screen
        for(std::vector<std::string> vec : dataList){
            geometry_msgs::Point wp;
            wp.x = std::stof(vec.at(0));
            wp.y = std::stof(vec.at(1));
            waypoints.push_back(wp);
        }
        // extract equally spaced points
        int curr = 0;
        int next =1;
        Point_ref p_start;
        p_start.x = waypoints.at(0).x; p_start.y = waypoints.at(0).y; p_start.theta = 0.0;
        centerline.push_back(p_start);
        float theta = 0.0;

        while(next < waypoints.size()){
            float dist = sqrt(pow(waypoints.at(next).x-waypoints.at(curr).x, 2)
                    +pow(waypoints.at(next).y-waypoints.at(curr).y, 2));
            float dist_to_start = sqrt(pow(waypoints.at(next).x-waypoints.at(0).x, 2)
                                       +pow(waypoints.at(next).y-waypoints.at(0).y, 2));
            if (dist>space){
                theta += dist;
                Point_ref p;
                p.x = waypoints.at(next).x; p.y = waypoints.at(next).y; p.theta = theta;
                centerline.push_back(p);
                curr = next;
            }
            next++;
            // terminate when finished a lap
            if (next > waypoints.size()/2 && dist_to_start<space){
                break;
            }
        }

        length = theta + sqrt(pow(centerline.back().x-waypoints.at(0).x, 2)
                              +pow(centerline.back().y-waypoints.at(0).y, 2));
        Point_ref p_final;
        p_final.x = waypoints.at(0).x; p_final.y = waypoints.at(0).y; p_final.theta = length;
        centerline.push_back(p_final);   //close the loop

        vector<double> X;
        vector<double> Y;
        vector<double> thetas;
        for (int i=0; i<centerline.size(); i++){
            X.push_back(centerline.at(i).x);
            Y.push_back(centerline.at(i).y);
            thetas.push_back(centerline.at(i).theta);
        }

        X_spline.set_points(thetas, X);
        Y_spline.set_points(thetas, Y);
    }

    double findTheta(double x, double y, double theta_guess, bool global_search= false){
        /* return: projected theta along centerline, theta is between [0, length]
         * */
             wrapTheta(theta_guess);

            int start, end;
            if(global_search){
                start = 0;
                end = centerline.size();
            }
            else{
                int mid = int(floor(theta_guess/space));
                start = mid - SEARCH_RANGE;
                end = mid + SEARCH_RANGE;
            }
            int min_ind, second_min_ind;
            double min_dist2 = 10000000.0;
            for (int i=start; i<end; i++){
                if (i>centerline.size()-1){ i=0;}
                if (i<0){ i=centerline.size()-1;}
                double dist2 = pow(x-centerline.at(i).x, 2) + pow(y-centerline.at(i).y, 2);
                if( dist2 < min_dist2){
                    min_dist2 = dist2;
                    min_ind = i;
                }
            }
            if (sqrt(min_dist2)>HALF_WIDTH_MAX && !global_search){
                return findTheta(x, y, theta_guess, true);
            }
            Eigen::Vector2d p, p0, p1;
            int min_ind_prev = min_ind-1;
            int min_ind_next = min_ind+1;
            if (min_ind_next>centerline.size()-1){ min_ind_next=0;}
            if (min_ind_prev<0){ min_ind_prev=centerline.size()-1;}

            //closest line segment: either [min_ind ,min_ind+1] or [min_ind,min_ind-1]
            if (pow(x-centerline.at(min_ind_next).x, 2) + pow(y-centerline.at(min_ind_next).y, 2) <
                    pow(x-centerline.at(min_ind_prev).x, 2) + pow(y-centerline.at(min_ind_prev).y, 2)){
                    second_min_ind = min_ind_next;
            }
            else{
                second_min_ind = min_ind_prev;
            }

            p(0) = x;  p(1) = y;
            p0(0) = centerline.at(min_ind).x;  p0(1) = centerline.at(min_ind).y;
            p1(0) = centerline.at(second_min_ind).x;  p1(1) = centerline.at(second_min_ind).y;

            double projection = abs((p - p0).dot(p1 - p0)/(p1 - p0).norm());
            double theta;

            if (min_ind > second_min_ind){
                theta = centerline.at(min_ind).theta - projection;
            }
            else{
                theta = centerline.at(min_ind).theta + projection;
            }

            if (theta>length){ theta -= length;}
            if (theta<0){ theta += length;}

            return theta;
    }

    void wrapTheta(double& theta){
        while(theta > length){ theta -= length;}
        while(theta < 0){theta += length; }
    }

    double x_eval(double theta){
        wrapTheta(theta);
        return X_spline(theta);
    }

    double y_eval(double theta){
        wrapTheta(theta);
        return Y_spline(theta);
    }

    double x_eval_d(double theta){
        wrapTheta(theta);
        return X_spline.eval_d(theta);
    }

    double y_eval_d(double theta){
        wrapTheta(theta);
        return Y_spline.eval_d(theta);
    }

    double x_eval_dd(double theta){
        wrapTheta(theta);
        return X_spline.eval_dd(theta);
    }

    double y_eval_dd(double theta){
        wrapTheta(theta);
        return Y_spline.eval_dd(theta);
    }

    double getPhi(double theta){
        wrapTheta(theta);

        double dx_dtheta = X_spline.eval_d(theta);
        double dy_dtheta = Y_spline.eval_d(theta);

        return atan2(dy_dtheta, dx_dtheta);
    }

    double getHalfWidth(double theta){
        wrapTheta(theta);
        int ind = int(floor(theta/space));
        ind = max(0, ind);
        ind= min(int(centerline.size()-1), ind);
        return centerline.at(ind).half_width;
    }

};