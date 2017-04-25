#include "kalman_filter.h"
#include <iostream>       //cout

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

#define PI 3.14159

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_laser_in, MatrixXd &H_radar_in,
                        MatrixXd &R_laser_in, MatrixXd &R_radar_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_laser_ = H_laser_in;
  H_radar_ = H_radar_in;
  R_laser_ = R_laser_in;
  R_radar_ = R_radar_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;

  //new state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  P_ = (I - K * H_laser_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
    * // y = z - h(x')
  */

  //1. The predicted measurement vector x​′​​ is a vector containing values in the form [​p​x​​,p​y​​,v​x​​,v​y​​​​]
  // In order to calculate y for the radar sensor, we need to convert x​′​​ to polar coordinates

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  float rho = sqrt(px * px + py * py); // range

  float phi;  //bearing
  if (px == 0 && py == 0)
    phi = 0;
  else{
    phi = atan2(py, px);
    // normalization check
    if (phi > PI){
      cout << "Caution: phi value is greater than PI: " << phi << endl;
      phi = fmod(phi,(2.*PI)) - 2*PI;

    }
    else if (phi < (-1*PI)){
      cout << "Caution: phi value is smaller than -PI: " << phi << endl;
      phi = fmod(phi,(2.*PI)) + 2*PI;
    }
  }
  float rho_dot; // radial rate/velocity
  // avoiding zero division
  if (rho == 0)
    rho_dot = 0;
  else
    rho_dot = (px * vx + py * vy) / rho;

  // 2. Update
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;

  VectorXd y = z - z_pred;
  MatrixXd Ht = H_radar_.transpose();
  MatrixXd S = H_radar_ * P_ * Ht + R_radar_;
  MatrixXd K = P_ * Ht * S.inverse();

  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_radar_) * P_;
}
