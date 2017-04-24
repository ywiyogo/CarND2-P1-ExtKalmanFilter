#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define PI 3.14159

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
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
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;

  //new state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
    * // y = z - h(x')
  */
  
  //The predicted measurement vector x​′​​ is a vector containing values in the form [​p​x​​,p​y​​,v​x​​,v​y​​​​]
  // In order to calculate y for the radar sensor, we need to convert x​′​​ to polar coordinates
  MatrixXd Ht = H_.transpose();
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  float ro = sqrt(px * px + py * py); // range
  float phi;  //bearing
  if (px == 0 && py == 0)
    phi = 0;
  else
    phi = atan2(py, px);

  float ro_dot; // radial velocity
  if (ro == 0)
    ro_dot = 0;
  else
    ro_dot = (px * vx + py * vy) / ro;
  
  VectorXd z_pred(3);
  z_pred << ro, phi, ro_dot;
  VectorXd y = z - z_pred;
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
    
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}
