#pragma once
#ifndef UTILITY_H
#define UTILITY_H

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>

#define M_PI 3.1415926535

typedef Eigen::Triplet<double> T;

void fillGlobal(std::vector<T>& _M, Eigen::MatrixXd& block, int index0, int index1);

int sign(double number);

#endif