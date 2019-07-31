#include "utility.h"

void fillGlobalInertia(std::vector<T>& _M, Eigen::MatrixXd& block, int index0, int index1)
{
	for (int i = 0; i < block.rows(); i++)
	{
		for (int j = 0; j < block.cols(); j++)
		{
			_M.push_back(T(index0 + i, index1 + j, block(i, j)));
		}
	}
}

void fillGlobalStiffness(std::vector<T>& _M, Eigen::MatrixXd& block, int index0, int index1, double h)
{
	for (int i = 0; i < block.rows(); i++)
	{
		for (int j = 0; j < block.cols(); j++)
		{
			_M.push_back(T(index0 + i, index1 + j, -h*h*block(i, j)));
		}
	}
}

int sign(double number)
{
	return number >= 0 ? 1 : -1;
}