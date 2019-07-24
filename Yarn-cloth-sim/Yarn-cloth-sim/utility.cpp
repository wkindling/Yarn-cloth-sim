#include "utility.h"

void fillBlock(std::vector<T>& _K, Eigen::MatrixXd& block, int index0, int index1)
{
	for (int i = 0; i < block.rows(); i++)
	{
		for (int j = 0; j < block.cols(); j++)
		{
			_K.push_back(T(index0 + i, index1 + j, block(i, j)));
		}
	}
}
