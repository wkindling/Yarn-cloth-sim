#include "node.h"
#include <iostream>

using namespace Eigen;
using namespace std;

Node::Node(Vector3d pos, double _u, double _v)
{
	position = pos;
	velocity.setZero();

	u = _u;
	v = _v;
	anchorU = _u;
	anchorV = _v;
	velocityUV.setZero();

	normal.setZero();
	//whichUp and index assigned externally
}

Node::~Node() {}


// Use SVD to get the vertex normal according to neighbors
Vector3d Node::getNormal()
{
	int neighbor_count = 0;
	Vector3d mean=Vector3d::Zero();
	mean += this->position;
	for (int i = 0; i < 4; i++)
	{
		if (neighbor[i] != NULL)
		{
			mean += neighbor[i]->position;
			neighbor_count++;
		}
	}
	mean = mean / (double)(1 + neighbor_count);
	MatrixXd A;
	A.resize(neighbor_count + 1, 3);

	A.row(0) = (this->position - mean).transpose();
	int row_index = 1;
	for (int i = 0; i < 4; i++)
	{
		if (neighbor[i] != NULL)
		{
			A.row(row_index) = (neighbor[i]->position - mean).transpose();
			row_index++;
		}
	}
	
	JacobiSVD<MatrixXd> svd(A, ComputeFullU | ComputeFullV);
	MatrixXd V = svd.matrixV();
	Vector3d nor = V.col(2);
	nor.normalize();
	normal = nor;

	return normal;
}
