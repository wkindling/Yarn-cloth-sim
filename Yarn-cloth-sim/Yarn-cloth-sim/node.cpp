#include "node.h"

using namespace Eigen;

Node::Node(Vector3d pos, double _u, double _v)
{
	position = pos;
	velocity.setZero();

	u = _u;
	v = _v;
	anchorU = _u;
	anchorV = _v;
	velocityUV.setZero();

	//whichUp and index assigned externally
}

Node::~Node() {}


