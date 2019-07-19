#include "stretch.h"

using namespace std;
using namespace Eigen;

StretchSpring::StretchSpring(Node* n0, Node* n1, double _Y, double _R, YarnType type)
{
	node0 = n0;
	node1 = n1;
	Y = _Y;
	R = _R;

	Ks = Y * M_PI*R*R;
	stretchEnergy = 0;
	springType = type;
}

StretchSpring::~StretchSpring() {}

void StretchSpring::solve()
{
	if (springType == Weft) solveV();
	
	else if (springType == Warp) solveU();

	return;
}

void StretchSpring::solveU()
{
	double x0 = node0->position.x();
	double y0 = node0->position.y();
	double z0 = node0->position.z();

	double x1 = node1->position.x();
	double y1 = node1->position.y();
	double z1 = node1->position.z();

	double u0 = node0->u;
	double u1 = node1->u;

	//Fill in the global matrix instead of calculating the force and updating locally




}

void StretchSpring::solveV()
{



}

