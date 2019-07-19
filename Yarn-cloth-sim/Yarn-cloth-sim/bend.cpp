#include "bend.h"

using namespace std;
using namespace Eigen;

BendSpring::BendSpring(Node* n0, Node* n1, Node* n2, double _B, double _R, YarnType type)
{
	node0 = n0;
	node1 = n1;
	node2 = n2;

	B = _B;
	R = _R;
	Kb = B * M_PI*R*R;

	bendEnergy = 0;

	springType = type;
}

BendSpring::~BendSpring() {}

void BendSpring::solve()
{
	if (springType == Warp) solveU();

	else if (springType == Weft) solveV();

	return;
}


void BendSpring::solveU()
{


}

void BendSpring::solveV()
{


}
