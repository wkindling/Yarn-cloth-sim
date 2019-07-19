#include "shear.h"

using namespace std;
using namespace Eigen;

ShearSpring::ShearSpring(Node* n0, Node* n1, Node* n2, double _S, double _R, double _L, YarnType type)
{
	node0 = n0;
	node1 = n1;
	node2 = n2;
	S = _S;
	R = _R;
	L = _L;
	
	Kx = S * R*R;
	shearEnergy = 0;
	springType = type;
}

ShearSpring::~ShearSpring() {}

void ShearSpring::solve()
{



}