#include "stretch.h"


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


}

void StretchSpring::solveU()
{

}

void StretchSpring::solveV()
{



}

