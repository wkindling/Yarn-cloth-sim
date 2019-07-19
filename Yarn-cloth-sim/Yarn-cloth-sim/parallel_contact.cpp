#include "parallel_contact.h"

using namespace std;
using namespace Eigen;

ParallelContactSpring::ParallelContactSpring(Node* n0, Node* n1, double _Kc, double _R, double _L, YarnType type)
{
	node0 = n0;
	node1 = n1;

	Kc = _Kc;
	L = _L;
	R = _R;

	parallelContactEnergy = 0;

	springType = type;
}

ParallelContactSpring::~ParallelContactSpring() {}
