#ifndef _CutOff_H_
#define _CutOff_H_

#include "UnitStruct.h"
#include "cavc/polyline.hpp"
#include "cavc/polylineintersects.hpp"

using namespace std;
using namespace cavc;

//≈≈–Ú
class indexcomp {
public:
	bool operator()(cavc::PlineIntersect<double>  l1, cavc::PlineIntersect<double> l2) {
		return l1.sIndex1 < l2.sIndex1;
	}
};


#endif // !_CutOff_H_