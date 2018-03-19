// BCInterpPt.cpp: implementation of the BCInterpPt class.
//
//////////////////////////////////////////////////////////////////////

#include "BCInterpPt.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BCInterpPt::BCInterpPt()
{

}

BCInterpPt::BCInterpPt(istream&is):InterpPt()
{
	int num = -990;
	double xcoord=0.0, ycoord=0.0;
	while(is) {
		if(!(is >> xcoord)){
			num = -992;
			break;
		}
		if(!(is >> ycoord)){
			num = -993;
			break;
		}
		break;
	}
	setn(num);
	assign(xcoord, ycoord);
}

