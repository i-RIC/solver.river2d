
//		InterpPt.cpp

//		InterpPt class function definitions

#include <iostream>
#include "Interppt.h"
using namespace std;

InterpPt::InterpPt(	int name,double xcoord,double ycoord)
			:Item(name)
{
	assign(xcoord, ycoord);
}

InterpPt::InterpPt(istream& is) : Item()
{
	int num = -990;
	double xcoord=0.0, ycoord=0.0;
	while(is) {
		if(!(is >> num)){
			num = -991;
			break;
		}
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

void InterpPt::assign(double xcoord,double ycoord)
{
	x = xcoord;
	y = ycoord;
}

ostream& operator <<(ostream& os, InterpPt* n)
{
	if(n != NULL){
		os.width(5); 	os 	<< n->n << " ";
		os.precision(6);
		os.setf(ios::fixed,ios::floatfield);
		os 	<< n->x << " ";
		os 	<< n->y << " ";
		os << "\n";

	}

	return os;
}

