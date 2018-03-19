
//			InterpPt.h
//
//	Defines an class for a pair of numbers which are to form a list for interpolation,

//	Subclassed from Item, so it can be link listed.
//
//	Initiated:	January 7, 1997
//
//	Last Revised: January 7, 1997

#ifndef INTERPPT_H
#define INTERPPT_H

#include "Item.h"
#include <fstream>
using namespace std;

class InterpPt : public Item
{
	protected:
		double x;			// 	abcissa							
		double y;			// 	ordinate							
		
	public:
		InterpPt(int name=1,double xcoord=10.0,double ycoord=10.0);
		InterpPt(istream& is);
		double getXc() const {return x;}	
		double getYc() const {return y;}
		void assign(double xcoord,double ycoord);
		friend ostream& operator << (ostream& os, InterpPt* n);
};

#endif
