// PhiElementMtx.h: interface for the PhiElementMtx class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PHIELEMENTMTX_H
#define PHIELEMENTMTX_H

#include "ElementMtx.h"

class PhiElementMtx : public ElementMtx

{

public:
	PhiElementMtx(int name, Element *elP, int size, double wsElevIn);
	double wseIn;
	virtual ~PhiElementMtx();
	virtual void setElMtx();		//this function allocates the memory for the element (depending of the size)
									//and it also calculates the entries of the matrix

};

#endif
