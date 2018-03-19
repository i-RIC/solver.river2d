// PsiElementMtx.h: interface for the PsiElementMtx class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PSIELEMENTMTX_H
#define PSIELEMENTMTX_H

#include "ElementMtx.h"

class PsiElementMtx : public ElementMtx

{

public:
	PsiElementMtx(int name, Element *elP, int size);
	virtual ~PsiElementMtx();
	virtual void setElMtx();		//this function allocates the memory for the element (depending of the size)
									//and it also calculates the entries of the matrix

};

#endif
