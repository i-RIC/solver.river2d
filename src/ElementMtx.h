// ElementMtx.h: interface for the ElementMtx class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ELEMENTMTX_H__C30CA8E0_DB75_44D1_A1DE_A5082E2E1832__INCLUDED_)
#define AFX_ELEMENTMTX_H__C30CA8E0_DB75_44D1_A1DE_A5082E2E1832__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Item.h"
#include "Element.h"

class ElementMtx : public Item  
{

protected:
	double *KEMtx;			//pointer to the element matrix (either 3X3 (triangular) 2X2 (linear))
	double *FEMtx;			//RHS contribution of the equation at the element level
	Element *elementP;		//A pointer to the element for which the element matrix is for
	int mtxSize;			//the size of the element matrix
public:
	ElementMtx(int name, Element* elP, int size); //contructor calls setElMtx()
	virtual ~ElementMtx();
	virtual void setElMtx();//this function allocates the memory for the element (depending of the size)
							//and it also calculates the entries of the matrix
	double KE(int i, int j);//this function serves as an easy way to access the different entries of the matrix
	double FE(int i) {return FEMtx[i];}//returns the RHS contibution
	Element* getElement() {return elementP;}
	int getMtxSize() {return mtxSize;}

};

#endif // !defined(AFX_ELEMENTMTX_H__C30CA8E0_DB75_44D1_A1DE_A5082E2E1832__INCLUDED_)
