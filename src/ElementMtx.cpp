// ElementMtx.cpp: implementation of the ElementMtx class.
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
//#include "River2D.h"
#include "ElementMtx.h"
#include "Triangle.h"
#include "Habitat.h"
#include <math.h>


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ElementMtx::ElementMtx(int name, Element *elP, int size):Item(name)
{
	elementP = elP;
	mtxSize = size;
}

ElementMtx::~ElementMtx()
{
	delete [] KEMtx;
	delete [] FEMtx;
}

void ElementMtx::setElMtx() 
{

}

double ElementMtx::KE(int i, int j)
{
	return KEMtx[i*mtxSize + j];
}

