// PsiElementMtx.cpp: implementation of the PsiElementMtx class.
//
//////////////////////////////////////////////////////////////////////


#include "PsiElementMtx.h"
#include "Triangle.h"
#include "Habitat.h"
#include <math.h>
#include "fe_PS.h"

extern "C" struct transient tvals;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


PsiElementMtx::PsiElementMtx(int name, Element *elP, int size):ElementMtx(name, elP, size)
{
	if(mtxSize == 2)
	{
		KEMtx = new double [4];
		FEMtx = new double [2];
	}
	if(mtxSize == 3)
	{
		KEMtx = new double [9];
		FEMtx = new double [3];
	}
	setElMtx();
}
PsiElementMtx::~PsiElementMtx()
{

}

void PsiElementMtx::setElMtx()
{
	if (mtxSize == 2)  //boundary element
	{
		RivBound *segP = (RivBound*) elementP;

		HabitatNode *hNP0 = (HabitatNode*)segP->getNode(0);
		HabitatNode *hNP1 = (HabitatNode*)segP->getNode(1);

		double X21 = hNP1->getXc() - hNP0->getXc();
		double Y21 = hNP1->getYc() - hNP0->getYc();

		double le = /*fabs(sqrt(X21*X21 + Y21*Y21))*/1.0; //length of the element

		double psi = (hNP0->getPsi() + hNP1->getPsi())/2.0;  // psi for b.c.

		//penalty added at the solution level
		
		KEMtx[0] = le/3.0;		// -   -
		KEMtx[1] = le/6.0;		// |0 1|
		KEMtx[2] = KEMtx[1];	// |2 3|
		KEMtx[3] = KEMtx[0];	// -   -

		FEMtx[0]= psi*le/2.0;
		FEMtx[1]= FEMtx[0];

	}

	if(mtxSize == 3)	//internal element
	{
		Triangle *triP = (Triangle*) elementP;
		
		HabitatNode *hNP0 = (HabitatNode*)triP->getNode(0);
		HabitatNode *hNP1 = (HabitatNode*)triP->getNode(1);
		HabitatNode *hNP2 = (HabitatNode*)triP->getNode(2);

		double qx1 = hNP0->getPar(4);
		double qy1 = hNP0->getPar(5);
		if(hNP0->getPar(3) < tvals.gwH) { qx1 = 0; qy1 = 0;};

		double qx2 = hNP1->getPar(4);
		double qy2 = hNP1->getPar(5);
		if(hNP1->getPar(3) < tvals.gwH) { qx2 = 0; qy2 = 0;};
		double qx3 = hNP2->getPar(4);
		double qy3 = hNP2->getPar(5);
		if(hNP2->getPar(3) < tvals.gwH) { qx3 = 0; qy3 = 0;};

		double X21 = hNP1->getXc() - hNP0->getXc();
		double X13 = hNP0->getXc() - hNP2->getXc();
		double X32 = hNP2->getXc() - hNP1->getXc();

		double Y12 = hNP0->getYc() - hNP1->getYc();
		double Y31 = hNP2->getYc() - hNP0->getYc();
		double Y23 = hNP1->getYc() - hNP2->getYc();

		double detJ = fabs(X13*Y23 - Y31*X32);

		KEMtx[0] = (Y23*Y23 + X32*X32)/(2*detJ);	// -	 -
		KEMtx[1] = (Y23*Y31 + X32*X13)/(2*detJ);	// |0 1 2|
		KEMtx[3] = KEMtx[1];						// |	 |
		KEMtx[4] = (Y31*Y31 + X13*X13)/(2*detJ);	// |3 4 5|
		KEMtx[2] = (Y23*Y12 + X32*X21)/(2*detJ);	// |	 |
		KEMtx[6] = KEMtx[2];						// |6 7 8|
		KEMtx[8] = (Y12*Y12 + X21*X21)/(2*detJ);	// -	 -
		KEMtx[5] = (Y31*Y12 + X13*X21)/(2*detJ);	
		KEMtx[7] = KEMtx[5];

		FEMtx[0] =(1.0/6.0)*(Y23*qy1 + Y31*qy2 + Y12*qy3 - X32*qx1 - X13*qx2 - X21*qx3);
		FEMtx[1] = FEMtx[0];
		FEMtx[2] = FEMtx[0];
	}
}