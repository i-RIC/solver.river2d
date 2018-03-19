// PhiElementMtx.cpp: implementation of the PhiElementMtx class.
//
//////////////////////////////////////////////////////////////////////


#include "PhiElementMtx.h"
#include "Triangle.h"
#include "Habitat.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


PhiElementMtx::PhiElementMtx(int name, Element *elP, int size, double wsElevIn):ElementMtx(name, elP, size)
{
	wseIn = wsElevIn;
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
PhiElementMtx::~PhiElementMtx()
{

}

void PhiElementMtx::setElMtx()
{
	if (mtxSize == 2)  //boundary element
	{
		RivBound *segP = (RivBound*) elementP;

		HabitatNode *hNP0 = (HabitatNode*)segP->getNode(0);
		HabitatNode *hNP1 = (HabitatNode*)segP->getNode(1);

		double X21 = hNP1->getXc() - hNP0->getXc();
		double Y21 = hNP1->getYc() - hNP0->getYc();

		double le = /*fabs(sqrt(X21*X21 + Y21*Y21))*/1.0; //length of the element

		double phi;
		
		if(segP->getBcCode() == 1)		//inflow
			phi = wseIn;
		else if(segP->getBcCode() == 3)	//outflow constant elevation
			phi = segP->getBcValue();
		else							//outflow rating curve
			phi = (hNP0->getPar(6) + hNP1->getPar(6))/2.0;

		//penalty added at the solution level
		
		KEMtx[0] = le/3.0;		// -   -
		KEMtx[1] = le/6.0;		// |0 1|
		KEMtx[2] = KEMtx[1];	// |2 3|
		KEMtx[3] = KEMtx[0];	// -   -

		FEMtx[0]= phi*le/2.0;
		FEMtx[1]= FEMtx[0];

	}

	if(mtxSize == 3)	//internal element
	{
		Triangle *triP = (Triangle*) elementP;
		
		HabitatNode *hNP0 = (HabitatNode*)triP->getNode(0);
		HabitatNode *hNP1 = (HabitatNode*)triP->getNode(1);
		HabitatNode *hNP2 = (HabitatNode*)triP->getNode(2);

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

		FEMtx[0] =0;
		FEMtx[1] = FEMtx[0];
		FEMtx[2] = FEMtx[0];
	}
}