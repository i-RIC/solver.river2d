// HabitatTIN.cpp: implementation of the HabitatTIN class.
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
//#include "River2D.h"
#include "HabitatTIN.h"
#include "Shallow.h"
#include "ElementMtx.h"
#include "PsiElementMtx.h"
#include "PhiElementMtx.h"
#include "Solver.h"
#include <math.h>


//#ifdef _DEBUG
//#undef THIS_FILE
//static char THIS_FILE[]=__FILE__;
//#define new DEBUG_NEW
//#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

HabitatTIN::HabitatTIN(Physics *theProb):TIN(theProb)
{
 
}

HabitatTIN::~HabitatTIN()
{

}

HabitatNode* HabitatTIN::insertOneFloatingNode(double x, double y, TIN& boundTIN, TIN& bedTIN, TIN& meshTIN)
//	Puts one Node into the TIN (from physics spec).
//	Requires a three previously triangulated TIN.
//	boundTIN is used when the new node is being added within a user defined boundary - eg region refine
//	Bed elevation and roughness of the new node are interpolated from the bedTIN
//	while the other parameters (qx, qy, water surface elevation) are interpolated from 
//	the meshTIN.
//	The esulting node is added to tnodeL and  can be later accepted into nodeL or rejected.

{
	int nGoodNodes = 0;
	int nodeNum = getNextN();
	Triangle *tBedP, *tBoundP, *tMeshP, *triP;
	HabitatNode *nP;
	ShallowNode sN;
	double dz;

	nP = (HabitatNode *)physics->makeNewNode(nodeNum,x,y);
	if( (tBedP = bedTIN.whichTriangle(nP)) != NULL){
		if( (tMeshP = meshTIN.whichTriangle(nP)) != NULL){
			if( (tBoundP = boundTIN.whichTriangle(nP)) != NULL){
				if(tBoundP->getStatus() == active){
					tMeshP->locateNodeAt(x,y,nP);
					sN.assign(nP->getXc(), nP->getYc(), nP->getZc());
					triP = bedTIN.whichTriangle(&sN);
					if(triP != NULL)
					{
						triP->locateNodeAt(&sN);
						nP->setKs(sN.getKs());
						dz = sN.getZc() - nP->getZc();
						nP->setZc(nP->getZc() + dz);
						nP->setD(nP->getD() - dz);
					}

					nGoodNodes += 1;
					nodeNum += 1;
					tnodeL.appendItem(nP);
				}
				else{
					delete nP;
					nP = NULL;
				}
			}
			else{
				delete nP;
				nP = NULL;
			}
		}
		else{
			delete nP;
			nP = NULL;
		}
	}
	else{
		delete nP;
		nP = NULL;
	}

	return nP;
}

double HabitatTIN::smoothHabitatMesh(int nTimes, TIN& dataTIN, double bias)
{
	HabitatNode *nP;
	ShallowNode sN;
	Triangle *tP, *atP, *dtP;
	Node *nP1;
	double a, sumx, sumy, suma, *xNew, *yNew, xN, yN, dx, dy, lx, ly, d, dz;
	int	j, index, iin=0;

   bias = sqrt(bias);

	xNew = new double[nodeL.numberOfItems()];
	yNew = new double[nodeL.numberOfItems()];

	for(int i=0;i<nTimes;i++){
		nP = (HabitatNode *) nodeL.firstItem();
	
		while(nP != NULL){
			xNew[iin] = nP->getXc();
			yNew[iin] = nP->getYc();
			if( (nP->getFixed() == floating) && (nP->getBound() == interior) ) {
				if((tP = whichTriHasNode(nP)) != NULL) {
					sumx = sumy = suma = 0.0;
					atP = tP;
					do {
						a = ((1.0-bias) + bias*fabs(dataTIN.triCenterDz(atP))) * sqrt(fabs(atP->area())); // atP->area();
						suma += 4.0*a;
						sumx += nP->getXc() * a;
						sumy += nP->getYc() * a;
						for(j=0;j<3;j++) {
							nP1 = atP->getNode(j);
							sumx += nP1->getXc() * a;
							sumy += nP1->getYc() * a;
							if(nP1 == nP)
								index = j;
						}
						index -= 1;
						if(index == -1)
							index = 2;
						atP = (Triangle *) atP->getAdjt(index);
					} while (atP != tP);
					xNew[iin] = (nP->getXc() + sumx/suma) *0.5;
					yNew[iin] = (nP->getYc() + sumy/suma) *0.5;
				}
			}
			nP = (HabitatNode *) nodeL.nextItem();
			iin += 1;
		}

		Segment *fsP = firstfSeg();
		Segment *nsP;
		while (fsP != NULL){
			nP = (HabitatNode *)fsP->getNode(1);
			if(nP->getFixed() == sliding) {
				tP = whichTriHasNode(nP);
				sumx = sumy = suma = 0.0;
				atP = tP;
				do {
					a = ((1.0-bias) + bias*fabs(dataTIN.triCenterDz(atP))) * sqrt(fabs(atP->area())); // atP->area();
					suma += 4.0*a;
					sumx += nP->getXc() * a;
					sumy += nP->getYc() * a;
					for(j=0;j<3;j++) {
						nP1 = atP->getNode(j);
						sumx += nP1->getXc() * a;
						sumy += nP1->getYc() * a;
						if(nP1 == nP)
							index = j;
					}
					index -= 1;
					if(index == -1)
						index = 2;
					atP = (Triangle *) atP->getAdjt(index);
				} while (atP != tP);
				xN = (nP->getXc() + sumx/suma) *0.5;
				yN = (nP->getYc() + sumy/suma) *0.5;
				if( (nsP = (Segment *) fsP->getNextOne()) != NULL){
					if(nsP->getNode(0) == nP){
						dx = xN - fsP->getNode(0)->getXc();
						dy = yN - fsP->getNode(0)->getYc();
						lx = nsP->getNode(1)->getXc() - fsP->getNode(0)->getXc();
						ly = nsP->getNode(1)->getYc() - fsP->getNode(0)->getYc();
						d = (dx*lx + dy*ly)/(lx*lx + ly*ly);
						if(d < 0.33333)
							d = 0.33333;
						else if(d > 0.66667)
							d = 0.66667;
						xN = fsP->getNode(0)->getXc() + d*lx;
						yN = fsP->getNode(0)->getYc() + d*ly;

						sN.assign(xN,yN,100.0);
						dtP = dataTIN.whichTriangle(&sN);
						if (dtP != NULL)
						{
							dtP->locateNodeAt(&sN);
							dz= sN.getZc() - nP->getZc();
							nP->assign(xN,yN,nP->getZc()+dz);
							nP->setKs(sN.getKs());
							nP->setD(nP->getD()-dz);
						}

					}
				}
			}
			fsP = nextfSeg();
		}

		Segment *bsP = firstSeg();
		Element *elP;
		while (bsP != NULL){
			nP = (HabitatNode *) bsP->getNode(1);
			if(nP->getFixed() == sliding) {
				tP = (Triangle *) bsP->getAdjt(0);
				sumx = sumy = suma = 0.0;
				atP = tP;
				do {
					a = sqrt(fabs(atP->area()));
					suma += 3.0*a;
					for(j=0;j<3;j++) {
						nP1 = atP->getNode(j);
						sumx += nP1->getXc() * a;
						sumy += nP1->getYc() * a;
						if(nP1 == nP)
							index = j;
					}
					index -= 1;
					if(index == -1)
						index = 2;
					elP = atP->getAdjt(index);
					if(elP->numNodes() == 3)
						atP = (Triangle *) elP;
					else
						break;
				} while (atP != tP);
				xN = (nP->getXc() + sumx/suma) *0.5;
				yN = (nP->getYc() + sumy/suma) *0.5;
				nsP = (Segment *) bsP->getNextOne();
				dx = xN - bsP->getNode(0)->getXc();
				dy = yN - bsP->getNode(0)->getYc();
				lx = nsP->getNode(1)->getXc() - bsP->getNode(0)->getXc();
				ly = nsP->getNode(1)->getYc() - bsP->getNode(0)->getYc();
				d = (dx*lx + dy*ly)/(lx*lx + ly*ly);
				if(d < 0.33333)
					d = 0.33333;
				else if(d > 0.66667)
					d = 0.66667;
				xN = bsP->getNode(0)->getXc() + d*lx;
				yN = bsP->getNode(0)->getYc() + d*ly;

				sN.assign(xN,yN,100.0);
				dtP = dataTIN.whichTriangle(&sN);
				if (dtP != NULL)
				{
					dtP->locateNodeAt(&sN);
					dz= sN.getZc() - nP->getZc();
					nP->assign(xN,yN,nP->getZc()+dz);
					nP->setKs(sN.getKs());
					nP->setD(nP->getD()-dz);
				}
				
			}
			bsP = nextSeg();
		}

		nP = (HabitatNode *) firstNode();
		iin = 0;
		while(nP != NULL){
			if( nP->getFixed() == floating ) {

				sN.assign(xNew[iin],yNew[iin],100.0);
				dtP = dataTIN.whichTriangle(&sN);
				if (dtP != NULL)
				{
					dtP->locateNodeAt(&sN);
					dz= sN.getZc() - nP->getZc();
					nP->assign(xNew[iin],yNew[iin],nP->getZc()+dz);
					nP->setKs(sN.getKs());
					nP->setD(nP->getD()-dz);
				}

			}
			nP = (HabitatNode *) nextNode();
			iin += 1;
		}
	}
	delete [] xNew;
	delete [] yNew;

	return 1.0;
}

int HabitatTIN::refineHabitatRegion(TIN &boundTIN, TIN &dataTIN)
{
//	Refines a HabitatTIN by placing a new node in the centre of each existing triangle
//	within the boundTIN region.
//	Requires a previously triangulated dataTIN, dataTIN,
//	which masks inactive areas and provides interpolation data.
//	Resulting set of nodes are in tnodeL which can be later accepted into nodeL or rejected.
//	The number of nodes created are returned.

{
	int nGoodNodes = 0;
	int nodeNum = getNextN();
	Triangle *tP, *tBP, *eTP;
	ShallowNode cNode;
	HabitatNode *nP;
	double dz;

	eTP = firstTri();  //meshP
	while(eTP != NULL){
		if(eTP->getStatus() == active){
			nP = (HabitatNode *)physics->makeNewNode(nodeNum);
			eTP->locateNodeAtCenter(nP); //interpolates a new node based on mesh
			cNode.assign(nP->getXc(),nP->getYc(), nP->getZc()); //assigns the shallownode the same x and y
			if( (tP = dataTIN.whichTriangle(nP)) != NULL){
				if( (tBP = boundTIN.whichTriangle(nP)) != NULL){
					if(tBP->getStatus() == active){
						tP->locateNodeAt(&cNode);
						nP->setKs(cNode.getKs());
						dz = cNode.getZc() - nP->getZc();
						nP->setZc(nP->getZc()+ dz);
						nP->setD(nP->getD() - dz);
						nGoodNodes += 1;
						nodeNum += 1;
						tnodeL.appendItem(nP);
					}
					else
						delete nP;
				}
				else
					delete nP;
			}
			else
				delete nP;
		}
      eTP = nextTri();
	}
	return nGoodNodes;
}
}

int HabitatTIN::checkNoOfInflowSegs()
{
	RivBound *bsP;
	int prevBcCode = 0, currentBcCode, Nseg = 0;
	
	
	bsP = (RivBound *)firstSeg();

	while(bsP != NULL)
	{
		currentBcCode = bsP->getBcCode();
		if(currentBcCode != prevBcCode)
		{
			if (currentBcCode == 1 )
			{
				Nseg += 1;
			}
			prevBcCode = currentBcCode;
		}
		bsP = (RivBound*)nextSeg();
	}
	return Nseg;
}

Segment* HabitatTIN::setBoundaryPsi(HabitatNode* startnP)
{
	RivBound *bsP, *startbsP;
	HabitatNode *nP1, *nP2, *firstnP;
	int prevBcCode, currentBcCode;
	double qx, qy, dx, dy, discharge, d, QinTotal = 0, QoutTotal= 0, erTotal, erPercent;
	
	// Calculate the error between inflow and outflow to be applied to the outflow
	// Also resetting the psi values along the boundary to 0
	
	bsP = (RivBound *)firstSeg();
	nP1 = (HabitatNode*)bsP->getNode(0);
	nP1->setPsi(0);
	
	while(bsP != NULL)
	{
		nP1 = (HabitatNode*)bsP->getNode(0);
		nP2 = (HabitatNode*)bsP->getNode(1);
		nP2->setPsi(0);

		if(bsP->getBcCode()==1)
		{
			dx = nP2->getXc()-nP1->getXc();
			dy = nP2->getYc()-nP1->getYc();
			d = fabs(sqrt(dx*dx+dy*dy));
			discharge = bsP->getBcValue()*d;
			QinTotal += discharge;
		}
		if(bsP->getBcCode()==3 || bsP->getBcCode()==5)
		{
			dx = nP2->getXc()-nP1->getXc();
			dy = nP2->getYc()-nP1->getYc();
			qx = (nP1->getPar(4)+nP2->getPar(4))/2;
			qy = (nP1->getPar(5)+nP2->getPar(5))/2;
			discharge = qx*dy - qy*dx;
			QoutTotal += discharge;
		}

		bsP = (RivBound *)nextSeg();
	}

	erTotal = QinTotal - QoutTotal;

	//check to indicate to the user if inflow/outflow error is too large
	
	erPercent = erTotal / QinTotal;

	if( erPercent > 0.10 )
		//AfxMessageBox("Note: The solution may not have yet reached steady state.\n\n"
		//			  "The net outflow is more than 10 percent of the total inflow.\n"
		//			  "As a result, the calculated cumulative discharge may\nnot be very accurate.");

	// Finding the starting node at a corner where an inflow boundary element 
	// meets a noflow boundary element

	bsP = (RivBound *)firstSeg();
	startbsP = bsP;

	firstnP = (HabitatNode*)bsP->getNode(0);  //firstnP is first node in the list defining the external boundary

	if(startnP == NULL) // for the case of 1 inflow
	{
		do
		{
			prevBcCode = bsP->getBcCode();
			bsP = (RivBound *)nextSeg();
			currentBcCode = bsP->getBcCode();
			if (prevBcCode == 1 && currentBcCode == 0)
				startbsP = bsP;
	
		}while (bsP->getNode(1)!=firstnP);
	}
	else // multiple inflows - starting node startnP input by user 
		 // this loop finds the starting boundary segment based on the starting node
	{
		while (bsP != NULL)
		{
			if(bsP->getNode(0) == startnP)
				startbsP = bsP;
			bsP = (RivBound *)nextSeg();
		}
	}

	//Using the new starting point cycle through the boundary
	//element list and set the stream function value - previously
	//initialized to zero throughout the entire domain

	bsP = startbsP;
	bsegL.setCurrentItem(bsP); // had to add this line so that bsP is the "current"
							   // item on the list, without this line it may not be

	startnP = (HabitatNode *)startbsP->getNode(0);	//startnP is the starting node for setting the boundary psi values
	do
	{
		nP1 = (HabitatNode*)bsP->getNode(0);
		nP2 = (HabitatNode*)bsP->getNode(1);
		
		if(bsP->getBcCode()==0)
		{
			nP2->setPsi(nP1->getPsi());
		}
		if(bsP->getBcCode()==1)
		{
			dx = nP2->getXc()-nP1->getXc();
			dy = nP2->getYc()-nP1->getYc();
			d = fabs(sqrt(dx*dx+dy*dy));
			discharge = bsP->getBcValue()*d;
			nP2->setPsi(nP1->getPsi()-discharge);
		}
		if(bsP->getBcCode()==3 || bsP->getBcCode()==5)
		{
			dx = nP2->getXc()-nP1->getXc();
			dy = nP2->getYc()-nP1->getYc();
			qx = (nP1->getPar(4)+nP2->getPar(4))/2;
			qy = (nP1->getPar(5)+nP2->getPar(5))/2;
			discharge = (qx*dy - qy*dx)*(erTotal/QoutTotal + 1);
			nP2->setPsi(nP1->getPsi()+discharge);
		}
		
		//this if else is here so that we skip the internal boundaries
		
		if( firstnP != startnP && bsP->getNode(1) == firstnP) 
			bsP = (RivBound *) firstSeg();  //need to check for firstnP == startnP or get caught in loop 

		else
		{
			bsP = (RivBound *) nextSeg();
			if(bsP == NULL)
				bsP = (RivBound *) firstSeg();
		}
		
	}while (bsP->getNode(1) != startnP);

	bsegL.setCurrentItem(startbsP);
	return startbsP;

}

void HabitatTIN::setAllPsi(HabitatNode* startnP)
{
	RivBound *startbsP, *bsP;
	HabitatNode *startbnP, *nP, *nP1;
	Element *aElP;
	Triangle *atP, *tP;
	ItemList tList, aList;
	double sum;
	int count;

	int j, index;
		
	double qx, qy, dx, dy, discharge;

	//aList is for nodes with psi set that have not 
	//been checked for adjacent nodes boundary or interior nodes
	
	//tList is for nodes that have psi set and have been checked for
	//adjacent nodes

	//before starting...create an array of pointers and store a pointer
	//to each node which will be used to restore the node ordering after
	//sweeping through the domain to set the psi values
	//also setting n for all nodes to 0 as this value is used as a flag
	//to see whether or not a node has been removed from nodeL
	//n = 0 -> still on nodeL
	//n = 1 -> no longer on nodeL

	Node **nodeList = new Node*[nodeL.numberOfItems()];
	Node *nodeP = firstNode();
	int inode = 0;

	while (nodeP != NULL)
	{
		nodeList[inode] = nodeP;
		nodeP->setn(0);
		inode++;
		nodeP = nextNode();
	}

	startbsP = (RivBound *)setBoundaryPsi(startnP);
	
	//cycle through all EXTERIOR boundary nodes and put on the aList.
	
	bsP = (RivBound *)firstSeg();
	startbnP = (HabitatNode*)bsP->getNode(0);

	while (-1)								
	{
		nP = (HabitatNode*)bsP->getNode(0);
		nodeL.deleteItem(nP);
		aList.appendItem(nP);
		nP->setn(1);
		if(bsP->getNode(1) == startbnP) break;
		bsP = (RivBound *)nextSeg();	
	}

	//need to reorder tList so that startbnP is at the top of the list

	startbnP = (HabitatNode*)startbsP->getNode(0);
	nP = (HabitatNode *)aList.firstItem();

	while(nP != startbnP)
	{
		aList.deleteItem(nP);			
		aList.appendItem(nP);
		nP = (HabitatNode*) aList.firstItem();
	}

	//cycle through aList which starts as only boundary nodes
	//start by popping the top item off the list and pushing the adjacent
	//nodes onto the aList once their psi value has been set
	//and pushing the node that has been popped of onto the tList.

	nP = (HabitatNode *)aList.firstItem();
	tList.emptyList();
	
	while (nP != NULL)
	{
		
		tList.appendItem(aList.pop());	//had to change to appendItem from 
										//push because push cannot
										//cannot set previous item which 
										//causes deleteItem to crash
		if((tP = whichTriHasNode(nP)) != NULL) 
		{
			atP = tP;
			do {
				for(j=0;j<3;j++) 
				{
					nP1 = (HabitatNode *)atP->getNode(j);
					if(nP1 == nP)
						index = j;
				}
				index += 1;
				if(index == 3)
					index = 0;
				aElP = atP->getAdjt(index);
				if(aElP == NULL) //added because in orderNodesRCM
					break;		 //added because in orderNodesRCM
				if(aElP->numNodes() == 2)
					break;
				atP = (Triangle *) aElP;
			} while (atP != tP);
			tP = atP;
			do {
				for(j=0;j<3;j++) 
				{
					nP1 = (HabitatNode *)atP->getNode(j);
					if(nP1 == nP)
						index = j;
					else 
					{
						if(nP1->getn() == 0)
						{
							dx = nP1->getXc()-nP->getXc();
							dy = nP1->getYc()-nP->getYc();
							qx = (nP1->getPar(4)+nP->getPar(4))/2;
							qy = (nP1->getPar(5)+nP->getPar(5))/2;
							discharge = (qx*dy - qy*dx);
						//	if (nP1->getPar(3) < 0)
						//		discharge = 0;
							nP1->setPsi(nP->getPsi()+discharge);
							nodeL.deleteItem(nP1);
							aList.push(nP1);
							nP1->setn(1);
						
						}
					}
				}
				index -= 1;
				if(index == -1)
					index = 2;
				aElP = atP->getAdjt(index);
				if(aElP == NULL) //added because in orderNodesRCM
					break;		 //added because in orderNodesRCM	
				if(aElP->numNodes() == 2)
					break;
				atP = (Triangle *) aElP;
			} while (atP != tP);
		}

		
	nP = (HabitatNode*) aList.firstItem();
	}

	//restore the node ordering that existed before
	//the sweep using the array of node pointers

	for(int ii=0; ii<inode; ii++)
	{
		tList.deleteItem(nodeList[ii]);
		nodeL.appendItem(nodeList[ii]);
	}

	delete [] nodeList;

	//reset n for each node according to
	//the node ordering

	nP = (HabitatNode *) nodeL.firstItem();
	int nn = 1;
	while (nP != NULL)
	{
		nP->setn(nn);
		nn++;
		nP = (HabitatNode *) nP->getNextOne();
	}
	
	//Now that psi has been calculated everywhere, go back, average the psi
	//values for each internal boundary and reset each value to the average

	bsP = (RivBound *)firstSeg();
	startbnP = (HabitatNode*)bsP->getNode(0);
	//loop to get to the first internal boundary
	do
		{
			bsP = (RivBound *)nextSeg();	
		}while (bsP->getNode(1)!=startbnP);
	
	bsegL.setCurrentItem(bsP);
	bsP = (RivBound *)bsegL.nextItem();
	if(bsP != NULL)
		nP = (HabitatNode*)bsP->getNode(0);
	
	while (bsP != NULL)
	{
		startbnP = (HabitatNode*)bsP->getNode(0);
		startbsP = bsP;
		sum = 0;
		count = 0;
		//loop through internal boundary first time
		//to sum up the psi values
		while(-1)
		{
			nP = (HabitatNode*)bsP->getNode(0);
			sum += nP->getPsi();
			count += 1;
			if(bsP->getNode(1) == startbnP) break;
			bsP = (RivBound *)nextSeg();
		}
		double psi = sum/count;
		nP = (HabitatNode*)startbsP->getNode(0);
		bsP = startbsP;
		bsegL.setCurrentItem(bsP);
		//loop through again to reset the psi values
		while(-1)
		{
			nP = (HabitatNode*)bsP->getNode(0);
			nP->setPsi(psi);
			if(bsP->getNode(1) == startbnP) break;
			bsP = (RivBound *)nextSeg();
		}
		//go to the next internal boundary if any...
		bsP = (RivBound *)nextSeg();
	}

}

int HabitatTIN::extractBoundary(TIN* meshTIN, HabitatNode* hNUSBound, HabitatNode* hNDSBound, double leftBound, double rightBound, double spacing, int leftBoundType, int rightBoundType)
{
	HabitatNode *hNP1, *hNP2, *hNP3, *hNP4;
	Triangle *triP, *checkTriP, *prevTriP;
	RivBound *boundSegDS, *boundSegUS;
	int segNum;
	double wse; //downstream water surface elevation
	double dischargeIntensity;

	if (rightBoundType == 2)  //right boundary chosen along computation boundary
	{
		mapBoundaryStreamLine(meshTIN, hNUSBound, hNDSBound, 0);
		hNP1 = (HabitatNode *)firstNode();
		hNP2 = (HabitatNode *)nodeL.lastItem();
	}
	else if (rightBoundType == 0) //right boundary chosen along a streamline
	{
		hNP1 = getBoundaryNode(meshTIN, hNUSBound, rightBound, 0, 15);
		appendNode(hNP1);
		
		hNP2 = getBoundaryNode(meshTIN, hNDSBound, rightBound, 0, 15);

		triP = meshTIN->whichTriangle(hNP1);
		prevTriP = NULL; //used to check if mapping is oscillating between triangles
		while(-1)
		{
			checkTriP = mapStreamLine(meshTIN, triP, hNP1, hNP2, rightBound, 0, 15);
			if(triP == checkTriP) break;	
			if(checkTriP == prevTriP) 
				checkTriP = findTriAndNodeWhenVelDirCannot(meshTIN, triP, rightBound, 0, 15);
			prevTriP = triP;
			triP = checkTriP;
		}

	}
	else	//right boundary chosen along a depth contour
	{
		hNP1 = getBoundaryNode(meshTIN, hNUSBound, rightBound, 0, 3);
		appendNode(hNP1);
		
		hNP2 = getBoundaryNode(meshTIN, hNDSBound, rightBound, 0, 3);

		triP = meshTIN->whichTriangle(hNP1);
		prevTriP = NULL; //used to check if mapping is oscillating between triangles
		while(-1)
		{
			checkTriP = mapStreamLine(meshTIN, triP, hNP1, hNP2, rightBound, 0, 3);
			if(triP == checkTriP) break;	
			if(checkTriP == prevTriP) 
				checkTriP = findTriAndNodeWhenVelDirCannot(meshTIN, triP, rightBound, 0, 3);
			prevTriP = triP;
			triP = checkTriP;
		}
	}

	if (leftBoundType == 2)	//left boundary chosen along computation boundary
	{
		
		mapBoundaryStreamLine(meshTIN, hNUSBound, hNDSBound, 1);
		hNP4 = (HabitatNode *)nodeL.lastItem();
	}
	
	else if (leftBoundType == 0) //left boundary chosen along a streamline
	{
		hNP3 = getBoundaryNode(meshTIN, hNDSBound, leftBound, 1, 15);
		appendNode(hNP3);
		wse = (hNP2->getPar(6)+hNP3->getPar(6))/2;
		segNum = bsegL.numberOfItems() + 1;
		boundSegDS = new RivBound(segNum, hNP2, hNP3, 3, wse);
		appendBSeg(boundSegDS);

	
		hNP4 = getBoundaryNode(meshTIN, hNUSBound, leftBound, 1, 15);

		triP = meshTIN->whichTriangle(hNP3);
		prevTriP = NULL;
		while(-1)
		{
			checkTriP = mapStreamLine(meshTIN, triP, hNP3, hNP4, leftBound, 1, 15);
			if(triP == checkTriP) break;
			if(checkTriP == prevTriP) 
				checkTriP = findTriAndNodeWhenVelDirCannot(meshTIN, triP, leftBound, 1, 15);
			prevTriP = triP;
			triP = checkTriP;
		}
	}
	else	//left boundary chosen along a depth contour
	{
		hNP3 = getBoundaryNode(meshTIN, hNDSBound, leftBound, 1, 3);
		appendNode(hNP3);
		wse = (hNP2->getPar(6)+hNP3->getPar(6))/2;
		segNum = bsegL.numberOfItems() + 1;
		boundSegDS = new RivBound(segNum, hNP2, hNP3, 3, wse);
		appendBSeg(boundSegDS);

	
		hNP4 = getBoundaryNode(meshTIN, hNUSBound, leftBound, 1, 3);

		triP = meshTIN->whichTriangle(hNP3);
		prevTriP = NULL;
		while(-1)
		{
			checkTriP = mapStreamLine(meshTIN, triP, hNP3, hNP4, leftBound, 1, 3);
			if(triP == checkTriP) break;
			if(checkTriP == prevTriP) 
				checkTriP = findTriAndNodeWhenVelDirCannot(meshTIN, triP, leftBound, 1, 3);
			prevTriP = triP;
			triP = checkTriP;
		}
	}

	dischargeIntensity = (leftBound - rightBound)/(hNP4->dist(hNP1));
	segNum = bsegL.numberOfItems() + 1;
	boundSegUS = new RivBound(segNum, hNP4, hNP1, 1, dischargeIntensity);
	appendBSeg(boundSegUS);

	resampleBoundary(spacing);  

	return 1;
}

Triangle* HabitatTIN::findTriAndNodeWhenVelDirCannot(TIN *meshTIN,Triangle *triP, double boundaryValue, int flag, int parChoice)
{
	HabitatNode *intersectNodes[3];
	HabitatNode *lastNode;
	RivBound *lastSeg;
	Segment *elementSegs[3];
	RivBound *boundSeg;
	int segNum;
	int i,j;

	//velocity vector has led us astray so last node and segment are not what we want
	//last node is essentially the same as the second last node and last segment is 
	//segment between these nodes
	//so first want to delete these from the lists

	lastNode = (HabitatNode*)nodeL.lastItem();
	nodeL.setCurrentItem(lastNode);
	nodeL.deleteCurrentItem();

	lastSeg = (RivBound*)bsegL.lastItem();
	bsegL.setCurrentItem(lastSeg);
	bsegL.deleteCurrentItem();

	//if not, get the next node and segment along the streamline to append to the list

	for (j = 0; j < 3; j++)
			intersectNodes[j] = NULL;

		for(int k = 0; k < 3 ; k ++)
		{
			int i = k + 1;
			if ( i == 3) i = 0;
			j = k + 2;
			if ( j == 3) j = 0;
			if ( j == 4) j = 1;
			elementSegs[k] = new Segment(1, triP->getNode(i), triP->getNode(j));
		}
		
		//find nodes where the streamline intersects the triangle triP

		for (int i = 0;i < 3; i++)
		{
			double parValue1, parValue2, r;
			
			parValue1 = elementSegs[i]->getNode(0)->getPar(parChoice);
			parValue2 = elementSegs[i]->getNode(1)->getPar(parChoice);

			//provisions for the case when the streamline intersects exactly
			//at a node - simply move the streamline to avoid the intersection at
			//a node
			
			if( parValue1 == boundaryValue || parValue2 == boundaryValue )
			{
				if (flag == 0) boundaryValue = boundaryValue + 0.0001;
				if (flag == 1) boundaryValue = boundaryValue - 0.0001;
			}

			if(((parValue1>boundaryValue) && (parValue2<boundaryValue)) || ((parValue1<boundaryValue) && (parValue2>boundaryValue)))
			{
				if(parValue2 > parValue1)
					r = (boundaryValue - parValue1)/(parValue2 - parValue1);
				else if (parValue2 < parValue1)
					r = 1.0 - (boundaryValue - parValue2)/(parValue1 - parValue2);
				else r = 0.0;

				intersectNodes[i] = (HabitatNode*)physics->makeNewNode(1,100,100);
				elementSegs[i]->locateNodeAt(r, intersectNodes[i]);	
			}	
		}

		//determine which of the two intersection nodes is the next node to appended
		//opposite to what it was when I called mapStreamLine

		for( i = 0; i < 3; i ++)
		{
			j = i+1;
			if (j == 3) j = 0;
			if((intersectNodes[i] != NULL) && (intersectNodes[j] != NULL))
			{
				double xDir = intersectNodes[j]->getXc() - intersectNodes[i]->getXc();
				double yDir = intersectNodes[j]->getYc() - intersectNodes[i]->getYc();
				double qxAvg = (intersectNodes[j]->getPar(4) + intersectNodes[i]->getPar(4))/2;
				double qyAvg = (intersectNodes[j]->getPar(5) + intersectNodes[i]->getPar(5))/2;
				double dotProd = xDir*qxAvg + yDir*qyAvg;
				if(dotProd<0)
				{
					if(flag == 0)
					{
						lastNode = (HabitatNode*)nodeL.lastItem();
						appendNode(intersectNodes[j]);
						if(bsegL.firstItem() == NULL) segNum = 1;
						else segNum = bsegL.numberOfItems() + 1;
						boundSeg = new RivBound(segNum, lastNode, intersectNodes[j], 0, 0);
						appendBSeg(boundSeg);
						delete intersectNodes[i];
						triP = (Triangle*)triP->getAdjt(j);
					}
					else
					{
						lastNode = (HabitatNode*)nodeL.lastItem();
						appendNode(intersectNodes[i]);
						if(bsegL.firstItem() == NULL) segNum = 1;
						else segNum = bsegL.numberOfItems() + 1;
						boundSeg = new RivBound(segNum, lastNode, intersectNodes[i], 0, 0);
						appendBSeg(boundSeg);
						delete intersectNodes[j];
						triP = (Triangle*)triP->getAdjt(i);

					}
				}
				else
				{
					if(flag == 0)
					{	
						lastNode = (HabitatNode*)nodeL.lastItem();
						appendNode(intersectNodes[i]);
						if(bsegL.firstItem() == NULL) segNum = 1;
						else segNum = bsegL.numberOfItems() + 1;
						boundSeg = new RivBound(segNum, lastNode, intersectNodes[i], 0, 0);
						appendBSeg(boundSeg);
						delete intersectNodes[j];
						triP = (Triangle*)triP->getAdjt(i);
					}
					else
					{
						
						lastNode = (HabitatNode*)nodeL.lastItem();
						appendNode(intersectNodes[j]);
						if(bsegL.firstItem() == NULL) segNum = 1;
						else segNum = bsegL.numberOfItems() + 1;
						boundSeg = new RivBound(segNum, lastNode, intersectNodes[j], 0, 0);
						appendBSeg(boundSeg);
						delete intersectNodes[i];
						triP = (Triangle*)triP->getAdjt(j);
					}
				}
				
			}
		}


	for (i = 0; i < 3; i ++)
	{
		delete elementSegs[i];
	}

	return triP;
}

HabitatNode* HabitatTIN::getBoundaryNode(TIN* meshTIN, HabitatNode* hNP, double boundaryValue, int flag, int parChoice)
{
	double dx, dy, d, dxLimit, dyLimit, dLimit;
	double xPath, yPath;
	double v1, v2, r;
	int i,j;

	HabitatNode *hNPath, *hNBoundNode;
	HabitatNode *intersectNodes[3];
	Segment* segPathP;
	Segment *elementSegs[3];
	Triangle *triP, *pathTriP, *checkTriP;

	//create path to follow to get the lower/upper boundary node

	box limits = meshTIN->getLimits();
	dx = hNP->getPar(4);  // uses direction of discharge intensity to get direction of path
	dy = hNP->getPar(5);
	d = sqrt(dx*dx + dy*dy);
	dxLimit = limits.x2 - limits.x1;
	dyLimit = limits.y2 - limits.y1;
	dLimit = sqrt(dxLimit*dxLimit + dyLimit*dyLimit);
	if(flag == 0)		//to get path to lower boundary node
	{
		xPath = (dLimit/d)*dy + hNP->getXc();
		yPath = -(dLimit/d)*dx + hNP->getYc();
	}
	else				// to get path to upper boundary node
	{
		xPath = -(dLimit/d)*dy + hNP->getXc();
		yPath = (dLimit/d)*dx + hNP->getYc();
	}

	hNPath = (HabitatNode *)physics->makeNewNode(1,xPath,yPath);
	segPathP = new Segment(1, hNP, hNPath);

	pathTriP = meshTIN->whichTriangle(hNPath);
	
	// check to see if boundary node is in the current triangle, if not get triangle with boundary node
	
	triP = meshTIN->whichTriangle(hNP);
	while(-1)
	{ 
		checkTriP = findTriWithBoundaryNode(triP, segPathP, boundaryValue, parChoice);
		if(triP != checkTriP) triP = checkTriP;
		else break;
	}

	//with correct triangle, get boundary node
	
		//find nodes where path intersects the correct triangle
	
	for (int j = 0; j < 3; j++)
			intersectNodes[j] = NULL;

	for(int k = 0; k < 3 ; k ++)
	{
		int i = k + 1;
		if ( i == 3) i = 0;
		j = k + 2;
		if ( j == 3) j = 0;
		if ( j == 4) j = 1;
		elementSegs[k] = new Segment(1, triP->getNode(i), triP->getNode(j));
	}
	for(int i = 0; i < 3; i++)
	{
		r = elementSegs[i]->intersectd(segPathP);
		if( (r > 0) && (r < 1.0) )
		{
			intersectNodes[i] = (HabitatNode *)physics->makeNewNode(1,100,100);
			elementSegs[i]->locateNodeAt(r, intersectNodes[i]);
		}	 
	}
	
		//reset path to extend only across the correct triangle
	
	for( i = 0; i < 3; i ++)
	{
		j = i+1;
		if (j == 3) j = 0;
		if((intersectNodes[i] != NULL) && (intersectNodes[j] != NULL))
		{
			segPathP->setNode(0, intersectNodes[i]);
			segPathP->setNode(1, intersectNodes[j]);
		}
	}

		//interpolate the boundary node within the correct triangle

	v1 = segPathP->getNode(0)->getPar(parChoice);
	v2 = segPathP->getNode(1)->getPar(parChoice);

	if(v2 > v1)
		r = (boundaryValue - v1)/(v2 - v1);
	else if (v2 < v1)
		r = 1.0 - (boundaryValue - v2)/(v1 - v2);
	else
		r = 0.0;
	
	hNBoundNode = (HabitatNode *)physics->makeNewNode(1,100,100);
	segPathP->locateNodeAt(r, hNBoundNode);

	for ( i = 0; i < 3; i ++)
	{
		delete elementSegs[i];
		delete intersectNodes[i];
	}

	return hNBoundNode;
}

Triangle* HabitatTIN::findTriWithBoundaryNode(Triangle* triP, Segment* segPathP, double boundaryValue, int parChoice)
{
	Segment *elementSegs[3];
	HabitatNode *intersectNodes[3];
	Node *movePath, *hNPath;
	double r , distance;
	int i,j,k;

	distance = 10e10;
		
	for (int j = 0; j < 3; j++)
		intersectNodes[j] = NULL;

	//find nodes where path intersects the triangle triP
	
		//first create the segments of the triangle

	for(int k = 0; k < 3 ; k ++)
	{
		int i = k + 1;
		if ( i == 3) i = 0;
		j = k + 2;
		if ( j == 3) j = 0;
		if ( j == 4) j = 1;
		elementSegs[k] = new Segment(1, triP->getNode(i), triP->getNode(j));
	}
	
		//check to see if intersection is at a node, if so move the path slightly to avoid this
		//this might not always work...need to think about

	for (int i = 0; i < 3; i++)
	{
		r = elementSegs[i]->intersectd(segPathP);
		if( (r == 0) || (r == 1.0) )
		{
			movePath = segPathP->getNode(0);
			movePath->assign(movePath->getXc() + 0.0001, movePath->getYc() + 0.0001, 100);
			movePath = segPathP->getNode(1);
			movePath->assign(movePath->getXc() + 0.0001, movePath->getYc() + 0.0001, 100);
		}
	}

		//find intersection points - assumes there are 2

	for(i = 0; i < 3; i++)
	{
		r = elementSegs[i]->intersectd(segPathP);
		if( (r > 0) && (r < 1.0) )
		{
			intersectNodes[i] = (HabitatNode *)physics->makeNewNode(1,100,100);
			elementSegs[i]->locateNodeAt(r, intersectNodes[i]);
		}	 
	}
	
	//find the next triangle in along the path 
	//based the distance of the point from hNPath which is the user input point 

	j = -1;
	hNPath = segPathP->getNode(1);
	for( i = 0; i < 3; i ++)	//find which node is closest to hNPath
	{
		if(intersectNodes[i] != NULL)
		{
			if(fabs(intersectNodes[i]->dist(hNPath)) < distance)
			{
				distance = fabs(intersectNodes[i]->dist(hNPath));
				j = i;
			}
		}
	}
	for( i = 0; i < 3; i++)	//find which nodes is farthest from hNPath
							//which will be the node that is not j
							//and not NULL
	{
		if ((intersectNodes[i] != NULL) && ( i != j))
			k = i;
	}

	double diffj, diffk, diffjk;

	diffj = fabs(intersectNodes[j]->getPar(parChoice) - boundaryValue);
	diffk = fabs(intersectNodes[k]->getPar(parChoice) - boundaryValue);
	diffjk = fabs(intersectNodes[j]->getPar(parChoice) - intersectNodes[k]->getPar(parChoice));

	if ((diffj <= diffjk) && (diffk <= diffjk)); //leave triP as is
	else
		triP = (Triangle*)triP->getAdjt(j);

	for ( i = 0; i < 3; i ++)
	{
		delete elementSegs[i];
		delete intersectNodes[i];
	}

	return triP;

}

Triangle* HabitatTIN::mapStreamLine(TIN* meshTIN, Triangle* triP, HabitatNode* hNPStart, HabitatNode* hNPEnd, double boundaryValue, int flag, int parChoice)
{
	HabitatNode *intersectNodes[3];
	HabitatNode *lastNode;
	Segment *elementSegs[3];
	RivBound *boundSeg;
	int segNum;
	int i,j,k;

	//if the current triangle triP contains the end point then append the end point
	//and return the triangle

	if(triP == meshTIN->whichTriangle(hNPEnd))
	{	
		lastNode = (HabitatNode*)nodeL.lastItem();
		appendNode(hNPEnd);
		if(bsegL.firstItem() == NULL) segNum = 1;
		else segNum = bsegL.numberOfItems() + 1;
		boundSeg = new RivBound(segNum, lastNode, hNPEnd, 0, 0);
		appendBSeg(boundSeg);
		return triP;
	}

	//if not, get the next node and segment along the streamline to append to the list

	else
	{
		for (int j = 0; j < 3; j++)
			intersectNodes[j] = NULL;

		for(int k = 0; k < 3 ; k ++)
		{
			int i = k + 1;
			if ( i == 3) i = 0;
			j = k + 2;
			if ( j == 3) j = 0;
			if ( j == 4) j = 1;
			elementSegs[k] = new Segment(1, triP->getNode(i), triP->getNode(j));
		}
		
		//find nodes where the streamline intersects the triangle triP

		for (int i = 0;i < 3; i++)
		{
			double parValue1, parValue2, r;
			
			parValue1 = elementSegs[i]->getNode(0)->getPar(parChoice);
			parValue2 = elementSegs[i]->getNode(1)->getPar(parChoice);

			//provisions for the case when the streamline intersects exactly
			//at a node - simply move the streamline to avoid the intersection at
			//a node
			//not sure how well this works...need to think about
			
			if( parValue1 == boundaryValue || parValue2 == boundaryValue )
			{
				if (flag == 0) boundaryValue = boundaryValue + 0.0001;
				if (flag == 1) boundaryValue = boundaryValue - 0.0001;
			}

			if(((parValue1>boundaryValue) && (parValue2<boundaryValue)) || ((parValue1<boundaryValue) && (parValue2>boundaryValue)))
			{
				if(parValue2 > parValue1)
					r = (boundaryValue - parValue1)/(parValue2 - parValue1);
				else if (parValue2 < parValue1)
					r = 1.0 - (boundaryValue - parValue2)/(parValue1 - parValue2);
				else r = 0.0;

				intersectNodes[i] = (HabitatNode*)physics->makeNewNode(1,100,100);
				elementSegs[i]->locateNodeAt(r, intersectNodes[i]);	
			}	
		}

		//determine which of the two intersection nodes is the next node to appended

		for( i = 0; i < 3; i ++)
		{
			j = i+1;
			if (j == 3) j = 0;
			if((intersectNodes[i] != NULL) && (intersectNodes[j] != NULL))
			{
				
				if(triP == meshTIN->whichTriangle(hNPStart))
				{
					double distance1, distance2;
					distance1 = intersectNodes[i]->dist(hNPEnd);
					distance2 = intersectNodes[j]->dist(hNPEnd);			
					if (distance1 > distance2)
					{	
						lastNode = (HabitatNode*)nodeL.lastItem();
						appendNode(intersectNodes[j]);
						if(bsegL.firstItem() == NULL) segNum = 1;
						else segNum = bsegL.numberOfItems() + 1;
						boundSeg = new RivBound(segNum, lastNode, intersectNodes[j], 0, 0);
						appendBSeg(boundSeg);
						delete intersectNodes[i];
						triP = (Triangle*)triP->getAdjt(j);
					}
					else
					{
						lastNode = (HabitatNode*)nodeL.lastItem();
						appendNode(intersectNodes[i]);
						if(bsegL.firstItem() == NULL) segNum = 1;
						else segNum = bsegL.numberOfItems() + 1;
						boundSeg = new RivBound(segNum, lastNode, intersectNodes[i], 0, 0);
						appendBSeg(boundSeg);
						delete intersectNodes[j];
						triP = (Triangle*)triP->getAdjt(i);
					}
				}

				else
				{
					double xDir = intersectNodes[j]->getXc() - intersectNodes[i]->getXc();
					double yDir = intersectNodes[j]->getYc() - intersectNodes[i]->getYc();
					double qxAvg = (intersectNodes[j]->getPar(4) + intersectNodes[i]->getPar(4))/2;
					double qyAvg = (intersectNodes[j]->getPar(5) + intersectNodes[i]->getPar(5))/2;
					double dotProd = xDir*qxAvg + yDir*qyAvg;
					
					if(dotProd>0)
					{
						if(flag == 0)
						{
							lastNode = (HabitatNode*)nodeL.lastItem();
							appendNode(intersectNodes[j]);
							if(bsegL.firstItem() == NULL) segNum = 1;
							else segNum = bsegL.numberOfItems() + 1;
							boundSeg = new RivBound(segNum, lastNode, intersectNodes[j], 0, 0);
							appendBSeg(boundSeg);
							delete intersectNodes[i];
							triP = (Triangle*)triP->getAdjt(j);
						}
						else
						{
							lastNode = (HabitatNode*)nodeL.lastItem();
							appendNode(intersectNodes[i]);
							if(bsegL.firstItem() == NULL) segNum = 1;
							else segNum = bsegL.numberOfItems() + 1;
							boundSeg = new RivBound(segNum, lastNode, intersectNodes[i], 0, 0);
							appendBSeg(boundSeg);
							delete intersectNodes[j];
							triP = (Triangle*)triP->getAdjt(i);

						}
					}
					else
					{
						if(flag == 0)
						{	
							lastNode = (HabitatNode*)nodeL.lastItem();
							appendNode(intersectNodes[i]);
							if(bsegL.firstItem() == NULL) segNum = 1;
							else segNum = bsegL.numberOfItems() + 1;
							boundSeg = new RivBound(segNum, lastNode, intersectNodes[i], 0, 0);
							appendBSeg(boundSeg);
							delete intersectNodes[j];
							triP = (Triangle*)triP->getAdjt(i);
						}
						else
						{
							
							lastNode = (HabitatNode*)nodeL.lastItem();
							appendNode(intersectNodes[j]);
							if(bsegL.firstItem() == NULL) segNum = 1;
							else segNum = bsegL.numberOfItems() + 1;
							boundSeg = new RivBound(segNum, lastNode, intersectNodes[j], 0, 0);
							appendBSeg(boundSeg);
							delete intersectNodes[i];
							triP = (Triangle*)triP->getAdjt(j);
						}
					}
				}		
			}
		}

	}

	for (int i = 0; i < 3; i ++)
	{
		delete elementSegs[i];
	}

	return triP;
}





void HabitatTIN::resampleBoundary(double spacing)
{
	HabitatNode *hNP1, *hNP2, *hNP3, *hNP4, *tempHNP1, *tempHNP2, *newHNP, *newStartHNP, *lastNodeP, *oldNodeP;
	RivBound *boundSeg;
	double lowerArc = 0, upperArc = 0, noOfNodesExact;
	Segment *tempSeg;
	double remainder , segLength, incrementLength, r, wse, dischargeIntensity;
	int noOfNodes, noOfNodesAppended, segNum;
	int i,j,k;

	//get the corner nodes of the boundary

	hNP1 = (HabitatNode*)nodeL.firstItem();
	hNP4 = (HabitatNode*)nodeL.lastItem();

	for(int i = 1; i <=bsegL.numberOfItems(); i++)
	{
		boundSeg = (RivBound*)bsegL.i(i);
		if (boundSeg->getBcCode() == 3)
		{
			hNP2 = (HabitatNode*)boundSeg->getNode(0);
			hNP3 = (HabitatNode*)boundSeg->getNode(1);
		}

	}

	//delete the boundary segments in the list, will create a new ones

	while(bsegL.numberOfItems()>0)
	{
		boundSeg = (RivBound*)bsegL.pop();
		delete boundSeg;		
	}
	
	//find the length of the upper and lower streamlines

	tempHNP1 = hNP1;
	nodeL.setCurrentItem(tempHNP1);
	tempHNP2 = (HabitatNode*)nodeL.nextItem();

	while(tempHNP1 != hNP2)
	{
		lowerArc += tempHNP1->dist(tempHNP2);
		tempHNP1 = tempHNP2;
		nodeL.setCurrentItem(tempHNP1);
		tempHNP2 = (HabitatNode*)nodeL.nextItem();	
	}

	tempHNP1 = hNP3;
	nodeL.setCurrentItem(tempHNP1);
	tempHNP2 = (HabitatNode*)nodeL.nextItem();

	while(tempHNP1 != hNP4)
	{
		upperArc += tempHNP1->dist(tempHNP2);
		tempHNP1 = tempHNP2;
		nodeL.setCurrentItem(tempHNP1);
		tempHNP2 = (HabitatNode*)nodeL.nextItem();
	}

	//finding the exact spacing between the nodes for the lower streamline
	//based on the user defined spacing

	noOfNodesExact = lowerArc/spacing;
	noOfNodes = int(noOfNodesExact);
	spacing = lowerArc/noOfNodes;

	//setting the first node at the beginning of the lower streamline 
	//identical to hNP1 except that it is fixed
	//and initiating some variables for the while loop

	tempHNP1 = hNP1;
	nodeL.setCurrentItem(tempHNP1);
	tempHNP2 = (HabitatNode*)nodeL.nextItem();
	tempSeg = new Segment(1, tempHNP1, tempHNP2);
	newStartHNP = (HabitatNode *)physics->makeNewNode(1,100,100);
	tempSeg->locateNodeAt(0, newStartHNP);
	newStartHNP->setFixed();
	appendNode(newStartHNP);

	incrementLength = spacing;
	remainder = 0;
	noOfNodesAppended = 0;

	//looping through all of the nodes in the lower streamline, appending 
	//new nodes at the new spacing and boundary segments between the nodes

	while(tempHNP1 != hNP2)
	{
		tempSeg->setNode(0, tempHNP1);
		tempSeg->setNode(1, tempHNP2);
		segLength= tempSeg->length();

		while ( segLength >= incrementLength)
		{
			if(noOfNodesAppended == (noOfNodes -1)) break;
			r = incrementLength / segLength;
			newHNP = (HabitatNode *)physics->makeNewNode(1,100,100);
			tempSeg->locateNodeAt(r, newHNP);
			newHNP->setFixed();
			lastNodeP = (HabitatNode*)nodeL.lastItem();
			appendNode(newHNP);
			noOfNodesAppended += 1;
			if(bsegL.firstItem() == NULL) segNum = 1;
			else segNum = bsegL.numberOfItems() + 1;
			boundSeg = new RivBound(segNum, lastNodeP, newHNP, 0, 0);
			appendBSeg(boundSeg);
			incrementLength += spacing;
			
		}
		remainder = incrementLength - segLength;
		incrementLength = remainder;
		
		tempHNP1 = tempHNP2;
		nodeL.setCurrentItem(tempHNP1);
		tempHNP2 = (HabitatNode*)nodeL.nextItem();	
	}
		
	//appending a node at the end of the lower streamline - same as hNP2 but fixed
	//appending the last boundary segment in the streamline

	lastNodeP = newHNP;
	tempSeg = new Segment(1, tempHNP1, tempHNP2);
	newHNP = (HabitatNode *)physics->makeNewNode(1,100,100);
	tempSeg->locateNodeAt(0, newHNP);  //may want to delete tempSeg when done with it
	newHNP->setFixed();
	appendNode(newHNP);
	segNum = bsegL.numberOfItems() + 1;
	boundSeg = new RivBound(segNum, lastNodeP, newHNP, 0, 0);
	appendBSeg(boundSeg);

	//appending the first node in the upper streamline (starting at the downstream end)
	//same as hNP3 and appending the outflow boundary segment

	lastNodeP = newHNP;
	newHNP = (HabitatNode *)physics->makeNewNode(1,100,100);
	tempSeg->locateNodeAt(1, newHNP);
	newHNP->setFixed();
	appendNode(newHNP);
	segNum = bsegL.numberOfItems() + 1;
	wse = (lastNodeP->getPar(6)+newHNP->getPar(6))/2;
	boundSeg = new RivBound(segNum, lastNodeP, newHNP, 3, wse);
	appendBSeg(boundSeg);

	//finding the exact spacing between the nodes for the upper streamline
	//based on the user defined spacing

	noOfNodesExact = upperArc/spacing;
	noOfNodes = int(noOfNodesExact);
	spacing = upperArc/noOfNodes;

	//initiating some variables for the while loop
	
	tempHNP1 = hNP3;
	nodeL.setCurrentItem(tempHNP1);
	tempHNP2 = (HabitatNode*)nodeL.nextItem();
	tempSeg = new Segment(1, tempHNP1, tempHNP2);
	
	incrementLength = spacing;
	remainder = 0;
	noOfNodesAppended = 0;
	
	while(tempHNP1 != hNP4)
	{
		tempSeg->setNode(0, tempHNP1);
		tempSeg->setNode(1, tempHNP2);
		segLength= tempSeg->length();

		while ( segLength >= incrementLength)
		{
			if(noOfNodesAppended == (noOfNodes -1)) break;
			r = incrementLength / segLength;
			newHNP = (HabitatNode *)physics->makeNewNode(1,100,100);
			tempSeg->locateNodeAt(r, newHNP);
			newHNP->setFixed();
			lastNodeP = (HabitatNode*)nodeL.lastItem();
			appendNode(newHNP);
			noOfNodesAppended += 1;
			if(bsegL.firstItem() == NULL) segNum = 1;
			else segNum = bsegL.numberOfItems() + 1;
			boundSeg = new RivBound(segNum, lastNodeP, newHNP, 0, 0);
			appendBSeg(boundSeg);	
			incrementLength += spacing;
	
		}
		remainder = incrementLength - segLength;
		incrementLength = remainder;
		
		tempHNP1 = tempHNP2;
		nodeL.setCurrentItem(tempHNP1);
		tempHNP2 = (HabitatNode*)nodeL.nextItem();
	}

	//appending a node at the end of the upper streamline - same as hNP4 but fixed
	//appending the last boundary segment in the streamline

	lastNodeP = newHNP;
	tempSeg = new Segment(1, tempHNP1, tempHNP2);
	newHNP = (HabitatNode *)physics->makeNewNode(1,100,100);
	tempSeg->locateNodeAt(0, newHNP);
	newHNP->setFixed();
	appendNode(newHNP);
	segNum = bsegL.numberOfItems() + 1;
	boundSeg = new RivBound(segNum, lastNodeP, newHNP, 0, 0);
	appendBSeg(boundSeg);

	//appending the outflow boundary segment

	dischargeIntensity = (newHNP->getPar(15) - newStartHNP->getPar(15))/(newHNP->dist(newStartHNP));
	segNum = bsegL.numberOfItems() + 1;
	boundSeg = new RivBound(segNum, newHNP, newStartHNP, 1, dischargeIntensity);
	appendBSeg(boundSeg);

	oldNodeP = (HabitatNode*)nodeL.firstItem();

	while(oldNodeP != NULL)
	{
		if(oldNodeP->getFixed() != fixednode)
		{
			nodeL.deleteItem(oldNodeP);
			oldNodeP = (HabitatNode*)nodeL.firstItem();
		}
		else
			oldNodeP = (HabitatNode*)nodeL.nextItem();
	}
}

void HabitatTIN::mapBoundaryStreamLine(TIN* meshTIN, HabitatNode* hNUSBound, HabitatNode* hNDSBound, int flag)
{
	double dx, dy, d, dxLimit, dyLimit, dLimit, r1, r2, wse, distance;
	double xPath, yPath;
	Segment *USsegPathP, *DSsegPathP;
	RivBound *bsP, *USbsP, *DSbsP, *nextbsP, *newbsP;
	HabitatNode *hNUSPath, *hNDSPath, *hNstartP, *hNendP, *nP, *lastNodeP;
	int segNum;

	box limits = meshTIN->getLimits();
	dxLimit = limits.x2 - limits.x1;
	dyLimit = limits.y2 - limits.y1;
	dLimit = sqrt(dxLimit*dxLimit + dyLimit*dyLimit);
	
	//get path for upstream end

	dx = hNUSBound->getPar(4);  // uses direction of discharge intensity to get direction of path
	dy = hNUSBound->getPar(5);
	d = sqrt(dx*dx + dy*dy);
	
	if(flag == 0)		//to get path to lower boundary node
	{
		xPath = (dLimit/d)*dy + hNUSBound->getXc();
		yPath = -(dLimit/d)*dx + hNUSBound->getYc();
	}
	else				// to get path to upper boundary node
	{
		xPath = -(dLimit/d)*dy + hNUSBound->getXc();
		yPath = (dLimit/d)*dx + hNUSBound->getYc();
	}

	hNUSPath = (HabitatNode *)physics->makeNewNode(1,xPath,yPath);
	USsegPathP = new Segment(1, hNUSBound, hNUSPath);

	//get path for downstream end

	dx = hNDSBound->getPar(4);  // uses direction of discharge intensity to get direction of path
	dy = hNDSBound->getPar(5);
	d = sqrt(dx*dx + dy*dy);
	
	if(flag == 0)		//to get path to lower boundary node
	{
		xPath = (dLimit/d)*dy + hNDSBound->getXc();
		yPath = -(dLimit/d)*dx + hNDSBound->getYc();
	}
	else				// to get path to upper boundary node
	{
		xPath = -(dLimit/d)*dy + hNDSBound->getXc();
		yPath = (dLimit/d)*dx + hNDSBound->getYc();
	}

	hNDSPath = (HabitatNode *)physics->makeNewNode(1,xPath,yPath);
	DSsegPathP = new Segment(1, hNDSBound, hNDSPath);

	//find the boundary element that is intersected by the upstream path
	//the one that is first intersected (closest)
	
	bsP = (RivBound *)meshTIN->firstSeg();
	USbsP = bsP;

	distance = 10e10;

	while (bsP != NULL)
	{
		r1 = bsP->intersectd(USsegPathP);
		r2 = USsegPathP->intersectd(bsP);

		if ( (r1 >= 0) && (r1 <= 1.0) && (r2 >= 0) && (r2 <= 1.0))
			if((bsP->getNode(0)->dist(USsegPathP->getNode(0))) < distance)
			{
				distance = bsP->getNode(0)->dist(USsegPathP->getNode(0));
				USbsP = bsP;
			}
		bsP = (RivBound *)meshTIN->nextSeg();

	}

	//find the boundary element that is intersected by the downstream path
	//the one that is first intersected (closest)

	bsP = (RivBound *)meshTIN->firstSeg();
	DSbsP = bsP;

	distance = 10e10;

	while (bsP != NULL)
	{
		r1 = bsP->intersectd(DSsegPathP);
		r2 = DSsegPathP->intersectd(bsP);

		if ( (r1 >= 0) && (r1 <= 1.0) && (r2 >= 0) && (r2 <= 1.0))
			if((bsP->getNode(0)->dist(DSsegPathP->getNode(0))) < distance)
			{
				distance = bsP->getNode(0)->dist(DSsegPathP->getNode(0));
				DSbsP = bsP;
			}
		bsP = (RivBound *)meshTIN->nextSeg();

	}
	
	//getting the starting and end nodes in the psi = 0 streamline
	
	if(flag == 0)
	{
		hNstartP = (HabitatNode*)physics->makeNewNode(1,100,100);
		r1 = USbsP->intersectd(USsegPathP);
		USbsP->locateNodeAt(r1, hNstartP);
		hNendP = (HabitatNode*)physics->makeNewNode(1,100,100);
		r1 = DSbsP->intersectd(DSsegPathP);
		DSbsP->locateNodeAt(r1, hNendP);
	}
	
	//getting the starting and end nodes int the psi = Qtotal streamline
	
	else
	{
		hNstartP = (HabitatNode*)physics->makeNewNode(1,100,100);
		r1 = DSbsP->intersectd(DSsegPathP);
		DSbsP->locateNodeAt(r1, hNstartP);
		hNendP = (HabitatNode*)physics->makeNewNode(1,100,100);
		r1 = USbsP->intersectd(USsegPathP);
		USbsP->locateNodeAt(r1, hNendP);	
	}

	//mapping the streamline that is along the right boundary
	
	if(flag == 0)
	{
		appendNode(hNstartP);
		bsP = USbsP;
	
		while (bsP != DSbsP)
		{
			lastNodeP = (HabitatNode*)nodeL.lastItem();
			nP = (HabitatNode*)physics->makeNewNode(1,100,100);
			bsP->locateNodeAt(1, nP);
			appendNode(nP);
			if(bsegL.firstItem() == NULL) segNum = 1;
			else segNum = bsegL.numberOfItems() + 1;
			newbsP = new RivBound(segNum, lastNodeP, nP, 0, 0);
			appendBSeg(newbsP);
			nextbsP = (RivBound *)meshTIN->nextSeg();
			if(nextbsP == NULL)
				nextbsP = (RivBound *)meshTIN->firstSeg();
			if(nextbsP->getNode(0) == bsP->getNode(1))
				bsP = nextbsP;
			else
			{
				while(nextbsP->getNode(0) != bsP->getNode(1))
				{
					nextbsP = (RivBound *)meshTIN->nextSeg();
					if(nextbsP == NULL)
						nextbsP = (RivBound *)meshTIN->firstSeg();
				}
				bsP = nextbsP;
			}
		}
	}
	
	//mapping the streamline that is along left boundary
	
	else
	{
		lastNodeP = (HabitatNode*)nodeL.lastItem();
		appendNode(hNstartP);
		wse = (lastNodeP->getPar(6)+ hNstartP->getPar(6))/2;
		if(bsegL.firstItem() == NULL) segNum = 1;
		else segNum = bsegL.numberOfItems() + 1;
		bsP = new RivBound(segNum, lastNodeP, hNstartP, 3, wse);
		appendBSeg(bsP);

		bsP = DSbsP;

		while (bsP != USbsP)
		{
			lastNodeP = (HabitatNode*)nodeL.lastItem();
			nP = (HabitatNode*)physics->makeNewNode(1,100,100);
			bsP->locateNodeAt(1, nP);
			appendNode(nP);
			if(bsegL.firstItem() == NULL) segNum = 1;
			else segNum = bsegL.numberOfItems() + 1;
			newbsP = new RivBound(segNum, lastNodeP, nP, 0, 0);
			appendBSeg(newbsP);
			nextbsP = (RivBound *)meshTIN->nextSeg();
			if(nextbsP == NULL)
				nextbsP = (RivBound *)meshTIN->firstSeg();
			if(nextbsP->getNode(0) == bsP->getNode(1))
				bsP = nextbsP;
			else
			{
				while(nextbsP->getNode(0) != bsP->getNode(1))
				{
					nextbsP = (RivBound *)meshTIN->nextSeg();
					if(nextbsP == NULL)
						nextbsP = (RivBound *)meshTIN->firstSeg();
				}
				bsP = nextbsP;
			}
		}
	}
	
	//appending last node and boundary element in streamline
	
	lastNodeP = (HabitatNode*)nodeL.lastItem();
	appendNode(hNendP);
	if(bsegL.firstItem() == NULL) segNum = 1;
	else segNum = bsegL.numberOfItems() + 1;
	bsP = new RivBound(segNum, lastNodeP, hNendP, 0, 0);
	appendBSeg(bsP);
}

int HabitatTIN::autoRefineHabitatRegion(TIN* boundTIN, TIN &dataTIN, double lowerlimit, double upperlimit, int parNum, int valueOrChange)
//	Refines a HabitatTIN by placing a new node in the centre of each triangle
//	that adheres to a particular critieria 
//	Requires a previously triangulated dataTIN (bed data)
//	which masks inactive areas and provides interpolation data.
//	Resulting set of nodes are in tnodeL which can be later accepted into nodeL or rejected.
//	The number of nodes created are returned.

{
	int nGoodNodes = 0;
	int nodeNum = getNextN();
	Triangle *tP, *tBP, *eTP, *eltP;
	Element *elP;
	ShallowNode cNode;
	HabitatNode *nP;
	double value0, value1, value2, average;
	double dz;

	eTP = firstTri();  //meshP
	while(eTP != NULL){
		if(eTP->getStatus() == active){
			value0 = eTP->getNode(0)->getPar(parNum);
			value1 = eTP->getNode(1)->getPar(parNum);
			value2 = eTP->getNode(2)->getPar(parNum);

			switch (valueOrChange) {

				case 0: 
				{
					average = (value0 + value1 + value2)/3.;
					if( (lowerlimit < average) && (average < upperlimit) )
					{
						eTP->setRefineOn();
						for(int i=0;i<3;i++){
							elP = eTP->getAdjt(i);
							if(elP != NULL && elP->numNodes() == 3){
								eltP = (Triangle *) elP;
								eltP->setRefineOn();
							}
						}
					}
					break;
				}

				case 1: //the change in
				{			
					if( ((lowerlimit<fabs(value0-value1))&&(fabs(value0-value1)<upperlimit))
						|| ((lowerlimit<fabs(value1-value2))&&(fabs(value1-value2)<upperlimit))
						|| ((lowerlimit<fabs(value2-value0))&&(fabs(value2-value0)<upperlimit)) ) 
					{
						eTP->setRefineOn();
						for(int i=0;i<3;i++){
							elP = eTP->getAdjt(i);
							if(elP != NULL && elP->numNodes() == 3){
								eltP = (Triangle *) elP;
								eltP->setRefineOn();
							}
						}
					}
					break;
				}

				case 2: //the change in
				{			
					if(    ((lowerlimit<value0) && (value0<upperlimit))
						|| ((lowerlimit<value1) && (value1<upperlimit))
						|| ((lowerlimit<value2) && (value2<upperlimit)) ) 
					{
						eTP->setRefineOn();
					}
					break;
				}
			}
		}
		eTP = nextTri();
	}
	eTP = firstTri();  //meshP
	while(eTP != NULL){
		if(eTP->isRefineOn() == 1){
			nP = (HabitatNode *)physics->makeNewNode(nodeNum);
			eTP->locateNodeAtCenter(nP); //interpolates a new node based on mesh
			cNode.assign(nP->getXc(),nP->getYc(), nP->getZc()); //assigns the shallownode the same x and y
			if( (tP = dataTIN.whichTriangle(nP)) != NULL){
				tBP = eTP;
				if(boundTIN != NULL){
					tBP = boundTIN->whichTriangle(nP);
				}
				if(tBP != NULL){
					if(tBP->getStatus() == active){
						tP->locateNodeAt(&cNode);
						nP->setKs(cNode.getKs());
						dz = cNode.getZc() - nP->getZc();
						nP->setZc(nP->getZc()+ dz);
						nP->setD(nP->getD() - dz);
						nGoodNodes += 1;
						nodeNum += 1;
						tnodeL.appendItem(nP);
					}
					else
						delete nP;
				}
				else
					delete nP;
			}
			else
				delete nP;
		}
		eTP->setRefineOff();
		eTP = nextTri();
	}
	return nGoodNodes;
}

int HabitatTIN::deletePrimaryMeshNode(Node * nP)
{
	ItemList tempfsegL;	//temporary brealine list
	Segment *segP;
	Node *bN, *eN;

	if((nP->getFixed() == fixednode) || (nP->getFixed() == sliding))
	{
		segP = firstfSeg();
		
		//remove all breakline segments that include node nP
		//and put on a temporary list
		
		while(segP != NULL)
		{
			if((segP->getNode(0) == nP) || (segP->getNode(1) == nP))
			{
				fsegL.deleteCurrentItem();
				tempfsegL.appendItem(segP);
				segP = (Segment*) fsegL.currentItem();
			}
			else
				segP = nextfSeg();
		}

		//fix all of the nodes that are connected to node nP
		//by breakline segments
		//then delete the breakline segments

		segP = (Segment*) tempfsegL.firstItem();
		while (segP != NULL)
		{
			bN = segP->getNode(0);
			eN = segP->getNode(1);
			if (bN == nP) eN->setFixed();
			else bN->setFixed();
			segP = (Segment*) tempfsegL.nextItem();
		}
		tempfsegL.emptyList();
	
		//Set node nP to sliding so that it can be deleted

		nP->setSliding();
	}
	if((nP->getFixed() == floating) || (nP->getFixed() == sliding))
	{
		nP->setDeleted();
		return 1;
	}
	else
		return 0;
}

HabitatNode* HabitatTIN::removeNode(Node*nP)
{
	HabitatNode *currentNode;
	currentNode = (HabitatNode*) nodeL.deleteItem(nP);
	return currentNode;
}

void HabitatTIN::initSolutionVector(double *X)
{
	if (FEMProb == psiSolution)
	{
		HabitatNode *hNP;
		int i = 0;
		hNP = (HabitatNode*)firstNode();

		while(hNP != NULL)
		{
			X[i] = hNP->getPsi();
			i++;
			hNP = (HabitatNode*)nextNode();	
		}
	}
	if (FEMProb == phiSolution)
	{
		HabitatNode *hNP;
		int i = 0;
		hNP = (HabitatNode*)firstNode();

		while (hNP != NULL)
		{
			X[i] = hNP->getPar(6);
			i++;
			hNP = (HabitatNode*)nextNode();
		}

	}
}

void HabitatTIN::updateSolToNodes(double *X)
{
	if (FEMProb == psiSolution)
	{
		HabitatNode *hNP;
		int i = 0;
		hNP = (HabitatNode*)firstNode();

		while (hNP != NULL)
		{
			hNP->setPsi(X[i]);
			i++;
			hNP = (HabitatNode*)nextNode();
		}
	}
	if(FEMProb == phiSolution)
	{
		HabitatNode *hNP;
		int i = 0;
		hNP = (HabitatNode*)firstNode();

		while (hNP != NULL)

		{
			hNP->setFlow(X[i]-hNP->getPar(1),0,0);
			hNP->setVandF();
			i++;
			hNP = (HabitatNode*)nextNode();
		}
	}
}

void HabitatTIN::destroyMatrices()
{

	ElementMtx *elMP, *BelMP;
	
	//delete the element matrices and vectors once solution is obtained

	//delete the triangular element matrices

	elMP = (ElementMtx *) ElMtxL.firstItem();
	while(elMP != NULL)
	{
		ElMtxL.deleteCurrentItem();
		delete elMP;
		elMP = (ElementMtx *) ElMtxL.firstItem();
	}

	//delete the linear element matrices (boundary elements)

	BelMP = (ElementMtx *) BElMtxL.firstItem();
	while(BelMP != NULL)
	{
		BElMtxL.deleteCurrentItem();
		delete BelMP;
		BelMP = (ElementMtx *) BElMtxL.firstItem();
	}
}

void HabitatTIN::buildPsiMatrices()
{

	Triangle *triP;
	RivBound *bSegP;
	PsiElementMtx *elMP, *BelMP;

	//create the triangular element matrices
	
	triP = firstTri();

	while (triP != NULL)
	{
		if (triP->getStatus() == active)
		{
			elMP = new PsiElementMtx(triP->getn(), triP, 3);
			ElMtxL.appendItem(elMP);
		}
		triP = nextTri();
	}

	//create the linear element matrices (boundary elements)

	bSegP = (RivBound*)firstSeg();

	while (bSegP != NULL)
	{
		if(bSegP->getBcCode() == 0)
		{
			BelMP = new PsiElementMtx(bSegP->getn(), bSegP, 2);
			BElMtxL.appendItem(BelMP);
		}
		bSegP = (RivBound*)nextSeg();
	}

}

void HabitatTIN::calculatePsi(HabitatNode *startnP)

{
	setAllPsi(startnP);
	buildPsiMatrices();
	FEMProb = psiSolution;
	Solver* psiSolver = new Solver(this);
	psiSolver->preconditionedCGs( &ElMtxL, &BElMtxL, 1);
	delete psiSolver;
	destroyMatrices();
}

void HabitatTIN::buildPhiMatrices(double wsElevIn)
{

	Triangle *triP;
	RivBound *bSegP;
	PhiElementMtx *elMP, *BelMP;

	//create the triangular element matrices
	
	triP = firstTri();

	while (triP != NULL)
	{
		if (triP->getStatus() == active)
		{
			elMP = new PhiElementMtx(triP->getn(), triP, 3, wsElevIn);
			ElMtxL.appendItem(elMP);
		}
		triP = nextTri();
	}

	//create the linear element matrices (boundary elements)

	bSegP = (RivBound*)firstSeg();

	while (bSegP != NULL)
	{
		if(bSegP->getBcCode() == 1 || bSegP->getBcCode() == 3 || bSegP->getBcCode() == 5)
		{
			BelMP = new PhiElementMtx(bSegP->getn(), bSegP, 2, wsElevIn);
			BElMtxL.appendItem(BelMP);
		}
		bSegP = (RivBound*)nextSeg();
	}

}

void HabitatTIN::calculatePhi(double wsElevIn)

{
	buildPhiMatrices(wsElevIn);
	FEMProb = phiSolution;
	Solver* phiSolver = new Solver(this);
	phiSolver->preconditionedCGs( &ElMtxL, &BElMtxL, 1);
	delete phiSolver;
	destroyMatrices();
}


