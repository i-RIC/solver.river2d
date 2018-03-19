
//		Shallow.cpp

#include "Shallow.h"
#include "math.h"
#include <string.h>

ShallowNode::ShallowNode(int nm, double xc, double yc, double zc, double rf,
									double depth, double xDis, double yDis)
				:Node(nm,xc,yc,zc)
{
//	setn(nm);
//	assign(xc, yc, zc);
	ks = rf;
	d = depth;
	qx = xDis;
	qy = yDis;
}

ShallowNode::ShallowNode(istream& is)
				:Node(is)
{
	int num = -990;
//	double xcoord=0.0, ycoord=0.0, zcoord=0.0;
	double rough = 0.05;
	while(is) {
/*		if(!(is >> num)){
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
		if(!(is >> zcoord)){
			num = -994;
			break;
		}
*/		if(!(is >> rough)){
			num = -995;
         setn(num);
			break;
		}
		break;
	}
/*
	is >> num;
	is >> xcoord;
	is >> ycoord;
	is >> zcoord;
	is >> rough;
*/
//	setn(num);
//	assign(xcoord, ycoord, zcoord);
	ks = rough;
	d = 1.0;
	qx = 0.0;
	qy = 0.0;
}

double ShallowNode::getPar(int i)
{
	switch(i){
		case 0 :
			return z;
		case 1 :
			return z;
		case 2 :
			return ks;
		case 3 :
			return d;
		case 4 :
			return qx;
		case 5 :
			return qy;
		case 6 :
			return z+d;     // water surface elevation
		case 7 :
			return sqrt(qx*qx+qy*qy)/d;  // Velocity
		case 8 :
			return sqrt(qx*qx+qy*qy)/(d*sqrt(fabs(d)*G));  //  Froude Number
	}
	return 0.0;
}

double ShallowNode::getVar(int i)
{
	switch(i){
		case 1 :
			return d;
		case 2 :
			return qx;
		case 3 :
			return qy;
	}
	return 0.0;
}

void ShallowNode::getParName(int i, char* parName)
{
	switch(i){
		case 0 :
			strcpy(parName,"Bed Elevation");
			break;
		case 1 :
			strcpy(parName,"Bed Elevation");
			break;
		case 2 :
			strcpy(parName,"Bed Roughness");
			break;
		case 3 :
			strcpy(parName,"Depth");
			break;
		case 4 :
			strcpy(parName,"qx");
			break;
		case 5 :
			strcpy(parName,"qy");
			break;
		case 6 :
			strcpy(parName,"Water Surface Elev");
			break;
		case 7 :
			strcpy(parName,"Velocity");
			break;
		case 8 :
			strcpy(parName,"Froude #");
			break;
	}
}

void ShallowNode::interp(int n, Node** nPtrs, double* wts)
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
	ks = 0.0;
	for(int i=0;i<n;i++) {
		if(nPtrs[i] != NULL) {
			x += wts[i] * nPtrs[i]->getXo();
			y += wts[i] * nPtrs[i]->getYo();
			z += wts[i] * nPtrs[i]->getZc();
			ks += wts[i] * nPtrs[i]->getPar(2);
		}
	}
	saveLocation();
}

RivBound::RivBound(int num, Node* nP1, Node* nP2,
						 int c, double v)
		: Segment(num,nP1,nP2)
{
	bcCode = c;
	bcValue = v;
	flowBound = NULL;
}

ostream& operator << (ostream& os, RivBound* seg)
{
	os.width(5);
	os 	<< seg->getn() ;
	
	if(seg != NULL){
	
		for(int i=0;i<seg->numNodes();i++){
			os.width(10);
			if(seg->getNode(i) != NULL)
				os 	<<  (seg->getNode(i))->getn() ;
			else
				os << -1;
		}
		os.width(10);
		os << seg->getBcCode() << " ";
		os.width(10);
		os << seg->getBcValue();
		os << "\n";
	}
	return os;
}

Node* Shallow::makeNewNode(int n, double x, double y, double z)
{
	return new ShallowNode(n,x,y,z,0.1);
}

Segment* Shallow::makeNewSegment(int n, Node *nP1, Node *nP2, Segment *segP)
{
	RivBound *bSeg, *obSeg;

	bSeg = new RivBound(n,nP1,nP2);
	if(segP != NULL){
		obSeg = (RivBound *) segP;
		bSeg->setBcCode(obSeg->getBcCode());
		bSeg->setBcValue(obSeg->getBcValue());
	}

	return bSeg;
}

Node* Shallow::readInNode(istream& is)
{
	return new ShallowNode(is);
}

Segment* Shallow::readInSeg(int n, Node* nP1, Node* nP2, istream& is)
{
	RivBound *nBSeg;
	int	code;
	double value;

	nBSeg = new RivBound(n,nP1,nP2);

	if ((is >> code))
		nBSeg->setBcCode(code);
	else
		nBSeg->setn(-890);

	if ((is >> value))
		nBSeg->setBcValue(value);
	else
		nBSeg->setn(-891);

	return(nBSeg);
}

void Shallow::writeBSeg(ostream& os, Segment *segP)
{
	os << (RivBound *) segP;
}

