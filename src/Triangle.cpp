
//		Triangle.cpp

//		Triangle Class method definitions

#include <math.h>
#include "Triangle.h"
#include <iostream>
using namespace std;

Triangle::Triangle(int name, Node* firstNodeP, Node* secondNodeP, Node* thirdNodeP)
		: Element(name)
{	
	nodes[0]=firstNodeP;
	nodes[1]=secondNodeP;
	nodes[2]=thirdNodeP;
	adj[0] = NULL;
	adj[1] = NULL;
	adj[2] = NULL;
	edges[0] = NULL;
	edges[1] = NULL;
	edges[2] = NULL;
	toBeRefined = 0;
}

int Triangle::nodeIndex(Node* nP)
{
	for(int i=0;i<numNodes();i++){
		if(nP == nodes[i])
			return i;
	}
	
	return -1;
}

double Triangle::area()
{
	if((nodes[0] != NULL) && (nodes[1] != NULL) && (nodes[2] != NULL)) {
		A =	((nodes[1]->getXc() - nodes[0]->getXc()) * (nodes[2]->getYc() - nodes[0]->getYc())
		  -	 (nodes[2]->getXc() - nodes[0]->getXc()) * (nodes[1]->getYc() - nodes[0]->getYc())  ) / 2.0;
		return A;
	}
	else 
		return -1.0;
}

double Triangle::tribAreas(double* an)
{
	if((nodes[0] != NULL) && (nodes[1] != NULL) && (nodes[2] != NULL)) {		
		double x1 = nodes[0]->getXc();
		double y1 = nodes[0]->getYc();
		double x2 = nodes[1]->getXc();
		double y2 = nodes[1]->getYc();
		double x3 = nodes[2]->getXc();
		double y3 = nodes[2]->getYc();
	 	double X21 = x2 - x1;
		double X32 = x3 - x2;
		double X13 = x1 - x3;
		double Y21 = y2 - y1;
		double Y32 = y3 - y2;
		double Y13 = y1 - y3;
		double d1 = - X13*X21 - Y13*Y21;
		double d2 = - X21*X32 - Y21*Y32;
		double d3 = - X32*X13 - Y32*Y13;
		double c1 = d2*d3;
		double c2 = d3*d1;
		double c3 = d1*d2;
		double c = c1 + c2 + c3;
		double xc = ((c2+c3)*x1 + (c3+c1)*x2 + (c1+c2)*x3)/(2.0*c);
		double yc = ((c2+c3)*y1 + (c3+c1)*y2 + (c1+c2)*y3)/(2.0*c);
		Node cNode(100,xc,yc);
		if(inside(&cNode) == this)
		{
			Triangle t1(1,nodes[0],&cNode,nodes[2]);
			Triangle t2(2,nodes[1],&cNode,nodes[0]);
			Triangle t3(3,nodes[2],&cNode,nodes[1]);
			an[0] = (t1.area() + t2.area())/2.0;
			an[1] = (t2.area() + t3.area())/2.0;
			an[2] = (t3.area() + t1.area())/2.0;
		}
		else
		{
			Segment l1(1,nodes[0],nodes[1]);
			Segment l2(2,nodes[1],nodes[2]);
			Segment l3(3,nodes[2],nodes[0]);
			if ((l1.length() > l2.length()) && (l1.length() > l3.length()))
			{
				Node bisectNode(100,100,100);
				l3.locateNodeAt(0.5,&bisectNode);
				Segment bisectSeg(1,&cNode,&bisectNode);
				double r = l1.intersectd(&bisectSeg);
				Node intersectNode(100,100,100);
				l1.locateNodeAt(r,&intersectNode);
				Triangle tri (1,&intersectNode,&bisectNode,nodes[0]);
				an[0] = tri.area();

				l2.locateNodeAt(0.5,&bisectNode);
				bisectSeg.setNode(1,&bisectNode);
				r = l2.intersectd(&bisectSeg);
				l1.locateNodeAt(r, &intersectNode);
				tri.setNode(1,nodes[1]);
				tri.setNode(2,&bisectNode);
				an[1] = tri.area();

				an[2] = this->area() - an[0] - an[1];	
			}
			else if(l2.length() > l3.length())
			{
				Node bisectNode(100,100,100);
				l1.locateNodeAt(0.5,&bisectNode);
				Segment bisectSeg(1,&cNode,&bisectNode);
				double r = l2.intersectd(&bisectSeg);
				Node intersectNode(100,100,100);
				l2.locateNodeAt(r,&intersectNode);
				Triangle tri (1,&intersectNode,&bisectNode,nodes[1]);
				an[1] = tri.area();

				l3.locateNodeAt(0.5,&bisectNode);
				bisectSeg.setNode(1,&bisectNode);
				r = l2.intersectd(&bisectSeg);
				l2.locateNodeAt(r, &intersectNode);
				tri.setNode(1,nodes[2]);
				tri.setNode(2,&bisectNode);
				an[2] = tri.area();

				an[0] = this->area() - an[1] - an[2];				
			}
			else
			{
				Node bisectNode(100,100,100);
				l2.locateNodeAt(0.5,&bisectNode);
				Segment bisectSeg(1,&cNode,&bisectNode);
				double r = l3.intersectd(&bisectSeg);
				Node intersectNode(100,100,100);
				l3.locateNodeAt(r,&intersectNode);
				Triangle tri (1,&intersectNode,&bisectNode,nodes[2]);
				an[2] = tri.area();

				l1.locateNodeAt(0.5,&bisectNode);
				bisectSeg.setNode(1,&bisectNode);
				r = l3.intersectd(&bisectSeg);
				l3.locateNodeAt(r, &intersectNode);
				tri.setNode(1,nodes[0]);
				tri.setNode(2,&bisectNode);
				an[0] = tri.area();

				an[1] = this->area() - an[0] - an[2];
			}
		}

		return (an[0] + an[1] + an[2]);
	}
	else 
		return -1.0;
}

double Triangle::quality()
{
	if((nodes[0] != NULL) && (nodes[1] != NULL) && (nodes[2] != NULL)) {
	 	double X21 = nodes[1]->getXc() - nodes[0]->getXc();
		double X32 = nodes[2]->getXc() - nodes[1]->getXc();
		double X13 = nodes[0]->getXc() - nodes[2]->getXc();
		double Y21 = nodes[1]->getYc() - nodes[0]->getYc();
		double Y32 = nodes[2]->getYc() - nodes[1]->getYc();
		double Y13 = nodes[0]->getYc() - nodes[2]->getYc();
		double d1 = - X13*X21 - Y13*Y21;
		double d2 = - X21*X32 - Y21*Y32;
		double d3 = - X32*X13 - Y32*Y13;
		double c = d2*d3 + d3*d1 + d1*d2;
		if(c == 0.0)
			return 0.0;
		double r = sqrt( (d1+d2)*(d2+d3)*(d3+d1)/c) / 2.;
		double q = 2.0*area()/(3.0*r*r);
		return q;
	}
	return 0.0;
}

Element* Triangle::inside(Node* otherNodeP)
{
	Segment testSeg1(1,nodes[0],nodes[1]);
	Segment testSeg2(1,nodes[1],nodes[2]);
	Segment testSeg3(1,nodes[2],nodes[0]);
	
	if(otherNodeP != NULL) {
		if(testSeg1.whichSide(otherNodeP) < -0.0000001) {
			return adj[2];
		}
		if(testSeg2.whichSide(otherNodeP) < -0.0000001) {
			return adj[0];
		}
		if(testSeg3.whichSide(otherNodeP) < -0.0000001)  {
			return adj[1];
		}
					
		return this;
	}
	else
		return NULL;
}

double Triangle::insideC(Node* otherNodeP)
{
	double xp, yp, x1, y1, x2, y2, x3, y3, numr, numi, denr, deni;
	
	if(otherNodeP != NULL) {
		xp = otherNodeP->getXc();
		yp = otherNodeP->getYc();
	}
	else
		return -1.0;
		
	x1 = nodes[0]->getXc();
	y1 = nodes[0]->getYc();
	x2 = nodes[1]->getXc();
	y2 = nodes[1]->getYc();
	x3 = nodes[2]->getXc();
	y3 = nodes[2]->getYc();
	
	numr = (x2 - x3)*(xp - x1) - (y2 - y3)*(yp - y1);
	numi = (x2 - x3)*(yp - y1) + (y2 - y3)*(xp - x1);
	denr = (x2 - x1)*(xp - x3) - (y2 - y1)*(yp - y3);
	deni = (x2 - x1)*(yp - y3) + (y2 - y1)*(xp - x3);
	
	return denr*numi - deni*numr;
}

Triangle* Triangle::contains(Node* nodeP)
{	
	if(nodeP != NULL) {
		if( (nodes[0] == nodeP) || (nodes[1] == nodeP) || (nodes[2] == nodeP) ) {
			return this;
		}
		Element *aP = inside(nodeP);
		while(aP->numNodes() != 3)
			aP = aP->inside(nodeP);
		return (Triangle *)aP;
	}
	return NULL;
}

Triangle* Triangle::actContains(Node* nodeP)
{
	if(nodeP != NULL) {
		int index = nodeIndex(nodeP);
		if(index < 0){
			Element *aP = inside(nodeP);
			if(aP != NULL){
				while(aP->numNodes() != 3){
					aP = aP->inside(nodeP);
					if((aP == NULL) || (aP == this))
						return NULL;
				}
				return (Triangle *)aP;
			}
			else
         	return NULL;
		}
		if(status == active)
			return this;
	}
	return NULL;
}

Triangle* Triangle::contains(Segment* segP)
{	
	if(segP != NULL) {
		if( (nodes[0] == segP->getNode(0)) &&
			(nodes[1] == segP->getNode(1)) ) {
			return this;
		}
		if( (nodes[1] == segP->getNode(0)) &&
			(nodes[2] == segP->getNode(1)) ) {
			return this;
		}
		if( (nodes[2] == segP->getNode(0)) &&
			(nodes[0] == segP->getNode(1)) ) {
			return this;
		}
	}
	
	return NULL;
}

Triangle* Triangle::adjContaining(Node* nP1)
{
	Element * pP;
	Triangle * tP;
	
	for(int i=0;i<3;i++) {
		pP = adj[i];
		if(pP != NULL)
			if(pP->numNodes() == 3) {
				tP = (Triangle *) pP;
				if(tP->nodeIndex(nP1) >= 0)
					return tP;
			}
	}
	return NULL;
}

Node* Triangle::locateNodeAt(double x, double y, Node* nP)
{
	double weights[3];
 	double X21 = nodes[1]->getXc() - nodes[0]->getXc();
	double X31 = nodes[2]->getXc() - nodes[0]->getXc();
	double Y21 = nodes[1]->getYc() - nodes[0]->getYc();
	double Y31 = nodes[2]->getYc() - nodes[0]->getYc();
	double X = x - nodes[0]->getXc();
	double Y = y - nodes[0]->getYc();
	double det = X21*Y31 - X31*Y21;
	double r = (Y31*X - X31*Y)/det;
	double s = (-Y21*X + X21*Y)/det;
	double t = 1.0 - r - s;

	weights[0] = t;
	weights[1] = r;
	weights[2] = s;

	nP->interp(3,nodes,weights);

	return nP;
}

Node* Triangle::locateNodeAtCenter(Node* nP)
{
	double weights[3];
	double r = 1.0/3.0;
	double s = r;
	double t = r;

	weights[0] = t;
	weights[1] = r;
	weights[2] = s;

	nP->interp(3,nodes,weights);

	return nP;
}

Node* Triangle::locateNodeAt(Node* nP)
{
	double weights[3];
 	double X21 = nodes[1]->getXc() - nodes[0]->getXc();
	double X31 = nodes[2]->getXc() - nodes[0]->getXc();
	double Y21 = nodes[1]->getYc() - nodes[0]->getYc();
	double Y31 = nodes[2]->getYc() - nodes[0]->getYc();
	double X = nP->getXc() - nodes[0]->getXc();
	double Y = nP->getYc() - nodes[0]->getYc();
	double det = X21*Y31 - X31*Y21;
	double r = (Y31*X - X31*Y)/det;
	double s = (-Y21*X + X21*Y)/det;
	double t = 1.0 - r - s;

	weights[0] = t;
	weights[1] = r;
	weights[2] = s;

	nP->interp(3,nodes,weights);

	return nP;
}

Node* Triangle::getClosestNode(Node* nP)
{
	Node *theNodeP = nodes[0];
	int	iNode;
	Element *adjtElP;
	double distance;
	double minDistance = nP->dist(theNodeP);
	for(int i=1;i<=2;i++){
		distance = nP->dist(nodes[i]);
		if(distance < minDistance){
			minDistance = distance;
			theNodeP = nodes[i];
		}
	}
/*
	// Check adjacent triangles as well
	for(int iadj=0;iadj<3;iadj++){
		adjtElP = getAdjt(iadj);
		if(adjtElP != NULL){
			iNode = adjtElP->reflectAdj(this);
			distance = nP->dist(adjtElP->getNode(iNode));
			if(distance < minDistance){
				minDistance = distance;
				theNodeP = adjtElP->getNode(iNode);
			}
		}
	}
*/
	return theNodeP;
}

