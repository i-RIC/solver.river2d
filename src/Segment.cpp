
//		Segment.cpp

//		Segment Class method definitions

#include <math.h>
#include "Segment.h"

Segment::Segment(int name, Node* firstNodeP, Node* secondNodeP)
		: Element(name)
{	
	nodes[0]=firstNodeP;
	nodes[1]=secondNodeP;
	adj[0] = NULL;
	adj[1] = NULL;
}

int Segment::nodeIndex(Node* nP)
{
	for(int i=0;i<numNodes();i++){
		if(nP == nodes[i])
			return i;
	}
	
	return -1;
}

double Segment::length()
{
	if((nodes[0] != NULL) && (nodes[1] != NULL))
		return nodes[0]->dist(nodes[1]);
		
	else 
		return -1.0;
}

double Segment::whichSide(Node* otherNodeP)
{
	if(otherNodeP != NULL)
		return (otherNodeP->getXc() - nodes[0]->getXc()) * (nodes[0]->getYc() - nodes[1]->getYc())
			 + (otherNodeP->getYc() - nodes[0]->getYc()) * (nodes[1]->getXc() - nodes[0]->getXc());
	else
		return -1.0;
}

double Segment::midDistance(Node* otherNodeP)
{
	if(otherNodeP != NULL){
		Node midPoint(0,(nodes[0]->getXc()+nodes[1]->getXc())/2.,
							 (nodes[0]->getYc()+nodes[1]->getYc())/2.,100.0);
		return midPoint.dist(otherNodeP);
	}
	else
		return -1.0;
}

Element* Segment::inside(Node* otherNodeP)
{
	if(otherNodeP != NULL) {
		if(whichSide(otherNodeP) > 0.0) 
			return adj[0];
		else 
			return adj[1];
	}
	else
		return NULL;
}


void Segment::locateNodeAt(double pos, Node* nP)
{
	double weights[2];
	
	if(pos <= 0.0) {
		weights[0] = 1.0;
		weights[1] = 0.0;
	}
	else if(pos >= 1.0) {
		weights[0] = 0.0;
		weights[1] = 1.0;
	}
	else {
		weights[0] = 1.0 - pos;
		weights[1] = pos;
	}
	nP->interp(2,nodes,weights);
	
}

int Segment::intersect(Segment* otherSegP)
{
	Node* nP1 = this->getNode(0);
	Node* nP2 = this->getNode(1);
	Node* oP1 = otherSegP->getNode(0);
	Node* oP2 = otherSegP->getNode(1);

//	if( ((nP1 == oP1) && (nP2 == oP2)) || ((nP1 == oP2) && (nP2 == oP1)) )
//		 return -1;
	if( (nP1 == oP1) || (nP1 == oP2) || (nP2 == oP1) || (nP2 == oP2) )
		 return 0;

	 double ny = oP2->getXc() - oP1->getXc();
	 double nx = oP1->getYc() - oP2->getYc();

	 double num = nx*(oP1->getXc() - nP1->getXc()) + ny*(oP1->getYc() - nP1->getYc());
	 double den = nx*(nP2->getXc() - nP1->getXc()) + ny*(nP2->getYc() - nP1->getYc());

	 if(den == 0.0)
		return 0;

	 double r = num/den;
	 if( (r >= 0.0) && (r <= 1.0) )
		return -1;

	 return 1;
}

double Segment::intersectd(Segment* otherSegP)
{
	Node* nP1 = this->getNode(0);
	Node* nP2 = this->getNode(1);
	Node* oP1 = otherSegP->getNode(0);
	Node* oP2 = otherSegP->getNode(1);

	if( ((nP1 == oP1) && (nP2 == oP2)) || ((nP1 == oP2) && (nP2 == oP1)) )
		 return -1.0;
	if( (nP1 == oP1) || (nP1 == oP2) || (nP2 == oP1) || (nP2 == oP2) )
		 return 0.0;

	 double ny = oP2->getXc() - oP1->getXc();
	 double nx = oP1->getYc() - oP2->getYc();

	 double num = nx*(oP1->getXc() - nP1->getXc()) + ny*(oP1->getYc() - nP1->getYc());
	 double den = nx*(nP2->getXc() - nP1->getXc()) + ny*(nP2->getYc() - nP1->getYc());

	 if(den == 0.0)
		return 0.0;

	 double r = num/den;
	 if( (r >= 0.0) && (r <= 1.0) )
		return r;

	 return -1.0;
}



