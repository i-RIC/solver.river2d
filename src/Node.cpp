
//		Node.cpp

//		Node class function definitions

#include <math.h>
#include <iostream>
#include <ctype.h>
#include "Node.h"
#include "Triangle.h"

Node::Node(	int name,double xcoord,double ycoord,double zcoord)
			:Item(name)
{
	bFlag = interior;
	fixFlag = floating;
	aTp = NULL;
	assign(xcoord, ycoord, zcoord);
	designation = NULL;
}

Node::Node(istream& is) : Item()
{
	int num = -990;
	char c = 'f';
	double xcoord=0.0, ycoord=0.0, zcoord=0.0;
	bFlag = interior;
	fixFlag = floating;
	designation = NULL;
	while(is) {
		if(!(is >> num)){
			num = -991;
			break;
		}
		while(is) {
			c = is.peek();
			if( (c == '-') || (c == '.') || (isdigit(c)) )
				break;
			else if (c == 'x')
				fixFlag = fixednode;
			else if (c == 's')
				fixFlag = sliding;
			c = is.get();
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
			num = -993;
			break;
		}
		break;
	}
	setn(num);
	assign(xcoord, ycoord, zcoord);
}

Node::~Node()
{
	if(designation != NULL)
		delete[] designation;
}

void Node::assign(double xcoord,double ycoord,double zcoord)
{
	x = xcoord;
	y = ycoord;
	z = zcoord;
	saveLocation();
}

void Node::assignt(double xcoord,double ycoord,double zcoord)
{
	x = xcoord;
	y = ycoord;
	z = zcoord;
}

void Node::setDesignation(const char* desString)
{
	if(designation != NULL)
		delete[] designation;
	designation = new char[strlen(desString)+1];
	strcpy(designation,desString);
}

ostream& operator <<(ostream& os, Node* n)
{
	if(n != NULL){
		os.width(5); 	os 	<< n->n << " ";
		if(n->getFixed() == fixednode)
			os << "x ";
		else if (n->getFixed() == sliding)
			os << "s ";
		else
      	os << "  ";
		os.precision(6);
		os.setf(ios::fixed,ios::floatfield);
		os.width(14);	os 	<< n->xo << " ";
		os.width(14);	os 	<< n->yo << " ";
		for(int i=1;i<=n->getNParams();i++){
			os.width(10);	os 	<< n->getPar(i) << " ";
		}
	}

	return os;
}

double Node::getPar(int i)
{
	if(i == 1){
		return z;
	}
	else
   	return 0.0;
}

double Node::getVar(int i)
{
	return 0.0;
}

double Node::dist(Node* otherNode)
{
	double dx;
	double dy;
	
	if(otherNode != NULL) {
		dx = otherNode->x - x;
		dy = otherNode->y - y;
	
		return sqrt(dx*dx + dy*dy);
	}
	else
		return (-1.0);
} 

void Node::interp(int n, Node** nPtrs, double* wts)
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
	for(int i=0;i<n;i++) {
		if(nPtrs[i] != NULL) {
			x += wts[i] * nPtrs[i]->getXo();
			y += wts[i] * nPtrs[i]->getYo();
			z += wts[i] * nPtrs[i]->getZc();
		}
	}
	saveLocation();
}

void Node::saveLocation()
{
	xo = x;
	yo = y;
	return;
}

void Node::restoreLocation()
{
	x = xo;
	y = yo;
	return;
}

