
//			Node.h
//
//	Defines an base node class for just about anything,
//	including finite element analysis. Subclass it for specific problems
//	and add variables and parameters. 

//	Subclassed from Item, so it can be link listed.
//
//	Initiated:	September 24, 1995
//
//	Last Revised: November 18, 1995

#ifndef NODE_H
#define NODE_H

#include "Item.h"
#include "string.h"
#include <fstream>
using namespace std;

enum boundFlag {interior, onBoundary};
enum fixedFlag {floating, fixednode, sliding, deleted};

class Triangle;

class Node : public Item
{
	private:
		double xo;
		double yo;

	protected:
		double x;			// 	x coordinate							
		double y;			// 	y coordinate							
		double z;			//	bed elevation
		boundFlag bFlag;	//	boundary flag for various purposes
		fixedFlag fixFlag;	//	can node be moved or not
		Triangle *aTp;		//  pointer to one of the triangles that contains this node
		char *designation;  // string containing user defined designation for node
	public:
		Node(int name=1,double xcoord=10.0,double ycoord=10.0,double zcoord=10.0);
		Node(istream& is);
		~Node();
		double getXc() const {return x;}	
		double getYc() const {return y;}
		double getZc() const {return z;}
		double getXo() const {return xo;}	
		double getYo() const {return yo;}
		char* getDesignation() const {return designation;}
		void setDesignation(const char* desString);
		virtual int getNParams() {return 0;}
		virtual int getNVars() {return 0;}
		virtual double getPar(int i);
		virtual double getVar(int i);
		boundFlag getBound() const {return bFlag;}
		void setOnBound() {bFlag = onBoundary;}
		void setInterior() {bFlag = interior;}
		fixedFlag getFixed() const {return fixFlag;}
		void setFixed() {fixFlag = fixednode;}
		void setFloating() {fixFlag = floating;}
		void setSliding() {fixFlag = sliding;}
		void setDeleted() {fixFlag = deleted;}
		void setATriangle(Triangle* theTp) {aTp = theTp;}
		Triangle* getATriangle() {return aTp;}
		void assign(double xcoord,double ycoord,double zcoord);
		void assignt(double xcoord,double ycoord,double zcoord);
		friend ostream& operator << (ostream& os, Node* n);
		double dist(Node* otherNode);
		virtual void interp(int n, Node** nPtrs, double* wts);
		void saveLocation();
		void restoreLocation();

};

#endif
