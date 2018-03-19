
//		Quadrant.h


//		Initiated:	June 2, 2001
//
//		Last Revised: June 2, 2001

#ifndef QUADRANT_H
#define QUADRANT_H

#include "Triangle.h"

struct box
{
	double	x1;
	double	y1;
	double	x2;
	double	y2;
};

#define MINDIST 0.0001

class Quadrant
{
	protected:
		Quadrant	*subQuadrant[4];  // 4 subquadrants ne, nw, sw, se
		Node		*firstNode, *lastNode;
		int			nNodes, nTriangles, nBoundSegs;
		Triangle	*firstTriangle, *lastTriangle;
		Triangle	*bArc;
		box			boundBox;
		
	public:
		Quadrant(box limits);
		void appendNode(Node *newNode);
		void appendTriangle(Triangle *newTri);
		void deleteTriangle(Triangle *tP);
		int triangulateNodes();
		int recoverSubQuadInfo(int i);
		int getNNodes() {return nNodes;}
		int getNTrinangles() {return nTriangles;}
		Node* getNodes() {return firstNode;}
		Triangle* getTriangles() {return firstTriangle;}
		Triangle* getBArcs() {return bArc;}
		Triangle* stitchLR(Triangle *bArcL, Triangle *bArcR);
		void delaunay(Triangle *tP);
		Triangle* pushTri(Triangle *stack, Triangle *tP);
		Triangle* popTri(Triangle *stack);
		void /*Quadrant::*/swapDiagonal(Triangle *tP1, Triangle *tP2);//commented May 14 
		int checkDuplicateNodes();
};

#endif
