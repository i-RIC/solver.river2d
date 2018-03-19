
//		Segment.h

//		Segment.h is a base class which defines a straight line segment 
//		between two points. The points are defined as from the Node class.
//		Basic geometric functionality is provided.
//		Subclass a Segment for one-dimensional elements.

//		Initiated:	September 24, 1995
//
//		Last Revised: September 24, 1995

#ifndef SEGMENT_H
#define SEGMENT_H

#include "Element.h"


class Segment : public Element
{
	protected:
		Node *	nodes[2];
		Element	*adj[2];			//    to left, then to right
		
	public:
		Segment(int num=1, Node* firstNodeP=NULL, Node* secondNodeP=NULL);
		Node** getNodes() {return nodes;}
		virtual Node* getNode(int i) {return nodes[i];}
		virtual void setNode(int i, Node* np) {if(i < 2) nodes[i] = np;}
		virtual int nodeIndex(Node *);
		virtual Element** getAdj() {return adj;}
		virtual Element* getAdjt(int i) {return adj[i];}
		virtual void setAdj(int i,Element* ap) {if(i < 2) adj[i] = ap;}
		virtual int numNodes() {return 2;}
		virtual int numAdj(){return 2;}
		virtual char* typeName() {return "Segment";}
		virtual int getBcCode() {return 0;}
		double length();
		double whichSide(Node* otherNodeP);   	//   positive is to the left of the segment
		double midDistance(Node* otherNodeP);   	//   distance from middle of segment to node
		virtual Element* inside(Node* otherNodeP);	//	just a pass through for searches
		void locateNodeAt(double pos, Node* nP);		//	creates a new node at pos along the segment and interpolates
      int intersect(Segment* otherSegP);			// returns 1 if no intersection, 0 if special, -1 if intersecting
      double intersectd(Segment* otherSegP);			// returns dist to intersection, -1.0 if no intersecting
};

#endif
