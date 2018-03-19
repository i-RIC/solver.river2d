
//		Segment.h

//		Segment.h is a base class which defines a straight line segment 
//		between two points. The points are defined as from the Node class.
//		Basic geometric functionality is provided.
//		Subclass a Segment for one-dimensional elements.

//		Initiated:	September 24, 1995
//
//		Last Revised: September 24, 1995

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Element.h"
#include "Segment.h"


class Triangle : public Element
{
	protected:
		Node *	nodes[3];
		Element	*adj[3];
		Segment *edges[3];
		int		toBeRefined;
		double  A;
		
	public:
		Triangle *ntP, *ptP;
		Triangle(int num=1, Node* firstNodeP=NULL, Node* secondNodeP=NULL, Node* thirdNodeP=NULL);
		Node** getNodes() {return nodes;}
		virtual Node* getNode(int i) {return nodes[i];}
		virtual void setNode(int i, Node* np) {if(i < 3) nodes[i] = np;}
		virtual Segment* getEdge(int i) {return edges[i];}
		virtual void setEdge(int i, Segment* ep) {if(i < 3) edges[i] = ep;}
		virtual int	nodeIndex(Node* nP);
		virtual Element** getAdj() {return adj;}
		virtual Element* getAdjt(int i) {return adj[i];}
		Triangle* tP0() {return (Triangle *)adj[0];}
		Triangle* tP1() {return (Triangle *)adj[1];}
		Triangle* tP2() {return (Triangle *)adj[2];}
		Node* node0() {return nodes[0];}
		Node* node1() {return nodes[1];}
		Node* node2() {return nodes[2];}
		void setRefineOn(){toBeRefined = 1;}
		void setRefineOff(){toBeRefined = 0;}
		int isRefineOn(){return toBeRefined;}
//		int tsn() {return reflectAdj(adj[2]);}
		int tsn() {return adj[2]->reflectAdj(this);}
		virtual void setAdj(int i,Element* ap) {if(i < 3) adj[i] = ap;}
		virtual int numNodes() {return 3;}
		virtual int numAdj(){return 3;}
		virtual char* typeName() {return "Triangle";}
		double area();
		double tribAreas(double* an);
		double quality();
		virtual Element* inside(Node* otherNodeP); //  pointer to self if inside, otherwise to next triangle
		Triangle* contains(Node* nodeP); //  pointer to self if nodeP is one the tri nodes, otherwise to next triangle
		Triangle* actContains(Node* nodeP); //  pointer to self if nodeP is one the tri nodes, otherwise to next triangle
		Triangle* contains(Segment* segP); //  pointer to self if segment inside, otherwise to next triangle
		Triangle* adjContaining(Node* nP1); // pointer to adjacent triangle which contains node
		double insideC(Node* otherNodeP);	//	positive if inside of circumcircle
		Node* locateNodeAt(double x, double y, Node* nP);	//	locates a node at x, y and interpolates
		Node* locateNodeAtCenter(Node* nP);	//	locates a node at triangle center and interpolates
		Node* locateNodeAt(Node* nP);	//	interpolates a node which already has x, y set
		Node* getClosestNode(Node* nP); // returns the triangle node closest to the given node
};

#endif
