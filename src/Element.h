
//		Element.h

//		The Element class is an abstract base class for various
//		objects comprised of connected sets of nodes. As such
//		it only sets out the basic form of the connectivity table
//		with the directed list of nodes and adjacent objects to be 
//		allocated by the specific objects.
//		This object declares some of the basic geometric operations
//		necessary in Finite Element analysis through pure functions.
//		Implementation is up to the specific subclasses.
//		Element is derived from Item, so it is listable with an
//		ItemList type controller.

#ifndef ELEMENT_H
#define ELEMENT_H

#include <iostream>
#include "Node.h"
using namespace std;

enum statusFlag {notActive, active};

class Element : public Item
{		
	protected:
		statusFlag status;
		 
	public:
		Element(int num=1):Item(num)
		{status = active;}
		virtual Node** getNodes() = 0;
		virtual Node* getNode(int i) = 0;
		virtual int nodeIndex(Node *) = 0;
		virtual void setNode(int i, Node* np) =0;
		virtual Element** getAdj() = 0;
		virtual Element* getAdjt(int i) = 0;
		virtual void setAdj(int i, Element* ap) =0;
		virtual int numNodes() = 0;
		virtual int numAdj() = 0;
		virtual char* typeName() = 0;
		virtual Element* inside(Node* otherNodeP) = 0;
		int reflectAdj(Element* ap);
      Node centroid();
		statusFlag getStatus() {return status;}
		void deactivate() {status = notActive;}
		void activate() {status = active;}
		friend ostream& operator << (ostream& os, Element* n);

};

#endif
