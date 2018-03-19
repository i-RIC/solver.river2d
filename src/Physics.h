
//		Physics.h

//		The Physics Class is an abstract base class which defines
//		basic functionality of a physical problem description.
//		Specific problem types should be derived from the physics
//		class. The derived class should define the necessary variables
//		and parameters at each node, create approproiate nodes, give
//		internal and boundary weak statements. Material properties
//		and behaviours should be defined there.
//		The TIN and derived mesh classes contain a pointer to a
//		physics class which is referred to for various operations.
//		This physics class must be specified when the TIN is
//		created and provision is made for changeing it.


#ifndef PHYSICS_H
#define PHYSICS_H


#include "Node.h"
#include "Segment.h"
#include <fstream>
using namespace std;


class Physics 
{		
	public:
		Physics()
		{}
		virtual Node* makeNewNode(int n=1,double x= 100.0,
				double y= 100.0,double z= 100.0) = 0;
		virtual Segment* makeNewSegment(int n=1,Node *nP1= NULL,
				Node *nP2= NULL, Segment *segP=NULL) = 0;
		virtual Node* readInNode(istream& is) = 0;
		virtual Segment* readInSeg(int n, Node *nP1, Node *nP2, istream& is) = 0;
		virtual void writeBSeg(ostream& os, Segment *segP) = 0;
		virtual int nVars() = 0;
		virtual int nParams() = 0;
		virtual char* typeName() = 0;
};

#endif
