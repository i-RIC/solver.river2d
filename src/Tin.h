
//			Tin.h

//			Basic Triangular Irregular Network Class
//			
//
//			Can be subclassed for various FEM meshes based on linear triangles.
//
#ifndef TIN_H
#define TIN_H

#include "Node.h"
#include "Element.h"
#include "Segment.h"
#include "Triangle.h"
#include "Itemlist.h"
#include "Physics.h"
#include "Quadrant.h"
#include <fstream>
using namespace std;


class TIN
{
	protected:
		ItemList	nodeL;		//		linked list of nodes
		ItemList	bsegL;		//		linked list of boundary segments
		ItemList	elmL;		//		linked list of elements (triangles)
		ItemList	edgeL;		//		linked list of triangle edge segments
		ItemList	fsegL;		//		linked list of feature segments (breaklines)
		ItemList	tnodeL;		//		linked list of temporary nodes
		ItemList	tbsegL;		//		linked list of temporary boundary segments
		ItemList	tfsegL;		//		linked list of temporary feature segments (breaklines)
		ItemList	deadNodeL;	//		linked list of nodes which have been deleted
		box			limits;		//		bounding box surrounding all nodes
		Physics*	physics;	//		pointer to physical problem description
		int	serialNum;

	public:
		TIN(Physics *theProb) {physics = theProb; serialNum = 0;}
		~TIN();
		int	getNextN();
		Physics* getPhysics() {return physics;}
		int Nnodes() {return nodeL.numberOfItems();}
		int Nbsegs() {return bsegL.numberOfItems();}
		int Nfsegs() {return fsegL.numberOfItems();}
		int Nelms() {return elmL.numberOfItems();}
		int Nedges() {return edgeL.numberOfItems();}
		int NTnodes() {return tnodeL.numberOfItems();}
		int NTfsegs() {return tfsegL.numberOfItems();}
		Segment* firstSeg() {return (Segment *)bsegL.firstItem();}
		Segment* nextSeg() {return (Segment *)bsegL.nextItem();}
		Segment* firstfSeg() {return (Segment *)fsegL.firstItem();}
		Segment* nextfSeg() {return (Segment *)fsegL.nextItem();}
		Segment* firstTfSeg() {return (Segment *)tfsegL.firstItem();}
		Segment* nextTfSeg() {return (Segment *)tfsegL.nextItem();}
		Node* firstNode() {return (Node *) nodeL.firstItem();}
		Node* nextNode() {return (Node *) nodeL.nextItem();}
		Node* firstTNode() {return (Node *) tnodeL.firstItem();}
		Node* nextTNode() {return (Node *) tnodeL.nextItem();}
		Triangle* firstTri() {return (Triangle *) elmL.firstItem();}
		Triangle* nextTri() {return (Triangle *) elmL.nextItem();}
		Segment* firstEdge() {return (Segment *)edgeL.firstItem();}
		Segment* nextEdge() {return (Segment *)edgeL.nextItem();}
		void appendNode(Node *nP) {nodeL.appendItem(nP);}
		void insertNode(Node *nP) {nodeL.insertItem(nP);}
		void setCurrentNode(Node *nP) {nodeL.setCurrentItem(nP);}
		Node* getNodeByN(int n) {return (Node *) nodeL.n(n);}
		void appendBSeg(Segment *bP) {bsegL.appendItem(bP);}
		int readNodes(istream& is);
		int readNodesDes(istream& is);
		int checkDupNode(Node *nP);
		int checkDupNodeTri(Node *nP, Triangle *tP);
		void checkAllNodes();
		void writeNodes(ostream& os);
		void writeNodesDes(ostream& os);
		int readBoundSegs(istream& is);
		void writeBoundSegs(ostream& os);
		int readFeatureSegs(istream& is);
		void checkAllFeatureSegs();
		Segment* addFeatureSeg(Segment* newfsegp, ItemList& fList);
		Segment* addTFeatureSeg(Segment* newfsegp);
		int checkFeatureSeg(Segment* newfsegp);
		void writeFeatureSegs(ostream& os);
		void writeTriangles(ostream& os);
		void dumpCSVFile(ostream& os,int index);
		void dumpCSVFile(ostream& os, int index1, int index2);
		double meshQuality(Triangle** thisOne);
		int boundaryRefine(double minQual);
		int triangulate();
		int triangulateCons();
		int	doEdges();
		int checkDupNodes();
		Triangle* whichTriangle(Node *nodeP);
		Triangle* currentTriangle() {return (Triangle *) elmL.currentItem();}
		Triangle* whichTriHasNode(Node *nodeP);
		Triangle* triangleListForNode(Node *nodeP, Triangle *stP);
		Triangle* whichTriHasSeg(Segment *segP, Triangle *startTP = NULL);
		Triangle* getFirstTriOnSeg(Segment *segP);
		Triangle* getSecondTriOnSeg(Triangle *fTP, Segment *segP);
		void insertBoundSeg(Segment *bP,Triangle *triP);
		void insertNewNode(Node *nodeP, Triangle *triP);
		void insertNewNodeCons(Node *nodeP, Triangle *triP, Segment *bP = NULL, Triangle *stP = NULL);
		void swapDiagonal(Triangle *tP1, Triangle *tP2);
		int doFeatures(TIN* dataTINP);
      Segment* threadFeature(Segment *fsegP, ItemList &segList);
		Segment* bisectSeg(Segment *segP, ItemList &segL, double dist = 0.5);
		Segment* bisectBSeg(Segment *segP, double dist = 0.5);
		double smoothMesh(int nTimes, TIN& dataTIN, double bias = 0.0);
		double triCenterDz(Triangle *triP);
		int updateMesh(TIN& dataTIN);
		box getLimits(double margin = 0.01);
		void scratchTriangle();
		void clearBsegs();
		void acceptNodes();
		void rejectNodes();
		void acceptBsegs();
		void acceptNewBsegs();
		void rejectBsegs();
		void acceptFsegs();
		void rejectFsegs();
		Node* insertOneNode(double x, double y, TIN& boundTIN, TIN &dataTIN);
		Node* insertOneNode(double x, double y);
		Segment* defNewBound(double x, double y,TIN &dataTIN);
		Segment* defNewBound(double x, double y);
		Segment* closeBound(TIN &dataTIN);
		Segment* closeBound();
		Segment* closeBoundLoop(Node *firstNP, TIN &dataTIN);
		int fillNodesUniform(double spacing, double theta, box dLimits,
									 TIN& boundTIN, TIN& dataTIN, TIN& bound2TIN);
		int fillNodesSource(double x0, double y0, double r0, int nRays,
								  double theta0, TIN& boundTIN, TIN &dataTIN);
		int refineRegion(TIN& boundTIN, TIN& dataTIN);
		int refineRegionYellow(TIN& boundTIN, TIN& dataTIN, double minDisplayDz);
		int defineBoundary(double spacing, TIN& boundTIN, TIN& dataTIN);
		int defineNewBoundLoop(double spacing, TIN& boundTIN, TIN& dataTIN, Segment* firstSP);
		int orderNodesRCM();
		int deleteNode(Node * nP);
		void removeFSeg(Segment *fSegP);
		void removeTFSeg(Segment *fSegP);
		double getMinPar(int nPar);
		double getMaxPar(int nPar);

};

#endif
