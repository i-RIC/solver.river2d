
//			Shallow.h


//		Shallow is a subclass of Physics which describes
//		river bed topography description.
//		The nodes produced are of the BedTopo
//		Node class with 3 coordinates plus a bed roughness height.

//		Initiated:	March 9, 1996
//
//		Last Revised: March 9, 1996

#ifndef SHALLOW_H
#define SHALLOW_H

#include "Physics.h"
#include "FlowBoundary.h"

const double G = 9.806;

class ShallowNode : public Node
{
	protected:
		double ks;			// 	roughness height
		double d;			//		flow depth
		double qx;			//		discharge intensity in x direction
		double qy;			//		discharge intensity in y direction

	public:
		ShallowNode(	int nm=1,
					double xc=1000.0,
					double yc=1000.0,
					double zc=100.0,
					double rf = 0.05,
					double depth = 1.0,
					double xDis = 0.0,
					double yDis = 0.0
					);
		ShallowNode(istream& is);
		void setZc(double elev) {z = elev;}
		double getKs() const {return ks;}
		void setKs(double rough) {ks = rough;}
		double getD() const {return d;}
		void setD(double depth) {d = depth;}
		void setFlow(double depth, double xDis, double yDis)
			{d = depth; qx = xDis; qy = yDis;}
		virtual int getNParams() {return 2;}
		virtual int getNVars() {return 3;}
		virtual double getPar(int i);
		virtual double getVar(int i);
		virtual void getParName(int i, char* parName);
		virtual void interp(int n, Node** nPtrs, double* wts);
};

class RivBound : public Segment
{
	protected:
		int	bcCode;
		double bcValue;
		double bcValue2;
		FlowBoundary *flowBound;

	public:
		RivBound(int num=1, Node* nP1=NULL, Node* nP2=NULL,
						 int c=0, double v=0.0);
		void setBcCode(int c) {bcCode = c;}
		virtual int getBcCode() {return bcCode;}
		void setBcValue(double d) {bcValue = d;}
		double getBcValue() {return bcValue;}
		void setBcValue2(double d) {bcValue2 = d;}
		double getBcValue2() {return bcValue2;}
		friend ostream& operator << (ostream& os, RivBound* seg);
		void setFlowBound(FlowBoundary* fBP) {flowBound = fBP;}
		FlowBoundary* getFlowBound() {return flowBound;}
};

class Shallow : public Physics
{
	protected:
		double minDepth;

	public:
		Shallow() : Physics()
		{minDepth = 0.02;}
		virtual Node* makeNewNode(int n=1,double x= 100.0,
				double y= 100.0,double z= 100.0);
		virtual Segment* makeNewSegment(int n=1,Node *nP1= NULL,
				Node *nP2= NULL, Segment *segP=NULL);
		virtual Node* readInNode(istream& is);
		virtual Segment* readInSeg(int n, Node* nP1, Node* nP2, istream& is);
		virtual void writeBSeg(ostream& os, Segment *segP);
		virtual int nVars() {return 3;}
		virtual int nParams() {return 2;}
		virtual char* typeName() {return "Shallow";}
		void setMinDepth(double d) {minDepth = d;}
		double getMinDepth() {return minDepth;}
};

#endif

