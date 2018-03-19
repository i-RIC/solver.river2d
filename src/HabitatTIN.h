// HabitatTIN.h: interface for the HabitatTIN class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_HABITATTIN_H__BFC1DD71_EB52_4B28_9E96_15EDABE22B76__INCLUDED_)
#define AFX_HABITATTIN_H__BFC1DD71_EB52_4B28_9E96_15EDABE22B76__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Tin.h"
#include "Habitat.h"
#include "Triangle.h"
#include "ItemList.h"

enum FEMProblem
{
	psiSolution,
	phiSolution
};


class HabitatTIN : public TIN  
{
private:
	ItemList BElMtxL;
	ItemList ElMtxL;
	FEMProblem FEMProb;

public:
	HabitatTIN(Physics *theProb);
	virtual ~HabitatTIN();
	HabitatNode* insertOneFloatingNode(double x, double y, TIN& boundTIN, TIN& bedTIN, TIN& meshTIN);
	double smoothHabitatMesh(int nTimes, TIN& dataTIN, double bias = 0.0);
	int refineHabitatRegion(TIN &boundTIN, TIN &dataTIN);
	Segment* setBoundaryPsi(HabitatNode* startnP);
	void setAllPsi(HabitatNode* startnP);
	int extractBoundary(TIN* meshTIN, HabitatNode* hNUSBound, HabitatNode* hNDSBound, double leftBound, double rightBound, double spacing, int leftBoundType, int rightBoundType);
	HabitatNode* getBoundaryNode(TIN* meshTIN, HabitatNode* hNP, double boundaryValue, int flag, int parChoice);
	Triangle* findTriWithBoundaryNode(Triangle* triP, Segment* segPathP, double boundaryValue, int parChoice);
	Triangle* mapStreamLine(TIN* meshTIN, Triangle* triP, HabitatNode* hNPStart, HabitatNode* hNPEnd, double psiBoundary, int flag, int parChoice);
	void resampleBoundary(double spacing);
	int checkNoOfInflowSegs();
	void mapBoundaryStreamLine(TIN* meshTIN, HabitatNode* hNUSBound, HabitatNode* hNDSBound, int flag);
	Triangle *findTriAndNodeWhenVelDirCannot(TIN* meshTIN, Triangle *triP, double boundaryValue, int flag, int parChoice);
	int autoRefineHabitatRegion(TIN* boundTIN, TIN &dataTIN, double lowerLimit, double upperLimit, int parNum, int valueOrChange);
	int deletePrimaryMeshNode(Node * nP);
	HabitatNode* removeNode(Node *nP);
	void buildPsiMatrices();
	void buildPhiMatrices(double wsElevIn);
	void destroyMatrices();
	void calculatePsi(HabitatNode* startnP);
	void calculatePhi(double wsElevIn);
	void initSolutionVector(double *X);
	void updateSolToNodes(double *X);
};

#endif // !defined(AFX_HABITATTIN_H__BFC1DD71_EB52_4B28_9E96_15EDABE22B76__INCLUDED_)
