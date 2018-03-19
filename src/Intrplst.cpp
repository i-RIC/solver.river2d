//
//		InterpLst.cpp
//

#include "Intrplst.h"
#include <iostream>
#include <ctype.h>
using namespace std;
	
InterpLst::InterpLst()
{
	numInterpPts = 0;
	theFirstInterpPt = NULL;
	theCurrentInterpPt = NULL;
	InterpPtIndex = NULL;
	indexState = lnotValid;
}

InterpLst::InterpLst(istream& is)
{
	numInterpPts = 0;
	theFirstInterpPt = NULL;
	theCurrentInterpPt = NULL;
	InterpPtIndex = NULL;
	indexState = lnotValid;

	char c = 'a';
	InterpPt *pt;

	while( is && (c != '(')){
		is.get(c);
		cout.put(c);
	}
	
	if(! is)
		return;
			
	while(is){
		pt = new InterpPt(is);
		if(is)
			appendInterpPt(pt);
	}
	is.clear();
	
	return;
}

InterpLst::~InterpLst()
{
	clearList();
}

void InterpLst::emptyList()
{
	firstInterpPt();
	while(numInterpPts>0)
		deleteCurrentInterpPt();
	delete [] InterpPtIndex;
	InterpPtIndex = NULL;
	indexState = lnotValid;
}

InterpPt* InterpLst::i(int n)
{
	if(indexState == lnotValid) {
		buildIndex();
	}
	if((n <= numInterpPts) && (n > 0))
		return (InterpPtIndex[n-1]);
	else
		return NULL;
}

int InterpLst::buildIndex()
{
	int index;

	delete [] InterpPtIndex;
	if(numInterpPts > 0) {
		InterpPtIndex = new InterpPt* [numInterpPts];
		InterpPtIndex[0] = firstInterpPt();
		theFirstInterpPt->setIndex(0);
		for(index=1;index<numInterpPts;index++){
			InterpPtIndex[index] = nextInterpPt();
			theCurrentInterpPt->setIndex(index);
		}
	}
	indexState = lvalid;
	return numInterpPts;
}

InterpPt* InterpLst::firstInterpPt()
{
	if(numInterpPts >0){
		theCurrentInterpPt = theFirstInterpPt;
		return(theFirstInterpPt);
	}
	else
		return(NULL);
}							
			
InterpPt* InterpLst::currentInterpPt()
{
	if(numInterpPts >0){
		return(theCurrentInterpPt);
	}
	else
		return(NULL);
}							
			
InterpPt* InterpLst::setCurrentInterpPt(InterpPt *InterpPtP)
{
	if(InterpPtP != NULL)
		theCurrentInterpPt = InterpPtP;
	else
		theCurrentInterpPt = theFirstInterpPt;
		
	return theCurrentInterpPt;
}							
			
InterpPt* InterpLst::nextInterpPt()
{
	if(numInterpPts > 0) {
		if(theCurrentInterpPt->getNextOne() != NULL){
			theCurrentInterpPt = (InterpPt *) theCurrentInterpPt->getNextOne();
			return (theCurrentInterpPt);
		}
		else
			return(NULL);
	}
	else
		return(NULL);
}
			
InterpPt* InterpLst::n(int name)
{
	if(numInterpPts > 0) {
		theCurrentInterpPt = theFirstInterpPt;
		while(theCurrentInterpPt->getn() != name){
			if(theCurrentInterpPt->getNextOne() != NULL){
				theCurrentInterpPt = (InterpPt *) theCurrentInterpPt->getNextOne();
			}
			else
				return(NULL);
		}
		return(theCurrentInterpPt);
	}
	else
		return(NULL);
}

InterpPt* InterpLst::appendInterpPt(InterpPt *theNewInterpPt)
{
	if(numInterpPts > 0) {
		while(theCurrentInterpPt->getNextOne() != NULL) {
			theCurrentInterpPt = (InterpPt *)theCurrentInterpPt->getNextOne();
		}
	}
	return (insertInterpPt(theNewInterpPt));
}
		
InterpPt* InterpLst::insertInterpPt(InterpPt *theNewInterpPt)
{	
	if(numInterpPts > 0) {
		theNewInterpPt->setNextOne(theCurrentInterpPt->getNextOne());
		theCurrentInterpPt->setNextOne(theNewInterpPt);
		theCurrentInterpPt = theNewInterpPt;
	}
	else {
		theNewInterpPt->setNextOne(NULL);
		theFirstInterpPt = theNewInterpPt;
		theCurrentInterpPt = theNewInterpPt;
	}
	numInterpPts += 1;
	indexState = lnotValid;
	
	return(theNewInterpPt);
}

void InterpLst::clearList()
{		
	InterpPt* nextInterpPtP;
	
	theCurrentInterpPt = theFirstInterpPt;
	
	while(numInterpPts > 0) {
		nextInterpPtP = (InterpPt *) theCurrentInterpPt->getNextOne();
		delete theCurrentInterpPt;
		theCurrentInterpPt = nextInterpPtP;
		numInterpPts -= 1;
	}
	theFirstInterpPt = NULL;
	delete [] InterpPtIndex;
	InterpPtIndex = NULL;
	indexState = lnotValid;
}

InterpPt* InterpLst::deleteCurrentInterpPt()
{
	InterpPt *prevInterpPt;

	prevInterpPt = theFirstInterpPt;

	if(theFirstInterpPt != NULL) {
		if(theCurrentInterpPt != theFirstInterpPt) {
			while(prevInterpPt->getNextOne() != theCurrentInterpPt) {
				prevInterpPt = (InterpPt *)prevInterpPt->getNextOne();
				if(prevInterpPt == NULL){
					return theCurrentInterpPt;
				}
			}
			prevInterpPt->setNextOne(theCurrentInterpPt->getNextOne());
			if(prevInterpPt->getNextOne() != NULL)
				theCurrentInterpPt = (InterpPt *) prevInterpPt->getNextOne();
			else
				theCurrentInterpPt = prevInterpPt;
		}
		else{
			theFirstInterpPt = (InterpPt *) theCurrentInterpPt->getNextOne();
			theCurrentInterpPt = theFirstInterpPt;
		}
		numInterpPts -= 1;
	}
	indexState = lnotValid;
	return (theCurrentInterpPt);
}

InterpPt* InterpLst::deleteInterpPt(InterpPt *iP)
{
	InterpPt *prevInterpPt = theFirstInterpPt;

	if(theFirstInterpPt == NULL)
		return NULL;

	if(iP == theFirstInterpPt) {
		theFirstInterpPt = (InterpPt *)iP->getNextOne();
		numInterpPts -= 1;
		indexState = lnotValid;
		if(iP == theCurrentInterpPt){
			theCurrentInterpPt = theFirstInterpPt;
		}
		return theCurrentInterpPt;
	}

	while(prevInterpPt != NULL) {
		if(iP == prevInterpPt->getNextOne()) {
			prevInterpPt->setNextOne(iP->getNextOne());
			numInterpPts -= 1;
			indexState = lnotValid;
			if(iP == theCurrentInterpPt){
				theCurrentInterpPt = prevInterpPt;
			}
			return theCurrentInterpPt;
		}
		prevInterpPt = (InterpPt *) prevInterpPt->getNextOne();
	}
	return theCurrentInterpPt;
}

double InterpLst::value(double xr)		
{
	InterpPt *pt1, *pt2;
	
	if((pt2 = firstInterpPt()) == NULL)
		return 0.0;
		
	if(xr <= pt2->getXc())
		return pt2->getYc();
		
	while(pt2 != NULL){
		pt1 = pt2;
		if( (pt2 = nextInterpPt()) == NULL)
			break;
		if(xr <= pt2->getXc()) {
			double r = (xr - pt1->getXc())/(pt2->getXc() - pt1->getXc());
			return (1-r) * pt1->getYc() + r * pt2->getYc();
		}
	}
	
	return pt1->getYc();
}

