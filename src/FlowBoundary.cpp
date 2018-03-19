//
//		FlowBoundary.cpp
//

#include "FlowBoundary.h"



FlowBoundary::FlowBoundary(int name, Node* start, Node* end, int code, double value, double value2, int tCode, string path):Item(name) //CString replaced by string May 14
{
	startNode = start;
	endNode = end;
	bcCode = code;
	bcValue = value;
	bcValue2 = value2;
	transCode = tCode;
	bcFilePathName = path;
	transLst = NULL;
	if (bcFilePathName.c_str() != "")
	{
		if (bcFilePathName[bcFilePathName.size()-1] == 'r')//GetLenght changed to size() May 14
			ratingCurve = true;
		else
			ratingCurve = false;
	}
	else 
		ratingCurve = false;
}	

void FlowBoundary::loadTransLst(ifstream &f)
{
	delete transLst;
	transLst = new BCInterpLst(f);
}
void FlowBoundary::deleteTransLst()
{
	delete transLst;
	transLst = NULL;
}

double FlowBoundary::getInterpValue(double timeOrQ)
{
	return (transLst->value(timeOrQ));
}

