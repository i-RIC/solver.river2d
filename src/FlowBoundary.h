
//		FlowBoundary.h


#ifndef FLOWBOUNDARY_H
#define FLOWBOUNDARY_H

#include "Item.h"
#include "BCInterpLst.h"
//#include "afx.h" check if this belongs to windows
#include "Node.h"

class FlowBoundary : public Item
{		
	protected:
		int bcCode;
		double bcValue;
		double bcValue2;
		Node* startNode;
		Node* endNode;
		int transCode;
		bool ratingCurve;
		string bcFilePathName;//changed CString by string May 14
		BCInterpLst *transLst;
		 
	public:
		FlowBoundary(int name = 1, Node* start = NULL, Node* end = NULL, int code = 0, double value = 0.0, double value2 = 0.0, int tCode = 0, string path="" );//CString changed by string May 14
		Node* getStartNode() {return startNode;}
		void setStartNode(Node* start) {startNode = start;}
		Node* getEndNode() {return endNode;}
		void setEndNode(Node* end) {endNode = end;}
		int getBcCode() {return bcCode;}
		void setBcCode(int code) {bcCode = code;}
		int getTransCode() {return transCode;}
		void setTransCode(int trans) {transCode = trans;}
		double getBcValue() {return bcValue;}
		void setBcValue(double value) {bcValue = value;}
		double getBcValue2() {return bcValue2;}
		void setBcValue2(double value2) {bcValue2 = value2;}
		string getFilePath() {return bcFilePathName;} //CString changer by string
		void setFilePath(string path) {bcFilePathName = path;}//CString changed by string
		bool getYesNoRatingCurve() {return ratingCurve;}
		void setYesNoRatingCurve(bool yesNoCurve) {ratingCurve = yesNoCurve;}
		BCInterpLst* gettransLst() {return transLst;}
		void loadTransLst(ifstream &f);
		void deleteTransLst();
		double getInterpValue(double timeOrQ);

			
};

#endif
