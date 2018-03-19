
//   InterpLst Class


#ifndef INTERPLST_H
#define INTERPLST_H

#include <fstream>
#include "Interppt.h"
using namespace std;

enum lstStatus {lvalid, lnotValid};
	
class InterpLst
{
	
	protected:
		InterpPt	*theFirstInterpPt;
		InterpPt	*theCurrentInterpPt;
		int		numInterpPts;
		InterpPt	**InterpPtIndex;
		lstStatus	indexState;
		
			
	public:
		InterpLst();								//	construcor
		InterpLst(istream& is);
		~InterpLst();							//	destructor
		void emptyList();						//	does not delete InterpPts
		int numberOfInterpPts() {return numInterpPts;}		
		InterpPt* i(int i);							//	return pointer to ith InterpPt
		int buildIndex();
		InterpPt* firstInterpPt();						//	get first InterpPt in list
		InterpPt* currentInterpPt();					//	get current InterpPt in list
		InterpPt* setCurrentInterpPt(InterpPt *InterpPtP);		//	set current InterpPt in list
		InterpPt* nextInterpPt();						//	get next InterpPt in list
		InterpPt* n(int name);						//	get InterpPt with n = name
		InterpPt* appendInterpPt(InterpPt *theNewInterpPt);		//	add new InterpPt to end of list
		InterpPt* insertInterpPt(InterpPt *theNewInterpPt);		//	add new InterpPt after current InterpPt
		void clearList();						//	delete (really delete) all InterpPts in list
		InterpPt* deleteCurrentInterpPt();				//	remove current InterpPt, next becomes current
		InterpPt* deleteInterpPt(InterpPt *iP);		// delete InterpPt pointed at, return currentInterpPt
		double value(double xr);		// return interpolated value
};

#endif
