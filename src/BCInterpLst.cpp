// BCInterpLst.cpp: implementation of the BCInterpLst class.
//
//////////////////////////////////////////////////////////////////////

#include "BCInterpLst.h"
#include "BCInterpPt.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
InterpPt* appendInterpPt(InterpPt *theNewInterpPt);//Added by Christian Frias May 14


BCInterpLst::BCInterpLst()
{

}

BCInterpLst::BCInterpLst(istream &is):InterpLst()
{
	char c = 'a';
	BCInterpPt *pt;
	int i = 0;

	while( is && (c != '(')){
		is.get(c);
//		cout.put(c);
	}
	
	if(! is)
		return;
			
	while(is){
		pt = new BCInterpPt(is);
		if(is)
		{
			i++;
			pt->setn(i);
			appendInterpPt(pt);
		}
	}
	is.clear();
	
	return;
}

BCInterpLst::~BCInterpLst()
{
	clearList();
}

