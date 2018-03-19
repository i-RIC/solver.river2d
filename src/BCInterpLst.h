// BCInterpLst.h: interface for the BCInterpLst class.
//
//////////////////////////////////////////////////////////////////////

#ifndef BCINTERPLST_H
#define BCINTERPLST_H


//#if _MSC_VER > 1000
//#pragma once
//#endif // _MSC_VER > 1000

#include "Intrplst.h"

class BCInterpLst : public InterpLst  
{
public:
	BCInterpLst();
	BCInterpLst(istream &is);
	~BCInterpLst();

};

#endif 

