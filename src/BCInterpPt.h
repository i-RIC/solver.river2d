// BCInterpPt.h: interface for the BCInterpPt class.
//
//////////////////////////////////////////////////////////////////////

#ifndef BCINTERPPT_H
#define BCINTERPPT_H

//#if _MSC_VER > 1000
//#pragma once
//#endif // _MSC_VER > 1000

#include "Interppt.h"

class BCInterpPt : public InterpPt  
{
public:
	BCInterpPt();
	BCInterpPt(istream &is);

};

#endif 
