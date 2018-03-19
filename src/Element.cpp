//
//		Element.cpp
//

#include "Element.h"

ostream& operator << (ostream& os, Element* seg)
{
	os.width(5);
	os 	<< seg->getn() ;
	
	if(seg != NULL){
	
		for(int i=0;i<seg->numNodes();i++){
			os.width(10);
			if(seg->getNode(i) != NULL)
				os 	<<  (seg->getNode(i))->getn() ;
			else
				os << -1;
		}
/*
		for(int j=0;j<seg->numAdj();j++){
			os.width(10);
			if(*(seg->getAdj()+j) != NULL)
				os 	<<  (*(seg->getAdj()+j))->getn() ;
			else
				os << -1;
		}
*/
		os << "\n";
	}
	return os;
}

int	Element::reflectAdj(Element* ap)
{
	for(int j=numAdj()-1;j>=0;j--){
		if(*(getAdj()+j) == ap)
			return j;
	}
	return -1;
}

Node Element::centroid()
{
	double x = 0.0, y = 0.0;

	for(int i=0;i<numNodes();i++){
		x += getNode(i)->getXc();
		y += getNode(i)->getYc();
	}
	Node nd(1,x/numNodes(),y/numNodes(),0.0);

	return nd;
}

