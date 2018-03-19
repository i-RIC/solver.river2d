
//		Quadrant.cpp

//		Quadrant Class method definitions

#include "Quadrant.h"

Quadrant::Quadrant(box limits)
{	
	subQuadrant[0]=NULL;
	subQuadrant[1]=NULL;
	subQuadrant[2]=NULL;
	subQuadrant[3]=NULL;
	nNodes = nTriangles = 0;
	firstNode = lastNode = NULL;
	firstTriangle = lastTriangle = NULL;
	bArc = NULL;
	boundBox = limits;
}

void Quadrant::appendNode(Node *newNode)
{
	if(newNode != NULL) {
		if(firstNode == NULL){
			firstNode = newNode;
		}
		else {
			lastNode->setNextOne(newNode);
		}
		lastNode = newNode;
		newNode->setNextOne(NULL);
		nNodes += 1;
	}
	return;
}

void Quadrant::appendTriangle(Triangle *newTri)
{
	if(newTri != NULL) {
		if(firstTriangle == NULL){
			firstTriangle = newTri;
			newTri->ptP = NULL;

		}
		else {
			lastTriangle->ntP = newTri;
			newTri->ptP = lastTriangle;
		}
		newTri->ntP = NULL;
		lastTriangle = newTri;
		nTriangles += 1;
	}
	return;
}

void Quadrant::deleteTriangle(Triangle *tP)
{
	if(tP->ntP != NULL){
		tP->ntP->ptP = tP->ptP;
	}
	else{
		lastTriangle = tP->ptP;
	}
	if(tP->ptP != NULL){
		tP->ptP->ntP = tP->ntP;
	}
	else{
		firstTriangle = tP->ntP;
	}
	tP->ntP = NULL;
	tP->ptP = NULL;
	nTriangles += 1;
	return;
}

int Quadrant::triangulateNodes()
{
	if(nNodes == 0){							// Make an arc to point nowhere
		bArc = NULL;
	}
	if(nNodes == 3){					
		Node *np1 = firstNode;
		Node *np2 = (Node *) np1->getNextOne();
		Node *np3 = (Node *) np2->getNextOne();
/*
		double minDist = 0.0001;		// Check for duplicate nodes
		if(np1->dist(np2) < minDist) {
			delete np2;
			np2 = NULL;
			np1->setNextOne(np3);
			nNodes--;
		}
		if(np1->dist(np3) < minDist) {
			delete np3;
			np3 = NULL;
			np2->setNextOne(NULL);
			lastNode = np2;
			nNodes--;
		}
		if((np2 != NULL) && (np1 != NULL)){
			if(np2->dist(np3) < minDist) {
				delete np3;
				np3 = NULL;
				np2->setNextOne(NULL);
				lastNode = np2;
				nNodes--;
			}
		}
		if(nNodes == 3){
*/			Triangle *tp = new Triangle(1,np1,np2,np3);	// Make a triangle surrounded by arcs
			firstTriangle = lastTriangle = tp;
			tp->ntP = tp->ptP = NULL;
			double a = tp->area();
			if(a < 0.0) {
				tp->setNode(1,np3);
				tp->setNode(2,np2);
			}
			Triangle *arc0 = new Triangle(1,tp->getNode(0),tp->getNode(1),NULL);
			Triangle *arc1 = new Triangle(1,tp->getNode(1),tp->getNode(2),NULL);
			Triangle *arc2 = new Triangle(1,tp->getNode(2),tp->getNode(0),NULL);
			tp->setAdj(0,arc1);
			tp->setAdj(1,arc2);
			tp->setAdj(2,arc0);
			arc0->setAdj(0,arc1);
			arc0->setAdj(1,arc2);
			arc0->setAdj(2,tp);
			arc1->setAdj(0,arc2);
			arc1->setAdj(1,arc0);
			arc1->setAdj(2,tp);
			arc2->setAdj(0,arc0);
			arc2->setAdj(1,arc1);
			arc2->setAdj(2,tp);
			bArc = arc1;
///		}
	}

	if(nNodes == 2){							
		Node *np1 = firstNode;
		Node *np2 = (Node *) np1->getNextOne();

/*		double minDist = 0.0001;		// Check for duplicate nodes
		if(np1->dist(np2) < minDist) {
			delete np2;
			np2 = NULL;
			np1->setNextOne(NULL);
			lastNode = np1;
			nNodes--;
		}
		if(nNodes == 2){
*/									// Make a double arc between 2 nodes
			Triangle *arc1 = new Triangle(1,np1,np2,NULL);
			Triangle *arc2 = new Triangle(1,np2,np1,NULL);
			arc1->setAdj(0,arc2);
			arc1->setAdj(1,arc2);
			arc1->setAdj(2,arc2);
			arc2->setAdj(0,arc1);
			arc2->setAdj(1,arc1);
			arc2->setAdj(2,arc1);
			bArc = arc1;
//		}
	}

	if(nNodes == 1){							// Make an arc to point to node
		Node *np1 = firstNode;
		Triangle *arc1 = new Triangle(1,np1,np1,NULL);
		arc1->setAdj(0,arc1);
		arc1->setAdj(1,arc1);
		bArc = arc1;
	}

		
	 
	if(nNodes > 3){								//  Subdivide
		double ym = (boundBox.y1 + boundBox.y2)/2.;
		double xm = (boundBox.x1 + boundBox.x2)/2.;

		box qBox;								// NE subquadrant
		qBox.x1 = xm; qBox.x2 = boundBox.x2; 
		qBox.y1 = ym; qBox.y2 = boundBox.y2;
		subQuadrant[0] = new Quadrant(qBox);
												// NW subquadrant
		qBox.x1 = boundBox.x1; qBox.x2 = xm; 
		qBox.y1 = ym; qBox.y2 = boundBox.y2;
		subQuadrant[1] = new Quadrant(qBox);
												// SW subquadrant
		qBox.x1 = boundBox.x1; qBox.x2 = xm; 
		qBox.y1 = boundBox.y1; qBox.y2 = ym;
		subQuadrant[2] = new Quadrant(qBox);
												// SE subquadrant
		qBox.x1 = xm; qBox.x2 = boundBox.x2; 
		qBox.y1 = boundBox.y1; qBox.y2 = ym;
		subQuadrant[3] = new Quadrant(qBox);
												// Distribute Nodes
		Node *nP = firstNode;
		Node *nextNode;
		while(nP != NULL) {
			nextNode = (Node*) nP->getNextOne();
			if(nP->getYc() >= ym) {
				if(nP->getXc() >= xm) {
					subQuadrant[0]->appendNode(nP);
				}
				else {
					subQuadrant[1]->appendNode(nP);
				}
			}
			else {
				if(nP->getXc() < xm)
					subQuadrant[2]->appendNode(nP);
				else
					subQuadrant[3]->appendNode(nP);
			}
			nP = nextNode;
		}
		nNodes = 0;
		firstNode = lastNode = NULL;

												// Recurse triangulation
		for(int i=0;i<4;i++) {
			subQuadrant[i]->triangulateNodes();
		}
		int NNE = recoverSubQuadInfo(0);
		Triangle *bArcE = subQuadrant[0]->getBArcs();
		int NNW = recoverSubQuadInfo(1);
		Triangle *bArcW = subQuadrant[1]->getBArcs();
		Triangle *bArcN;
		int NNN = NNE + NNW;
		if(NNE == 0){
			bArcN = bArcW;
		}
		else if(NNW == 0){
			bArcN = bArcE;
		}
		else{
			bArcN = stitchLR(bArcW,bArcE);
		}

		NNW = recoverSubQuadInfo(2);
		bArcW = subQuadrant[2]->getBArcs();
		NNE = recoverSubQuadInfo(3);
		bArcE = subQuadrant[3]->getBArcs();
		Triangle *bArcS;
		int NNS = NNE + NNW;
		if(NNE == 0){
			bArcS = bArcW;
		}
		else if(NNW == 0){
			bArcS = bArcE;
		}
		else{
			bArcS = stitchLR(bArcW,bArcE);
		}

		int NN = NNN + NNS;
		if(NNN == 0){
			bArc = bArcS;
		}
		else if(NNS == 0){
			bArc = bArcN;
		}
		else{
			bArc = stitchLR(bArcN,bArcS);
		}
	}
	for(int i=0;i<4;i++) {
		delete subQuadrant[i];
	}

	return nNodes;
}

int Quadrant::recoverSubQuadInfo(int i)
{
	Node *nP = subQuadrant[i]->getNodes();
	Node *nextNP;
	while(nP != NULL) {
		nextNP = (Node *)nP->getNextOne();
		appendNode(nP);
		nP = nextNP;
	}

	Triangle *tP = subQuadrant[i]->getTriangles();
	Triangle *nextTP;
	while(tP != NULL) {
		nextTP = tP->ntP;
		appendTriangle(tP);
		tP = nextTP;
	}
	return subQuadrant[i]->getNNodes();
}

Triangle* Quadrant::stitchLR(Triangle *bArcL, Triangle *bArcR)
{
	Triangle *topArcL = bArcL;
	Triangle *topArcR = bArcR;
	Triangle *topArc = new Triangle(1,topArcR->node0(),topArcL->node1(),NULL);
	int doneFlag = -1;
	while(doneFlag <1){
		doneFlag = 1;
		if(topArcL->node0() != topArcL->node1()){
			topArc->setNode(2,topArcL->node0());
			if(topArc->area() <= 0.0){
				topArcL = topArcL->tP1();
				topArc->setNode(1,topArcL->node1());
				doneFlag = -1;
			}
			else {
				topArc->setNode(2,topArcL->tP0()->node1());		
				if(topArc->area() < 0.0){
					topArcL = topArcL->tP0();
					topArc->setNode(1,topArcL->node1());
					doneFlag = -1;
				}
			}
		}
		if(topArcR->node0() != topArcR->node1()){
			topArc->setNode(2,topArcR->tP1()->node0());
			if(topArc->area() < 0.0){
				topArcR = topArcR->tP1();
				topArc->setNode(0,topArcR->node0());
				doneFlag = -1;
			}
			else{
				topArc->setNode(2,topArcR->node1());
				if(topArc->area() <= 0.0){
					topArcR = topArcR->tP0();
					topArc->setNode(0,topArcR->node0());
					doneFlag = -1;
				}
			}
		}
	}

	topArc->setNode(2,NULL);
	topArc->setAdj(0,topArcL->tP0());
	topArc->setAdj(1,topArcR->tP1());
	topArcL->tP0()->setAdj(1,topArc);
	topArcR->tP1()->setAdj(0,topArc);
	Triangle *botArc = new Triangle(1,topArcL->node1(),topArcR->node0(),NULL);
	topArc->setAdj(2,botArc);
	botArc->setAdj(2,topArc);
	botArc->setAdj(0,topArcR);
	botArc->setAdj(1,topArcL);
	topArcL->setAdj(0,botArc);
	topArcR->setAdj(1,botArc);
	Triangle *cwArc = topArcL;
	Triangle *ccwArc = topArcR;

	int snode = 0;
	if(topArcR->node0() == topArcR->node1()){
		topArc->setAdj(1,botArc);
		botArc->setAdj(0,topArc);
		ccwArc = topArc;
		delete topArcR;
		snode += 1;
	}
	if(topArcL->node0() == topArcL->node1()){
		topArc->setAdj(0,botArc);
		botArc->setAdj(1,topArc);
		cwArc = topArc;
		delete topArcL;
		snode += 2;
	}
	if(snode == 3)
		return topArc;

	Node *nP1 = ccwArc->node1();
	Node *nP2 = botArc->node1();
	Node *nP3 = botArc->node0();
	Node *nP4 = cwArc->node0();

	double areaL, areaR;
	int	goLeft;
	Triangle *tP;

	while(-1){
		areaL = areaR = -1.0;
		goLeft = 0;
		tP = new Triangle(1,nP4,nP2,nP3);
		if(snode != 2){
			areaL = tP->area();
		}
		tP->setNode(0,nP1);
		if(snode != 1){
			areaR = tP->area();
		}
		if((areaL <= 0.0) && (areaR <= 0.0)){
			delete tP;
			break;
		}
		if((areaL > 0.0) && (areaR <= 0.0)){
			goLeft = 1;
		}
		if((areaL <= 0.0) && (areaR > 0.0)){
			goLeft = -1;
		}
		if((areaL > 0.0) && (areaR > 0.0)){
			if(tP->insideC(nP4) >= 0.0){
				goLeft = 1;
			}
			else{
				goLeft = -1;
			}
		}
		if(goLeft > 0) {
			tP->setNode(0,nP4);
			tP->setAdj(0,botArc->getAdjt(2));
			botArc->getAdjt(2)->setAdj(botArc->tsn(),tP);
			tP->setAdj(1,cwArc->getAdjt(2));
			if(cwArc->getAdjt(2) != cwArc){
				cwArc->getAdjt(2)->setAdj(cwArc->tsn(),tP);
			}
			else{
				cwArc->tP1()->setAdj(2,tP);
			}
			tP->setAdj(2,botArc);
			botArc->setAdj(2,tP);
			botArc->setAdj(1,cwArc->tP1());
			botArc->setNode(0,cwArc->node0());
			cwArc->tP1()->setAdj(0,botArc);
			delete cwArc;
			delaunay(tP);
			cwArc = botArc->tP1();
			nP3 = botArc->node0();
			nP4 = cwArc->node0();
		}
		else{
			tP->setAdj(0,botArc->getAdjt(2));
			botArc->getAdjt(2)->setAdj(botArc->tsn(),tP);
			tP->setAdj(2,ccwArc->getAdjt(2));
			if(ccwArc->getAdjt(2) != ccwArc){
				ccwArc->getAdjt(2)->setAdj(ccwArc->tsn(),tP);
			}
			else{
				ccwArc->tP0()->setAdj(2,tP);
			}
			tP->setAdj(1,botArc);
			botArc->setAdj(2,tP);
			botArc->setAdj(0,ccwArc->tP0());
			botArc->setNode(1,ccwArc->getNode(1));
			ccwArc->tP0()->setAdj(1,botArc);
			delete ccwArc;
			delaunay(tP);
			ccwArc = botArc->tP0();
			nP2 = botArc->node1();
			nP1 = ccwArc->node1();
		}
	}

	return topArc;
}

void Quadrant::delaunay(Triangle *tP)
{
	int	index, changed;
	Triangle *atP, *ntP;
	Triangle *cStack = NULL;
	tP->deactivate();
	cStack = pushTri(cStack,tP);
	while(cStack != NULL){
		for(int i=0;i<3;i++){
			changed = 0;
			atP = (Triangle *) cStack->getAdjt(i);
			if(atP != NULL){
				index = atP->reflectAdj(cStack);
				if(cStack->insideC(atP->getNode(index)) > 0.0){
					swapDiagonal(cStack,atP);
					if(atP->getStatus() != notActive){
						atP->deactivate();
						deleteTriangle(atP);
						cStack = pushTri(cStack,atP);
					}
					changed = 1;
					break;
				}
			}
		}
		if(changed == 0){
			ntP = cStack->ntP;
			appendTriangle(cStack);
			cStack->activate();
			cStack = ntP;
		}
	}

	return;
}

Triangle* Quadrant::pushTri(Triangle *stack, Triangle *tP)
{
	if(tP != NULL){
		tP->ntP = stack;
		if(stack != NULL)
			stack->ptP = tP;
		return tP;
	}
	else{
		return stack;
	}
}

Triangle* Quadrant::popTri(Triangle *stack)
{

	return stack->ntP;
}

void Quadrant::swapDiagonal(Triangle *tP1, Triangle *tP2)
{
	Node *nP0, *nP1, *nP2, *nP3;
	Element *pP1, *pP2, *pP3;

	int index = tP2->reflectAdj(tP1);		// Find number of common edge in second triangle
	int indexp = index + 1;
	if(indexp == 3)							// Next edge number
		indexp = 0;
	int indexm = index - 1;
	if(indexm == -1)						// Previous edge number
		indexm = 2;

	int index1 = tP1->reflectAdj(tP2);		// Find number of common edge in first triangle
	int index1p = index1 + 1;
	if(index1p == 3)							// Next edge number
		index1p = 0;
	int index1m = index1 - 1;
	if(index1m == -1)						// Previous edge number
		index1m = 2;

	nP0 = tP1->getNode(index1);					// The node that started it all
	nP1 = tP1->getNode(index1p);
	nP2 = tP1->getNode(index1m);					// Save 4 nodes in quadrilateral
	nP3 = tP2->getNode(index);

	pP1 = tP1->getAdjt(index1m);
	if( pP1 != NULL) {
		index = pP1->reflectAdj(tP1);		// Adjust pointers from adjacent triangles
		pP1->setAdj(index, tP2);
	}
 
	pP2 = tP2->getAdjt(indexm);
	if( pP2 != NULL) {
		index = pP2->reflectAdj(tP2);
		pP2->setAdj(index, tP1);
	}

	pP3 = tP1->getAdjt(index1p);

	tP1->setNode(0,nP0);
	tP1->setNode(1,nP3);
	tP1->setNode(2,nP2);
	tP1->setAdj(0,pP2);	// Setup new triangle 1
	tP1->setAdj(1,pP3);
	tP1->setAdj(2,tP2);

	tP2->setNode(0,nP0);
	tP2->setNode(1,nP1);
	tP2->setNode(2,nP3);
	tP2->setAdj(0,tP2->getAdjt(indexp));	// Set up new triangle 2
	tP2->setAdj(1,tP1);
	tP2->setAdj(2,pP1);
}

int Quadrant::checkDuplicateNodes()
{
	double minDist = 0.0001;
	
	if(nNodes == 0 || nNodes == 1){	
		return nNodes;
	}

	if(nNodes == 2){							
		Node *np1 = firstNode;
		Node *np2 = (Node *) np1->getNextOne();

		if(np1->dist(np2) < MINDIST) {
			np2->setDeleted();
		}
		return nNodes;
	}
	 
	if(nNodes > 2){								
		if(boundBox.y2 - boundBox.y1 < minDist) {
			if(boundBox.x2 - boundBox.x1 < minDist){
				Node *np = firstNode;
				while(np != NULL){
					np->setDeleted();
					np = (Node *) np->getNextOne();
				}
			return nNodes;
			}
		}

		double ym = (boundBox.y1 + boundBox.y2)/2.;
		double xm = (boundBox.x1 + boundBox.x2)/2.;

		box qBox;								// NE subquadrant
		qBox.x1 = xm; qBox.x2 = boundBox.x2; 
		qBox.y1 = ym; qBox.y2 = boundBox.y2;
		subQuadrant[0] = new Quadrant(qBox);
												// NW subquadrant
		qBox.x1 = boundBox.x1; qBox.x2 = xm; 
		qBox.y1 = ym; qBox.y2 = boundBox.y2;
		subQuadrant[1] = new Quadrant(qBox);
												// SW subquadrant
		qBox.x1 = boundBox.x1; qBox.x2 = xm; 
		qBox.y1 = boundBox.y1; qBox.y2 = ym;
		subQuadrant[2] = new Quadrant(qBox);
												// SE subquadrant
		qBox.x1 = xm; qBox.x2 = boundBox.x2; 
		qBox.y1 = boundBox.y1; qBox.y2 = ym;
		subQuadrant[3] = new Quadrant(qBox);
												// Distribute Nodes
		Node *nP = firstNode;
		Node *nextNode;
		while(nP != NULL) {
			nextNode = (Node*) nP->getNextOne();
			if(nP->getYc() >= ym) {
				if(nP->getXc() >= xm) {
					subQuadrant[0]->appendNode(nP);
				}
				else {
					subQuadrant[1]->appendNode(nP);
				}
			}
			else {
				if(nP->getXc() < xm)
					subQuadrant[2]->appendNode(nP);
				else
					subQuadrant[3]->appendNode(nP);
			}
			nP = nextNode;
		}
		nNodes = 0;
		firstNode = lastNode = NULL;

												// Recursion
		for(int i=0;i<4;i++) {
			subQuadrant[i]->checkDuplicateNodes();
		}
		for(int i=0;i<4;i++) { // Added int keyword before i
			recoverSubQuadInfo(i);
			delete subQuadrant[i];
		}
	}

	return nNodes;
}

