//
//		TIN.cpp
//

#include "Tin.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <ctype.h>
using namespace std;

TIN::~TIN()
{
	Node *nP = (Node *) nodeL.firstItem();
	while(nP != NULL){
		nodeL.deleteCurrentItem();
		delete nP;
		nP = (Node *) nodeL.firstItem();
	}

	nP = (Node *) tnodeL.firstItem();
	while(nP != NULL){
		tnodeL.deleteCurrentItem();
		delete nP;
		nP = (Node *) tnodeL.firstItem();
	}

	nP = (Node *) deadNodeL.firstItem();
	while(nP != NULL){
		deadNodeL.deleteCurrentItem();
		delete nP;
		nP = (Node *) deadNodeL.firstItem();
	}

	Triangle *tP = (Triangle *) elmL.firstItem();
	while(tP != NULL){
		elmL.deleteCurrentItem();
		delete tP;
		tP = (Triangle *) elmL.firstItem();
	}

	Segment *eP = (Segment *) edgeL.firstItem();
	while(eP != NULL){
		edgeL.deleteCurrentItem();
		delete eP;
		eP = (Segment *) edgeL.firstItem();
	}

	Segment *sP = (Segment *) bsegL.firstItem();
	while(sP != NULL){
		bsegL.deleteCurrentItem();
		delete sP;
		sP = (Segment *) bsegL.firstItem();
	}

	sP = (Segment *) fsegL.firstItem();
	while(sP != NULL){
		fsegL.deleteCurrentItem();
		delete sP;
		sP = (Segment *) fsegL.firstItem();
	}

	sP = (Segment *) tfsegL.firstItem();
	while(sP != NULL){
		tfsegL.deleteCurrentItem();
		delete sP;
		sP = (Segment *) tfsegL.firstItem();
	}

	sP = (Segment *) tbsegL.firstItem();
	while(sP != NULL){
		tbsegL.deleteCurrentItem();
		delete sP;
		sP = (Segment *) tbsegL.firstItem();
	}

}

int TIN::getNextN()
{
	serialNum += 1;
	return serialNum;
}

int TIN::readNodes(istream& is)
{
	Node *nodeP=NULL, *pnodeP=NULL, *tnodeP=NULL, *snodeP=NULL;// *bNodeP;
	Segment *segP=NULL;
	char c = '\t';
	int nf = -1, bFlag = 0;
	int nMax = 0;

	while((c != '.')){
		c = is.peek();
		if(isdigit(c)){
			tnodeP = nodeP;
			nodeP = physics->readInNode(is);
			is.get(c);
			while(c != '\n'){
				is.get(c);
			}
			if((nodeP->getn() >= 0)/* && (checkDupNode(nodeP) == 0)*/){
				pnodeP = tnodeP;
				nodeL.appendItem(nodeP);
				if(nodeP->getn() > serialNum)
					serialNum = nodeP->getn();
				if(nf == 0){
					snodeP = nodeP;
					if(bFlag == 1){
//						bNodeP = new Node(nf+1,nodeP->getXc(),nodeP->getYc());
//						tnodeL.appendItem(bNodeP);
					}
					nf++;
				}
				else if(nf > 0){
					segP = new Segment(fsegL.numberOfItems()+1,pnodeP,nodeP);
					addFeatureSeg(segP,fsegL);
					if(bFlag == 1){
//						bNodeP = new Node(nf+1,nodeP->getXc(),nodeP->getYc());
//						tnodeL.appendItem(bNodeP);
						segP = physics->makeNewSegment(bsegL.numberOfItems()+1,pnodeP,nodeP);
						bsegL.appendItem(segP);
					}
					nf++;
				}
			}
			else {
				is.clear();
				delete nodeP;
				nodeP = tnodeP;
			}
		}
		else {
			is.get(c);
			if( (c == '(') || (c == '[') ) {
				nf = 0;
			}
			if( c == '{' ) {
				nf = 0;
				bFlag = 1;
			}
			if( c == ')' ){
				nf = -1;
			}
			if( (c == ']') || (c == '}') ) {
				segP = new Segment(fsegL.numberOfItems()+1,nodeP,snodeP);
				addFeatureSeg(segP,fsegL);
				nf = -1;
				if(c == '}'){
//					bNodeP = new Node(nf+1,snodeP->getXc(),snodeP->getYc());
//					tnodeL.appendItem(bNodeP);
					segP = physics->makeNewSegment(bsegL.numberOfItems()+1,nodeP,snodeP);
					bsegL.appendItem(segP);
					bFlag = 0;
				}
			}
		}
	}
//	is.clear();
//	while((c != '.') && is) {
//		is.get(c);
//		cout << c;
//	}
//	serialNum = Nnodes();
	return nodeL.numberOfItems();
}

int TIN::readNodesDes(istream& is)
{
	Node *nodeP=NULL, *pnodeP=NULL, *tnodeP=NULL, *snodeP=NULL;// *bNodeP;
	Segment *segP=NULL;
	char c = '\t';
	int nf = -1, bFlag = 0;
	int nMax = 0;
	char strbuf[10] = "";

	while((c != '.')){
		c = is.peek();
		if(isdigit(c)){
			tnodeP = nodeP;
			nodeP = physics->readInNode(is);
			is.get(c);
			if(c != '\n'){
				is.get(strbuf,21,'\n');
				nodeP->setDesignation(strbuf);
				is.get(c);
				while(c != '\n'){
					is.get(c);
				}
			}
			if((nodeP->getn() >= 0)/* && (checkDupNode(nodeP) == 0)*/){
				pnodeP = tnodeP;
				nodeL.appendItem(nodeP);
				if(nodeP->getn() > serialNum)
					serialNum = nodeP->getn();
				if(nf == 0){
					snodeP = nodeP;
					if(bFlag == 1){
					}
					nf++;
				}
				else if(nf > 0){
					segP = new Segment(fsegL.numberOfItems()+1,pnodeP,nodeP);
					addFeatureSeg(segP,fsegL);
					if(bFlag == 1){
						segP = physics->makeNewSegment(bsegL.numberOfItems()+1,pnodeP,nodeP);
						bsegL.appendItem(segP);
					}
					nf++;
				}
			}
			else {
				is.clear();
				delete nodeP;
				nodeP = tnodeP;
			}
		}
		else {
			is.get(c);
			if( (c == '(') || (c == '[') ) {
				nf = 0;
			}
			if( c == '{' ) {
				nf = 0;
				bFlag = 1;
			}
			if( c == ')' ){
				nf = -1;
			}
			if( (c == ']') || (c == '}') ) {
				segP = new Segment(fsegL.numberOfItems()+1,nodeP,snodeP);
				addFeatureSeg(segP,fsegL);
				nf = -1;
				if(c == '}'){
					segP = physics->makeNewSegment(bsegL.numberOfItems()+1,nodeP,snodeP);
					bsegL.appendItem(segP);
					bFlag = 0;
				}
			}
		}
	}
	return nodeL.numberOfItems();
}

int TIN::checkDupNode(Node* pnP)
{
	int dupFlag = 0;

	Node *nP = firstNode();
	while(nP != NULL) {
		if( fabs(nP->getXc() - pnP->getXc()) < 0.00001){
			if( fabs(nP->getYc() - pnP->getYc()) < 0.00001){
				dupFlag = -1;
				break;
			}
		}
      nP = nextNode();
	}

	return dupFlag;
}

//	New function added 5/2001, JDS
//	This function checks all nodes to see if they still belong in the 
//	dataset after a boundary redefine.  I originally tried the boundary
//	redefine letting the system handle the nodes, but there was a problem 
//	removing the first node from the list that resulted in a ghost node.

void TIN::checkAllNodes()
{
	Triangle *tP;
	Node *nP = (Node *)nodeL.firstItem();
	nodeL.setCurrentItem(nP);
	while(nP!=NULL){
		tP = whichTriHasNode(nP);
		if(tP==NULL || (tP->getStatus()==notActive)){
			nP->setFloating();
			deleteNode(nP);
			nodeL.deleteCurrentItem();
			nP = (Node *)nodeL.currentItem();
		}
		else{
			nP = (Node *)nodeL.nextItem();
		}
	}
}

int TIN::checkDupNodeTri(Node* pnP, Triangle* tP)
{
	int dupFlag = 0;

	Node *nP;
	for(int i=0;i<3;i++) {
		nP = tP->getNode(i);
		if( fabs(nP->getXc() - pnP->getXc()) < 0.00001){
			if( fabs(nP->getYc() - pnP->getYc()) < 0.00001){
				return (-1);
			}
		}
	}

	return (0);
}
	
void TIN::checkAllFeatureSegs()
{
	Triangle *triP1;
	Triangle *triP2;
	Segment *segP = firstfSeg();
	while(segP != NULL){
		triP1 = whichTriangle(segP->getNode(0));
		triP2 = whichTriangle(segP->getNode(1));
		if( (triP1 == NULL ) || (triP1->getStatus() == notActive) ||
			 (triP2 == NULL ) || (triP2->getStatus() == notActive) ||
			(checkFeatureSeg(segP) == -1) ){
				fsegL.deleteCurrentItem();
				delete segP;
				segP = (Segment *) fsegL.currentItem();
		}
		else
			segP = nextfSeg();
	}
	segP = firstTfSeg();
	while(segP != NULL){
		triP1 = whichTriangle(segP->getNode(0));
		triP2 = whichTriangle(segP->getNode(1));
		if( (triP1 == NULL ) || (triP1->getStatus() == notActive) ||
			 (triP2 == NULL ) || (triP2->getStatus() == notActive) ||
			 (checkFeatureSeg(segP) == -1) ){
				tfsegL.deleteCurrentItem();
				delete segP;
				segP = (Segment *) tfsegL.currentItem();
		}
		else
			segP = nextTfSeg();
	}
	return;
}

Segment* TIN::addFeatureSeg(Segment* newfsegp, ItemList& fList)
{
	int result = checkFeatureSeg(newfsegp);
	if(result != -1){
		fList.appendItem(newfsegp);
		return newfsegp;
	}
	else{
		delete newfsegp;
		return NULL;
	}
}

Segment* TIN::addTFeatureSeg(Segment* newfsegp)
{
	return addFeatureSeg(newfsegp,tfsegL);
}

int TIN::checkFeatureSeg(Segment* newfsegp)
{
	int i1 = 1, i2 = 1;

	if((newfsegp->getNode(0)->getFixed() == deleted) ||
		(newfsegp->getNode(1)->getFixed() == deleted)) {
		return -1;
	}

	Segment* fsegp = firstfSeg();
	while( (fsegp != NULL) && (fsegp != newfsegp)){
		i1 = newfsegp->intersect(fsegp);
		i2 = fsegp->intersect(newfsegp);
		if( (i1 == -1) && (i2 == -1) )
			return -1;
		fsegp = nextfSeg();
	}
	fsegp = firstTfSeg();
	while( (fsegp != NULL) && (fsegp != newfsegp)){
		i1 = newfsegp->intersect(fsegp);
		i2 = fsegp->intersect(newfsegp);
		if( (i1 == -1) && (i2 == -1) )
			return -1;
		fsegp = nextTfSeg();
	}
	fsegp = firstSeg();
	while( (fsegp != NULL) && (fsegp != newfsegp)){
		i1 = newfsegp->intersect(fsegp);
		i2 = fsegp->intersect(newfsegp);
		if( (i1 == -1) && (i2 == -1) )
			return -1;
		fsegp = nextSeg();
	}
	return 1;
}

int TIN::readBoundSegs(istream& is)
{
	Segment *segP;
	Node *nodeP1, *nodeP2;
	char c = 'a';
	int name, n1, n2;

	while(is){
		if ( !(is >> name))
			break;
		if ( !(is >> n1))
			break;
		nodeP1 = (Node *) nodeL.n(n1);
		if ( !(is >> n2))
			break;
		nodeP2 = (Node *) nodeL.n(n2);

		segP = physics->readInSeg(name,nodeP1,nodeP2,is);
		if(segP->getn() >= 0)
			bsegL.appendItem(segP);
	}
	is.clear();
	while((c != '.') && is) {
		is.get(c);
		cout << c;
	}
	return bsegL.numberOfItems();
}

int TIN::readFeatureSegs(istream& is)
{
	Segment *segP;
	Node *nodeP1, *nodeP2;
	char c = 'a';
	int name, n1, n2;

	while(is && (c != '.')){
		c = is.peek();
		if(isdigit(c)){
			if ( !(is >> name))
				break;
			if ( !(is >> n1))
				break;
			nodeP1 = (Node *) nodeL.n(n1);
			if ( !(is >> n2))
				break;
			nodeP2 = (Node *) nodeL.n(n2);
			if(nodeP1 == NULL || (nodeP2 == NULL))
				return -1;
			segP = new Segment(name,nodeP1,nodeP2);
			addFeatureSeg(segP,fsegL);
		}
		else
			is.get(c);
	}
	is.clear();
	return fsegL.numberOfItems();
}

void TIN::writeNodes(ostream& os)
{
	Node *nP;
	for(int i=1;i<=nodeL.numberOfItems();i++) {
		nP = (Node *) nodeL.i(i);
		nP->setn(i);
		os << nP << "\n";
	}
	os << "\nno more nodes.\n\n";
}

void TIN::writeNodesDes(ostream& os)
{
	Node *nP;
	for(int i=1;i<=nodeL.numberOfItems();i++) {
		nP = (Node *) nodeL.i(i);
		nP->setn(i);
		os << nP;
		if (nP->getDesignation() != NULL)
			os << nP->getDesignation();
		os << "\n";
	}
	os << "\nno more nodes.\n\n";
}

void TIN::writeBoundSegs(ostream& os)
{
	for(int i=1;i<=bsegL.numberOfItems();i++)
		physics->writeBSeg(os,(Segment *)bsegL.i(i));
	os << "\nno more boundary segments.\n\n";
}

void TIN::writeFeatureSegs(ostream& os)
{
	for(int i=1;i<=fsegL.numberOfItems();i++)
		os << (Segment *)fsegL.i(i);
	os << "\nno more breakline segments.\n\n";
}

void TIN::writeTriangles(ostream& os)
{
	for(int i=1;i<=elmL.numberOfItems();i++)
		os << (Triangle *)elmL.i(i);
}

void TIN::dumpCSVFile(ostream& os,int index)
{
	Node *nodeP;

	os.precision(4);
	os.setf(ios::fixed,ios::floatfield);

	for(int i=1;i<=nodeL.numberOfItems();i++){
		nodeP = (Node *)nodeL.i(i);
		os << nodeP->getn() << "," ;
		os << nodeP->getXo() << "," ;
		os << nodeP->getYo() << "," ;
		os << nodeP->getPar(index) << "\n" ;
	}
}

void TIN::dumpCSVFile(ostream& os, int index1, int index2)
{
	Node *nodeP;

	os.precision(4);
	os.setf(ios::fixed,ios::floatfield);

	for(int i=1;i<=nodeL.numberOfItems();i++){
		nodeP = (Node *)nodeL.i(i);
		os << nodeP->getn() << "," ;
		os << nodeP->getXo() << "," ;
		os << nodeP->getYo() << "," ;
		os << nodeP->getPar(index1) << ",";
		os << nodeP->getPar(index2) << "\n" ;
	}
}

double TIN::meshQuality(Triangle** thisOne)
{
	double minQual = 1.0;
   *thisOne = NULL;
	double elmQual;
	for(int i=1;i<=elmL.numberOfItems();i++){
		if( ((Triangle *)elmL.i(i))->getStatus() == active){
			elmQual = ((Triangle *)elmL.i(i))->quality();
			if (elmQual < minQual){
				minQual = elmQual;
				*thisOne = (Triangle *)elmL.i(i);
			}
		}
	}
	return minQual;
}

int TIN::boundaryRefine(double minQual)
{
	int index[3] = {0,0,0}, count = 0, ind;
	double length[3];
	Triangle *tP;

	tP = firstTri();
	while(tP != NULL){
   	index[0] = index[1] = index[2] = 0;
		if( (tP->getStatus() == active) && (tP->quality() < minQual)){
			if ( ((tP->getAdjt(0))->numNodes()) == 2)
				index[0] = 1;
			if ( ((tP->getAdjt(1))->numNodes()) == 2)
				index[1] = 1;
			if ( ((tP->getAdjt(2))->numNodes()) == 2)
				index[2] = 1;
			if ((index[0] + index[1] + index[2]) == 1){
				length[0] = (tP->getNode(1))->dist(tP->getNode(2));
				length[1] = (tP->getNode(2))->dist(tP->getNode(0));
				length[2] = (tP->getNode(0))->dist(tP->getNode(1));
				if((length[0] > length[1]) && (length[0] > length[2]))
					ind = 0;
				else if((length[1] > length[0]) && (length[1] > length[2]))
					ind = 1;
				else
					ind = 2;
				if(index[ind] == 1){
					bisectSeg((Segment *)tP->getAdjt(ind),bsegL);
					count++;
				}
			}
		}
		tP = nextTri();
	}
	return count;
}

int TIN::triangulate()
{	
	if(nodeL.numberOfItems() < 3)
		return(0);
											// delete old triangles
	Triangle *tP = (Triangle *) elmL.firstItem();
	while(tP != NULL){
		elmL.deleteCurrentItem();
		delete tP;
		tP = (Triangle *) elmL.firstItem();
	}

	checkDupNodes();
	getLimits();
	Quadrant *qtp = new Quadrant(limits);

	Node	**nodeList = new Node*[nodeL.numberOfItems()];
	Node	*nP = firstNode();
	int		inode = 0;

	double xn, yn, zn;
	while(nP != NULL){
		nodeL.deleteItem(nP);
		if(nP->getFixed() != deleted) {
			nodeList[inode] = nP;
			inode++;
			nP->restoreLocation();   // resets coordinates to prevent accumulated drift
			xn = nP->getXc() + (rand()%500 - 250)/100000000.0;  //  joggle
			yn = nP->getYc() + (rand()%500 - 250)/100000000.0;
			zn = nP->getZc();
			nP->assignt(xn,yn,zn);
			qtp->appendNode(nP);
		}
		else{
			deadNodeL.appendItem(nP);
		}
		nP = firstNode();
	}

	int numNodes =  qtp->triangulateNodes();

	for(int ii=0;ii<inode;ii++){
		nodeL.appendItem(nodeList[ii]);
	}
	delete [] nodeList;

	Triangle *nextTriangle;
	tP = qtp->getTriangles();
	while(tP != NULL){
		nextTriangle = tP->ntP;
		elmL.appendItem(tP);
		tP = nextTriangle;
	}

	Triangle *arcP = qtp->getBArcs();
	Triangle *nextArcP;
	arcP->tP1()->setAdj(0,NULL);
	while(arcP != NULL) {
		nextArcP = arcP->tP0();
		arcP->tP2()->setAdj(arcP->tsn(),NULL);
		delete arcP;
		arcP = nextArcP;
	}
	delete qtp;

	for(int i =1;i<= elmL.numberOfItems();i++){			//	Renumber triangles
		tP = (Triangle *) elmL.i(i);
		tP->setn((tP->getIndex())+1);
		for(int j=0;j<3;j++){
			tP->getNode(j)->setATriangle(tP);
		}
	}
	
	return elmL.numberOfItems();
}

int TIN::checkDupNodes()
{	
	getLimits();
	Quadrant *qtp = new Quadrant(limits);

	Node	**nodeList = new Node*[nodeL.numberOfItems()];
	Node	*nP = firstNode();
	int		inode = 0;

	while(nP != NULL){
		nodeL.deleteItem(nP);
		nodeList[inode] = nP;
		inode++;
		qtp->appendNode(nP);
		nP = firstNode();
	}

	int numNodes =  qtp->checkDuplicateNodes();

	for(int ii=0;ii<inode;ii++){
		nodeL.appendItem(nodeList[ii]);
	}
	delete [] nodeList;

	delete qtp;
	
	return elmL.numberOfItems();
}

void TIN::insertNewNode(Node *nodeP, Triangle *triP)
{
	Node *nP1, *nP2, *nP3;
	Triangle *tP1, *tP2, *tP3;
	Element *atP;
	ItemList triStack;
	int index;
	
	nP1 = triP->getNode(0);							//	Save original nodes
	nP2 = triP->getNode(1);
	nP3 = triP->getNode(2);
	
	tP1 = triP;
	elmL.setCurrentItem(triP);						//	Take existing element out of mesh
	elmL.deleteCurrentItem();						//  Start working on three new ones
	tP1->setNode(0,nodeP);							//	Put them on the stack for checking later
	triStack.push(tP1);
	tP2 = new Triangle(elmL.numberOfItems()+2,nodeP,nP3,nP1);
	triStack.push(tP2);
	tP3 = new Triangle(elmL.numberOfItems()+3,nodeP,nP1,nP2);
	triStack.push(tP3);
	
	tP2->setAdj(0,tP1->getAdjt(1));
	tP2->setAdj(1,tP3);
	tP2->setAdj(2,tP1);

	tP3->setAdj(0,tP1->getAdjt(2));
	tP3->setAdj(1,tP1);
	tP3->setAdj(2,tP2);

	tP1->setAdj(1,tP2);
	tP1->setAdj(2,tP3);

	if( (atP = tP2->getAdjt(0)) != NULL) {
		index = atP->reflectAdj(tP1);
		atP->setAdj(index, tP2);
	}
	
	if( (atP = tP3->getAdjt(0)) != NULL) {
		index = atP->reflectAdj(tP1);
		atP->setAdj(index, tP3);
	}
	
	while( (tP1 = (Triangle *) triStack.pop()) != NULL) {	// while any tringles on stack
		if( (atP = tP1->getAdjt(0)) != NULL) {				// if not adjacent to nothing
			if(atP->numNodes() == 3) {						// if adjacent to a triangle (could be a boundary segment)
				tP2 = (Triangle *) atP;
				if(tP2->insideC(nodeP) >= 0) {			//	If not Delauney
					elmL.setCurrentItem(tP2);
					elmL.deleteCurrentItem();			// Take second triangle out for now
					swapDiagonal(tP1,tP2);				//	We need to swap
					triStack.push(tP1);					//  New triangles need to be checked
					triStack.push(tP2);					//	afterward so put both on stack
				}
				else{
					elmL.insertItem(tP1);
					for(int i=0;i<3;i++){
						tP1->getNode(i)->setATriangle(tP1);
					}
				}
			}
			else{
				elmL.insertItem(tP1);			// Triangle good so accept back into mesh
				for(int i=0;i<3;i++){
					tP1->getNode(i)->setATriangle(tP1);
				}
			}
		}
		else{
			elmL.insertItem(tP1);
			for(int i=0;i<3;i++){
				tP1->getNode(i)->setATriangle(tP1);
			}
		}
	}
}

int TIN::triangulateCons()
{
	int nelms;
	if(nodeL.numberOfItems() < 3)
		return(0);
	nelms = triangulate();
	if(bsegL.numberOfItems() < 3)
		return nelms;

	Triangle *tP1, *tP2, *btP, *tP, *stP;
	Element *pP;
	Segment *nbP;
	Segment *bP = firstSeg();
	stP = NULL;
	while(bP != NULL) {
		bP->getNode(0)->setOnBound();
		int notGood = 0;						//	Now thread the boundary segments
		while( (btP = whichTriHasSeg(bP,stP)) == NULL) {
			tP1 = getFirstTriOnSeg(bP);
			tP2 = getSecondTriOnSeg(tP1,bP);
			swapDiagonal(tP1,tP2);
		}
		insertBoundSeg(bP,btP);
		for(int side=0;side<2;side++){
			tP1 = (Triangle *) bP->getAdjt(side);
			if(tP1 != NULL){
				for(int adj=0;adj<3;adj++){
					pP = tP1->getAdjt(adj);
					if(pP != NULL && pP != bP){
						if(pP->numNodes() == 3) {
							tP2 = (Triangle *) pP;
							if(tP1->area()<=0.0 || tP2->area()<=0.0){
								swapDiagonal(tP1,tP2);
							}
						}
					}
				}
			}
		}
		nbP = nextSeg();
		stP = NULL;
		if(nbP != NULL){
			if(nbP->getNode(0) == bP->getNode(1)){
				stP = (Triangle *) bP->getAdjt(0);
			}
		}
		bP = nbP;
	}
/*
	bP = firstSeg();
	while(bP != NULL) {
		btP = whichTriHasSeg(bP);
		insertBoundSeg(bP,btP);
		for(int side=0;side<2;side++){
			tP1 = (Triangle *) bP->getAdjt(side);
			if(tP1 != NULL){
				for(int adj=0;adj<3;adj++){
					pP = tP1->getAdjt(adj);
					if(pP != NULL && pP != bP){
						if(pP->numNodes() == 3) {
							tP2 = (Triangle *) pP;
							if(tP1->area()<=0.0 || tP2->area()<=0.0){
								swapDiagonal(tP1,tP2);
							}
						}
					}
				}
			}
		}
		bP = nextSeg();
	}
*/
	bP = (Segment *) bsegL.firstItem();				//	Deactivate (but retain) triangles
	while(bP != NULL){								//	which are outside the boundary
		if((pP = bP->getAdjt(1)) != NULL)				
			pP->deactivate();
		bP = (Segment *) bsegL.nextItem();
	}
	tP = (Triangle *) elmL.firstItem();				//	Catch the outside triangles with no
	while(tP != NULL){								//	edges on boundary
		if( (tP->getAdjt(0) == NULL) ||
			 (tP->getAdjt(1) == NULL) ||
			 (tP->getAdjt(2) == NULL)   )
			 tP->deactivate();
		else if((tP->getAdjt(0)->getStatus() == notActive) ||
				  (tP->getAdjt(1)->getStatus() == notActive) ||
				  (tP->getAdjt(2)->getStatus() == notActive) ) {
				tP->deactivate();
//				cout << "Caught'im " << tP;
		}
		if(tP->getStatus() == active){
			for(int j=0;j<3;j++){
				tP->getNode(j)->setATriangle(tP);
			}
		}
		tP = (Triangle *) elmL.nextItem();
	}
	for(int i =1;i<= elmL.numberOfItems();i++){			//	Renumber triangles
		tP = (Triangle *) elmL.i(i);
		tP->setn((tP->getIndex())+1);
	}

	Node *nP;
	for(int i =1;i<= nodeL.numberOfItems();i++){			//	Renumber nodes //the keyword int was added before the variable i May 14
		nP = (Node *) nodeL.i(i);
		nP->setn((nP->getIndex())+1);
	}

	return elmL.numberOfItems();

}

int TIN::doEdges()
{
	Segment *eP = (Segment *) edgeL.firstItem();		// Delete old edges
	while(eP != NULL){
		edgeL.deleteCurrentItem();
		delete eP;
		eP = (Segment *) edgeL.firstItem();
	}

	Triangle *tP = firstTri();							// Clear edge pointers in triangles
	while(tP != NULL){
		tP->setEdge(0,NULL);
		tP->setEdge(1,NULL);
		tP->setEdge(2,NULL);
		tP = nextTri();
	}

	int n1, n2, ia, ne = 0;
	Triangle *atP;
	Element  *aeP;
	tP = firstTri();
	while(tP != NULL){
		if(tP->getStatus() == active){
			for(int i=0;i<3;i++){
				if(tP->getEdge(i) == NULL){
					ne++;
					n1 = (i+1)%3;
					n2 = (i+2)%3;
					eP = new Segment(ne,tP->getNode(n1),tP->getNode(n2));
					eP->setAdj(0,tP);
					tP->setEdge(i,eP);
					aeP = tP->getAdjt(i);
					if(aeP != NULL){
						if(aeP->numNodes() == 3){
							atP = (Triangle *) aeP;
							ia = atP->reflectAdj(tP);
							eP->setAdj(1,atP);
							atP->setEdge(ia,eP);
						}
					}
					edgeL.appendItem(eP);
				}
			}
		}
		tP = nextTri();
	}
	return ne;
}

Triangle* TIN::getFirstTriOnSeg(Segment *segP)
{
	Node n(1);
	double ds = 0.01;
	segP->locateNodeAt(ds,&n);
	Triangle *triP = whichTriangle(&n);
	while( (triP->contains(segP->getNode(0))) != triP){
		ds *= 0.1;
		segP->locateNodeAt(ds,&n);
		triP = whichTriangle(&n);
	}
	return triP;
}

Triangle* TIN::getSecondTriOnSeg(Triangle *fTP, Segment *segP)
{
	int index = fTP->nodeIndex(segP->getNode(0));
   return (Triangle *) fTP->getAdjt(index);
}

void TIN::insertBoundSeg(Segment *bP,Triangle *triP)
{
	Element *atP;

	Node *bnP = bP->getNode(0);				//	Rotate the element indices so that
	if(bnP == triP->getNode(1)) {				//	the tail of the boundary segment is at 0
		triP->setNode(1,triP->getNode(2));
		triP->setNode(2,triP->getNode(0));
		triP->setNode(0,bnP);
		atP = triP->getAdjt(1);
		triP->setAdj(1,triP->getAdjt(2));
		triP->setAdj(2,triP->getAdjt(0));
		triP->setAdj(0,atP);
	}
	else if(bnP == triP->getNode(2)) {
		triP->setNode(2,triP->getNode(1));
		triP->setNode(1,triP->getNode(0));
		triP->setNode(0,bnP);
		atP = triP->getAdjt(2);
		triP->setAdj(2,triP->getAdjt(1));
		triP->setAdj(1,triP->getAdjt(0));
		triP->setAdj(0,atP);
	}

	atP = triP->getAdjt(2);
	bP->setAdj(0,triP);
	bP->setAdj(1,atP);
	triP->setAdj(2,bP);							//	triangles, 0 is to left (inside) of segment
	if(atP != NULL)
		atP->setAdj(atP->reflectAdj(triP),bP);

	return;
}

void TIN::swapDiagonal(Triangle *tP1, Triangle *tP2)
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
//	for(int i=0;i<3;i++){
//		tP1->getNode(i)->setATriangle(tP1);
//		tP2->getNode(i)->setATriangle(tP2);
//	}
}

int TIN::doFeatures(TIN* dataTINP)
{
	int startNum, nTimes = 0, errcode = 1;
	Segment *sP, *lastSP = NULL;
	Segment *fP;

//	checkAllFeatureSegs();            //		check feature segs for consistentency

	do {
		startNum = nodeL.numberOfItems();            // Thread feature lines through
		sP = firstfSeg();						//	the mesh, inserting additional
		fP = sP;
		lastSP = sP;
		while(sP != NULL) {
				sP = threadFeature(sP,fsegL);
			if(sP == lastSP)
				nTimes += 1;
			else
				nTimes = 0;
			if(nTimes > 4){
				errcode = - sP->getn();
				sP = (Segment*) fsegL.nextItem();
				break;
			}
			if(sP == fP)
				break;
			lastSP = sP;
		}
		if(errcode < 0)
			break;
	}while(nodeL.numberOfItems() != startNum);

	Triangle* tP;
	int i;
	for(i =1;i<= elmL.numberOfItems();i++){			//	Renumber triangles
		tP = (Triangle *) elmL.i(i);
		tP->setn((tP->getIndex())+1);
	}

/*	Node* nP;
	for(i =1;i<= nodeL.numberOfItems();i++){			//	Renumber nodes
		nP = (Node *) nodeL.i(i);
		nP->setn((nP->getIndex())+1);
	}
*/
	return errcode;
}

Segment* TIN::threadFeature(Segment *fsegp, ItemList &segList)
{
	Node n(1);
	fsegp->locateNodeAt(0.001,&n);
	Triangle *triP;
	if( (triP = whichTriangle(&n)) == NULL){
		return (Segment*) segList.nextItem();
	}
	if ((triP->contains(fsegp)) == triP){
		return (Segment*) segList.nextItem();
	}
	Segment revfseg(1,fsegp->getNode(1),fsegp->getNode(0));
	if ((triP->contains(&revfseg)) == triP){
		return (Segment*) segList.nextItem();
	}
	int nindex = triP->nodeIndex(fsegp->getNode(0));
	if(nindex < 0)
		return (Segment*) segList.nextItem();
	nindex += 1;
	if(nindex > 2)
		nindex = 0;
	Node *on1 = triP->getNode(nindex);
	nindex += 1;
	if(nindex > 2)
		nindex = 0;
	Node *on2 = triP->getNode(nindex);
	Segment opSeg(1,on1,on2);
	double dist = fsegp->intersectd(&opSeg);
	if((dist > 0.0) && (dist < 1.0)){
		bisectSeg(fsegp,segList,dist);
		if( (triP = whichTriangle(fsegp->getNode(1))) != NULL){
			insertNewNode(fsegp->getNode(1),triP);
		}
	}
	return (Segment *) 	segList.nextItem();
}

double TIN::smoothMesh(int nTimes, TIN& dataTIN, double bias)
{
	Node *nP;
	Triangle *tP, *atP, *dtP;
	Node *nP1;
	double a, sumx, sumy, suma, *xNew, *yNew, xN, yN, dx, dy, lx, ly, d;
	int	j, index, iin=0;

   bias = sqrt(bias);

	xNew = new double[nodeL.numberOfItems()];
	yNew = new double[nodeL.numberOfItems()];

	for(int i=0;i<nTimes;i++){
		nP = (Node *) nodeL.firstItem();
	
		while(nP != NULL){
			xNew[iin] = nP->getXc();
			yNew[iin] = nP->getYc();
			if( (nP->getFixed() == floating) && (nP->getBound() == interior) ) {
				if((tP = whichTriHasNode(nP)) != NULL) {
					if(tP->getStatus() == notActive) { // node has slipped outside boundary
						nP->setDeleted();
						xNew[iin] = nP->getXc();
						yNew[iin] = nP->getYc();
					}
					else{
						sumx = sumy = suma = 0.0;
						atP = tP;
						do {
							a = ((1.0-bias) + bias*fabs(dataTIN.triCenterDz(atP))) * sqrt(fabs(atP->area())); // atP->area();
							suma += 4.0*a;
							sumx += nP->getXc() * a;
							sumy += nP->getYc() * a;
							for(j=0;j<3;j++) {
								nP1 = atP->getNode(j);
								sumx += nP1->getXc() * a;
								sumy += nP1->getYc() * a;
								if(nP1 == nP)
									index = j;
							}
							index -= 1;
							if(index == -1)
								index = 2;
							atP = (Triangle *) atP->getAdjt(index);
						} while (atP != tP);
						xNew[iin] = (nP->getXc() + sumx/suma) * 0.5;
						yNew[iin] = (nP->getYc() + sumy/suma) * 0.5;
					}
				}
			}
			nP = (Node *) nodeL.nextItem();
			iin += 1;
		}

		Segment *fsP = firstfSeg();
		Segment *nsP;
		while (fsP != NULL){
			nP = fsP->getNode(1);
			if(nP->getFixed() == sliding) {
				tP = whichTriHasNode(nP);
				sumx = sumy = suma = 0.0;
				atP = tP;
				do {
					a = ((1.0-bias) + bias*fabs(dataTIN.triCenterDz(atP))) * sqrt(fabs(atP->area())); // atP->area();
					suma += 4.0*a;
					sumx += nP->getXc() * a;
					sumy += nP->getYc() * a;
					for(j=0;j<3;j++) {
						nP1 = atP->getNode(j);
						sumx += nP1->getXc() * a;
						sumy += nP1->getYc() * a;
						if(nP1 == nP)
							index = j;
					}
					index -= 1;
					if(index == -1)
						index = 2;
					atP = (Triangle *) atP->getAdjt(index);
				} while (atP != tP);
				xN = (nP->getXc() + sumx/suma) *0.5;
				yN = (nP->getYc() + sumy/suma) *0.5;
				if( (nsP = (Segment *) fsP->getNextOne()) != NULL){
					if(nsP->getNode(0) == nP){
						dx = xN - fsP->getNode(0)->getXc();
						dy = yN - fsP->getNode(0)->getYc();
						lx = nsP->getNode(1)->getXc() - fsP->getNode(0)->getXc();
						ly = nsP->getNode(1)->getYc() - fsP->getNode(0)->getYc();
						d = (dx*lx + dy*ly)/(lx*lx + ly*ly);
						if(d < 0.33333)
							d = 0.33333;
						else if(d > 0.66667)
							d = 0.66667;
						xN = fsP->getNode(0)->getXc() + d*lx;
						yN = fsP->getNode(0)->getYc() + d*ly;
						nP->assign(xN,yN,100.0);
						if( (dtP = dataTIN.whichTriangle(nP)) != NULL)
							dtP->locateNodeAt(xN,yN,nP);
					}
				}
			}
			fsP = nextfSeg();
		}

		Segment *bsP = firstSeg();
		Element *elP;
		while (bsP != NULL){
			nP = bsP->getNode(1);
			if(nP->getFixed() == sliding) {
				tP = (Triangle *) bsP->getAdjt(0);
				sumx = sumy = suma = 0.0;
				atP = tP;
				do {
					a = sqrt(fabs(atP->area()));
					suma += 3.0*a;
					for(j=0;j<3;j++) {
						nP1 = atP->getNode(j);
						sumx += nP1->getXc() * a;
						sumy += nP1->getYc() * a;
						if(nP1 == nP)
							index = j;
					}
					index -= 1;
					if(index == -1)
						index = 2;
					elP = atP->getAdjt(index);
					if(elP->numNodes() == 3)
						atP = (Triangle *) elP;
					else
						break;
				} while (atP != tP);
				xN = (nP->getXc() + sumx/suma) *0.5;
				yN = (nP->getYc() + sumy/suma) *0.5;
				nsP = (Segment *) bsP->getNextOne();
				dx = xN - bsP->getNode(0)->getXc();
				dy = yN - bsP->getNode(0)->getYc();
				lx = nsP->getNode(1)->getXc() - bsP->getNode(0)->getXc();
				ly = nsP->getNode(1)->getYc() - bsP->getNode(0)->getYc();
				d = (dx*lx + dy*ly)/(lx*lx + ly*ly);
				if(d < 0.33333)
					d = 0.33333;
				else if(d > 0.66667)
					d = 0.66667;
				xN = bsP->getNode(0)->getXc() + d*lx;
				yN = bsP->getNode(0)->getYc() + d*ly;
				nP->assign(xN,yN,100.0);
				if( (dtP = dataTIN.whichTriangle(nP)) != NULL)
					dtP->locateNodeAt(xN,yN,nP);
			}
			bsP = nextSeg();
		}

		nP = firstNode();
		iin = 0;
		while(nP != NULL){
			if( nP->getFixed() == floating ) {
				nP->assign(xNew[iin],yNew[iin],100.0);
				if( (dtP = dataTIN.whichTriangle(nP)) != NULL)
					dtP->locateNodeAt(xNew[iin],yNew[iin],nP);
			}
			nP = nextNode();
			iin += 1;
		}
	}
	delete [] xNew;
	delete [] yNew;

	return 1.0;
}

double TIN::triCenterDz(Triangle *triP)
{
	Node tCNode(100), tinNode(101);
	triP->locateNodeAtCenter(&tCNode);
	Triangle* tinTri;
	double dz = 0.0;
	if( (tinTri = whichTriangle(&tCNode)) != NULL) {
		tinTri->locateNodeAt(tCNode.getXc(),tCNode.getYc(),&tinNode);
		dz = tCNode.getZc() - tinNode.getZc();
	}
	return (dz);
}

int TIN::updateMesh(TIN& dataTIN)
{
	int err = 0;

	Triangle* dtP;
	Node* nP = firstNode();
	while(nP != NULL){
		if( (dtP = dataTIN.whichTriangle(nP)) != NULL)
			dtP->locateNodeAt(nP->getXc(),nP->getYc(),nP);
		else {
			err = -1;
			break;
		}
		nP = nextNode();
	}
	return err;
}

box TIN::getLimits(double margin)
{
	Node *nP;
	
	nP = (Node *) nodeL.firstItem();
	
	if(nP != NULL) {
		limits.x1 = nP->getXc();
		limits.y1 = nP->getYc();
		limits.x2 = nP->getXc();
		limits.y2 = nP->getYc();
		nP = (Node *) nodeL.nextItem();
	}

	while(nP != NULL){

		if(nP->getXc() < limits.x1)
			limits.x1 = nP->getXc();
		if(nP->getXc() > limits.x2)
			limits.x2 = nP->getXc();
		if(nP->getYc() < limits.y1)
			limits.y1 = nP->getYc();
		if(nP->getYc() > limits.y2)
			limits.y2 = nP->getYc();

		nP = (Node *) nodeL.nextItem();
	}
	
	double dx = limits.x2 - limits.x1;
	double dy = limits.y2 - limits.y1;
	limits.x1 -= margin * dx;
	limits.x2 += margin * dx;
	limits.y1 -= margin * dy;
	limits.y2 += margin * dy;
	return limits;
}
	  
double TIN::getMinPar(int nPar)
{
	Node *nP;
	double min;
	
	nP = (Node *) nodeL.firstItem();
	
	if(nP != NULL) 
		min = nP->getPar(nPar);
	else
		return 0.0;

	while(nP != NULL){
		if(nP->getPar(nPar) < min)
			min = nP->getPar(nPar);
		nP = (Node *) nodeL.nextItem();
	}

	return min;
}
	  
double TIN::getMaxPar(int nPar)
{
	Node *nP;
	double max;
	
	nP = (Node *) nodeL.firstItem();
	
	if(nP != NULL) 
		max = nP->getPar(nPar);
	else
		return 0.0;

	while(nP != NULL){
		if(nP->getPar(nPar) > max)
			max = nP->getPar(nPar);
		nP = (Node *) nodeL.nextItem();
	}

	return max;
}


void TIN::scratchTriangle()
{
	Triangle* tP;
	Element* pP;
	int ind;

	tP = (Triangle *) elmL.currentItem();
	for(int i=0;i<3;i++) {
		pP = tP->getAdjt(i);
		if(pP != NULL) {
			ind = pP->reflectAdj(tP);
			if(ind >= 0)
				pP->setAdj(ind,NULL);
		}
	}
	elmL.deleteCurrentItem();
	delete tP;
}
		
Triangle* TIN::whichTriangle(Node * nodeP)
{
	Element *retP = NULL, *triP;

	triP = (Element *) elmL.currentItem();
	if(triP == NULL)
		return NULL;
	while(-1) {
		retP = triP->inside(nodeP);
		if ((retP == triP) || (retP == NULL))
			break;
		triP = retP;
		if(retP->numNodes() == 3)
			elmL.setCurrentItem(retP);
	}
	if(retP != NULL){
		elmL.setCurrentItem(retP);
	}
/*	else {
		if( (nodeP->dist(triP->getNode(0)) < 0.0001) ||
			(nodeP->dist(triP->getNode(1)) < 0.0001) ||
			(nodeP->dist(triP->getNode(2)) < 0.0001)    )
				retP = triP;
	}
*/
	return (Triangle *) retP;
}

Triangle* TIN::whichTriHasNode(Node *nodeP)
{
	Triangle *triP;
	triP = nodeP->getATriangle();
	if(triP != NULL){
		elmL.setCurrentItem(triP);
	}
	return triP;
}

Triangle* TIN::triangleListForNode(Node *nP, Triangle *sTP)
{
	Triangle *tPf;
	if(sTP == NULL)
		tPf = nP->getATriangle();
	else
		tPf = sTP;

	Triangle *tPs = tPf;
	Triangle *tP1;
	Node *nP1;
	int	index;
	int hitEdge = 0;
	Element *aEP;
	if(tPs != NULL){
		tPs->ptP = tPs->ntP = NULL;
		do{
			for(int i=0;i<3;i++){
				nP1 = tPs->getNode(i);
				if(nP1 == nP){
					index = i-1;
					if(index == -1)
						index = 2;
					break;
				}
			}
			aEP = tPs->getAdjt(index);
			if(aEP != NULL){
				if(aEP->numNodes() == 3){
					tP1 = (Triangle *)aEP;
					if(tP1 != tPf){
						tPs->ntP = tP1;
						tP1->ptP = tPs;
						tP1->ntP = NULL;
					}
					tPs = tP1;
				}
				else 
					hitEdge = 1;
			}
			else
				hitEdge = 1;
		} while((tPs != tPf) && (hitEdge == 0));
		if(hitEdge == 1){
			tPs = tPf;
			hitEdge = 0;
			do{
				for(int i=0;i<3;i++){
					nP1 = tPs->getNode(i);
					if(nP1 == nP){
						index = i+1;
						if(index == 3)
							index = 0;
						break;
					}
				}
				aEP = tPs->getAdjt(index);
				if(aEP != NULL){
					if(aEP->numNodes() == 3){
						tP1 = (Triangle *)aEP;
						tPs->ptP = tP1;
						tP1->ntP = tPs;
						tP1->ptP = NULL;
						tPs = tP1;
					}
					else 
						hitEdge = 1;
				}
				else
					hitEdge = 1;
			} while(hitEdge == 0);
			tPf = tPs;
		}
	}

	return tPf;
}

Triangle* TIN::whichTriHasSeg(Segment *segP, Triangle *sTP)
{
	Triangle *triP;

	triP = triangleListForNode(segP->getNode(0),sTP);
	while(triP != NULL) {
		if ((triP->contains(segP)) == triP)
			break;
		triP = triP->ntP; 
	}

	return triP;
}

Segment* TIN::bisectSeg(Segment *segP, ItemList &segL, double dist)
{
	Node *nodeP;
	Segment *nSegP;

	if(dist > 0.99)
		dist = 0.99;
	if(dist < 0.01)
		dist = 0.01;
	nodeP = physics->makeNewNode(getNextN());
	if((segP->getNode(0)->getBound() == onBoundary) &&
			(segP->getNode(1)->getBound() == onBoundary))
		nodeP->setOnBound();
	if(		(segP->getNode(0)->getFixed() == fixednode)
		 || 	(segP->getNode(0)->getFixed() == sliding) )
		nodeP->setSliding();
	nodeL.appendItem(nodeP);
	segP->locateNodeAt(dist,nodeP);
	nSegP = physics->makeNewSegment(segL.numberOfItems()+1,
			nodeP,segP->getNode(1),segP);

//	int iin;
//	cout <<"Bisecting segment " << segP;
//	cout <<"New Node " << nodeP;
//	cout <<"New Segment " << nSegP;
//	writeTriangles(cout);
//	cin >> iin;

	segP->setNode(1,nodeP);
	segL.setCurrentItem(segP);
	segL.insertItem(nSegP);

	Segment *sp = (Segment *) segL.firstItem();
	if(sp == segP)
		return sp;
	while(sp->getNextOne() != segP)
		sp = (Segment *) segL.nextItem();
	return sp;
}

//	New function added 4/2001, JDS
//	Same as bisectSeg but specific to boundary segments.  This function
//	handles the operation when the user chooses to manually bisect
//	a boundary segment.

Segment* TIN::bisectBSeg(Segment *segP, double dist)
{
	Node *nodeP;
	Segment *nSegP;

	if(dist > 0.99) dist = 0.99;
	if(dist < 0.01) dist = 0.01;
	nodeP = physics->makeNewNode(getNextN());
	nodeP->setOnBound();
	nodeP->setFixed();
	nodeL.appendItem(nodeP);
	segP->locateNodeAt(dist,nodeP);
	nSegP = physics->makeNewSegment(bsegL.numberOfItems()+1,
			nodeP,segP->getNode(1),segP);

	segP->setNode(1,nodeP);
	bsegL.setCurrentItem(segP);
	bsegL.insertItem(nSegP);

	Segment *sp = (Segment *) bsegL.firstItem();
	if(sp == segP)
		return sp;
	while(sp->getNextOne() != segP)
		sp = (Segment *) bsegL.nextItem();
	return sp;
}

void TIN::acceptNodes()
{
	nodeL.catList(&tnodeL);
	tnodeL.emptyList();

	Node *nP = firstNode();
	int i = 1;
	while(nP != NULL){
		nP->setn(i++);
		nP = nextNode();
	}
}

void TIN::clearBsegs()
{
	Segment *sP = (Segment *) bsegL.firstItem();
	while(sP != NULL){
		bsegL.deleteCurrentItem();
		delete sP;
		sP = (Segment *) bsegL.firstItem();
	}
}

void TIN::rejectNodes()
{
	Node *nP = (Node *) tnodeL.firstItem();
	while(nP != NULL){
		tnodeL.deleteCurrentItem();
		delete nP;
		nP = (Node *) tnodeL.firstItem();
	}
}

void TIN::acceptBsegs()
{
	Segment *sP = (Segment *) bsegL.firstItem();
	while(sP != NULL){
		bsegL.deleteCurrentItem();
		delete sP;
		sP = (Segment *) bsegL.firstItem();
	}
	bsegL.catList(&tbsegL);
    tbsegL.emptyList();
}

void TIN::acceptNewBsegs()
{
	bsegL.catList(&tbsegL);
	tbsegL.emptyList();
}


void TIN::rejectBsegs()
{
	Node *nP = (Node *) tnodeL.firstItem();
	while(nP != NULL){
		tnodeL.deleteCurrentItem();
		delete nP;
		nP = (Node *) tnodeL.firstItem();
	}
	Segment *sP = (Segment *) tbsegL.firstItem();
	while(sP != NULL){
		tbsegL.deleteCurrentItem();
		delete sP;
		sP = (Segment *) tbsegL.firstItem();
	}
}

void TIN::acceptFsegs()
{
	fsegL.catList(&tfsegL);
	tfsegL.emptyList();
}

void TIN::rejectFsegs()
{
	Segment *sP = (Segment *) tfsegL.firstItem();
	while(sP != NULL){
		tfsegL.deleteCurrentItem();
		delete sP;
		sP = (Segment *) tfsegL.firstItem();
	}

}

Segment* TIN::defNewBound(double x, double y,TIN &dataTIN)

//	Puts one Node into the TIN (from physics spec).
//  Then joins that node to the current node to form a new boundary segment.
//	Requires a previously triangulated data TIN, dataTIN,
//	which masks inactive areas and provides interpolation data.
//	Resulting node is in tnodeL which can be later accepted into nodeL or rejected.
//	The segment created is returned.
//  If something went wrong, NULL is returned.

{
	Segment *nSegP = NULL;
	Node *pNP = (Node *) tnodeL.lastItem();
	Node *newNP = insertOneNode(x,y,dataTIN,dataTIN);
	if(newNP != NULL){
		nSegP = physics->makeNewSegment(bsegL.numberOfItems()+1,
													pNP,newNP,NULL);
		if(checkFeatureSeg(nSegP) == 1)
			bsegL.appendItem(nSegP);
		else {
			delete nSegP;
			tnodeL.deleteCurrentItem();
			delete newNP;
			tnodeL.setCurrentItem(pNP);
			nSegP = NULL;
		}
	}

	return nSegP;
}

Segment* TIN::defNewBound(double x, double y)
{
	Segment *nSegP = NULL;
	Node *pNP = (Node *) tnodeL.lastItem();
	Node *newNP = insertOneNode(x,y);
	if(newNP != NULL){
		nSegP = physics->makeNewSegment(bsegL.numberOfItems()+1,
													pNP,newNP,NULL);
		if(checkFeatureSeg(nSegP) == 1)
			bsegL.appendItem(nSegP);
		else {
			delete nSegP;
			tnodeL.deleteCurrentItem();
			delete newNP;
			tnodeL.setCurrentItem(pNP);
			nSegP = NULL;
		}
	}

	return nSegP;
}

Segment* TIN::closeBound(TIN &dataTIN)

// Joins the last node to the first node to close the boundary loop.
//	Requires a previously triangulated data TIN, dataTIN,
//	which masks inactive areas and provides interpolation data.
//	The segment created is returned.
// If something went wrong, NULL is returned.

{
	Segment *nSegP = NULL;
	Node *pNP = (Node *) tnodeL.lastItem();
	Node *newNP = firstSeg()->getNode(0);
	if(newNP != NULL){
		nSegP = physics->makeNewSegment(bsegL.numberOfItems()+1,
													pNP,newNP,NULL);
		if(checkFeatureSeg(nSegP) == 1)
			bsegL.appendItem(nSegP);
		else {
			delete nSegP;
			nSegP = NULL;
		}
	}

	return nSegP;
}

Segment* TIN::closeBound()
{
	Segment *nSegP = NULL;
	Node *pNP = (Node *) tnodeL.lastItem();
	Node *newNP = firstSeg()->getNode(0);
	if(newNP != NULL){
		nSegP = physics->makeNewSegment(bsegL.numberOfItems()+1,
													pNP,newNP,NULL);
		if(checkFeatureSeg(nSegP) == 1)
			bsegL.appendItem(nSegP);
		else {
			delete nSegP;
			nSegP = NULL;
		}
	}

	return nSegP;
}

Segment* TIN::closeBoundLoop(Node *firstNP, TIN &dataTIN)

// Joins the last node to the first node to close the boundary loop.
//	Requires a previously triangulated data TIN, dataTIN,
//	which masks inactive areas and provides interpolation data.
//	The segment created is returned.
// If something went wrong, NULL is returned.

{
	Segment *nSegP = NULL;
	Node *pNP = (Node *) tnodeL.lastItem();
	if(firstNP != NULL){
		nSegP = physics->makeNewSegment(bsegL.numberOfItems()+1,
													pNP,firstNP,NULL);
		bsegL.appendItem(nSegP);
	}

	return nSegP;
}

Node* TIN::insertOneNode(double x, double y, TIN &boundTIN, TIN &dataTIN)

//	Puts one Node into the TIN (from physics spec).
//	Requires a previously triangulated data TIN, dataTIN,
//	which masks inactive areas and provides interpolation data.
//	Resulting nodeis in tnodeL which can be later accepted into nodeL or rejected.
//	The number of nodes created are returned.

{
	int nGoodNodes = 0;
	int nodeNum = getNextN();
	Triangle *tP, *tBP;
	Node *nP;

	nP = physics->makeNewNode(nodeNum,x,y);
	if( (tP = dataTIN.whichTriangle(nP)) != NULL){
		if( (tBP = boundTIN.whichTriangle(nP)) != NULL){
			if(tBP->getStatus() == active){
				tP->locateNodeAt(x,y,nP);
				nGoodNodes += 1;
				nodeNum += 1;
				tnodeL.appendItem(nP);
			}
			else{
				delete nP;
				nP = NULL;
			}
		}
		else{
			delete nP;
			nP = NULL;
		}
	}
	else{
		delete nP;
		nP = NULL;
	}

	return nP;
}

Node* TIN::insertOneNode(double x, double y)
{
	int nodeNum = getNextN();
	Node *nP;

	nP = physics->makeNewNode(nodeNum,x,y);
	tnodeL.appendItem(nP);

	return nP;
}

int TIN::fillNodesUniform(double spacing, double theta,
									box dLimits, TIN &boundTIN, TIN &dataTIN, TIN &bound2TIN)

//	Fills TIN with uniformly spaced, equilaterally arranged Nodes (from physics spec).
//	Spacing is approximate (goal) distance between nodes.
//	Theta (in degrees, 0 to 90) is the angle of the pattern to the x direction.
//	dLimits is a box defining the area to be filled.
//	Requires a previously triangulated data TIN, dataTIN,
//	which masks inactive areas and provides interpolation data.
//	Resulting set of nodes are in tnodeL which can be later accepted into nodeL or rejected.
//	The number of nodes created are returned.

{
	double r, s, x, y, r0;
	int nGoodNodes = 0;
	int nodeNum = getNextN();
	Triangle *tP, *tBP;
	Node *nP;

	if( (theta < 0.0) || (theta > 90.0) )
		return 0;
		
	double ct = cos(theta/45. * atan(1.0));
	double st = sin(theta/45. * atan(1.0));
	double rmax = (dLimits.x2-dLimits.x1)*ct + (dLimits.y2-dLimits.y1)*st;
	double smax = (dLimits.x2-dLimits.x1)*st + (dLimits.y2-dLimits.y1)*ct;
	double x0 = dLimits.x1 + (dLimits.x2-dLimits.x1)*(1.0 - ct*ct);
	double y0 = dLimits.y1 - (dLimits.x2-dLimits.x1)*ct*st;
	double dr = spacing;
	double ds = spacing * sqrt(3.0) / 2.0;
	int nr = (int) (rmax/dr + 1.5);
	dr = rmax/nr;
	int ns = (int) (smax/ds +1.5);
	ds = smax/ns;
	
	for(int j=1;j<ns;j++){
		s = j * ds;
		r0 = ((j%2) - 1) * dr/2.;
		for(int i=1;i<nr;i++){
			r = r0 + i*dr;
			x = x0 + r*ct - s*st + (rand()%500)/100000.0;
			y = y0 + r*st + s*ct + (rand()%500)/100000.0;
			nP = physics->makeNewNode(nodeNum,x,y);
			
			if( (tP = dataTIN.whichTriangle(nP)) != NULL){
				if( (tBP = boundTIN.whichTriangle(nP)) != NULL){
					if(tBP->getStatus() == active){
						if((tBP = bound2TIN.whichTriangle(nP)) != NULL){
							if(tBP->getStatus() == active){
								tP->locateNodeAt(x,y,nP);
								nGoodNodes += 1;
								nodeNum += 1;
								tnodeL.appendItem(nP);
							}
							else
								delete nP;
						}
						else
							delete nP;
					}
					else
						delete nP;
				}
				else
					delete nP;
			}
			else
				delete nP;
		}
	}
	return nGoodNodes;
}


int TIN::fillNodesSource(double x0, double y0, double r0, int nRays,
								 double theta0, TIN &boundTIN, TIN &dataTIN)

//	Fills TIN with equilaterally arranged Nodes (from physics spec) in an expanding circular pattern.
//	The spacing is set analagous to to a point source potential flow net.
//	Spacing is approximate (goal) distance between nodes.
//	Theta1 (in degrees, 0 to 360) is the start angle of the pattern from the x direction.
//	Theta2 (in degrees, 0 to 360) is the finish angle of the pattern from the x direction.
//	Requires a previously triangulated data TIN, dataTIN,
//	which masks inactive areas and provides interpolation data.
//	Resulting set of nodes are in tnodeL which can be later accepted into nodeL or rejected.
//	The number of nodes created are returned.

{
	double r, rmax, theta, deltheta, fac, x, y;
	int j;
	int nGoodNodes = 0;
	int nodeNum = getNextN();
	box dLimits;
	Triangle *tP, *tBP;
	Node *nP;


	deltheta = 8.0*atan(1.0) / nRays;
	theta0 *= 8.0*atan(1.0) /360.0;
//	fac = (nRays + sqrt(3.0)/2.)/(nRays - sqrt(3.0)/2.);

	fac = exp(deltheta/sqrt(2.0));

	dLimits = dataTIN.getLimits();
	if(fabs(x0 - dLimits.x1) > fabs(x0 - dLimits.x2))
		x = fabs(x0 - dLimits.x1);
	else
		x = fabs(x0 - dLimits.x2);
	if(fabs(y0 - dLimits.y1) > fabs(y0 - dLimits.y2))
		y = fabs(y0 - dLimits.y1);
	else
		y = fabs(y0 - dLimits.y2);
	rmax = sqrt(x*x + y*y);
	r = r0;
	j = 0;
	
	while(r < rmax){
		for(int i=0;i<nRays;i++){
			theta = theta0  + (i + (j%2)*0.5)*deltheta + 0.001;
			x = x0 + r*cos(theta);
			y = y0 + r*sin(theta);
			nP = physics->makeNewNode(nodeNum,x,y);
			
			if( (tP = dataTIN.whichTriangle(nP)) != NULL){
				if( (tBP = boundTIN.whichTriangle(nP)) != NULL){
					if(tBP->getStatus() == active){
						tP->locateNodeAt(x,y,nP);
						nGoodNodes += 1;
						nodeNum += 1;
						tnodeL.appendItem(nP);
					}
					else
               	delete nP;
				}
				else
					delete nP;
			}
			else
				delete nP;
		}
		r *= fac;
		j++;
	}
	return nGoodNodes;
}

int TIN::refineRegion(TIN &boundTIN, TIN &dataTIN)

//	Refines a TIN by placing a new node in the centre of each existing triangle
//	within the boundTIN region.
//	Requires a previously triangulated data TIN, dataTIN,
//	which masks inactive areas and provides interpolation data.
//	Resulting set of nodes are in tnodeL which can be later accepted into nodeL or rejected.
//	The number of nodes created are returned.

{
	int nGoodNodes = 0;
	int nodeNum = getNextN();
	Triangle *tP, *tBP, *eTP;
	Node cNode, *nP;

	eTP = firstTri();
	while(eTP != NULL){
		if(eTP->getStatus() == active){
			eTP->locateNodeAtCenter(&cNode);
			nP = physics->makeNewNode(nodeNum,cNode.getXc(),cNode.getYc());
			if( (tP = dataTIN.whichTriangle(nP)) != NULL){
				if( (tBP = boundTIN.whichTriangle(nP)) != NULL){
					if(tBP->getStatus() == active){
						tP->locateNodeAt(nP);
						nGoodNodes += 1;
						nodeNum += 1;
						tnodeL.appendItem(nP);
					}
					else
						delete nP;
				}
				else
					delete nP;
			}
			else
				delete nP;
		}
      eTP = nextTri();
	}
	return nGoodNodes;
}

int TIN::refineRegionYellow(TIN &boundTIN, TIN &dataTIN, double minDz)

//	New function added 3/2001, JDS
//	Same as refineRegion but refines only large elevation difference triangles.
//
//	Refines a TIN by placing a new node in the centre of each existing triangle
//	within the boundTIN region.
//	Requires a previously triangulated data TIN, dataTIN,
//	which masks inactive areas and provides interpolation data.
//	Resulting set of nodes are in tnodeL which can be later accepted into nodeL or rejected.
//	The number of nodes created are returned.

{
	int nGoodNodes = 0;
	int nodeNum = getNextN();
	Triangle *tP, *tBP, *eTP;
	Node cNode, *nP;

	eTP = firstTri();
	while(eTP != NULL){
		if((eTP->getStatus() == active) && fabs(dataTIN.triCenterDz(eTP))>minDz){
			eTP->locateNodeAtCenter(&cNode);
			nP = physics->makeNewNode(nodeNum,cNode.getXc(),cNode.getYc());
			if( (tP = dataTIN.whichTriangle(nP)) != NULL){
				if( (tBP = boundTIN.whichTriangle(nP)) != NULL){
					if(tBP->getStatus() == active){
						tP->locateNodeAt(nP);
						nGoodNodes += 1;
						nodeNum += 1;
						tnodeL.appendItem(nP);
					}
					else
						delete nP;
				}
				else
					delete nP;
			}
			else
				delete nP;
		}
      eTP = nextTri();
	}
	return nGoodNodes;
}

int TIN::defineBoundary(double spacing, TIN &boundTIN, TIN &dataTIN)

//	Generates new boundary nodes and segments from physics spec.
//	Outline is provided by boundTIN.
//	Spacing is approximate (goal) distance between new nodes.
//	Resulting set of nodes are in tnodeL which can be later accepted into nodeL or rejected.
//	Resulting boundary segments are in tbsegL. Accepting segements automatically accepts nodes.
//	The number of segments (= total number of nodes on boundary) created are returned.

{
	int	i, nNewSegs;
	Segment *bP = (Segment *) boundTIN.firstSeg();
	Node *startLoopP;
	Segment *nSegP;
	Node *nodeP1, *nodeP2, *startNodeP;
	Triangle *tP;
	
	while(bP != NULL) {
		startLoopP = bP->getNode(0);
		nodeP1 = physics->makeNewNode(getNextN());
		startNodeP = nodeP1;
		nodeP1->setFixed();
		tnodeL.appendItem(nodeP1);
		bP->locateNodeAt(0.0,nodeP1);
		tP = dataTIN.whichTriangle(nodeP1);
		if(tP != NULL)
			tP->locateNodeAt(nodeP1);
		while(bP != NULL) {
			nNewSegs = 1 + (int) (bP->length() / spacing);
			for(i=1;i<=nNewSegs;i++){
				nodeP2 = physics->makeNewNode(getNextN());
				nodeP2->setSliding();
				tnodeL.appendItem(nodeP2);
				bP->locateNodeAt((1.0*i/nNewSegs),nodeP2);
				tP = dataTIN.whichTriangle(nodeP2);
				if(tP != NULL)
					tP->locateNodeAt(nodeP2);
				nSegP = physics->makeNewSegment(tbsegL.numberOfItems()+1,nodeP1,nodeP2,bP);
				tbsegL.appendItem(nSegP);
				nodeP1 = nodeP2;
			}
			nodeP1->setFixed();
			if(bP->getNode(1) == startLoopP){
				bP = (Segment *) boundTIN.nextSeg();
				break;
			}
			bP = (Segment *) boundTIN.nextSeg();
		}
		nSegP->setNode(1,startNodeP);
		tnodeL.deleteCurrentItem();
		delete nodeP2;
	}
	
	return tbsegL.numberOfItems();
}
int TIN::defineNewBoundLoop(double spacing, TIN &boundTIN, TIN &dataTIN, Segment* firstSP)

//	Generates new boundary nodes and segments from physics spec.
//	Outline is provided by boundTIN, starting with firstSP, which is presumably a new loop.
//	Spacing is approximate (goal) distance between new nodes.
//	Resulting set of nodes are in tnodeL which can be later accepted into nodeL or rejected.
//	Resulting boundary segments are in tbsegL. Accepting segements automatically accepts nodes.
//	The number of segments (= total number of nodes on boundary) created are returned.

{
	int	i, nNewSegs;
	Segment *bP = firstSP;
	Node *startLoopP;
	Segment *nSegP;
	Node *nodeP1, *nodeP2, *startNodeP;
	Triangle *tP;
	
	while(bP != NULL) {
		startLoopP = bP->getNode(0);
		nodeP1 = physics->makeNewNode(getNextN());
		startNodeP = nodeP1;
		nodeP1->setFixed();
		tnodeL.appendItem(nodeP1);
		bP->locateNodeAt(0.0,nodeP1);
		tP = dataTIN.whichTriangle(nodeP1);
		if(tP != NULL)
			tP->locateNodeAt(nodeP1);
		while(bP != NULL) {
			nNewSegs = 1 + (int) (bP->length() / spacing);
			for(i=1;i<=nNewSegs;i++){
				nodeP2 = physics->makeNewNode(getNextN());
				nodeP2->setSliding();
				tnodeL.appendItem(nodeP2);
				bP->locateNodeAt((1.0*i/nNewSegs),nodeP2);
				tP = dataTIN.whichTriangle(nodeP2);
				if(tP != NULL)
					tP->locateNodeAt(nodeP2);
				nSegP = physics->makeNewSegment(tbsegL.numberOfItems()+1,nodeP1,nodeP2,bP);
				tbsegL.appendItem(nSegP);
				nodeP1 = nodeP2;
			}
			nodeP1->setFixed();
			if(bP->getNode(1) == startLoopP){
				bP = (Segment *) bP->getNextOne();
				break;
			}
			bP = (Segment *) bP->getNextOne();
		}
		nSegP->setNode(1,startNodeP);
		tnodeL.deleteCurrentItem();
		delete nodeP2;
	}
	
	return tbsegL.numberOfItems();
}

//	cout <<"Bisecting segment " << segP; 
//	cout <<"New Node " << nodeP;
//	cout <<"New Segment " << nSegP;
//	writeTriangles(cout);
//	cin >> iin;
	
int TIN::orderNodesRCM()
{
	int ii, i, j, index, count = 0;
	Element *aElP;
	Triangle *atP, *tP;
	Node *nP1, *mnP, *tnP, *nP;

	for(ii =1;ii<= elmL.numberOfItems();ii++){			//	Renumber triangles
		tP = (Triangle *) elmL.i(ii);
		tP->setn((tP->getIndex())+1);
	}
/*	
	for(ii =1;ii<= nodeL.numberOfItems();ii++){			//	Renumber nodes
		nP = (Node *) nodeL.i(ii);
		nP->setn((nP->getIndex())+1);
	}
*/
//	Set Node indices to zero

	nP = firstNode();
	while(nP != NULL){
		nP->setIndex(0);
		nP->setn(0);
		nP = nextNode();
	}

//	Set index for each node to number of connecting elements

	tP = firstTri();
	while(tP != NULL){
		if(tP->getStatus() == active){
			for(i=0;i<3;i++){
				nP = tP->getNode(i);
				nP->setIndex(nP->getIndex() - 1);
			}
		}
		tP = nextTri();
	}
//	Initialize process
	ItemList tList, aList;

	nP = firstNode();
	if(nP != NULL){
		nodeL.deleteItem(nP);
		tList.appendItem(nP);
		nP->setn(1);
		count++;
	}

//	Main ordering loop, nodelist to tempList

	while(nP != NULL){

//	Get List of adjacent nodes

		if((tP = whichTriHasNode(nP)) != NULL) {
			atP = tP;
			do {
				for(j=0;j<3;j++) {
					nP1 = atP->getNode(j);
					if(nP1 == nP)
						index = j;
				}
				index += 1;
				if(index == 3)
					index = 0;
				aElP = atP->getAdjt(index);
				if(aElP==NULL)
					break;
				if(aElP->numNodes() == 2)
					break;
				atP = (Triangle *) aElP;
			} while (atP != tP);
			tP = atP;
			do {
				for(j=0;j<3;j++) {
					nP1 = atP->getNode(j);
					if(nP1 == nP)
						index = j;
					else {
						if(nP1->getn() == 0){
							nodeL.deleteItem(nP1);
							aList.appendItem(nP1);
							nP1->setn(1);
						}
					}
				}
				index -= 1;
				if(index == -1)
					index = 2;
				aElP = atP->getAdjt(index);
				if(aElP==NULL)
					break;
				if(aElP->numNodes() == 2)
					break;
				atP = (Triangle *) aElP;
			} while (atP != tP);
		}

// Put in tempList in order of increasing degree

		while(aList.numberOfItems() > 0){
			tnP = (Node *) aList.firstItem();
			mnP = tnP;
			while(tnP != NULL){
				if(tnP->getIndex() > mnP->getIndex())
					mnP = tnP;
				tnP = (Node *) aList.nextItem();
			}
			aList.deleteItem(mnP);
			tList.appendItem(mnP);
			count++;
		}

		nP = (Node *) nP->getNextOne();
	}

//	transfer back to nodeList in reverse order

	while(tList.numberOfItems() > 0){
		nodeL.push(tList.pop());
	}

	nP = (Node *) nodeL.firstItem();
	int nn = 1;
	while(nP != NULL){
		nP->setn(nn);
		nn++;
		nP = (Node *) nP->getNextOne();
	}

	return count;

}

//	Added a line to check for a NULL nsegP.  Returns 2 if so.
//	Add 4/2001 JDS

int	TIN::deleteNode(Node *nP)
{
	if(nP->getFixed() == sliding){
		Segment *segP = firstSeg();
		Segment *nsegP;
		while(segP != NULL) {
			if(segP->getNode(1) == nP){
				nsegP = (Segment *) segP->getNextOne();
				if(nsegP==NULL) return 2;	//	Added line.  JDS
				nsegP->setNode(0,segP->getNode(0));
				bsegL.deleteItem(segP);
				delete segP;
				nP->setInterior();
				break;
			}
			segP = nextSeg();
		}
		segP = firstfSeg();
		while(segP != NULL) {
			if(segP->getNode(1) == nP){
				nsegP = (Segment *) segP->getNextOne();
				nsegP->setNode(0,segP->getNode(0));
				fsegL.deleteItem(segP);
				delete segP;
            break;
			}
			segP = nextfSeg();
		}
	}
	if( (nP->getFixed() == floating) || (nP->getFixed() == sliding) ){
		nP->setDeleted();
		return 1;
	}
	else
		return 0;
}

void	TIN::removeFSeg(Segment *fSegP)
{
	fsegL.deleteItem(fSegP);
}

void	TIN::removeTFSeg(Segment *fSegP)
{
	tfsegL.deleteItem(fSegP);
}
