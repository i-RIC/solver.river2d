
//		outcdg2d.cpp

//		This file is for the purpose of outputing a CDG2D
//		data file from a generated mesh. Just instantiate
//		with a pointer to the generated mesh and an estimate of
//		the inflow elevation for the initial guess.
//		Call output with the ostream to generate the file.

#include <iostream>

#include <stdio.h>
#include <math.h>
#include <fstream>
#include <ctype.h>
#include "Tin.h"
#include "Outcdg2d.h"
#include "FlowBoundary.h"
using namespace std;

OutCDG2D::OutCDG2D(TIN *meshP, double elIn)
{
	theMeshP = meshP;
	wsElevIn = elIn;
	con.trans = 1;
	con.meshtype = 1;
	tvals.nsteps = 1;
    tvals.dtfac = 1.0;
    tvals.t = 0.0;
    tvals.dt = 10.0;
    tvals.theta = 1.0;
    tvals.UW = 0.5;
    tvals.DIF = 0.5;
	tvals.epsilon1 = 0;
	tvals.epsilon3 = 0;
    tvals.latitude = 0.0;
    tvals.S = 1.0;
    tvals.diffusivewave = 0;
    tvals.uwJ = 0.0;
    tvals.plotcode = 2;
    tvals.transbcs = 0;
    tvals.maxitnum = 9;
    tvals.smallH = 1;
    tvals.JE = 1;
    tvals.minH = 0.01;
    tvals.gwH = 0.01;
    tvals.GWD = 0.1;
    tvals.T = 0.1;
	con.dims = 2;
	con.vars = 3;
	for(int ii=0;ii<((con.vars+1)*con.vars);ii++)
		con.Keqns[ii] = 1;
    con.params = 3;
    con.bparams = 7;
    con.nodes = 0;
    con.elms = 0;
    con.belms = 0;
    con.boundarysegs = 0;
}

OutCDG2D::~OutCDG2D()
{
	FlowBoundary *flowBP = (FlowBoundary *) flowBoundL.firstItem();
	while(flowBP != NULL)
	{
		flowBoundL.deleteCurrentItem();
		delete flowBP;
		flowBP = (FlowBoundary *) flowBoundL.firstItem();
	}
}

void OutCDG2D::output(ofstream& f)
{
	ShallowNode	*np ;
	Triangle	*elp ;
	RivBound *belp ;

	put_control(f);
	f << "\n Node Information \n";
	f <<"\n Node #, Coordinates, Parameters, Variables\n\n";
	wsElevOut = getWsElev() ;
	np = (ShallowNode *) theMeshP->firstNode() ;
	while (np != NULL) {
		put_node(f,np) ;
		np = (ShallowNode *) theMeshP->nextNode() ;
	}

	f << "\n Element Information \n";
	f << "\n Element #, vtype, gtype, nodes\n\n";
	elp = theMeshP->firstTri() ;
	while(elp != NULL) {
		put_elm(f,elp) ;
		elp = theMeshP->nextTri() ;
	}

	f<<"\n Boundary Element #, vtype, gtype, nodes, boundary condition codes\n";
	belp = (RivBound *) theMeshP->firstSeg() ;
	while(belp != NULL) {
		put_belm(f,belp) ;
		belp = (RivBound *) theMeshP->nextSeg() ;
	}

	f<<"\n Boundary Seg #,Boundary type,stage,QT,start node #,end node #\n\n";
	put_bsegs(f);

	theMeshP->writeFeatureSegs(f);

	return;
}

double OutCDG2D::getWsElev()
{
	int countOut = 0;
	double elOut = 0.0;

	RivBound *bP = (RivBound *) theMeshP->firstSeg();
	while(bP != NULL){
		if((bP->getBcCode() == 3) || (bP->getBcCode() == 5)){
			outSegP = bP;
			countOut += 1;
			elOut += bP->getBcValue();
		}
		else if(bP->getBcCode() == 1){
			inSegP = bP;
		}
		bP = (RivBound *) theMeshP->nextSeg();
	}
	if(countOut > 0)
		return elOut/countOut;
	else
		return 0.0;
}

void OutCDG2D::put_control(ofstream& f)
{
	f << " Transient analysis = " << con.trans << "\n" ;
	f << " Mesh type = " << con.meshtype << "\n" ;
	f << " Number of Time Steps = " << tvals.nsteps << "\n" ;
	f << " Delta t Acceleration Factor = " << tvals.dtfac << "\n";
	f << " Time = " << tvals.t << "\n" ;
	f << " Delta t = " << tvals.dt << "\n";
	f << " Theta = " << tvals.theta << "\n";
	f << " UW = " << tvals.UW << "\n";
	f << " Eddy Viscosity Bed Shear Parameter = " << tvals.DIF << "\n";
	f << " Latitude \t\t = " << tvals.latitude << "\t degrees\n";
	f << " Groundwater Storativity\t\t = " << tvals.S << "\n";
	f << " Diffusive wave Solution = " << tvals.diffusivewave << "\t zero for fully dynamic only\n";
	f << " UW Jacobian terms included = " << tvals.uwJ << "\t zero for not included\n";
	f << " Plot Code \t\t\t\t\t= " << tvals.plotcode << "\t zero for xsec one for contour two for velocity and three for threeD\n";
	f << " Transient Boundary Condition = " << tvals.transbcs << "\t zero for Steady BCs\n";
	f << " Maximum Number of Iterations = " << tvals.maxitnum << "\n";
	f << " Small Depths Occur \t= " << tvals.smallH << "\t zero for no small depth calculations\n";
	f << " Jacobian Terms included = " << tvals.JE << "\t zero for not included\n";
	f << " Eddy Viscosity Constant = " << tvals.epsilon1 << "\n";
	f << " Minimum Depth for Groundwater Flow Calculation = " << tvals.gwH << "\n";
	f << " Eddy Viscosity Horizontal Shear Paramter = " << tvals.epsilon3 << "\n";
	f << " Groundwater Transmissivity\t\t = " << tvals.T << "\n";

	f << " Dimensions = " << con.dims << "\n";
	f << " Number of Variables = " << con.vars << "\n";
	f << " [K] governing equation numbers \n" ;
		for(int i=0;i<3;i++) {
			for(int ii=0;ii<3;ii++) {
				f << "\t 1 ";
			}
			f << "\t\t 1 \n";
		}
	f << " Number of Parameters = " << con.params << "\n";
	f << " Number of Boundary Parameters = " << con.bparams << "\n";
	f << " Number of Nodes = " << theMeshP->Nnodes() << "\n";
	Triangle *elp = theMeshP->firstTri() ;
	int eCount = 0;
	while(elp != NULL) {
		if(elp->getStatus() == active)
			eCount += 1; ;
		elp = theMeshP->nextTri() ;
	}
	f << " Number of Elements = " << eCount << "\n";
	f << " Number of  Boundary Elements = " << theMeshP->Nbsegs() << "\n";
	f << " Number of  Boundary Segments = " << con.boundarysegs << "\n";
	f <<  " \n" ;
	
	return;
}

void OutCDG2D::put_node(ofstream& f, ShallowNode *nP)
{
	f << nP->getn() << "\t";
	if(nP->getFixed() == fixednode)
		f << "x ";
	else if (nP->getFixed() == sliding)
		f << "s ";
	else
      	f << "  ";
	f.setf(ios::fixed,ios::floatfield);
	f << nP->getXo() << "\t";
	f << nP->getYo() << "\t";
	f.setf(ios::fmtflags(0),ios::floatfield); //modified by C.Frias June 10 of 2010
	f << nP->getZc() << "\t";
	f << nP->getKs() << "\t";
	f << "10 \t";
	if( wsElevOut <= 0.0){
		f << wsElevIn << "\t";
		f << 0.0 << "\t";
		f << 0.0 << "\n";
	}
	else if( wsElevIn > 0.0 ){
		double dout = fabs(outSegP->whichSide(nP))/outSegP->length();
		double din =  fabs(inSegP->whichSide(nP))/inSegP->length();
		double r = dout/(din + dout);
		f << (r*wsElevIn + (1-r)*wsElevOut - (nP->getZc()) ) << "\t";
		f << "0.0\t0.0\n";
	}
	else{
		f << nP->getD() << "\t";
		f << nP->getPar(4) << "\t";
		f << nP->getPar(5) << "\n";
	}
}

void OutCDG2D::put_elm(ofstream& f,Triangle *tP)
{
	if(tP->getStatus() == active){
		f << tP->getn() << "\t";
		f << "210\t210\t";
		f << tP->getNode(0)->getn() << "\t";
		f << tP->getNode(1)->getn() << "\t";
		f << tP->getNode(2)->getn() << "\t";
		f << "0.0\t0.0\t0.0\n";
	}
}

void OutCDG2D::put_belm(ofstream& f,RivBound *bP)
{
	f << bP->getn() << "\t";
	f << "111\t111\t";
	f << bP->getNode(0)->getn() << "\t";
	f << bP->getNode(1)->getn() << "\t";
	if(bP->getBcCode() == 0)
		f << "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
	else if(bP->getBcCode() == 1)
		f << "0\t" << bP->getBcValue()<< "\t0\t0\t0\t0\t0\t1\t0\t0\n";
	else if(bP->getBcCode() == 3)
		f << bP->getBcValue()<< "\t0\t0\t0\t0\t0\t0\t3\t0\t0\n";
	else if(bP->getBcCode() == 5)
		f << bP->getBcValue()<< "\t" <<bP->getBcValue2()<< "\t0\t0\t0\t0\t0\t5\t0\t0\n";
}

void OutCDG2D::put_bsegs(ofstream& f)
{
	FlowBoundary *fBP = firstFBound();
	while(fBP != NULL)
	{
		f << fBP->getn() << "\t";
		f << fBP->getBcCode() << "\t";
		f << fBP->getBcValue() << "\t";
		f << fBP->getBcValue2() << "\t";
		f << ((Node*)fBP->getStartNode())->getn() << "\t";
		f << ((Node*)fBP->getEndNode())->getn();
		if(fBP->getTransCode() == 0)
			f << "\n";
		else
			f << "\t" << fBP->getFilePath().c_str() << "\n";
		fBP = nextFBound();
	}
}
void OutCDG2D::buildFlowBoundaries()
{
	int Nseg = 0, oldBcCode=0, count = 0;
	double sum = 0.0;
	RivBound *bP, *startBP, *oldBP;
	FlowBoundary *fBP;

	bP = (RivBound *) theMeshP->firstSeg();
	while(bP != NULL){
		if(bP->getBcCode() != oldBcCode){
			if(oldBcCode == 1){
				Nseg += 1;
				fBP = new FlowBoundary(Nseg,startBP->getNode(0), oldBP->getNode(1), 1, 0, sum);
				flowBoundL.appendItem(fBP);
				con.boundarysegs++;
			}
			if(oldBcCode == 3){
				Nseg += 1;
				fBP = new FlowBoundary(Nseg,startBP->getNode(0), oldBP->getNode(1), 3, startBP->getBcValue(), 0);
				flowBoundL.appendItem(fBP);
				con.boundarysegs++;
			}
			if(oldBcCode == 5){
				Nseg += 1;
				fBP = new FlowBoundary(Nseg,startBP->getNode(0), oldBP->getNode(1), 5, startBP->getBcValue(),startBP->getBcValue2() );
				flowBoundL.appendItem(fBP);
				con.boundarysegs++;
			}
			oldBcCode = bP->getBcCode();
			startBP = bP;
			count = 0;
			sum = 0.0;
		}
		count += 1;
		sum += bP->getBcValue() * bP->length();
		oldBP = bP;
		bP = (RivBound *) theMeshP->nextSeg();
	}
	if(oldBcCode == 1){
		Nseg += 1;
		fBP = new FlowBoundary(Nseg,startBP->getNode(0), oldBP->getNode(1), 1, 0, sum);
		flowBoundL.appendItem(fBP);
		con.boundarysegs++;
		}
	if(oldBcCode == 3){
		Nseg += 1;
		fBP = new FlowBoundary(Nseg,startBP->getNode(0), oldBP->getNode(1), 3, startBP->getBcValue(), 0);
		flowBoundL.appendItem(fBP);
		con.boundarysegs++;
		}
	if(oldBcCode == 5){
		Nseg += 1;
		fBP = new FlowBoundary(Nseg,startBP->getNode(0), oldBP->getNode(1), 5, startBP->getBcValue(),startBP->getBcValue2() );
		flowBoundL.appendItem(fBP);
		con.boundarysegs++;
	}
}

int OutCDG2D::num_bsegs()
{
	int Nseg = 0, oldBcCode=0, count = 0;
	double sum = 0.0;

	RivBound *bP = (RivBound *) theMeshP->firstSeg();
	while(bP != NULL){
		if(bP->getBcCode() != oldBcCode){
			if(oldBcCode > 0){
				Nseg += 1;
			}
			oldBcCode = bP->getBcCode();
		}
		bP = (RivBound *) theMeshP->nextSeg();
	}
	if(oldBcCode > 0){
		Nseg += 1;
	}
	return Nseg;
}


void OutCDG2D::input(ifstream& f)
{
	getControl(f);
	for(int i=0;i<con.nodes;i++)
		getNode(f);
	for(int j=0;j<con.elms;j++)
		getElm(f);
	for(int k=0;k<con.belms;k++)
		getBElm(f);
	for(int l=0;l<con.boundarysegs;l++)
		getBSeg(f);
	theMeshP->readFeatureSegs(f);
}

int OutCDG2D::getInt(ifstream& f)
{
	char c;
	int i;

	while(!isdigit( c = f.peek())){
		if(c == '-')
			break;
		f.get(c);
	}
	f >> i;
	return i;
}

double OutCDG2D::getDouble(ifstream& f)
{
	char c;
	double d;

	while(!isdigit(c = f.peek())){
		if(c == '-')
			break;
		f.get(c);
	}
	f >> d;
	return d;
}

void OutCDG2D::getControl(ifstream& f)
{
	con.trans = getInt(f);
	con.meshtype = getInt(f);
	if(con.trans == 1){
		tvals.nsteps = getInt(f) ;
        tvals.dtfac = getDouble(f) ;
        tvals.t = getDouble(f) ;
        tvals.dt = getDouble(f) ;
        tvals.theta = getDouble(f) ;
        tvals.UW = getDouble(f) ;
        tvals.DIF     = getDouble(f);
        tvals.latitude = getDouble(f);
        tvals.S = getDouble(f);
        tvals.diffusivewave = getInt(f);
        tvals.uwJ = getDouble(f);
        tvals.plotcode = getInt(f);
        tvals.transbcs = getInt(f);
        tvals.maxitnum = getInt(f);
        tvals.smallH = getInt(f);
        tvals.JE =getInt(f);
//      tvals.minH = getDouble(f);
		tvals.epsilon1 = getDouble(f);
        tvals.gwH = getDouble(f);
//      tvals.GWD = getDouble(f);
		tvals.epsilon3 = getDouble(f);
        tvals.T = getDouble(f);
	}
	con.dims = getInt(f);
	con.vars = getInt(f) ;
	for(int ii=0;ii<((con.vars+1)*con.vars);ii++)
		con.Keqns[ii] = getInt(f) ;
    con.params = getInt(f) ;
    con.bparams = getInt(f) ;
    con.nodes = getInt(f) ;
    con.elms = getInt(f) ;
    con.belms = getInt(f) ;
    con.boundarysegs = getInt(f) ;
}

void OutCDG2D::getNode(ifstream& f)
{
	int n = getInt(f);
	char c;
	fixedFlag fx = floating;
	while(f) {
		c = f.peek();
		if( (c == '-') || (c == '.') || (isdigit(c)) )
			break;
		else if (c == 'x')
			fx = fixednode;
		else if (c == 's')
			fx = sliding;
		c = f.get();
	}
	double x = getDouble(f);
	double y = getDouble(f);
	double z = getDouble(f);
	double k = getDouble(f);
	double w = getDouble(f);
	double d = getDouble(f);
	double qx = getDouble(f);
	double qy = getDouble(f);

	HabitatNode* nP = new HabitatNode((Habitat*)theMeshP->getPhysics(),n,x,y,z,k,5.0,d,qx,qy);
	if(fx == fixednode)
		nP->setFixed();
	if(fx == sliding)
		nP->setSliding();
	theMeshP->appendNode(nP);
	return;
}

void OutCDG2D::getElm(ifstream& f)
{
	int n = getInt(f);
	int vt = getInt(f);
	int gt = getInt(f);
	int nn1 = getInt(f);
	int nn2 = getInt(f);
	int nn3 = getInt(f);
	double p1 = getDouble(f);
	double p2 = getDouble(f);
	double p3 = getDouble(f);
	return;
}

void OutCDG2D::getBElm(ifstream& f)
{
	int n = getInt(f);
	int vt = getInt(f);
	int gt = getInt(f);
	int nn1 = getInt(f);
	int nn2 = getInt(f);
	double p1 = getDouble(f);
	double p2 = getDouble(f);
	double p3 = getDouble(f);
	double p4 = getDouble(f);
	double p5 = getDouble(f);
	double p6 = getDouble(f);
	double p7 = getDouble(f);
	int c1 = getInt(f);
	int c2 = getInt(f);
	int c3 = getInt(f);
	Node *n1P = theMeshP->getNodeByN(nn1);
	Node *n2P = theMeshP->getNodeByN(nn2);
	double bValue = 0.0;
	if(c1 == 1)
		bValue = p2;
	else if(c1 == 3)
		bValue = p1;
	RivBound* bSeg = new RivBound(n,n1P,n2P,c1,bValue);
	theMeshP->appendBSeg(bSeg);
	return;
}

void OutCDG2D::getBSeg(ifstream& f)
{
	int n = getInt(f);
	int code = getInt(f);
	double p0 = getDouble(f);
	double p1 = getDouble(f);
	int startnode = getInt(f);
	int endnode = getInt(f);

	char path[200] = "";
	char c;
	f.get(c);
	if(c != '\n')
	{
		f.get(path,200,'\n');
		f.get(c);
		while(c != '\n')
			f.get(c);
	}
	string filePath = path; //changed CString to string May 14
	int transCode;
	if(filePath.c_str() == "")
		transCode = 0;
	else
		transCode = 1;

	Node* start = theMeshP->getNodeByN(startnode);
	Node* end	= theMeshP->getNodeByN(endnode);

	FlowBoundary *flowBP = new FlowBoundary(n,start,end,code,p0,p1,transCode,filePath);
	flowBoundL.appendItem(flowBP);
	setBElementFlowBound(flowBP, start, end);
	ifstream bcfile(filePath.c_str(), ios::in);//modified by C.Frias June 10 2010
	if(!bcfile.fail())
		flowBP->loadTransLst(bcfile);
	return;
}

void OutCDG2D::setBElementFlowBound(FlowBoundary* fBP, Node* startNP,Node* endNP)
{

	RivBound *segP;

	segP = (RivBound*)theMeshP->firstSeg();
	if (segP == NULL)
		return;
	while (segP->getNode(0) != startNP)
	{
		segP = (RivBound*)theMeshP->nextSeg();
		if (segP == NULL)
			segP = (RivBound*)theMeshP->firstSeg();
	}

	segP->setFlowBound(fBP);

	while (segP->getNode(1) != endNP)
	{
		segP = (RivBound*)theMeshP->nextSeg();
		if (segP == NULL)
			segP = (RivBound*)theMeshP->firstSeg();
		segP->setFlowBound(fBP);
	}
	return;
}

void OutCDG2D::appendFBound(FlowBoundary *fBP)
{
	flowBoundL.appendItem(fBP);
	flowBoundL.buildIndex();
	con.boundarysegs++;
}
void OutCDG2D::removeFBound(FlowBoundary *fBP)
{
	flowBoundL.deleteItem(fBP);
	flowBoundL.buildIndex();
	con.boundarysegs--;
	fBP = firstFBound();
	int i = 1;
	while (fBP != NULL)
	{
		fBP->setn(i);
		i++;
		fBP = nextFBound();
	}	
}	

void OutCDG2D::setUW(double upwinding)
{
	tvals.UW = upwinding;
}

void OutCDG2D::setDIF(double diffusivity)
{
	tvals.DIF = diffusivity;
}

void OutCDG2D::setgwH(double gwHeight)
{
	tvals.gwH = gwHeight;
}

void OutCDG2D::setT(double transmissivity)
{
	tvals.T = transmissivity;
}

void OutCDG2D ::setS(double storativity)
{
	tvals.S = storativity;
}

void OutCDG2D::setEpsilon1(double Epsilon1)
{
	tvals.epsilon1 = Epsilon1;
}

void OutCDG2D::setEpsilon3(double Epsilon3)
{
	tvals.epsilon3 = Epsilon3;
}

void OutCDG2D::sett(double time)
{
	tvals.t = time;
}

void OutCDG2D::setdt(double deltat)
{
	tvals.dt = deltat;
}

void OutCDG2D::settheta(double theta)
{
	tvals.theta = theta;
}

void OutCDG2D::setmaxitnum(int maxitnum)
{
	tvals.maxitnum = maxitnum;
}
