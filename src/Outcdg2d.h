
//		Outcdg2d.h
#ifndef OUTCDG2D_H
#define OUTCDG2D_H

#include "Tin.h"
#include "Habitat.h"
#include "FlowBoundary.h"

struct Control {
	int	trans	;
	int	meshtype;
	int	dims	;
	int	vars	;
	int	params	;
	int	bparams ;
	int	Keqns[12]	;
	int	nodes	;
	int	elms	;
	int	belms	;
	int	boundarysegs;
	} ;

struct Transient {
	int	nsteps ;
	double	dtfac ;
	double	t ;
	double	dt ;
	double	theta ;
	double	UW;
	double	DIF;
	double	minH;
	double	gwH;
	double	GWD;
	double	T;
	double	S;
	double	latitude ;
	int	diffusivewave ;
	int	plotcode;	
	int	transbcs;
	int	maxitnum;
	int	smallH;
	int	JE;
	double	uwJ;
	double epsilon1;
	double epsilon3;
	} ;

class OutCDG2D
{
	protected:
		TIN *theMeshP;
		RivBound *outSegP, *inSegP;
		Control con;
		Transient tvals;
		ItemList flowBoundL;
		double wsElevOut, wsElevIn;
		double getWsElev();
		void put_control(ofstream& f);
		void put_node(ofstream& f, ShallowNode *Np);
		void put_elm(ofstream& f, Triangle *Tp);
		void put_belm(ofstream& f, RivBound *Bp);
		void put_bsegs(ofstream& f);
		int num_bsegs();
		void getControl(ifstream& f);
		void getElm(ifstream& f);
		void getNode(ifstream& f);
		void getBSeg(ifstream& f);
		void getBElm(ifstream& f);
		int getInt(ifstream& f);
		double getDouble(ifstream& f);

	public:
		OutCDG2D(TIN *meshP, double elIn);
		~OutCDG2D();
		void output(ofstream& f);
		void input(ifstream& f);
		double getUW() const { return tvals.UW; }
		double getDIF() const { return tvals.DIF; }
		double getgwH() const { return tvals.gwH; }
		double getT() const { return tvals.T; }
		double getS() const { return tvals.S; }
		double getEpsilon1() const { return tvals.epsilon1;}
		double getEpsilon3() const { return tvals.epsilon3;}
		void setUW(double upwinding);
		void setDIF(double diffusivity);
		void setgwH(double gwHeight);
		void setT(double transmissivity);
		void setS(double storativity);
		void setEpsilon1(double Epsilon1);
		void setEpsilon3(double Epsilon3);
		void appendFBound(FlowBoundary *fBP);
		void removeFBound(FlowBoundary *fBP);
		int NFlowBounds() {return con.boundarysegs;}
		FlowBoundary* firstFBound() {return (FlowBoundary*) flowBoundL.firstItem();}
		FlowBoundary* lastFBound() {return (FlowBoundary*) flowBoundL.lastItem();}
		FlowBoundary* nextFBound() {return (FlowBoundary*) flowBoundL.nextItem();}
		void setBElementFlowBound(FlowBoundary* fBP, Node* startNP, Node* endNP);
		void buildFlowBoundaries();
		void sett(double time);
		void setdt(double deltat);
		void setmaxitnum(int maxitnum);
		void settheta(double theta);
		double gettheta() const { return tvals.theta;}
};

#endif
