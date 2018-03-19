#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#define  MNSF    6
#define  MRNSF    3
#define  MDOF    36
#define  MGPTS   12
#define  NDIMS    2
#define  NPARAMS  3
#define  NBPARAMS 7
#define  NVARS    3
#define  NKEQNS  20
#define	MNEIGHBOURS	51

struct shapefuncs {
	int	dof ;
	double	f[MNSF] ;
	double 	dfdr[MNSF] ;
	double	dfds[MNSF] ;
	double	dfdt[MNSF] ;
	double	d2fdr[MNSF] ;
	} ;

struct govfuncs {
	int	dof ;
	double	detJ ;
	double	w ;
	double	p[NPARAMS] ;
	double	f[MNSF] ;
	double	dfdx[MNSF] ;
	double	dfdy[MNSF] ;
	double	dfdz[MNSF] ;
	double	d2fdx[MNSF] ;
	} ;

struct gausspts {
	double	x ;
	double	y ;
	double	z ;
	double	w ;
	} ;

struct	node  {
	int		n ;
	int		i ;
	int		ui ;
	int		fxc	;
	double	x[NDIMS]	;
	double	u[NVARS]	;
	double	p[NPARAMS]	;
	double	ice[NPARAMS];
	double	uo[NVARS]	;
	double	uoo[NVARS]	;
	double	ud[NVARS];
	float	f[NDIMS]	;
	double	*blockDiagKp	;
	struct node	*nextnp	;

	double	depth;		// current depth of flow (m)
	double	deptho;		// previous depth of flow (m)
	double	concentration;	// current concentration (m^3/m^3)
	double	concentrationo;	// previous concentraton (m^3/m^3)
	int		dry;		// hysteresis flag for wet and dry
	double	DGW;
	double	DGWo;

	double K1, K2, K3, K4;

	//solution parameter for transport system
	double ucon[1];			// ucon[0] = concentration (m^3/m^3)
	double ucono[1];		

	} ;

struct	unode  {
	int		n ;
	short	int		fxc	;
	short	int		loc;
	int		nneighbours;
	int		nbelms;
	short	int		reject;
	short	int		incnp;
	double	x[NDIMS]	;
	double	p[2];
	double	u[NVARS];
	struct	tempelement	*neighbour[MNEIGHBOURS];
	struct	tempbelement	*belp[2];
	struct	node	*rnp;
	struct	unode	*nextunp;
	struct	unode	*previousunp;
	} ;

struct	closenode {
	int		nlamda;
	double	dist;
	double	alpha;
	struct	lamda		*firstl;
	struct	unode	*node;
	struct	closenode	*nextcnp;
	struct	closenode	*previouscnp;
}	;
	
struct	lamda {
	double	l;
	struct	lamda	*nextl;
}	;
 	
struct	element  {
	int	n    	;
	int	vtype	;
	int	gtype	;
	int	nnds	;
	struct node	*nps[MNSF];
	double	p[NPARAMS]	;
	double	*matrices ;
	float	elinfo ;
	int	gpcode[9];
	double FEn[9];
	struct element	*nextelp ;
	} ;

struct	belement  {
	int	n    	;
	int	vtype	;
	int	gtype	;
	int	nnds	;
	struct node	*nps[MRNSF];
	struct	element	*elp;
	struct	boundaryseg	*bseg;
	double	p[NBPARAMS]	;
	double	pcon[1];	// boundary conditions for transport processes
	double	*matrices ;
	double FEn[6];
	int	bcs[NVARS];
	struct belement	*nextbelp ;
	} ;

struct	boundaryseg	{
	int		n;
	int		bcs[NVARS];
	int		nstart;
	int		nend;
	double	p[NBPARAMS];
	struct	boundaryseg	*nextbseg;
	}	;

struct	tempelement  {
	int	n ;
	int	loc;
	struct unode	*nps[3];
	struct tempelement	*nexttelp ;
	} ;

struct	tempbelement  {
	int	n    	;
	int	bcs;
	double	p[NVARS];
	struct unode	*nps[2];
	struct tempbelement	*nexttbp;
	} ;


struct control {
	int	trans	;
	int	meshtype;
	int	dims	;
	int	vars	;
	int	params	;
	int	bparams ;
	int	Keqns[NKEQNS]	;
	int	nodes	;
	int	unodes	;
	int	uinodes;
	int	closest;
	int	elms	;
	int	telms	;
	int	ielms;
	int	belms	;
	int	boundarysegs;
	int	bsegs	;
	int	bnodes	;
	int	ibnodes;
	int	frontsegs;
	int	ordersegs;
	int	xsecs;
	int	datapts;
	int	frnds	;
	int	ukns	;
	int	sym	;
	int	nonlin	;

	// control variables for river transport modelling
	int transportNodes;			// nubmer of transport nodes

	} ;

struct eqnset {
	double	*Kp ;
	double	*Lp ;
	double	*Fp ;
	double	*Bp ;
	double	*Qp ;
	double	*Dp ;
	double	*Rp ;
	double	*Hp ;
	double	*Vp ;
	double	*Mp ;
	double  *Yp ;
	long	*diagp ;
	int	neqns ;
	int m;
	int	NelemInK;
	int assembleFlag;
	} ;
	
struct pointers {
	struct	node    *N ;
	struct	unode	*UN ;
	struct	unode	*UIN;
	struct	closenode	*C;
	struct	element *El ;
	struct	belement *B ;
	struct	boundaryseg	*BSEG;
	struct	tempbelement	*TB;
	struct	tempelement		*TE;
	struct	tempelement		*IE;
	struct	bsegment	*BS ;
	struct	frontseg	*FS;
	struct	orderseg	*OR;
	struct	xsec	*X;
	struct 	node **iptrs ;
	double	*S ;
	double	*M ;
	}   ;
	
struct elmpointers {
	double	*K ;
	double	*F ;
	double	*S ;
	double	*M ;
	double	*J ;
	double	*DM ;
	int	n  ;
	}  ;

struct fxnodes {
	int	i ;
	int	j[MDOF] ;
	double	k[MDOF] ;
	double	f ;
	} ;

struct transient {
	int	nsteps ;
	int	method ;
	int	iter ;
	int itnum ;
	int	crdnum;
	double	dtfac ;
	double	t ;
	double	dt ;
	double	theta ;
	double	tfinal;
	double	uc;
	double	crdval;
	double	minH;
	double	gwH;
	double	GWD;
	double	T;
	double	S;
	double	sg_ice;
	double itnumTransport;
	} ;

struct settings {
	double	latitude ;
	int	oldeqnew ;
	int	diffusivewave ;
	int	plotcode;	
	int	transbcs;
	int	maxitnum;
	int	smallH;
	int	JE;
	double	uwJ;
	int	plotn[4];
	double	plotv[4];
	
	} ;


struct	RegMesh {
	int	nx ;
	int	ny ;
	int	nbx ;
	int	nby ;
	int	elsinbx[MNSF] ;
	int	elsinby[MNSF] ;
	int	eltype ;
	int	maptype ;
	struct control	con ;
	struct blockmap	*firstmap ;
	} ;

struct	UnstrMesh {
	int	eltype;
	int	nloops;
	int	nd;
	int	nnodes;
	int	ordernodes[2*MDOF];
	int	nordersegs;
	int	nmax;
	double	ymin;
	double	ymax;
	double	dx;
	double	dy;
	struct	control	con;
	};

struct	bsegment	{
	int	n;
	int	loop;
	struct	unode	*nps[MRNSF];
	double	p[NVARS];
	int	bcs;
	struct	bsegment	*nextbsp;
	};
	
struct	frontseg  {
	int	n;
	struct unode	*nps[2];
	struct frontseg	*nextfsp ;
	struct frontseg	*previousfsp;
	} ;

struct	orderseg  {
	struct unode	*nps[2];
	struct orderseg	*nextosp ;
	} ;
	
struct	xsec	{
	int		n;
	int		ndatapts;
	int		nwet;
	int		ifirst;
	int		ilast;
	double	loc;
	double	xcoord;
	double	weight;
	double	sozf;
	double	wsl;
	double	Q;
	double	sigmahp;
	double	sigmaHp;
	double	qx;
	struct	datapoint	*dp[MDOF*2];
	struct	xsec	*nextxsp;
	struct	xsec	*previousxsp;
} ;

struct	datapoint {
	int		n;
	int		nwet;
	double	b;
	double	y;
	double	man;
	double	xcoord;
	double	ycoord;
	double	u;
	double	v;
	double	weight;
	double	h;
	double	hb;
	double	qxb;
	double	qx;
	double	qy;
} ;

struct blockmap {
	float	xnodes[MNSF] ;
	float	ynodes[MNSF] ;
	} ;

struct transportSettings {	// river transport modelling settings

	double	D;		// diffusion coefficient

	double theta;				// implicitness
	double tolerance;			// Newton Raphson solution tolerance;
	int	   itnum;				// number of Newton Raphson iterations
	int	   maxitnum;			// maximum number of Newton Raphson iterations
	

};
	
