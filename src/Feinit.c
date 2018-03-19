//#pragma options cpluscmt	
#include "Fe_PS.h"

extern struct control  N ;
extern struct pointers gp;
extern int	BandWidth ;
extern struct transient	tvals ;
extern struct RegMesh		Mesh ;

int		Nxb, Nyb ;
double	Utop, UW ,DIF, epsilon1, epsilon3;


int		localinit() 

{
	BandWidth = 69;
	Nxb = 5 ;
	Nyb = 5 ;
	Utop = 1.0 ;
	UW = 0.5 ;
	DIF = 0.21;
	epsilon1 = 0;
	epsilon3 = 0;
	
	return(0) ;
}

int Nodal_Values(nin)
int nin;
{
	struct node *np;
	int	i;
	
	np=gp.N;
	for(i=0;i<N.nodes;i++) {
		setnodevalue(np);
		np=np->nextnp;
	}
	return(0);
}

int	setnodevalue(np)
struct	node	*np;
{
//circular dambreak
	if(sqrt((np->x[0]-25.0)*(np->x[0]-25.0)+(np->x[1]-25.0)*(np->x[1]-25.0))<=11.0) 
		np->u[0]=10.0;
	else
		np->u[0]=5.0;
		np->u[1]=0.0;
		np->u[2]=0.0;
		np->p[0]=100.00;
		np->p[1]=0.0;
		np->p[2]=10.0;
	return(0);
}
int		nodevalues(np)
struct node	*np ;

{	double	sigma = 20.0, xo = 120.0, yo = 40.0;
	double	width = 60.0, depth = 3.0;
	double	So,So1,elevation,freeboard,Yo,sslope;
	double	topwidth,height,height1,trans1,trans2,trans3,outlet,reachlength,Qtotal,sslope1;

	
/* convert dimensions from centimeter to meter for small domains 
	np->x[0]/=100.0;
	np->x[1]/=100.0; 	*/
	/* phi values - A and Q */

/*	if(np->x[0]>=125.&&np->x[0]<=375.&&np->x[1]>=125.&&np->x[1]<=375.)
		np->u[0] = 10.0;
	else	*/
//	if(np->x[0]<=100)
	np->u[0] = (3.0-np->x[0]/100.0);
	if(np->u[0]<0.01)
		np->u[0]=0.01;
//	else
//		np->u[0] = 1.5;
	if(np->u[0]>0.01)
		np->u[1] =  0.1;
	else
		np->u[1]=0.01;
	np->u[2] =  0.0;
	
	/* parameters */
	/* roughness height */
	np->p[0] = 100.00-.0005*np->x[0]+1.8*exp(-((np->x[0]-xo)*(np->x[0]-xo)/2./sigma/sigma
				+ (np->x[1]-yo)*(np->x[1]-yo)/2./sigma/sigma));
	np->p[1] = 0.035 ;
	/* b */
	np->p[2] = 10.0 ;	

//	if ((np->n)>((N.nodes-11)/2.))
//		np->u[0] = 1.0;
//	else
	np->u[0] = 1.9 -1.8*exp(-((np->x[0]-xo)*(np->x[0]-xo)/2./sigma/sigma
				+ (np->x[1]-yo)*(np->x[1]-yo)/2./sigma/sigma));
	np->u[1] = 3.0;
	np->u[2] = 0.0;
/*
		np->u[0] =1.0+0.5*exp(-((np->x[0]-xo)*(np->x[0]-xo)/2./sigma/sigma
				+ (np->x[1]-yo)*(np->x[1]-yo)/2./sigma/sigma));		np->u[1] = 6.26/sqrt(2.0)  ;  np->u[2] = 6.26/sqrt(2.0) ;
*/
	np->u[0]=3.5;
	np->u[1]=10.0;
	np->u[2]=0.0;
	np->p[0]=100.00-0.0013345*np->x[0];
	np->p[1]=0.2;
	np->p[2]=10.0;
	
//trapezoidal section
/*
	So=0.0013345;
	elevation=100.00;
	freeboard=0.5;
	Yo=1.0;
	sslope=5.0;
	width = 10;
	if(np->x[1]<=(Yo*sslope))
		np->p[0]=elevation+Yo-np->x[1]/sslope-So*np->x[0];
	if((np->x[1]>(Yo*sslope))&&(np->x[1]<=(Yo*sslope+width)))
		np->p[0]=elevation-So*np->x[0];
	if(np->x[1]>(Yo*sslope+width))
		np->p[0]=elevation+(np->x[1]-(Yo*sslope+width))/sslope-So*np->x[0];
	
	np->u[0]=elevation+freeboard+Yo-So*np->x[0]-np->p[0];
	if((np->x[1]<(Yo*sslope-1))||(np->x[1]>(Yo*sslope+width+1)))
		np->u[1]=1.0;
	else
		np->u[1]=1.0;
	np->u[2]=0.0;
*/
// trapezoidal section with inlet and outlet transitions

	So=0.0013345;
	elevation=100.00;
	Yo=2.5;
	freeboard=0.5;
	width = 5.0;
	topwidth=35.0;
	trans1=50.0;
	trans2=50.0;
	trans3=50.0;
	outlet=50.0;
	reachlength=100.0;
	Qtotal=80;
	height = Yo-freeboard;
	sslope=(topwidth-width)/2.0/height;
	
// 	inlet rectangular reach
	if(np->x[0]<=trans1)
		np->p[0]=elevation-So*np->x[0];
//		transition from rectangular to trapezoidal
	if((np->x[0]>trans1)&&(np->x[0]<(trans1+trans2))){
		height1 = (np->x[0]-trans1)/trans2*height;
		sslope1=	(topwidth-width)/2.0/height1;
		if(np->x[1]<=((topwidth-width)/2.0))
			np->p[0]=elevation-So*np->x[0]+height1-np->x[1]/sslope1;
		if(np->x[1]>((topwidth-width)/2.0)&&(np->x[1]<=((topwidth+width)/2.0)))
			np->p[0]=elevation-So*np->x[0];
		if(np->x[1]>((topwidth+width)/2.0))
			np->p[0]=elevation-So*np->x[0]+(np->x[1]-(topwidth+width)/2.0)/sslope1;
	}
	
	if((np->x[0]>=(trans1+trans2))&&(np->x[0]<(trans1+trans2+reachlength))) {
		height1=height;
		sslope1=sslope;
		if(np->x[1]<=((topwidth-width)/2.0))
			np->p[0]=elevation-So*np->x[0]+height1-np->x[1]/sslope1;
		if(np->x[1]>((topwidth-width)/2.0)&&(np->x[1]<=((topwidth+width)/2.0)))
			np->p[0]=elevation-So*np->x[0];
		if(np->x[1]>((topwidth+width)/2.0))
			np->p[0]=elevation-So*np->x[0]+(np->x[1]-(topwidth+width)/2.0)/sslope1;
	}
	
	if((np->x[0]>=(trans1+trans2+reachlength))&&(np->x[0]<(trans1+trans2+reachlength+trans3))) {
		height1=height-(np->x[0]-(trans1+trans2+reachlength))/trans3*height;
		sslope1=	(topwidth-width)/2.0/height1;
		if(np->x[1]<=((topwidth-width)/2.0))
			np->p[0]=elevation-So*np->x[0]+height1-np->x[1]/sslope1;
		if(np->x[1]>((topwidth-width)/2.0)&&(np->x[1]<=((topwidth+width)/2.0)))
			np->p[0]=elevation-So*np->x[0];
		if(np->x[1]>((topwidth+width)/2.0))
			np->p[0]=elevation-So*np->x[0]+(np->x[1]-(topwidth+width)/2.0)/sslope1;
	}
	
	if(np->x[0]>=(trans1+trans2+reachlength+trans3)) {
		np->p[0]=elevation-So*np->x[0];
	}
	
	np->u[0]=elevation+Yo-So*np->x[0]-np->p[0];
	np->u[1]=Qtotal/topwidth*(np->u[0]/Yo);
	np->u[2]=0.0;

//circular dambreak
	if(sqrt((np->x[0]-25.0)*(np->x[0]-25.0)+(np->x[1]-25.0)*(np->x[1]-25.0))<=11.0) 
		np->u[0]=10.0;
	else
		np->u[0]=5.0;
		np->u[1]=0.0;
		np->u[2]=0.0;
		np->p[0]=100.00;
		np->p[1]=0.0;
		np->p[2]=10.0;
//rectangular bay with quadraticall varying bathymetry and periodic forcing function in x dir
	np->u[0]=3.048*(1.5*(np->x[0]/(6.0*15240)) + 1.0)*(1.5*(np->x[0]/(6.0*15240)) + 1.0);
	np->u[1]=0.0;
	np->u[2]=0.0;
	np->p[0]=100.0-np->u[0];
	np->p[1]=0.0;
	np->p[2]=10.0;

//rectangular bay with quadraticall varying bathymetry and periodic forcing function in y dir
	np->u[0]=3.048*(1.5*(np->x[1]/(6.0*15240)) + 1.0)*(1.5*(np->x[1]/(6.0*15240)) + 1.0);
	np->u[1]=0.0;
	np->u[2]=0.0;
	np->p[0]=100.0-np->u[0];
	np->p[1]=0.0;
	np->p[2]=10.0;
//two-dimensional hydraulic jump (horizontal frictionless channel)
	So=0.0013345;
	So1 = 0.1;
	if(np->x[0]<20.0) {
		np->u[0]=1.0;
		np->u[1]=10.0;
		np->p[0]=3.0-So1*np->x[0];
	}
	else {
		np->u[0]=2.0;
		np->u[1]=4.0;
		np->p[0]=3.0-So1*20.0-So*(np->x[0]-20.0);
	}
	np->u[2]=0.0;
	np->p[1]=0.1;
	np->p[2]=10.0;
// one-dimensional hj
	if(np->x[0]>75.0)
		np->u[0]=1.4052559;
	else
		np->u[0]=1.0;
	np->u[1]=4.0717195;
	np->u[2]=0.0;
	np->p[0]=0.0;
	np->p[1]=0.0;
	np->p[2]=10.0;
// one-dimensional db 0.5
	if(np->x[0]>1000.0)
		np->u[0]=5.0;
	else
		np->u[0]=10.0;
	np->u[1]=0.0;
	np->u[2]=0.0;
	np->p[0]=0.0;
	np->p[1]=0.0;
	np->p[2]=10.0;
// one-dimensional db 0.5 on 45 degrees
	if(sqrt(np->x[0]*np->x[0]+np->x[1]*np->x[1])>1000.0)
		np->u[0]=5.0;
	else
		np->u[0]=10.0;
	np->u[1]=0.0;
	np->u[2]=0.0;
	np->p[0]=0.0;
	np->p[1]=0.0;
	np->p[2]=10.0;
// one-dimensional hj on 45 degrees
	if(sqrt(np->x[0]*np->x[0]+np->x[1]*np->x[1])>75)
		np->u[0]=1.4052559;
	else
		np->u[0]=1.0;
	np->u[1]=2.87914047;
	np->u[2]=2.87914047;
	np->p[0]=0.0;
	np->p[1]=0.0;
	np->p[2]=10.0;
// Schoklitsch dambreak experiment
	if(np->x[0]>110.0)
		np->u[0]=0.0;
	else
		np->u[0]=0.074;
	np->u[1]=0.0;
	np->u[2]=0.0;
	np->p[0]=0.0;
	np->p[1]=0.00148;
	np->p[2]=10.0;

//supercritical contraction (Ippen)

	np->u[0]=0.03058;
	np->u[1]=0.06686;
	np->u[2]=0.0;
	np->p[0]=100.0;
	np->p[1]=0.0;
	np->p[2]=10.0;

	return(0) ;
  
}

int		elvalues(elp)
struct element	*elp ;

{
	return(0) ;
}

/*FAYE: 88-08-28 boundary parameters*/

int		leftb(belp)
struct belement	*belp ;

{
	
	belp->bcs[0] = 2 ;
	bndparm (belp) ;
	
	return(0) ;
}

int		rightb(belp)
struct belement	*belp ;

{

	belp->bcs[0] = 4 ;
	bndparm (belp) ;
	
	return(0) ;
}

int		bottomb(belp)
struct belement	*belp ;

{
	belp->bcs[0] = 0 ;
	bndparm(belp);

	return(0) ;
}

int		topb(belp)
struct belement	*belp ;

{
	belp->bcs[0] = 0 ;
	bndparm(belp);

	return(0) ;
}


int		bndparm (belp)
struct 	belement	*belp ;

{
 
	if (belp->bcs[0] == 0){ /*solid boundary: set Qx and Qy equal to zero */
		if(belp->nps[0]->x[0] <= 1.0 || belp->nps[0]->x[0] >= 499) {
			if(belp->nps[1]->x[0] <= 125 || belp->nps[1]->x[0] >= 375)
				belp->p[0] = 0.0;
			else
				belp->p[0] = 0.0 ;
		}
		else
		belp->p[0] = 0.0;
		belp->p[1] = 0.0 ;
		belp->p[2] = 0.0 ;
	}
	if (belp->bcs[0] == 1){ /*inflow boundary (subcritical)*/
		belp->p[0] = 3.5 ;
		belp->p[1] = 10 ;
		belp->p[2] = 0.0 ;
	}
	if (belp->bcs[0] == 2){ /*inflow boundary (supercritical)*/
		belp->p[0] = 0.03058 ;
		belp->p[1] = 0.06686 ;
		belp->p[2] = 0.0 ;
	}
	if (belp->bcs[0] == 3){ /*outflow boundary (subcritical)*/
		belp->p[0] = 1.4052559 ;
		belp->p[1] = 4.0717195 ;
		belp->p[2] = 0.0 ;
	}
	if (belp->bcs[0] == 4){ /*outflow boundary (supercritical)*/
		belp->p[0] =1.5 ;
		belp->p[1] =10.0 ;
		belp->p[2] = 0.0 ;
	}
	return(0) ;
}

double		uexact(nvar,x,y,z,p)
int		nvar ;
double		x, y, z, p ;

{
	double	t, s0, x0, xb, s2 ;
	s0 = 400.0 ;
	x0 = 2000.0 ;
	xb = x0+0.5*tvals.t ;
	s2 = s0*s0 + 2.0*p*tvals.t ;
	t = s0/sqrt(s2)*exp(-(x-xb)*(x-xb)/2.0/s2) ;
	return(t) ;
}

int		update_eta()

{
	int	i;
	
	for(i=0;i<N.nodes;i++) {
		gp.iptrs[i]->p[3] = gp.iptrs[i]->u[1] ;
	}
	return(0) ;
}

int		set_etabc() 

{
	int	i, j ;
	double	dn, eta_wall() ;
	struct belement		*belp ;
	struct node		*nw, *ni ;
	
	belp = gp.B ;
	for(i=0;i<Nxb;i++) {
		for(j=1;j<belp->nnds;j++) {
			if((i == Nxb-1) && (j == belp->nnds-1))
				break ;
			nw = belp->nps[j] ;
			ni = gp.iptrs[(nw->i)+1] ;
			dn = ni->x[1] - nw->x[1] ;
			nw->u[1] = eta_wall(nw->u[0],ni->u[0],dn,0.0,nw->u[1],ni->u[1]);
		}
		belp = belp->nextbelp ;
	}
	nw = belp->nps[0] ;
	ni = gp.iptrs[(nw->i)-Nyb] ;
	dn = sqrt((nw->x[0]-ni->x[0])*(nw->x[0]-ni->x[0])+(nw->x[1]-ni->x[1])*(nw->x[1]-ni->x[1])) ;
	nw->u[1] = eta_wall(nw->u[0],ni->u[0],dn,0.0,nw->u[1],ni->u[1]) ;
	for(i=0;i<Nyb;i++) {
		for(j=1;j<belp->nnds;j++) {
			if((i == Nyb-1) && (j == belp->nnds-1))
				break ;
			nw = belp->nps[j] ;
			ni = gp.iptrs[(nw->i)-Nyb-1] ;
			dn = nw->x[0] - ni->x[0] ;
			nw->u[1] = eta_wall(nw->u[0],ni->u[0],dn,0.0,nw->u[1],ni->u[1]) ;
		}
		belp = belp->nextbelp ;
	}
	nw = belp->nps[0] ;
	ni = gp.iptrs[(nw->i)-1] ;
	dn = sqrt((nw->x[0]-ni->x[0])*(nw->x[0]-ni->x[0])+(nw->x[1]-ni->x[1])*(nw->x[1]-ni->x[1])) ;
	nw->u[1] = eta_wall(nw->u[0],ni->u[0],dn,Utop/2.0,nw->u[1],ni->u[1]) ;
/*
	printf(" %d %d %f %f\n",nw->n,ni->n,dn,nw->u[1]) ;
*/
	for(i=0;i<Nxb;i++) {
		for(j=1;j<belp->nnds;j++) {
			if((i == Nxb-1) && (j == belp->nnds-1))
				break ;
			nw = belp->nps[j] ;
			ni = gp.iptrs[(nw->i)-1] ;
			dn = nw->x[1] - ni->x[1] ;
			nw->u[1] = eta_wall(nw->u[0],ni->u[0],dn,Utop,nw->u[1],ni->u[1]) ;
		}
		belp = belp->nextbelp ;
	}
	nw = belp->nps[0] ;
	ni = gp.iptrs[(nw->i)-1] ;
	dn = sqrt((nw->x[0]-ni->x[0])*(nw->x[0]-ni->x[0])+(nw->x[1]-ni->x[1])*(nw->x[1]-ni->x[1])) ;
	nw->u[1] = eta_wall(nw->u[0],ni->u[0],dn,0.0,nw->u[1],ni->u[1]) ;
	for(i=0;i<Nyb;i++) {
		for(j=1;j<belp->nnds;j++) {
			if((i == Nyb-1) && (j == belp->nnds-1))
				break ;
			nw = belp->nps[j] ;
			ni = gp.iptrs[(nw->i)+Nyb+1] ;
			dn = ni->x[0] - nw->x[0] ;
			nw->u[1] = eta_wall(nw->u[0],ni->u[0],dn,0.0,nw->u[1],ni->u[1]) ;
		}
		belp = belp->nextbelp ;
	}
	belp = gp.B ;
	nw = belp->nps[0] ;
	ni = gp.iptrs[(nw->i)+Nyb+2] ;
	dn = sqrt((nw->x[0]-ni->x[0])*(nw->x[0]-ni->x[0])+(nw->x[1]-ni->x[1])*(nw->x[1]-ni->x[1])) ;
	nw->u[1] = eta_wall(nw->u[0],ni->u[0],dn,0.0,nw->u[1],ni->u[1]);
	return(0) ;
}

double		eta_wall(pw,pi,h,u,ew,ei)
double		pw,pi,h,u,ew,ei ;

{
	double	t, r ;
	
	r = 0.005 ;
/*
	if(tvals.t < 1.0)
		u *= tvals.t ;
*/
	t = (1-r)*ew + r*(-3.0*(pi-pw)/h/h - ei/2.0 - 3.0*u/h) ;
	return(t) ;
}


int		update_p(np)
struct node	*np ;

{
	np->uoo[0] = np->uo[0] ;
	np->uo[0] = np->u[0] ;
	np->uoo[1] = np->uo[1] ;
	np->uo[1] = np->u[1] ;
	np->uo[2] = np->u[2] ;
	np->uoo[2] = np->uo[2];


	return(0) ;
}
