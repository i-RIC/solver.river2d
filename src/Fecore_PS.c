#include "Fe_PS.h"
#include <time.h>

#define KE(A,B) ( *(ep.K + ep.n*(A) + (B)) )
#define SE(A,B) ( *(ep.S + ep.n*(A) + (B)) )
#define ME(A,B) ( *(ep.M + ep.n*(A) + (B)) )
#define JE(A,B) ( *(ep.J + ep.n*(A) + (B)) )
#define DM(A,B,C) ( *(ep.DM + ep.n*((A) + (B)) + (C)) )
#define FE(B)   ( *(ep.F + B) )
#define GRAV    9.806
#define E	2.718

struct control          N ;
struct pointers         gp;
struct RegMesh          Mesh ;
struct fxnodes          *first_fnp, *next_fnp ;
struct transient        tvals ;
int                     Nukns, Nrfixed,  Fac, *BFnodes, BandWidth, assembleFlag, NelemInK ;
struct elmpointers      ep ;
struct eqnset           eqnsets[4] ; 
struct settings sets;
int						lumpS;
double                  uchange, maxChange,*fprime, memJacobian;
int						maxNode,maxVar;
int						gpIndex;
extern double           UW ,DIF, epsilon1, epsilon3;
extern double           defndp , defelp ;
double                  wx[3][3], wy[3][3],  I[3][3];
double                  um[3][3];
double                  dxx[3][3],dyx[3][3],dxy[3][3],dyy[3][3];
double                  dnf[3][3];
extern  int report(int itnum, double uchange);
extern  int localinit();
double					goalt;
int						outputts, outputtsnum;
double					resL2, resLinf;
int						resmaxi;
int						jacobianType;
int				preconditionerType;
int		put_BFKe(), UpDateMes(), StartMes(), isvarfixed(), get_bKeNumJ(), get_bKeAnalJ();
int		put_Ke(), reset_u(), get_PGKeAnalJ(), get_PGKeNumJ(), get_PGSe(), probNodes();
void	updateVelocities();
extern int	partID;
extern int	transportData;
extern int	uchangeTransport, maxChangeTransportGlobal;
struct transportSettings	transportSet;
void	findMaxchange();//added by mostafa Nov 01,09
void	setGWFtoZero();// added by mostafa Nov 01,09
double avgCelerity(); // added by mostafa Nov 02,09
struct  node   *maxnodeP;//added by mostafa Nov 01,09

int             init_gp()
{
        int     i,j ;

/*      Ashraf Aug 15, 1992     */

        for(i=0;i<3;i++) {
                for(j=0;j<3;j++) {
                        wx[i][j] = 0.0;
                        wy[i][j] = 0.0;
                        um[i][j] = 0.0;
                        dxx[i][j] = 0.0;
                        dxy[i][j] = 0.0;
                        dyx[i][j] = 0.0;
                        dyy[i][j] = 0.0;
                        dnf[i][j] = 0.0;

                }
        }               
        for(i=0;i<3;i++) {
                for(j=0;j<3;j++) {
                        if(i==j)
                                I[i][j] = 1.0;
                        else
                                I[i][j] = 0.0;
                }                               
        }               
                        
        gp.N = NULL ;
        gp.El = NULL ;
        gp.B = NULL ;
        gp.S = NULL ;
        gp.M = NULL ;
        for(i=0;i<4;i++) {
                eqnsets[i].Kp = NULL ;
                eqnsets[i].Fp = NULL ;
                eqnsets[i].diagp = NULL ;
                eqnsets[i].neqns = 0 ;
        }
        ep.K = NULL ;
        ep.F = NULL ;
        ep.S = NULL ;
        ep.M = NULL ;
        ep.J = NULL ;
        ep.DM = NULL ;
        fprime = NULL ;
        first_fnp = NULL ;
        Mesh.nx = 10 ;
        Mesh.ny = 1 ;
        Mesh.nbx = 1 ;
        Mesh.nby = 1 ;
        Mesh.eltype = 111 ;
        Mesh.maptype = 111 ;
        Mesh.firstmap = (struct blockmap *) calloc(24,sizeof(struct blockmap)) ;
        tvals.nsteps = 1 ;
        tvals.t = 654.321 ;
        tvals.dt = 1.0 ;        
        tvals.theta = 0.5 ;
        tvals.dtfac = 1.0 ;
        tvals.tfinal = 100.0;
        tvals.uc = 0.10;
        tvals.method = 0 ;
        tvals.iter = 1 ;
        tvals.gwH = 0.05 ;
        tvals.minH = 0.05;
        tvals.GWD = 0.0;
        tvals.T = 1.0;
		tvals.S = 1.0;
        sets.latitude = 0.0;
        sets.oldeqnew=1;
        sets.diffusivewave=0;
        sets.uwJ=0.0;
        sets.plotcode=2;
        sets.transbcs=0;
        sets.maxitnum=7;
        sets.smallH=1;
        sets.JE=1;
        N.sym = 1 ;
        N.nonlin = 1 ;
		N.trans = 1;
        NelemInK = 0;
        eqnsets[0].assembleFlag = 0;
		eqnsets[1].assembleFlag = 0;
		memJacobian = 0;
		jacobianType = 1;
		preconditionerType = 0;
        lumpS = 0;
        localinit() ;
		return (0);
}

int     MkLapMesh(theMesh)
struct  RegMesh *theMesh ;

{
		extern int MakeMesh();
        N.trans = theMesh->con.trans = 1 ;
        N.dims = theMesh->con.dims = 2 ;
        N.vars = theMesh->con.vars = 3 ;
        N.Keqns[0] = theMesh->con.Keqns[0] = 1 ;
        N.Keqns[1] = theMesh->con.Keqns[1] = 1 ;
        N.params = theMesh->con.params = 3 ;
        N.bparams = theMesh->con.bparams = 7 ;
        defndp = 0.00 ;
        defelp = 0.00 ;
        MakeMesh(theMesh) ;
        return(0) ;
}
   



double     update_PGKe(gbfp,gtfp,nbf,ntf,np,vcode,theElp,wt) 
int             nbf, ntf, vcode  ;
struct element  *theElp ;
struct node     *np[] ;
struct govfuncs *gbfp, *gtfp ;
double  wt;
{
        int     i, j, l, m, n ;
		int		update_gwKe();

        double  K00(),K01(),K02(),K10(),K11(),K12(),K20(),K21(),K22();
		double  J00(),J01(),J02(),J10(),J11(),J12(),J20(),J21(),J22();
        
        double  S(),FF(),dotm();
        double  UV, cx, cy,cs,us,vs;
        double  Cstar, Qxnew, Qynew, Hnew, Unew, Vnew, Sox, Soy, ff, root, nu, Htemp;
        double	Qxold, Qyold, Hold, Uold, Vold;
        double	dHdx, dHdy, dQXdx, dQXdy, dQYdx, dQYdy;
	    double	dnudH, dnudQX, dnudQY, SXi, SYi;
        double	dSXdHj, dSXdQXj, dSXdQYj, dSYdHj, dSYdQXj, dSYdQYj;
        double	dMXdHj, JMXdH, dMYdHj, JMYdH;
        double	dMXdQXj, JMXdQX, dMYdQXj, JMYdQX;
        double	dMXdQYj, JMXdQY, dMYdQYj, JMYdQY;
		double	dUdx, dUdy, dVdx, dVdy, Ph, dec;
		double  verticalflux;
		double	Dnew, Dold;


        double  F=0.0,c,c2,c22,U2,V2, e2;
//      double	Un[MNSF], Vn[MNSF];
		double	Uo[MNSF], Vo[MNSF];
        double	Unmax = 0.0, Unmin = 0.0;
        double	Vnmax = 0.0, Vnmin = 0.0;
		double	Ud, Vd, Hd;
        
        double  se1,se2,se3,Iu,IIu,IIIu,Ic,IIc,IIIc,div;
        double  Ax[3][3], Ay[3][3], A2[3][3], UI[3][3];
        
        double uw();
        double dotm();

        double  cs2, dCdH; 
	
		double kice;			//underice roughness	
		double tice;			//submerged ice thickness
		double fac;				//fac = 1.0 when open water, fac = 2.0 when ice covered
		double rdetJ, rSWE, rGWE;
		double UWTemp;
		double fCoriolis, dtNatural;
		int	   nTrans;
		double ff1;
	
		kice = tice = 0.0;
                        
		if(vcode<0)
                n = N.vars ;
        else
                n = 1 ;
        nTrans = n;

        Qxnew = Qynew = Hnew = Unew = Vnew = Sox = Soy = ff = Cstar = 0.0 ;
        Qxold = Qyold = Hold = Uold = Vold = 0.0;
		Ud = Vd = Hd = 0.0;
        dHdx = dHdy = dQXdx = dQXdy = dQYdx = dQYdy = 0.0;
	    dUdx = dUdy = dVdx = dVdy = 0.0;
		rdetJ = rSWE = rGWE = 0.0;

        for(i=0;i<nbf;i++) {
                Hnew += gbfp->f[i]*np[i]->u[0] ;
				tice += gbfp->f[i]*np[i]->ice[0] ;
        }
		tice = tvals.sg_ice * tice;

		if (Hnew <= 0.0){
			update_gwKe(gbfp,gtfp,nbf,ntf,np,vcode,theElp,wt);
			dtNatural = gbfp->detJ/gbfp->w*wt/4.0/tvals.T * tvals.S;
        	return (dtNatural);
		}
		
        for(i=0;i<nbf;i++) {
        	if((np[i]->uo[0]-(np[i]->ice[0]*tvals.sg_ice)) > tvals.gwH) {
        		Uo[i] = np[i]->uo[1] / (np[i]->uo[0]-(np[i]->ice[0]*tvals.sg_ice));
        		Vo[i] = np[i]->uo[2] / (np[i]->uo[0]-(np[i]->ice[0]*tvals.sg_ice));
        	}
        	else {
        		Uo[i] = 0.0;
        		Vo[i] = 0.0;
        	}
/*        	Unmax = (Un[i] > Unmax) ? Un[i] : Unmax;
        	Unmin = (Un[i] < Unmin) ? Un[i] : Unmin;
        	Vnmax = (Vn[i] > Vnmax) ? Vn[i] : Vnmax;
        	Vnmin = (Vn[i] < Vnmin) ? Vn[i] : Vnmin;
*/       }

        for(i=0;i<nbf;i++) {
                Qxnew += gbfp->f[i]*np[i]->u[1] ;
                Qynew += gbfp->f[i]*np[i]->u[2] ;
                Hold += gbfp->f[i]*np[i]->uo[0] ;
                Qxold += gbfp->f[i]*np[i]->uo[1] ;
                Qyold += gbfp->f[i]*np[i]->uo[2] ;
                ff += gbfp->f[i]*np[i]->p[1] ;
                Sox += - gbfp->dfdx[i]*np[i]->p[0] ;
                Soy += - gbfp->dfdy[i]*np[i]->p[0] ;
                dHdx += gbfp->dfdx[i]*np[i]->u[0];
                dHdy += gbfp->dfdy[i]*np[i]->u[0];
                dQXdx += gbfp->dfdx[i]*np[i]->u[1];
                dQXdy += gbfp->dfdy[i]*np[i]->u[1];
                dQYdx += gbfp->dfdx[i]*np[i]->u[2];
                dQYdy += gbfp->dfdy[i]*np[i]->u[2];
                dUdx += gbfp->dfdx[i]*Uo[i];
                dUdy += gbfp->dfdy[i]*Uo[i];
                dVdx += gbfp->dfdx[i]*Vo[i];
                dVdy += gbfp->dfdy[i]*Vo[i];
				kice += gbfp->f[i]*np[i]->ice[1];
				Hd += gbfp->f[i]*np[i]->ud[0];
				Ud += gbfp->f[i]*np[i]->ud[1];
				Vd += gbfp->f[i]*np[i]->ud[2];


        }

	    for(i=0;i<N.vars;i++) {
                for(j=0;j<N.vars;j++) {
                        wx[i][j] = 0.0;
                        wy[i][j] = 0.0;
                        dxx[i][j] = 0.0;
                        dxy[i][j] = 0.0;
                        dyx[i][j] = 0.0;
                        dyy[i][j] = 0.0;
                        um[i][j] = 0.0;
                        Ax[i][j] = 0.0;
                        Ay[i][j] = 0.0;
                        A2[i][j] = 0.0;
                        UI[i][j] = 0.0;
                        
                }
        }

		Dnew = Hnew - tice;
		if (Dnew < tvals.gwH){
			Dnew = tvals.gwH;
			tice = 0.0;
		}
		
		Dold = Hold - tice;
		if (Dold < tvals.gwH){
			Dold = tvals.gwH;
			tice = 0.0;
		}

		UWTemp = UW;
		if (Dnew == tvals.gwH)
			UW = 0.0;
		
        for(i=0;i<nbf;i++) {              
//                switch(sets.buw) {

//                        case 0: //no buw
                                gtfp->dfdx[i] *= (UW * sqrt(gbfp->detJ/gbfp->w*wt/4.0)) ;
                                gtfp->dfdy[i] *= (UW * sqrt(gbfp->detJ/gbfp->w*wt/4.0)) ;
//                                break;
/*
                        case 1: //no xbuw
                                gtfp->dfdx[i] = UW * gbfp->dfdx[i] * sqrt(gbfp->detJ/gbfp->w*wt/4.0) *(1-np[i]->fxc);
                                gtfp->dfdy[i] = UW * gbfp->dfdy[i] * sqrt(gbfp->detJ/gbfp->w*wt/4.0);
                                break;
                        case 2: //no ybuw
                                gtfp->dfdx[i] = UW * gbfp->dfdx[i] * sqrt(gbfp->detJ/gbfp->w*wt/4.0);
                                gtfp->dfdy[i] = UW * gbfp->dfdy[i] * sqrt(gbfp->detJ/gbfp->w*wt/4.0) *(1-np[i]->fxc);
                                break;
                        case 3: //with all buw
                                gtfp->dfdx[i] = UW * gbfp->dfdx[i] * sqrt(gbfp->detJ/gbfp->w*wt/4.0);
                                gtfp->dfdy[i] = UW * gbfp->dfdy[i] * sqrt(gbfp->detJ/gbfp->w*wt/4.0);
                                break;
                }
*/        }

		if((np[0]->ice[0]>0)&&(np[1]->ice[0]>0)&&(np[2]->ice[0]>0))
		{
			ff = pow((pow(kice,0.25)+pow(ff,0.25))/2,4);
			fac = 2.0;
		}
		else
		{
			tice = 0.0;
			kice = 0.0;
			fac = 1.0;
		}
               
        Unew = Qxnew/Dnew; 
        Vnew = Qynew/Dnew ;
        Uold = Qxold/Dold; 
        Vold = Qyold/Dold;

		if (Dnew == tvals.gwH){
			Unew = Vnew = 0.0;
		}

		if (Dold == tvals.gwH){
			Uold = Vold = 0.0;
		}

		if(jacobianType == 1)
		{
			Uold = Unew; 
			Vold = Vnew;
			Hold = Hnew;
			Dold = Dnew;
		}

	    UV = Unew*Vnew;
	    cx = GRAV*Dnew-Unew*Unew;
	    cy = GRAV*Dnew-Vnew*Vnew;
	    U2 = 2.*Unew;
	    V2 = 2.*Vnew;
	    
	    c22 = GRAV*Hnew/2.;

        Htemp = 12.0*Dnew/(fac*ff); 
		e2 = E*E;
        if(Htemp > e2) {
            Cstar = 2.5 * log(Htemp);
            dCdH = 2.5/Dnew;
        }
        else {
        	Cstar = 2.5 + 2.5/e2 * Htemp;
        	dCdH = 30.0/e2/ff;
        }
 
// Alternate friction formulation from Graeme Smart via Craig Addley
//  Assumes z0 = k/30
/*		Htemp = 11.0*Dnew/(fac*ff); 
		e2 = E/30.;
        if(Htemp > e2) {
            Cstar = 2.5 * (log(Htemp)+ 0.5/Htemp);
            dCdH = 2.5/Dnew*(1 - 0.5/Htemp);
        }
        else {
        	Cstar = 1.255 * Htemp;
        	dCdH = 13.8/ff;
        }
*/                
        root = sqrt( Qxnew*Qxnew + Qynew*Qynew );
        cs2 = Cstar*Cstar;               
        ff = fac*root/Dnew/Dnew/cs2 ;

		if(Dnew == tvals.gwH)
			ff = 100.;

		fCoriolis = 0.000146*sin(sets.latitude*0.0174533);
 
        dnf[0][0] = 0.0;
        dnf[0][1] = 0.0;
        dnf[0][2] = 0.0;
        dnf[1][0] = -2.0*Qxnew*fac*root/Dnew/Dnew/Dnew/cs2 * (1.0 + Dnew*dCdH / Cstar);
        dnf[2][0] = -2.0*Qynew*fac*root/Dnew/Dnew/Dnew/cs2 * (1.0 + Dnew*dCdH / Cstar);
        if(root == 0.0) {
                dnf[1][1] = 0.0;
                dnf[1][2] = 0.0;
                dnf[2][1] = 0.0;
                dnf[2][2] = 0.0;
        }
        else {    
                dnf[1][1] = Qxnew*Qxnew*fac/root/Dnew/Dnew/cs2;
                dnf[1][2] = Qxnew*Qynew*fac/root/Dnew/Dnew/cs2;                        
                dnf[2][1] = Qxnew*Qynew*fac/root/Dnew/Dnew/cs2;
                dnf[2][2] = Qynew*Qynew*fac/root/Dnew/Dnew/cs2;
        }
		
		if(UW > 0.0) {

		        c = sqrt(fabs(Dold*GRAV));
		        cs= c*c;        
		        us= Uold*Uold;  
		        vs= Vold*Vold;  
		        c2 = 2.*c;

		        Ax[0][0] = 0.0;
		        Ax[0][1] = 1.0;
		        Ax[0][2] = 0.0;

		        Ax[1][0] = cs-us;
		        Ax[1][1] = 2. *Uold;
		        Ax[1][2] = 0.0;

		        Ax[2][0] = - Uold * Vold;
		        Ax[2][1] = Vold;
		        Ax[2][2] = Uold;

		        
		        Ay[0][0] = 0.0;
		        Ay[0][1] = 0.0;
		        Ay[0][2] = 1.0;

		        Ay[1][0] = - Uold * Vold;
		        Ay[1][1] = Vold;
		        Ay[1][2] = Uold;

		        Ay[2][0] = cs-vs;
		        Ay[2][1] = 0.0;
		        Ay[2][2] = 2. * Vold;
		        

		        for(l=0;l<3;l++) {
		                for(m=0;m<3;m++)
		                        A2[l][m] = dotm(Ax,Ax,l,m,n) + dotm(Ay,Ay,l,m,n)
		                                                /* + dotm(nf,nf,l,m,n) shakib */;
		        }
		        
		        se1 = sqrt(us + vs + cs);
		        se2 = sqrt(3.*cs + 2.* us + 2. * vs + c*sqrt( cs + 16.*us + 16.*vs))/sqrt(2.);
		        se3 = sqrt(3.*cs + 2.* us + 2. * vs - c*sqrt( cs + 16.*us + 16.*vs))/sqrt(2.);

		        Iu = se1 + se2 + se3;
		        IIu = se1*se2 + se2*se3 + se3*se1;
		        IIIu = se1*se2*se3;

		        Ic = se1*se1 + se2*se2 + se3*se3;
		        IIc = se1*se1*se2*se2 + se2*se2*se3*se3 + se3*se3*se1*se1;
		        IIIc = IIIu*IIIu;

		        div = IIIu*IIIu*(IIIu+Iu*Ic)+Iu*Iu*(Iu*IIIc+IIIu*IIc);

		        for(l=0;l<3;l++) {
		                for(m=0;m<3;m++)

		                        UI[l][m] = (Iu*(Iu*IIu-IIIu)*dotm(A2,A2,l,m,n)-
		                                        (Iu*IIu-IIIu)*(IIIu+Iu*Ic)*A2[l][m]+
		                                        (IIu*IIIu*(IIIu+Iu*Ic)+Iu*Iu*(IIu*IIc+IIIc))*I[l][m])/div;
		        }
		        
		        for(l=0;l<3;l++) {
		                for(m=0;m<3;m++) {
		                                wx[l][m] = dotm(UI,Ax,l,m,n);
		                                wy[l][m] = dotm(UI,Ay,l,m,n);
		                }
		        }       

		        for(i=0;i<ntf;i++) {
		                for(l=0;l<n;l++) {
		                        for(m=0;m<n;m++)                        
		                                um[l][m] = uw(l,m,i,gtfp);
		                }
		        }
		}

		UW = UWTemp;

/* diffusion terms */

		Ph = sqrt(2.*dUdx*dUdx + (dUdy+dVdx)*(dUdy+dVdx) + 2.*dVdy*dVdy);
		dec = dHdx*Qxnew + dHdy*Qynew;

        if(Cstar>1.0)
                nu = epsilon1 + DIF*root/Cstar + epsilon3*epsilon3*Hold*Hold*Ph;
        else 
                nu = epsilon1 + DIF*root/1.0 + epsilon3*epsilon3*Hold*Hold*Ph;
		
		if(Dnew == tvals.gwH)
			nu = 0.0;

        dxx[1][0] = -Unew*(2.0*nu);
        dxx[1][1] = 2.0*nu;
        dxx[2][0] = - Vnew*nu;
        dxx[2][2] = 1.0*nu;

        dyx[2][0] = -Unew*nu;
        dyx[2][1] = nu;

        dxy[1][0] = -Vnew*nu;
        dxy[1][2] = nu;

        dyy[1][0] = -Unew*nu;
        dyy[1][1] = nu;
        dyy[2][0] = -Vnew*(2.0*nu);
        dyy[2][2] = 2.0*nu;

		if (jacobianType == 0)
		{
			dnudH = -nu/Cstar*dCdH;
			if(root > 0.0){
        		dnudQX = nu/root/root*Qxnew;
        		dnudQY = nu/root/root*Qynew;
			}
			else{
        		dnudQX = 0.0;
        		dnudQY = 0.0;
			}
        
    		JMXdH = ( (GRAV + 2*Unew*Unew/Dnew)*dHdx 
    				  - 2*Unew/Dnew*dQXdx + 2*Unew*Vnew/Dnew*dHdy
    				  - Vnew/Dnew*dQXdy - Unew/Dnew*dQYdy ) + dnf[1][0];
    		JMYdH = ( (GRAV + 2*Vnew*Vnew/Dnew)*dHdy 
    				  - 2*Vnew/Dnew*dQYdy + 2*Unew*Vnew/Dnew*dHdx
    				  - Unew/Dnew*dQYdx - Vnew/Dnew*dQXdx ) + dnf[2][0];
    		JMXdQX = -2*Unew/Dnew*dHdx + 2/Dnew*dQXdx - Vnew/Dnew*dHdy + dQYdy/Dnew + dnf[1][1] ;
    		JMYdQX = -Vnew/Dnew*dHdx + dQYdx/Dnew + dnf[2][1] ;
    		JMXdQY = -Unew/Dnew*dHdy + dQXdy/Dnew + dnf[1][2] ;
    		JMYdQY = -2*Vnew/Dnew*dHdy + 2/Dnew*dQYdy - Unew/Dnew*dHdx + dQXdx/Dnew + dnf[2][2] ;
		}
       
		//factor for dampening at D=tvals.gwH gausspoints
		if(tice > 0.0)
			ff1 = 100.0;
		else
			ff1 = 10.0;
                
        for(i=0;i<ntf;i++) {

			if (jacobianType == 0)
			{
				SXi = gbfp->dfdx[i]*2*(dQXdx-dHdx*Unew) 
        			+ gbfp->dfdy[i]*(-Unew*dHdy - Vnew*dHdx + dQXdy + dQYdx);
        		SYi = gbfp->dfdy[i]*2*(dQYdy-dHdy*Vnew) 
        			+ gbfp->dfdx[i]*(-Vnew*dHdx - Unew*dHdy + dQYdx + dQXdy);
			}
       			
                for(j=0;j<nbf;j++) {
                        KE(i*n,j*n) += K00(gtfp,gbfp,i,j,cx,cy,UV,Sox,Soy,ff,tice);
                        KE(i*n,j*n+1) += K01(gtfp,gbfp,i,j,U2,Vnew,Sox,Soy,ff,fCoriolis);
                        KE(i*n,j*n+2) += K02(gtfp,gbfp,i,j,Unew,V2,Sox,Soy,ff,fCoriolis);
                        KE(i*n+1,j*n) += K10(gtfp,gbfp,i,j,c22,cx,cy,UV,Sox,Soy,ff,tice);
                        KE(i*n+1,j*n+1) += K11(gtfp,gbfp,i,j,Unew,Vnew,U2,Sox,Soy,ff,fCoriolis);
                        KE(i*n+1,j*n+2) += K12(gtfp,gbfp,i,j,Unew,V2,Sox,Soy,ff,fCoriolis);
                        KE(i*n+2,j*n) += K20(gtfp,gbfp,i,j,c22,cx,cy,UV,Sox,Soy,ff,tice);
                        KE(i*n+2,j*n+1) += K21(gtfp,gbfp,i,j,U2,Vnew,Sox,Soy,ff,fCoriolis);
                        KE(i*n+2,j*n+2) += K22(gtfp,gbfp,i,j,Unew,Vnew,V2,Sox,Soy,ff,fCoriolis);

						if((Dnew == tvals.gwH)){
							KE(i*n+1,j*n+1) += ff1 * gbfp->f[i] * gbfp->f[j] * gbfp->detJ;
							KE(i*n+2,j*n+2) += ff1 * gbfp->f[i] * gbfp->f[j] * gbfp->detJ;
						}

						if (jacobianType == 0)
						{

							if(N.nonlin > 0) {
                        		dMXdHj = JMXdH *gbfp->f[j];
                        		dMYdHj = JMYdH *gbfp->f[j];
                        		dMXdQXj = JMXdQX *gbfp->f[j];
                        		dMYdQXj = JMYdQX *gbfp->f[j];
                        		dMXdQYj = JMXdQY *gbfp->f[j];
                        		dMYdQYj = JMYdQY *gbfp->f[j];
                        		
                        		dSXdHj = (SXi*dnudH + ((2*gbfp->dfdx[i]*dHdx + gbfp->dfdy[i]*dHdy)*Unew/Hnew 
                        								+ gbfp->dfdy[i]*dHdx*Vnew/Hnew)*nu) *gbfp->f[j];
                        		dSYdHj = (SYi*dnudH + ((2*gbfp->dfdy[i]*dHdy + gbfp->dfdx[i]*dHdx)*Vnew/Hnew 
                        								+ gbfp->dfdx[i]*dHdy*Unew/Hnew)*nu) *gbfp->f[j];
                        		dSXdQXj = (SXi*dnudQX -(2*gbfp->dfdx[i]*dHdx + gbfp->dfdy[i]*dHdy)/Hnew*nu) *gbfp->f[j];
                        		dSYdQXj = (SYi*dnudQX -(gbfp->dfdx[i]*dHdy)/Hnew*nu) *gbfp->f[j];
                        		dSXdQYj = (SXi*dnudQY -(gbfp->dfdy[i]*dHdx)/Hnew*nu) *gbfp->f[j];
                        		dSYdQYj = (SYi*dnudQY -(2*gbfp->dfdy[i]*dHdy + gbfp->dfdx[i]*dHdx)/Hnew*nu) *gbfp->f[j];
                        		
								JE(i*n,j*n) += J00(gtfp, gbfp, i, j, dMXdHj,dMYdHj) ;
								JE(i*n,j*n+1) += J01(gtfp, gbfp, i, j,Hnew,Unew,Vnew, ff,dMXdQXj,dMYdQXj) ;
								JE(i*n,j*n+2) += J02(gtfp, gbfp, i, j, Hnew, Unew,Vnew, ff,dMXdQYj,dMYdQYj) ;
								JE(i*n+1,j*n) += J10(gtfp, gbfp, i, j,Hnew, Unew,Vnew, ff,nu,dMXdHj,dMYdHj,dSXdHj) ;
								JE(i*n+1,j*n+1) += J11(gtfp, gbfp, i, j,Hnew, Unew,Vnew, ff,nu,dMXdQXj,dMYdQXj,dSXdQXj) ;
								JE(i*n+1,j*n+2) += J12(gtfp, gbfp, i, j,Hnew, Unew,Vnew, ff,nu,dMXdQYj,dMYdQYj,dSXdQYj) ;
								JE(i*n+2,j*n) += J20(gtfp, gbfp, i, j,Hnew, Unew,Vnew, ff,nu,dMXdHj,dMYdHj,dSYdHj) ;
								JE(i*n+2,j*n+1) += J21(gtfp, gbfp, i, j,Hnew, Unew,Vnew, ff,nu,dMXdQXj,dMYdQXj,dSYdQXj) ;
								JE(i*n+2,j*n+2) += J22(gtfp, gbfp, i, j,Hnew, Unew,Vnew, ff,nu,dMXdQYj,dMYdQYj,dSYdQYj) ;
							}
						}
		        }
                FE(i*n) += tvals.T * (gbfp->dfdx[i]*Sox + gbfp->dfdy[i]*Soy) * gbfp->detJ -(uw(0,1,i,gtfp)*GRAV*tice*Sox + uw(0,2,i,gtfp)*GRAV*tice*Soy)* gbfp->detJ;
				verticalflux = tvals.T * (-gbfp->dfdx[i]*dHdx - gbfp->dfdy[i]*dHdy + gbfp->dfdx[i]*Sox + gbfp->dfdy[i]*Soy) * gbfp->detJ;
				if(verticalflux > 0.0) verticalflux = 0.0;
				FE(i*n+1) += Unew*verticalflux -(gbfp->f[i]*GRAV*tice*Sox + uw(1,1,i,gtfp)*GRAV*tice*Sox + uw(1,2,i,gtfp)*GRAV*tice*Soy)* gbfp->detJ;
				FE(i*n+2) += Vnew*verticalflux -(gbfp->f[i]*GRAV*tice*Soy + uw(2,1,i,gtfp)*GRAV*tice*Sox + uw(2,2,i,gtfp)*GRAV*tice*Soy)* gbfp->detJ;
        }

        for(i=0;i<ntf;i++) {
                for(j=0;j<nbf;j++) {
                        for(l=0;l<nTrans;l++){       
                                for(m=0;m<nTrans;m++){
                                        SE(i*n+l,j*n+m) += S(gtfp,gbfp,i,j,l,m) ;
                                }
                        }
                }
        }
		dtNatural = sqrt(gbfp->detJ/gbfp->w*wt/4.0)/(fabs(Vnew)+sqrt(fabs(GRAV*Hnew)));
        return(dtNatural) ;
}

int	update_gwKe(gbfp,gtfp,nbf,ntf,np,vcode,theElp,wt) 
int             nbf, ntf, vcode  ;
struct element  *theElp ;
struct node     *np[] ;
struct govfuncs *gbfp, *gtfp ;
double  wt;
{
        int     i, j, n ;
        double  ff1, ff2, Sox = 0.0, Soy = 0.0, tice = 0.0;
        
        if(vcode<0)
                n = N.vars ;
        else
                n = 1 ;
        
        for(i=0;i<nbf;i++) {
                Sox += - gbfp->dfdx[i]*np[i]->p[0] ;
                Soy += - gbfp->dfdy[i]*np[i]->p[0] ;
				tice += gbfp->f[i]*np[i]->ice[0] ;
        }

		//factor for dampening
		if(tice > 0.0)
			ff1 = 100.0;
		else
			ff1 = 10.0;

		ff2 = 1.0e12;
                
        for(i=0;i<ntf;i++) {
                for(j=0;j<nbf;j++) {
			        KE(i*n,j*n) += tvals.T * (gbfp->dfdx[i]*gbfp->dfdx[j] + gbfp->dfdy[i]*gbfp->dfdy[j]) * gbfp->detJ;
                	SE(i*n,j*n) += tvals.S * (gbfp->f[i] * gbfp->f[j]) * gbfp->detJ;
					KE(i*n+1,j*n+1) += ff1 * gbfp->f[i] * gbfp->f[j] * sqrt(gbfp->detJ /6.0);
                    KE(i*n+2,j*n+2) += ff1 * gbfp->f[i] * gbfp->f[j] * sqrt(gbfp->detJ /6.0);
					if((i == j) && ((np[i]->u[0]) <= 0.0)){
                 		KE(i*n+1,j*n+1) += ff2 ;
						KE(i*n+2,j*n+2) += ff2 ;
					}
				}
                                     
                FE(i*n) += tvals.T * (gbfp->dfdx[i]*Sox + gbfp->dfdy[i]*Soy) * gbfp->detJ;
        }

        return(0) ;
}


double          S(v,f,i,j,l,m)                                              
int                     i, j, l, m;
struct govfuncs         *v, *f ;

{
        double  t ;

		if(lumpS == 0) {
        	if(l==m)
                t = (f->f[i]*f->f[j]+ (wx[l][m]*v->dfdx[i]*f->f[j]+wy[l][m]*v->dfdy[i]*f->f[j])) * f->detJ;
        	else
                t =  (wx[l][m]*v->dfdx[i]*f->f[j]+wy[l][m]*v->dfdy[i]*f->f[j])*f->detJ ;
		}
		else {
//			if( (l == 0) && (m == 0))
	       	if(l==m)
					t = (f->f[i]*f->f[j]) * f->detJ; 
        	else
        			t = 0.0;
        }

        return(t);
}


double          K00(v,f,i,j,cx,cy,UV,Sox,Soy,ff,ts)
int                     i,j;
struct govfuncs         *v, *f;
double                  cx,cy,UV,Sox,Soy,ff,ts;

{
        register        double  t3,t4,t5,t6,t;
        double          uw(),dotm();

        t3 = (cx*uw(0,1,i,v)-UV*uw(0,2,i,v))*f->dfdx[j];
        t4 = (cy*uw(0,2,i,v)-UV*uw(0,1,i,v))*f->dfdy[j];
        t5 = (-uw(0,1,i,v)*GRAV*Sox-uw(0,2,i,v)*GRAV*Soy);
        t6 = tvals.T * (f->dfdx[i]*f->dfdx[j] + f->dfdy[i]*f->dfdy[j]) ;
        
        t = (t3+t4+t5*f->f[j]+t6)*f->detJ;
        return(t);
}
        

double          K01(v,f,i,j,U2,V,Sox,Soy,ff,fCoriolis)
int                     i,j;
struct govfuncs         *v, *f;
double                  U2,V,Sox,Soy,ff,fCoriolis;

{
        register        double  t1,t3,t4,t6,t7,t;
        double          uw(),dotm();

        t1 = -f->dfdx[i];
        t3 = (uw(0,0,i,v)+U2*uw(0,1,i,v)+V*uw(0,2,i,v))*f->dfdx[j];
        t4 = V*uw(0,1,i,v)*f->dfdy[j];
        t6 = /*dotm(um,nf,0,1,3)*/ uw(0,1,i,v)*ff;
		t7 = uw(0,2,i,v)*fCoriolis;

        t = ((t1+t6+t7)*f->f[j]+t3+t4)*f->detJ;
        return(t);
}


double          K02(v,f,i,j,U,V2,Sox,Soy,ff,fCoriolis)
int                     i,j;
struct govfuncs         *v, *f;
double                  U,V2,Sox,Soy,ff,fCoriolis;

{
        register        double  t2,t3,t4,t6,t7,t;
        double          uw(),dotm();

        t2 = -f->dfdy[i];

        t3 = U*uw(0,2,i,v)*f->dfdx[j];
        t4 = (uw(0,0,i,v)+U*uw(0,1,i,v)+V2*uw(0,2,i,v))*f->dfdy[j];
        t6 = uw(0,2,i,v)*ff;
		t7 = -fCoriolis*uw(0,1,i,v);

        t = ((t2+t6+t7)*f->f[j]+t3+t4)*f->detJ;
        return(t);
}


double          K10(v,f,i,j,c22,cx,cy,UV,Sox,Soy,F,ts)
int                     i,j;
struct govfuncs         *v, *f;
double                  c22,cx,cy,UV,Sox,Soy,F,ts;

{
        register        double  t1,t2,t3,t4,t5,t6,t7,t8,t;
        double          uw(),dotm();

        t1 = -c22*f->dfdx[i];
		t2 = -GRAV*ts*f->f[i]*f->dfdx[j];
        t3 = (cx*uw(1,1,i,v)-UV*uw(1,2,i,v))*f->dfdx[j];
        t4 = (cy*uw(1,2,i,v)-UV*uw(1,1,i,v))*f->dfdy[j];
        t5 = -GRAV*Sox*f->f[i];
        t6 = (-uw(1,1,i,v)*GRAV*Sox-uw(1,2,i,v)*GRAV*Soy);
        t7 = dxx[1][0]*f->dfdx[i]*f->dfdx[j]+dxy[1][0]*f->dfdy[i]*f->dfdx[j];
        t8 = dyy[1][0]*f->dfdy[i]*f->dfdy[j];

        t = ((t1+t5+t6)*f->f[j]+t2+t3+t4+t7+t8)*f->detJ;
        return(t);
}
        

double          K11(v,f,i,j,U,V,U2,Sox,Soy,ff,fCoriolis)
int                     i,j;
struct govfuncs         *v, *f;
double                  U,U2,V,Sox,Soy,ff,fCoriolis;

{
        register        double  t1,t2,t3,t4,t5,t6,t8,t9,t10,t;
        double          uw(),dotm();

        t1 = -U*f->dfdx[i];
        t2 = -V*f->dfdy[i];

        t3 = (uw(1,0,i,v)+U2*uw(1,1,i,v)+V*uw(1,2,i,v))*f->dfdx[j];
        t4 = V*uw(1,1,i,v)*f->dfdy[j];
        t5 = ff*f->f[i];
        t6 = uw(1,1,i,v)*ff;
        t8 = dxx[1][1]*f->dfdx[i]*f->dfdx[j];
        t9 = dyy[1][1]*f->dfdy[i]*f->dfdy[j];
		t10 = fCoriolis*uw(1,2,i,v);

        t = ((t1+t2+t5+t6+t10)*f->f[j]+t3+t4+t8+t9)*f->detJ;
        return(t);
}


double          K12(v,f,i,j,U,V2,Sox,Soy,ff,fCoriolis)
int                     i,j;
struct govfuncs         *v, *f;
double                  U,V2,Sox,Soy,ff,fCoriolis;

{
        register        double  t3,t4,t6,t8,t9,t10,t;
        double          uw();

        t3 = U*uw(1,2,i,v)*f->dfdx[j];
        t4 = (uw(1,0,i,v)+U*uw(1,1,i,v)+V2*uw(1,2,i,v))*f->dfdy[j];
        t6 = uw(1,2,i,v)*ff;
        t8 = dxy[1][2]*f->dfdy[i]*f->dfdx[j];
		t9 = -fCoriolis*f->f[i]*f->f[j];
		t10 = -fCoriolis*uw(1,1,i,v);

        t = ((t6+t10)*f->f[j]+t3+t4+t8+t9)*f->detJ;
        return(t);
}



double          K20(v,f,i,j,c22,cx,cy,UV,Sox,Soy,ff,ts)
int                     i,j;
struct govfuncs         *v, *f;
double                  c22,cx,cy,UV,Sox,Soy,ff,ts;

{
        register        double  t1,t2,t3,t4,t5,t6,t7,t8,t;
        double          uw();

        t1 = -c22*f->dfdy[i];
		t2 = -GRAV*ts*f->f[i]*f->dfdy[j];
        t3 = (cx*uw(2,1,i,v)-UV*uw(2,2,i,v))*f->dfdx[j];
        t4 = (cy*uw(2,2,i,v)-UV*uw(2,1,i,v))*f->dfdy[j];
        t5 = -GRAV*Soy*f->f[i];
        t6 = (-uw(2,1,i,v)*GRAV*Sox-uw(2,2,i,v)*GRAV*Soy);
        t7 = dxx[2][0]*f->dfdx[i]*f->dfdx[j]+dyx[2][0]*f->dfdx[i]*f->dfdy[j];
        t8 = dyy[2][0]*f->dfdy[i]*f->dfdy[j];

        t = ((t1+t5+t6)*f->f[j]+t2+t3+t4+t7+t8)*f->detJ;
        return(t);
}
        

double          K21(v,f,i,j,U2,V,Sox,Soy,ff,fCoriolis)
int                     i,j;
struct govfuncs         *v, *f;
double                  U2,V,Sox,Soy,ff,fCoriolis;

{
        register        double  t3,t4,t6,t8,t9,t10,t;
        double          uw();


        t3 = (uw(2,0,i,v)+U2*uw(2,1,i,v)+V*uw(2,2,i,v))*f->dfdx[j];
        t4 = V*uw(2,1,i,v)*f->dfdy[j];
        t6 = uw(2,1,i,v)*ff;
        t8 = dyx[2][1]*f->dfdx[i]*f->dfdy[j];
		t9 = fCoriolis*f->f[i]*f->f[j];
		t10 = fCoriolis*uw(2,2,i,v);

        t = ((t6+t10)*f->f[j]+t3+t4+t8+t9)*f->detJ;
        return(t);
}


double          K22(v,f,i,j,U,V,V2,Sox,Soy,ff,fCoriolis)
int                     i,j;
struct govfuncs         *v, *f;
double                  U,V,V2,Sox,Soy,ff,fCoriolis;

{
        register        double  t1,t2,t3,t4,t5,t6,t8,t9,t10,t;
        double          uw(),dotm();


        t1 = -U*f->dfdx[i];
        t2 = -V*f->dfdy[i];
        t3 = U*uw(2,2,i,v)*f->dfdx[j];
        t4 = (uw(2,0,i,v)+U*uw(2,1,i,v)+V2*uw(2,2,i,v))*f->dfdy[j];
        t5 = ff*f->f[i];
        t6 = uw(2,2,i,v)*ff;
        t8 = dxx[2][2]*f->dfdx[i]*f->dfdx[j];
        t9 = dyy[2][2]*f->dfdy[i]*f->dfdy[j];
		t10 = -fCoriolis*uw(2,1,i,v);

        t = ((t1+t2+t5+t6+t10)*f->f[j]+t3+t4+t8+t9)*f->detJ;
        return(t);
}



double          uw(l,m,i,v)
struct  govfuncs        *v;
int             l,m,i;
{
        double  z;
/*      printf("i=%d\tl=%d\tm=%d\t",i,l,m);*/

        z = wx[l][m]*v->dfdx[i]+wy[l][m]*v->dfdy[i];/*printf("%g\n",z);*/
        return(z);

}


double           J00(v,f,i,j,dMXdHj,dMYdHj)
int                     i, j ;
struct govfuncs         *v, *f ;
double                  dMXdHj,dMYdHj;

{
        register        double  t3,t4,t;
        double          uw();
        
        t3 = uw(0,1,i,v)*dMXdHj ;
        t4 = uw(0,2,i,v)*dMYdHj;
        
        t = (t3+t4)* f->detJ  ;
        return(t) ;
}


double          J01(v,f,i,j,H,U,V,ff,dMXdQX,dMYdQX)               
int                     i, j ;
struct govfuncs         *v, *f ;
double                  H, U,V , ff,dMXdQX,dMYdQX;
{
        register        double  t3, t4,t;
        double          uw();

        t3 = uw(0,1,i,v)*dMXdQX ;
        t4 = uw(0,2,i,v)*dMYdQX;

        
        t = (t3 + t4)* f->detJ  ;
        return(t) ;
}

double          J02(v,f,i,j,H,U,V,ff,dMXdQY,dMYdQY)               
int                     i, j ;
struct govfuncs         *v, *f ;
double                  H, U,V , ff,dMXdQY,dMYdQY ;
{
        register        double  t3, t4,t;
        double          uw();

        t3 = uw(0,1,i,v)*dMXdQY ;
        t4 = uw(0,2,i,v)*dMYdQY;

        t = (t3 + t4) * f->detJ  ;

        return(t) ;
}

double          J10(v,f,i,j,H,U,V,ff,nu,dMXdHj,dMYdHj,dSXdHj)            
int                     i, j ;
struct govfuncs         *v, *f ;
double                  H, U,V , ff,nu, dMXdHj,dMYdHj,dSXdHj;
{
        register        double  t1,t2,t3,t4,t5,d8,t;
        double          uw(),dotm();

        t1 = -(GRAV * H / 2.0 - U * U)* f->dfdx[i] ;
        t2 =  U*V*f->dfdy[i];
        
        t3 = uw(1,1,i,v)*dMXdHj;
        t4 = uw(1,2,i,v)*dMYdHj;
        
        t5 = dnf[1][0] *f->f[i];
        
        d8 = dSXdHj;
        
        t = ((t1+t2+t5) * f->f[j] +(t3+t4)+ d8) * f->detJ  ;
        return(t) ;
}

double          J11(v,f,i,j,H,U,V,ff,nu,dMXdQX,dMYdQX,dSXdQXj)            
int                     i, j ;
struct govfuncs         *v, *f ;
double                  H, U ,V, ff,nu,dMXdQX,dMYdQX,dSXdQXj ;
{
        register        double  t1,t3,t4,t5,d8,t;
        double          uw();

        t1 = -U* f->dfdx[i];
        t3 = uw(1,1,i,v)*dMXdQX ;
        t4 = uw(1,2,i,v)*dMYdQX;

        t5 = dnf[1][1] *f->f[i];
        
        d8 = dSXdQXj;
        
        t = ((t1+t5) * f->f[j]+(t3+t4)+d8) * f->detJ  ;
        return(t) ;
}

double          J12(v,f,i,j,H,U,V,ff,nu,dMXdQY,dMYdQY,dSXdQYj)            
int                     i, j ;
struct govfuncs         *v, *f ;
double                  H, U ,V, ff,nu,dMXdQY,dMYdQY,dSXdQYj;
{
        register        double  t2,t3,t4,t5,d8,t;
        double          uw();

        t2 = -U*f->dfdy[i];
        
        t3 = uw(1,1,i,v)*dMXdQY ;
        t4 = uw(1,2,i,v)*dMYdQY;
        t5 = dnf[1][2] *f->f[i];
        
        d8 = dSXdQYj;
        
        t = ((t2+t5) * f->f[j]+(t3+t4)+d8) * f->detJ  ;
        return(t) ;
}

double          J20(v,f,i,j,H,U,V,ff,nu,dMXdHj,dMYdHj,dSYdHj)            
int                     i, j ;
struct govfuncs         *v, *f ;
double                  H, U ,V, ff,nu, dMXdHj,dMYdHj,dSYdHj ;
{
        register        double  t1,t2,t3,t4,t5,d8,t;
        double          uw();

        t1 = U*V* f->dfdx[i] ;
        t2 = -(GRAV * H / 2.0 - V*V)*f->dfdy[i];
        
        t3 = uw(2,1,i,v)*dMXdHj;
        t4 = uw(2,2,i,v)*dMYdHj;
        
        t5 = dnf[2][0] *f->f[i];
        
        d8 = dSYdHj;
        
        t = ((t1+t2+t5) * f->f[j] +(t3+t4)+d8)* f->detJ  ;
        return(t) ;
}

double          J21(v,f,i,j,H,U,V,ff,nu,dMXdQX,dMYdQX,dSYdQXj)            
int                     i, j ;
struct govfuncs         *v, *f ;
double                  H, U,V , ff ,nu,dMXdQX,dMYdQX,dSYdQXj;
{
        register        double  t1,t3,t4,t5,d8,t;
        double          uw();

        t1 = -V* f->dfdx[i];
        t3 = uw(2,1,i,v)*dMXdQX ;
        t4 = uw(2,2,i,v)*dMYdQX;

        t5 = dnf[2][1] *f->f[i];
        
        d8 = dSYdQXj;
        
        t = ((t1+t5) * f->f[j] +(t3+t4) + d8) * f->detJ  ;
        return(t) ;
}

double          J22(v,f,i,j,H,U,V,ff,nu,dMXdQY,dMYdQY,dSYdQYj)            
int                     i, j ;
struct govfuncs         *v, *f ;
double                  H, U ,V, ff,nu,dMXdQY,dMYdQY,dSYdQYj ;
{
        register        double  t2,t3,t4,t5,d8,t;
        double          uw();

        t2 = -V*f->dfdy[i];
        
        t3 = uw(2,1,i,v)*dMXdQY ;
        t4 = uw(2,2,i,v)*dMYdQY;

        t5 = dnf[2][2] *f->f[i];
        
        d8 = dSYdQYj;
        
        t = ((t2+t5) * f->f[j] +(t3+t4) + d8) * f->detJ  ;
        return(t) ;
}

double          FF(v,f,i,m)                                     
int                     i ;
struct govfuncs         *v, *f ;

{
        double  t ;
        t = 0.;
        return(t) ;
}

int     average_stage()

{
        struct  boundaryseg     *bseg;
        struct  belement                *belp,*bstart,*bend;
        double  twse,awse;
        int             i,j,ndata;
        
        bseg=gp.BSEG;
        belp=gp.B;
        
        for(i=0;i<N.boundarysegs;i++) {
                if((bseg->bcs[0]==1)||(bseg->bcs[0]==2)) {
                        twse=0.0;
                        while(belp->nps[0]->n!=bseg->nstart)
                                belp=belp->nextbelp;
                        bstart=belp;
                        ndata=0;
                        while(belp->nps[0]->n!=bseg->nend) {
                                for(j=0;j<belp->nnds;j++) {
                                        twse+=(belp->nps[j]->p[0]+belp->nps[j]->u[0]);
                                        ndata++;
                                }
                                belp=belp->nextbelp;
                        }
                        bend=belp;
                        awse=twse/ndata;
printf("awse=%g\n",awse);
//now calculate u[0] for all nodes on this bseg
                        belp=bstart;
                        while(belp!=bend) {
                                for(j=0;j<belp->nnds;j++)
                                        belp->nps[j]->u[0]=awse-belp->nps[j]->p[0];
                                belp=belp->nextbelp;
                        }
                }
                bseg=bseg->nextbseg;
        }
        return(0);
}       

double	getQin()
{
	struct boundaryseg	*bseg;

	bseg = gp.BSEG;
	while(bseg->bcs[0] != 1){
		bseg = bseg->nextbseg;
		if(bseg == NULL)
			return(0.0);
	}
	return(bseg->p[1]);
}

int	setQin(newQ)
double newQ;
{
	struct boundaryseg	*bseg;

	bseg = gp.BSEG;
	while(bseg->bcs[0] != 1){
		bseg = bseg->nextbseg;
		if(bseg == NULL)
			return(0);
	}
	bseg->p[1] = newQ;
	return(0);
}

struct boundaryseg* getObseg()
{
	struct boundaryseg	*bseg;

	bseg = gp.BSEG;
	while((bseg->bcs[0] != 3) && (bseg->bcs[0] != 5) ){
		bseg = bseg->nextbseg;
		if(bseg == NULL)
			return(NULL);
	}
	return(bseg);
}


int     update_BCs()

{

        struct  boundaryseg     *bseg;
        struct  belement                *belp,*bstart,*bend,*prevBelp;
        struct  gausspts                g[MGPTS];
        struct  shapefuncs      *fv,tfv; /*check*/
        int             i,j,ig,nbf,ni,dryFlag;
		int		nsf(), get_gspts(), get_shape();
		double  tice, fac;
		double wsesum, wse, awse;
        
        bseg    = gp.BSEG;
        belp    = gp.B;
        for(i=0;i<N.belms;i++) {
                for(j=0;j<NVARS;j++)
                        belp->bcs[j]=0;
                for(j=0;j<NBPARAMS;j++)
                        belp->p[j]=0.0;
                belp=belp->nextbelp;
        }
        belp    = gp.B;
        
        for(i=0;i<N.boundarysegs;i++) {
                bseg->p[2] = 0.0;
                bseg->p[5] = 0.0; 
				// sigma di5/3 * wi 
                while(belp->nps[0]->n!=bseg->nstart)
                        belp=belp->nextbelp;
                bstart = belp;
				wsesum = 0.0;

//	Loop through bseg to find average water surface elevation
				while(1){
                    belp->p[5]=sqrt((belp->nps[1]->x[0]-belp->nps[0]->x[0])*
                                    (belp->nps[1]->x[0]-belp->nps[0]->x[0])
                                   +(belp->nps[1]->x[1]-belp->nps[0]->x[1])*
                                    (belp->nps[1]->x[1]-belp->nps[0]->x[1]));
					bseg->p[5] += belp->p[5];
					wse =  (  belp->nps[0]->p[0]+belp->nps[0]->u[0]
							+ belp->nps[0]->p[0]+belp->nps[0]->u[0]) * 0.5;
					wsesum += wse*belp->p[5];
					if(belp->nps[1]->n!=bseg->nend)
						belp=belp->nextbelp;
					else
						break;
				}
				awse = wsesum/bseg->p[5];
				belp = bstart;

//	Loop through bseg again , estimating relative capacity of each belm based on awse
                while(1) {
						dryFlag = 0;
                        belp->bcs[0] = bseg->bcs[0];
						tice = 0.0;
                        belp->p[3]=0.0; // Hnew 
                        belp->p[4]=0.0; // BL      
                        belp->bseg=bseg;

                        nbf=nsf(belp->vtype);
                        ni=get_gspts(belp->vtype,g);
                        for(ig=0;ig<ni;ig++) {
                                get_shape(belp->vtype,&g[ig],&tfv);
                                fv=&tfv;
                                for(j=0;j<nbf;j++) {
										if((belp->nps[j]->u[0]-(belp->nps[j]->ice[0]*tvals.sg_ice))<tvals.gwH)
											dryFlag = 1;
										tice += fv->f[j]*tvals.sg_ice*belp->nps[j]->ice[0];
                                        belp->p[3] += fv->f[j]*belp->nps[j]->u[0];
                                        belp->p[4] += fv->f[j]*belp->nps[j]->p[0];
                                }
                        }
                        belp->p[3]/=ni;
                        belp->p[4]/=ni;
						belp->p[3] = awse - belp->p[4];
						if ((belp->nps[0]->ice[0]>0)&&(belp->nps[1]->ice[0]>0)) {
							belp->p[3] = belp->p[3]-tice;
							fac = 1.0/1.5874;
						}
						else
							fac = 1.0;
                        belp->p[0]=bseg->p[0]; 
                        belp->p[1]=bseg->p[1]; // spec H 
							if( ((belp->p[3]) <= tvals.gwH) || dryFlag == 1) {
                                belp->p[3] = 0.0; // the belp is dry -> no flow bound 
	                        }      
							if(belp->nps[1]->n!=bseg->nend)
							belp=belp->nextbelp;
						else
							break;
                }
                bend=belp;

//   Loop through the bseg again, reducing capacity of belms which are adjacent to zero flow belms.
				prevBelp = NULL; 
				belp = bstart;
                if((bseg->bcs[0]==1)||(bseg->bcs[0]==2)) {
					while(1){
						if(prevBelp != NULL){
							if (prevBelp->p[3] < 0.001)
								belp->p[3] /= 3.0;
						}
						if(belp != bend){
							if(belp->nextbelp->p[3] < 0.001)
								belp->p[3] /= 3.0;
						}
						if ((belp->nps[0]->ice[0]>0)&&(belp->nps[1]->ice[0]>0))
							fac = 1.0/1.5874;
						else
							fac = 1.0;
						bseg->p[2] += fac*pow(belp->p[3],(5.0/3.0)) * belp->p[5];
						if(belp->nps[1]->n!=bseg->nend){
							prevBelp = belp;
							belp=belp->nextbelp;
						}
						else
							break;
					}
				}
// Finally, calculate the q's for all belps of this bseg

                if((bseg->bcs[0]==1)||(bseg->bcs[0]==2)) {
					belp=bstart;
					while(1){
						if((belp->bcs[0]==1)||(belp->bcs[0]==2)){
							if ((belp->nps[0]->ice[0]>0)&&(belp->nps[1]->ice[0]>0))
								fac = 1.0/1.5874;
							else
								fac = 1.0;
							if(bseg->p[2] > 0.0)
								belp->p[1]= fac*pow(belp->p[3],(5.0/3.0))*bseg->p[1]/bseg->p[2];
							else
								belp->p[1] = bseg->p[1] / bseg->p[5];
							belp->p[2]=0.0;
						}
						if(belp!=bend)
							belp=belp->nextbelp;
						else
							break;
					}
				}
				bseg=bseg->nextbseg;
		}
        

        return(0);
}

int     update_BKe(belp,gbf,gtf,nbf,ntf,np,nvar,dx,dy) 
int                     nbf, ntf ;
double                  dx,dy;
struct node             *np[] ;
struct belement         *belp ;
struct shapefuncs       *gbf, *gtf ;

{
    int     i, j, n , l, m, bcCode;
    double  Qx, Qy, H, U, V, Qn, Qt, dS, Ha, Qxa, Qya, qn, zb, K;
    double  AK[3][3],cbcc[3], AJK[3][3],A[3][3],AF[3],bcc[3],phi[3],phin[3];
    double	Un[MRNSF], Vn[MRNSF];
    double  Qxnew, Qynew, Hnew, Unew, Vnew, tice;
	double	D0, D1, qf, dQn, qn0, qn1;
    double  ex;

    if(nvar<0)
            n = N.vars ;
    else
            n = 1 ;

/*  Compute Parameters Here      */     
/* Note: Ha, Qxa and Qya are variables which are used to calculate entries in the matrix A only */

    Qxnew = Qynew = Hnew = Unew = Vnew = Ha = Qxa = Qya = zb = tice = 0.0 ;
    dS = sqrt( dx*dx + dy*dy );

	for(i=0;i<nbf;i++)
		tice += gbf->f[i]*np[i]->ice[0] ;

	if((np[0]->ice[0]>0)&&(np[1]->ice[0]>0))
		tice = tvals.sg_ice * tice;	//submerged depth
	else
		tice = 0.0;


 	for(i=0;i<nbf;i++) {
		if((np[i]->u[0]-(np[i]->ice[0]*tvals.sg_ice)) > tvals.gwH) {
			Un[i] = np[i]->u[1] / (np[i]->u[0]-tice);
			Vn[i] = np[i]->u[2] / (np[i]->u[0]-tice);
		}
		else {
			Un[i] = 0.0;
			Vn[i] = 0.0;
		}
	}
	
 for(i=0;i<nbf;i++) {
        zb += gbf->f[i]*np[i]->p[0];
        Hnew += gbf->f[i]*np[i]->u[0] ;
        Qxnew += gbf->f[i]*np[i]->u[1] ;
        Qynew += gbf->f[i]*np[i]->u[2] ;
        Unew += gbf->f[i]*Un[i];
        Vnew += gbf->f[i]*Vn[i];
    }
    qn = -Qxnew*dy/dS + Qynew*dx/dS;
	qn0 = -np[0]->u[1]*dy/dS + np[0]->u[2]*dx/dS;
 	qn1 = -np[1]->u[1]*dy/dS + np[1]->u[2]*dx/dS;
               
    H = belp->p[0] - zb;
    Qn = belp->p[1];
	if(belp->bcs[0] == 5)
		Qn = 0.0;
	if(belp->bcs[0] == 1){
		D0 = np[0]->u[0]-(np[0]->ice[0]*tvals.sg_ice);
		D1 = np[1]->u[0]-(np[1]->ice[0]*tvals.sg_ice);
		if((D0 > 0.0) && (D1 > 0.0)){
			qf = pow(D0,1.67)/pow(D1,1.67);
			dQn = (qf - 1.0)/(qf + 1.0) * Qn;
			Qn = ((Qn+dQn)*gbf->f[0] + (Qn-dQn)*gbf->f[1]);
		}
		else if((D0 > 0.0) && (D1 <= 0.0)){
			Qn *= 2.0*gbf->f[1];
		}
		else if((D0 <= 0.0) && (D1 > 0.0)){
			Qn *= 2.0*gbf->f[0];
		}
		else {
			Qn = 0.0;
		}
	}
    Qt = belp->p[2];

    Qx = (-Qn*dy-Qt*dx)/dS;
    Qy = ( Qn*dx-Qt*dy)/dS;

    phi[0] = H;
    phi[1] = Qx;
    phi[2] = Qy;

    if((Hnew-tice)>tvals.gwH) {
        Unew = Qxnew/(Hnew-tice);
        Vnew = Qynew/(Hnew-tice);
    }
    else {
        Unew=0.0;
        Vnew=0.0;
    }
    
    phin[0] = Hnew;
    phin[1] = Qxnew;
    phin[2] = Qynew;
    
    bcCode = belp->bcs[0];
	if(bcCode == 3 && (H < 0.0)){
		bcCode = 0;
	}
    if((Hnew-tice) > tvals.gwH){
        if((qn/(Hnew-tice)) > sqrt((Hnew-tice)*GRAV)){
				if(bcCode == 1){
                        H = Hnew;
						bcCode = 1;
				}
                if(bcCode == 3){
					if((H-tice) > tvals.gwH){
                        Qx = -(H-tice)*sqrt((H-tice)*GRAV) * dy/dS; 
                        Qy =  (H-tice)*sqrt((H-tice)*GRAV) * dx/dS; 
					}
					else{
                        Qx = 0.0; 
                        Qy = 0.0;
					}
						for(i=0;i<ntf;i++){
							Ha = belp->p[0] - np[i]->p[0];
                            KE(i*n,i*n) += 1.0E12;
                            FE(i*n) += 1.0E12 * Ha;
                            KE(i*n+1,i*n+1) += 1.0E12;
                            KE(i*n+2,i*n+2) += 1.0E12;
						    if((Ha-tice) > tvals.gwH){
                                FE(i*n+1) += 1.0E12 * (-(Ha-tice)*sqrt((Ha-tice)*GRAV) * dy/dS);
                                FE(i*n+2) += 1.0E12 * ((Ha-tice)*sqrt((Ha-tice)*GRAV) * dx/dS);
						    }
						    else{
                                FE(i*n+1) += 0.0;
                                FE(i*n+2) += 0.0;
						    }
					}

//					bcCode = 2;
                }
 //               bcCode = 2;                                             //      Supercritical Inflow
        }
        if(qn/(Hnew-tice) < -sqrt((Hnew-tice)*GRAV))
                bcCode = 0;                                             //      Supercritical Outflow
    }

    switch(bcCode) {

        case 0:
            bcc[0] = 0.0;
            bcc[1] = 1.0;
            bcc[2] = 1.0;
            break;

        case 1:
            bcc[0] = 0.0;
            bcc[1] = 1.0;
            bcc[2] = 1.0;
            break;

        case 2:
            bcc[0] = 1.0;
            bcc[1] = 1.0;
            bcc[2] = 1.0;
            break;

        case 3:
            bcc[0] = 1.0;
            bcc[1] = 0.0;
            bcc[2] = 0.0;
            break;

        case 4:
            bcc[0] = 0.0;
            bcc[1] = 0.0;
            bcc[2] = 0.0;   
            break;
            
        case 5:
            bcc[0] = 0.0;
            bcc[1] = 0.0;
            bcc[2] = 0.0;   
            break;
    }

    for(i=0;i<n;i++)
        cbcc[i] = fabs(bcc[i]-1.);
                
    Ha  = bcc[0] * H  + cbcc[0] * Hnew;
    Qxa = bcc[1] * Qx + cbcc[1] * Qxnew;
    Qya = bcc[2] * Qy + cbcc[2] * Qynew;

    if((Ha < 0.0) ) {
        if(belp->bcs[0] == 3){
                for(i=0;i<ntf;i++){
                        Ha = belp->p[0] - np[i]->p[0];
                        if((Ha < 0.0)) {
                                KE(i*n,i*n) += 1.0E12;
                                FE(i*n) += 1.0E12 * Ha;
                        }
                }
        }
        return (0);
    }
                
    if((Ha-tice) >= tvals.gwH){
    	if(belp->bcs[0] > 5) {
    		U = Unew;
    		V = Vnew;
    	}
    	else {
        	U = Qxa/(Ha-tice);
        	V = Qya/(Ha-tice);
        }        
    }
    else {
        U=0.0;
        V=0.0;
    }

	AF[0] = 0.0;
	AF[1] = 0.0;
	AF[2] = 0.0;
    
// A are associated with known boundary values
        
    A[0][0] = 0.0;
    A[0][1] =  dy;
    A[0][2] = -dx;        
    A[1][0] =  GRAV*Ha/2.*dy;
    A[2][0] = -GRAV*Ha/2.*dx;
    A[1][1] =  (U*dy - V*dx);
    A[1][2] =  0.0;
    A[2][1] =  0.0;
    A[2][2] =  (U*dy - V*dx);

// AK are associated with unknown boundary values

    AK[0][0] = 0.0;
    AK[0][1] =  dy;
    AK[0][2] = -dx;
    AK[1][0] =  GRAV*Ha/2.*dy;
    AK[2][0] = -GRAV*Ha/2.*dx;            
    if(qn < 0.0){
        AK[1][1] = (U*dy - V*dx);
        AK[1][2] = 0.0;
        AK[2][1] = 0.0;
        AK[2][2] = (U*dy - V*dx);
    }        
    else {
        AK[1][1] = (U*dy);
        AK[1][2] = (V*dy);
        AK[2][1] = (-U*dx);
        AK[2][2] = (-V*dx);
    }
        
// AJK are the contribution of the boundary terms to the Jacobian  

    AJK[0][0] = 0.0;
    AJK[0][1] = 0.0;
    AJK[0][2] = 0.0;
    AJK[1][0] =  (GRAV*Ha/2. - U*U)*dy + (U*V)*dx;
    AJK[1][1] =  U*dy;
    AJK[2][0] = -(GRAV*Ha/2. - V*V)*dx - (U*V)*dy;
    AJK[2][2] = -V*dx;
    if(qn < 0.0){
        AJK[1][2] = -U*dx;
        AJK[2][1] =  V*dy;
    }        
    else {
        AJK[1][2] =  V*dy;
        AJK[2][1] = -U*dx;
    }
    
    if((belp->bcs[0] == 5) && ((Ha-tice) > tvals.gwH) ) {
    	ex = belp->p[1];
		K = belp->p[0]*pow((Ha-tice),ex-1.0);
    	AK[0][0] =  K*dS;
    	AK[0][1] = 0.0;
    	AK[0][2] = 0.0;
    	AJK[0][0] =  K *(ex-1.0) * dS;
    	AJK[0][1] = 0.0;
    	AJK[0][2] = 0.0;
    	AK[1][0] = (GRAV/2*Ha) * dy;
    	AJK[1][0] = (GRAV/2*Ha) * dy + (ex-1.0)*K*U*dS;
    	AK[2][0] = -(GRAV/2*Ha) * dx;
    	AJK[2][0] = -(GRAV/2*Ha) * dx + (ex-1.0)*K*V*dS;
    	AK[1][1] = K*dS;
    	AJK[1][1] = 0.0;
    	AK[1][2] = 0.0;
    	AJK[1][2] = 0.0;
    	AK[2][1] = 0.0;
    	AJK[2][1] = 0.0;
    	AK[2][2] = K*dS;
    	AJK[2][2] = 0.0;
		AF[0] = K*tice*dS;

    	
    }

    for(i=0;i<ntf;i++) {
//              if((qn > 0.0) && (belp->bcs[0] == 3)) {
//                      KE(i*n,i*n) += 1.0E12;
//                      FE(i*n) += 1.0E12 * (belp->p[0] - np[i]->p[0]);
//              }
		for(l=0;l<n;l++)
				FE(i*n+l) += AF[l]*gbf->f[i];
        for(j=0;j<nbf;j++) {
            for(l=0;l<n;l++) {
                for(m=0;m<n;m++) {
                    KE(i*n+l,j*n+m) += AK[l][m] * cbcc[m] * gbf->f[i] * gbf->f[j];
                    FE(i*n+l) -= A[l][m] * bcc[m] * phi[m]*gbf->f[i]*gbf->f[j];
					if (jacobianType == 0)
						JE(i*n+l,j*n+m) += AJK[l][m] * gbf->f[i] * gbf->f[j]*cbcc[m];
                }
            }                                       
        }
    }
    return (0);    
}

int     find_neighbour()
{
        struct  element *elp;
        struct  belement        *belp;
        int             i,j,k,l,found,n;

//function to find the neighbouring element for each belement

        belp=gp.B;
        
        for(l=0;l<N.belms;l++) {
                found=0;
                elp=gp.El;
                for(i=0;i<N.elms;i++) {
                        for(k=0;k<elp->nnds;k++) {
                                if(belp->nps[0]==elp->nps[k]) {
                                        found=1;
                                        for(j=1;j<belp->nnds;j++) {
                                                if((k+j)>(elp->nnds-1))
                                                        n=k+j-elp->nnds;
                                                else
                                                        n=k+j;
                                                if(belp->nps[j]!=elp->nps[n]){
                                                        found=0;
                                                        break;
                                                }
                                        }
                                        if(found==1) {
                                                belp->elp=elp;
                                                break;
                                        }
                                }
                        }
                        if(found==1)
                                break;
                        else
                                elp=elp->nextelp;
                }
                if(found==0) {
                        printf("Error, no neighbor was found for belement # %d\n",belp->n);
                        return(-1);
                }
                belp=belp->nextbelp;
        }
/*
belp=gp.B;
for(i=0;i<N.belms;i++){
printf("belement=%d\telement=%d\n",belp->n,belp->elp->n);
belp=belp->nextbelp;
}
*/
        return(0);
}

int     transient(nvar,elcode,tol, goaldt, fp,i)
int     nvar,elcode,i ;
double  tol, goaldt ;
FILE* fp;

{
        int     j, itnum, err ;
		double	olddt, oldt;
        FILE    *fptr;
        double  fac,test_outflow();
		int		update_BCs(), set_ptou(), assemble(), solve(), initialize(), assembleHydroDynamics();
		void	setGWFtoZero(), findMaxchange(), useHalfdp();
		double	maxChangeOld;
        
        N.trans = 1 ;
        lumpS = 0;  //should be 0 for transient problems but could try 1 to see if any difference

		outputts = 0;
		set_ptou(0);	//moved here so that thermal has access to old and new nodal values
        
        if((fptr = fp ) != NULL)
		{ 
			olddt = 0.0;
			oldt  = tvals.t;
			tvals.t +=tvals.dt;
			if (tvals.t > goalt)
			{
				olddt = tvals.dt;
				tvals.dt = goalt - oldt;
				tvals.t = goalt;
			}
			update_BCs();
			maxChange = 1.0;
			tvals.itnum = 0;
			itnum = 0;
			for(j=0;j<sets.maxitnum;j++) 
			{  
				if(fabs(maxChange)<tol)
					break ;
			//	if((err = assemble(nvar,0,0,elcode,itnum,fptr)) != 0) return (err) ;
				if((err = initialize(nvar,0,elcode,itnum,fptr)) != 0) return (err) ;
				assembleHydroDynamics(nvar,0,0,elcode,itnum,fptr);
				eqnsets[0].assembleFlag = 1;
				solve(nvar,2,0) ;
				tvals.itnum += 1 ;
				itnum += 1;
				setGWFtoZero();
				maxChangeOld = maxChange;
				findMaxchange();
			//	fprintf(fptr,"itnum = %d\t maxChange = %f\t maxNode = %d\n", tvals.itnum, maxChange,maxnodeP->n);

				if(maxChange == 0.0){
					maxChange = 1.0;				
					break;
				}
				if(fabs(maxChange) > tvals.uc*avgCelerity())
					break;

				if(fabs(maxChange) > fabs(0.9*maxChangeOld)){
					useHalfdp(0);
					findMaxchange();
				}
			}

			if(fabs(maxChange) > tol) //reject the time step
			{
				fprintf(fptr,"\n\tdt = %f didn't work\n",tvals.dt);
				fprintf(fptr,"\tMaximum change at node %d \tx = %f \ty = %f\n\t\t from h = %f qx = %f qy = %f \n\t\t to   h = %f qx = %f qy = %f\n",
                      maxnodeP->n,maxnodeP->x[0],maxnodeP->x[1],maxnodeP->uo[0],
					  maxnodeP->uo[1],maxnodeP->uo[2],
                      maxnodeP->u[0],maxnodeP->u[1],maxnodeP->u[2]);

				tvals.t -= tvals.dt;
				fac = 0.5;
				reset_u(0); //reset the values at the nodes to the values 
							//at the previous time step
				tvals.dt *= fac;
			}
			else
			{
                fprintf(fptr,"time step # %d\tt = %fs\tdt = %fs\ttotal inflow = %g\ttotal outflow = %g\t# of N.R. iterations = %d Maxchange = %g\n",i,tvals.t,tvals.dt,test_outflow(1),test_outflow(2),tvals.itnum,maxChange);
        //       set_ptou(0) ;
				updateVelocities();
				fac = 1.5;
				tvals.dt *= fac;
				if (olddt != 0.0)
					tvals.dt = olddt;
				if (tvals.t == goalt)
				{
					goalt += goaldt;
					outputts = 1;
					outputtsnum += 1;
				}
			}
			if(tvals.dt > goaldt)
				tvals.dt = goaldt;
        }
        return(0) ;
}

int     steady(nvar,elcode,tol,dtmax,fp,i)
int     nvar,elcode,i ;
double  tol,dtmax ;
FILE* fp;

{
        int     itnum, err ;
        FILE    *fptr;
		double	fac, allowChg=0.0, avgCel;
		int		update_BCs(), assemble(), solve(), reset_u(), set_ptou(), initialize(), assembleHydroDynamics();
		double	avgCelerity(), test_outflow();
		void	setGWFtoZero(), findMaxChange();
		
        N.trans = 1 ;
        lumpS = 1;
		
		if((fptr = fp ) != NULL)
		{
            tvals.t += tvals.dt ;
			update_BCs();
			itnum = 0 ;
		
		//	if( (err = assemble(nvar,0,0,elcode,itnum,fptr)) != 0) return(err);
			if((err = initialize(nvar,0,elcode,itnum,fptr)) != 0) return (err) ;
			assembleHydroDynamics(nvar,0,0,elcode,itnum,fptr);
			eqnsets[0].assembleFlag = 1;
			solve(nvar,2,0);
			
			// Start of new code // Sept 11, 2009

			setGWFtoZero();
			findMaxchange();
			avgCel = avgCelerity();
			uchange=fabs(maxChange)/avgCel;
			allowChg=tvals.uc*avgCel;
			
			if(fabs(maxChange) > 1.5 * allowChg) //if maxChg is too big reject iteration and write info to logfile							
			{				
				fprintf(fptr,"\n\tmaxChange > 1.5*allowChange Time Step Rejected\n");
				fprintf(fptr,"# %d  t= %f  dt= %f allowablechange=%E maxchange=%E solution change= %E\n",
								i,tvals.t,tvals.dt,allowChg, maxChange,uchange);
				fprintf(fptr,"\tmax change occur at node %d \tx = %f \ty = %f\n\t\t from h = %f qx = %f qy = %f \n\t\t to h = %f qx = %f qy = %f\n\n",
					maxnodeP->n,maxnodeP->x[0],maxnodeP->x[1],maxnodeP->uo[0],
					maxnodeP->uo[1],maxnodeP->uo[2],
					maxnodeP->u[0],maxnodeP->u[1],maxnodeP->u[2]);

				tvals.t -= tvals.dt;
				fac = 0.5 * allowChg/fabs(maxChange);  //lower the time step
				reset_u(0); //reset the values at the nodes to the values 
							//at the previous time step
			}
			else if(fabs(maxChange) > allowChg) // if true, then reduce the time step but accept iteration and write info to logfile
			{
				fprintf(fptr,"\n\tmaxChg > allowChg Time Step Accepted but deltat reduced\n");
				fprintf(fptr,"# %d  t= %f  dt= %f  allowablechange= %E  maxchange= %E solution change= %E\n",
								i,tvals.t,tvals.dt,allowChg,maxChange,uchange);
				fprintf(fptr,"\tmax change occur at node %d \tx = %f \ty = %f\n\t\t from h = %f qx = %f qy = %f \n\t\t to h = %f qx = %f qy = %f\n\n",
					maxnodeP->n,maxnodeP->x[0],maxnodeP->x[1],maxnodeP->uo[0],
					maxnodeP->uo[1],maxnodeP->uo[2],
					maxnodeP->u[0],maxnodeP->u[1],maxnodeP->u[2]);

				
				set_ptou(0); //copies the values at the nodes this time step to the space 
							 //allocated for the those at the previous time step
				updateVelocities();
				fac = 0.5*allowChg/fabs(maxChange);
			}
			else //amount of change is acceptable; write results to logfile and increase time step
			{
				fprintf(fptr,"\n\tmaxChg < allowChg Time Step Accepted\n");
				fprintf(fptr,"# %d  t= %f  dt= %f  allowablechange= %E  maxchange= %E uchange= %E\n",
								i,tvals.t,tvals.dt,allowChg,maxChange,uchange);
					fprintf(fptr,"\tmax change occur at node %d \tx = %f \ty = %f\n\t\t from h = %f qx = %f qy = %f \n\t\t to h = %f qx = %f qy = %f\n\n",
					maxnodeP->n,maxnodeP->x[0],maxnodeP->x[1],maxnodeP->uo[0],
					maxnodeP->uo[1],maxnodeP->uo[2],
					maxnodeP->u[0],maxnodeP->u[1],maxnodeP->u[2]);


				updateVelocities();
				set_ptou(0); //copies the values at the nodes this time step to the space 
							 //allocated for the those at the previous time step
				fac = allowChg/fabs(maxChange);
				if(fac > 1.5) fac = 1.5;
			}

			// end of new Sept 11, 2009 code

			tvals.dt *=fac;
			if(tvals.dt > dtmax)
				tvals.dt = dtmax;
		}
        return(0);
}


/****  this version of steady produces a modified logfile used for the purpose of testing the speed of the model solvers
int     steady(nvar,elcode,tol,dtmax,fp,i)
int     nvar,elcode,i ;
double  tol,dtmax ;
FILE* fp;

{
        int     j, jj, itnum, day, hr, min, sec, hnum1, hnum2, err ;
        long    time ;
        struct  node    *nodep, *maxnodep,*minnodep,*zeronp , *maxFrnodep;
        FILE    *fptr, *log_fptr() ;
        double  minh,maxcourant,courant,dx,u,c,zeroh, froude, maxFr, fac, test_outflow();
        double	uchange1, chgInH, maxChgInH;
        void            trans_bcs();
		clock_t startass, finishass, startsolve, finishsolve, start, finish;
		double durationass, durationsolve, durationtotal, duration;
        
        N.trans = 1 ;
        lumpS = 1;
		start = clock();
		
		if((fptr = fp ) != NULL)
		{
            tvals.t += tvals.dt ;
			update_BCs();
			uchange = 1.0;
			itnum = 0 ;
			for(j=0;j<1;j++) //for statement here so that break statement can be used
			{
				if(uchange<tol) break;
				startass = clock();
                if( (err = assemble(nvar,0,0,elcode,itnum,fptr)) != 0) return(err);
				finishass = clock();
				durationass = (double)(finishass - startass) / CLOCKS_PER_SEC;
                assembleFlag = 1;
				startsolve = clock();
                solve(nvar,2,0);
				finishsolve = clock();
				durationsolve = (double)(finishsolve - startsolve) / CLOCKS_PER_SEC;
				durationtotal = (double)(finishsolve - startass) / CLOCKS_PER_SEC;
				
                if(itnum == 0) uchange1 = uchange;//this will always be true because itnum
												 //is set to 0 every time steady is called
												 //itnum (# of iterations) is incremented
											     //in front end
				itnum += 1;

			}
            if(uchange1 > 1.25 * tvals.uc)  //if uchange is too big reject iteration 
											//and write info to logfile
			{
				nodep = gp.N;
				maxnodep = gp.N;
				maxChgInH = 0.0;
				for(jj=0;jj<N.nodes;jj++) //figuring out at which node 
										  //the max depth change occurred
				{
					chgInH = fabs(nodep->u[0] - nodep->uo[0]);
					if(chgInH > maxChgInH)
					{
						maxChgInH = chgInH;
						maxnodep=nodep;
					}
					nodep=nodep->nextnp;
				}
				nodep = maxnodep;
				
				fprintf(fptr,"\n\tdt = %f didn't work\n",tvals.dt);
				fprintf(fptr,"\tmax depth change at node %d \tx = %f \ty = %f\n\t\t from h = %f qx = %f qy = %f \n\t\t to   h = %f qx = %f qy = %f\n",
                      nodep->n,nodep->x[0],nodep->x[1],nodep->uo[0],
					  nodep->uo[1],nodep->uo[2],
                      nodep->u[0],nodep->u[1],nodep->u[2]);

				tvals.t -= tvals.dt;
				fac = 0.5 * tvals.uc/uchange1;  //lower the time step
				reset_u(0); //reset the values at the nodes to the values 
							//at the previous time step
			}
			else // if uchange is o.k. then write results to logfile and increase time step
			{
				fprintf(fptr,"%d, t = %fs\tdt = %fs\tnet outflow = %g\tuchange = %f\t",
                                			i,tvals.t,tvals.dt,test_outflow(),uchange);
				set_ptou(0); //copies the values at the nodes this time step to the space 
							 //allocated for the those at the previous time step
				updateVelocities();
				fac = tvals.uc/uchange1; //factor = goalchange (default 0.05) / uchange 
				if(fac >1.5) fac = 1.5;
			}
			tvals.dt *=fac;
			if(tvals.dt > dtmax)
				tvals.dt = dtmax;
		}
		finish = clock();
		duration = duration = (double)(finish - start) / CLOCKS_PER_SEC;
		fprintf(fptr,"%f\t %f\t %f\t %f\n", durationass, durationsolve, durationtotal, duration);
        return(0);
}
****/

int             assemble(varnum,code,nset,elcode,itnum,logfptr)
int             varnum, code, nset, elcode, itnum  ;
FILE			*logfptr;

{
        double					*K_init(), *L_init(),*F_init() ;
        struct fxnodes          *fn_init();
        struct element          *elp, tel ;
        struct belement         *belp, tbel ;
        double                  res ;
        register double			*dp;
        int						i,ns,ntf,istart,iend,ind ;
        FILE					*put_fptr() ;
		int						EndMes(), ep_init();
        
        istart = 0 ;
        iend = N.elms ;
        if(N.trans == 0)
                tvals.dt = 1.0 ;
        
/* note: nset is the global matrix number.  up to 4 different global matrices possible. only one is used. nset is set to zero. */       
        if(code< 2) 
		{
            if(eqnsets[nset].assembleFlag == 0)
			{
                if((eqnsets[nset].Kp != NULL) && (eqnsets[nset].assembleFlag == 0) )
                        free(eqnsets[nset].Kp) ;
        
                if( (eqnsets[nset].Kp = K_init(varnum,&eqnsets[nset].diagp,itnum)) == NULL )
                        return(-1) ;
        
                if(N.sym != 0) 
				{
					if((eqnsets[nset].Lp != NULL) && (eqnsets[nset].assembleFlag == 0) ) 
					   free(eqnsets[nset].Lp) ;

                    if( (eqnsets[nset].Lp = L_init(varnum,&eqnsets[nset].diagp,itnum)) == NULL )
                       return(-1) ;
                }
                        
                eqnsets[nset].neqns = Nukns ;
				eqnsets[nset].NelemInK = NelemInK;
     
                if((eqnsets[nset].Fp != NULL) && (eqnsets[nset].assembleFlag == 0))
                        free(eqnsets[nset].Fp) ;
        
                if( ( eqnsets[nset].Fp = F_init(Nukns) )  == NULL  )
                        return(-1) ;

				if((eqnsets[nset].Bp != NULL) && (eqnsets[nset].assembleFlag == 0))
					free(eqnsets[nset].Bp) ;
				if((eqnsets[nset].Qp != NULL) && (eqnsets[nset].assembleFlag == 0))
					free(eqnsets[nset].Qp) ;
				if((eqnsets[nset].Rp != NULL) && (eqnsets[nset].assembleFlag == 0))
					free(eqnsets[nset].Rp) ;
				if((eqnsets[nset].Dp != NULL) && (eqnsets[nset].assembleFlag == 0))
					free(eqnsets[nset].Dp) ;
				if((eqnsets[nset].Hp != NULL) && (eqnsets[nset].assembleFlag == 0))
					free(eqnsets[nset].Hp) ;
				if((eqnsets[nset].Vp != NULL) && (eqnsets[nset].assembleFlag == 0))
					free(eqnsets[nset].Vp) ;
				if((eqnsets[nset].Mp != NULL) && (eqnsets[nset].assembleFlag == 0))
					free(eqnsets[nset].Mp) ;
				if((eqnsets[nset].Yp != NULL) && (eqnsets[nset].assembleFlag == 0))
					free(eqnsets[nset].Yp) ;

            }
            else if(eqnsets[nset].assembleFlag == 1)
			{
                        dp = eqnsets[nset].Kp;
                for(i=0;i<eqnsets[nset].NelemInK;i++)
                       *(dp++) = 0.0 ;
                        dp = eqnsets[nset].Lp;
                for(i=0;i<eqnsets[nset].NelemInK;i++)
                       *(dp++) = 0.0 ;
            }
        }
		else if (code == 2) 
		{
                ind = 0 ;
                for(i=0;i<N.nodes;i++) 
				{
                        gp.iptrs[i]->ui = ind ;
                        ind++ ;
                }
        }
        if((code < 1) && (itnum == 0)) 
		{       
			 if( (ep.n = ep_init() ) == 0 )
				 return(-1) ;
        }
        for(i=0;i<eqnsets[nset].neqns;i++)
                *(eqnsets[nset].Fp + i) = 0.0 ;
        
		elp = gp.El ;
        
		for(i=istart;i<iend;i++)
		{
			if(jacobianType == 0)
				ns = get_PGKeAnalJ(elp,&tel,varnum,elcode,&ntf,itnum) ;
			else
 				ns = get_PGKeNumJ(elp,&tel,varnum,elcode,&ntf,itnum) ;
            put_Ke(ns,ntf,elp,varnum,eqnsets[nset].Kp,eqnsets[nset].Lp,eqnsets[nset].Fp,eqnsets[nset].diagp,code) ;
            elp = elp->nextelp ;
        }

        belp = gp.B ;

        for(i=0;i<N.belms;i++)
		{ 
			if(jacobianType == 0)
				ns = get_bKeAnalJ(belp,&tbel,varnum,elcode,&ntf,itnum) ;
			else
				ns = get_bKeNumJ(belp,&tbel,varnum,elcode,&ntf,itnum) ;
            if(ns > 0 )
				put_Ke(ns,ntf,belp,varnum,eqnsets[nset].Kp,eqnsets[nset].Lp,eqnsets[nset].Fp,eqnsets[nset].diagp,code) ;
            belp = belp->nextbelp ;
        }
	
//	Node withdrawals
/*		np = gp.N;
		for(i=0;i<N.nodes;i++) 
		{
			if(np->p[2] <0.0)
				*(eqnsets[nset].Fp + 3*i) += np->p[2]*tvals.dt;
			np = np->nextnp;
		}
*/
		resL2 = 0.0;
		resLinf = 00;
		resmaxi = -1;
        for(i=0;i<eqnsets[nset].neqns;i++){
                res = *(eqnsets[nset].Fp + i);
				resL2 += res*res;
				if(fabs(res) > fabs(resLinf)){
					resLinf = res;
					resmaxi = i;
				}
		}
		resL2 = sqrt(resL2);

//		if(resmaxi >= 0)
//			probNodes(0.5*resL2,varnum,code,elcode,nset);

        if(code < 2)
		{       
			EndMes() ;
        }
   
        return(0) ;
}

int             initialize(varnum,nset,elcode,itnum,logfptr)
int             varnum, nset, elcode, itnum  ;
FILE			*logfptr;

{
        double					*K_init(), *L_init(),*F_init() ;
        register double			*dp;
        int						i ;
        FILE					*put_fptr() ;
		int						ep_init();
        
/* note: nset is the global matrix number.  up to 4 different global matrices possible. only one is used. nset is set to zero. */       
        
            if(eqnsets[nset].assembleFlag == 0)
			{
                if((eqnsets[nset].Kp != NULL) && (eqnsets[nset].assembleFlag == 0) )
                        free(eqnsets[nset].Kp) ;
        
                if( (eqnsets[nset].Kp = K_init(varnum,&eqnsets[nset].diagp,itnum)) == NULL )
                        return(-1) ;
        
                if(N.sym != 0) 
				{
					if((eqnsets[nset].Lp != NULL) && (eqnsets[nset].assembleFlag == 0) ) 
					   free(eqnsets[nset].Lp) ;

                    if( (eqnsets[nset].Lp = L_init(varnum,&eqnsets[nset].diagp,itnum)) == NULL )
                       return(-1) ;
                }
                        
                eqnsets[nset].neqns = Nukns ;
				eqnsets[nset].NelemInK = NelemInK;
     
                if((eqnsets[nset].Fp != NULL) && (eqnsets[nset].assembleFlag == 0))
                        free(eqnsets[nset].Fp) ;
        
                if( ( eqnsets[nset].Fp = F_init(Nukns) )  == NULL  )
                        return(-1) ;

				if((eqnsets[nset].Bp != NULL) && (eqnsets[nset].assembleFlag == 0)){
					free(eqnsets[nset].Bp) ;
					eqnsets[nset].Bp = NULL;
				}
				if((eqnsets[nset].Qp != NULL) && (eqnsets[nset].assembleFlag == 0)){
					free(eqnsets[nset].Qp) ;
					eqnsets[nset].Qp = NULL;
				}
				if((eqnsets[nset].Rp != NULL) && (eqnsets[nset].assembleFlag == 0)){
					free(eqnsets[nset].Rp) ;
					eqnsets[nset].Rp = NULL;
				}
				if((eqnsets[nset].Dp != NULL) && (eqnsets[nset].assembleFlag == 0)){
					free(eqnsets[nset].Dp) ;
					eqnsets[nset].Dp = NULL;
				}
				if((eqnsets[nset].Hp != NULL) && (eqnsets[nset].assembleFlag == 0)){
					free(eqnsets[nset].Hp) ;
					eqnsets[nset].Hp = NULL;
				}
				if((eqnsets[nset].Vp != NULL) && (eqnsets[nset].assembleFlag == 0)){
					free(eqnsets[nset].Vp) ;
					eqnsets[nset].Vp = NULL;
				}
				if((eqnsets[nset].Mp != NULL) && (eqnsets[nset].assembleFlag == 0)){
					free(eqnsets[nset].Mp) ;
					eqnsets[nset].Mp = NULL;
				}

            }
            else if(eqnsets[nset].assembleFlag == 1)
			{
                        dp = eqnsets[nset].Kp;
                for(i=0;i<eqnsets[nset].NelemInK;i++)
                       *(dp++) = 0.0 ;
                        dp = eqnsets[nset].Lp;
                for(i=0;i<eqnsets[nset].NelemInK;i++)
                       *(dp++) = 0.0 ;
            }
        
        if((itnum == 0)) 
		{       
			 if( (ep.n = ep_init() ) == 0 )
				 return(-1) ;
        }
        for(i=0;i<eqnsets[nset].neqns;i++)
                *(eqnsets[nset].Fp + i) = 0.0 ;
   
        return(0) ;
}

int             assembleHydroDynamics(varnum,code,nset,elcode,itnum,logfptr)
int             varnum, code, nset, elcode, itnum  ;
FILE			*logfptr;

{
        double					*K_init(), *L_init(),*F_init() ;
        struct element          *elp, tel ;
        struct belement         *belp, tbel ;
        double                  res ;
        register double			*dp;
        int						i,ns,ntf,ind ;
        FILE					*put_fptr() ;
		int						EndMes(), ep_init();
        
		elp = gp.El ;
        
		for(i=0;i<N.elms;i++)
		{
			if(jacobianType == 0)
				ns = get_PGKeAnalJ(elp,&tel,varnum,elcode,&ntf,itnum) ;
			else
 				ns = get_PGKeNumJ(elp,&tel,varnum,elcode,&ntf,itnum) ;
            put_Ke(ns,ntf,elp,varnum,eqnsets[nset].Kp,eqnsets[nset].Lp,eqnsets[nset].Fp,eqnsets[nset].diagp,code) ;
            elp = elp->nextelp ;
        }

        belp = gp.B ;

        for(i=0;i<N.belms;i++)
		{ 
			if(jacobianType == 0)
				ns = get_bKeAnalJ(belp,&tbel,varnum,elcode,&ntf,itnum) ;
			else
				ns = get_bKeNumJ(belp,&tbel,varnum,elcode,&ntf,itnum) ;
            if(ns > 0 )
				put_Ke(ns,ntf,belp,varnum,eqnsets[nset].Kp,eqnsets[nset].Lp,eqnsets[nset].Fp,eqnsets[nset].diagp,code) ;
            belp = belp->nextbelp ;
        }
   
        return(0) ;
}

int	probNodes(resLine,varnum,code,elcode,nset)
double resLine;
int varnum, code, elcode, nset;

{
	struct element	*elp;
	struct node		*np;
	int				iN, iV, ntf, ns;
	int				nProbElms = 0;

	elp = gp.El;
	while(elp != NULL){
		for(iN=0;iN<elp->nnds;iN++){
			np = elp->nps[iN];
			for(iV=0;iV<N.vars;iV++){
				if(*(eqnsets[nset].Fp + (np->i)*N.vars + iV) > resLine){
					ns = get_PGSe(elp,varnum,elcode,&ntf);
		            put_Ke(ns,ntf,elp,varnum,eqnsets[nset].Kp,eqnsets[nset].Lp,eqnsets[nset].Fp,eqnsets[nset].diagp,code);
					iV = N.vars;
					iN = elp->nnds;
					nProbElms++;
				}
			}
		}
		elp = elp->nextelp;
	}
	return(nProbElms);
}

int             get_PGSe(theElp,nvar,code,ntf)
struct element	*theElp ;
int             nvar, code, *ntf;

{
        struct gausspts         g[MGPTS] ;
        struct shapefuncs       bf, tf;
        struct node				*np[MNSF] ;
        struct govfuncs			gbf, gtf ;
        int		i, j, nbf, ni, n, jj ;
        double	dt, dtmin, sdiag ,wt;
        double	get_mlt();
		int		initm(), get_PGgf(), check_mixed(), get_mixedgps();

        if(nvar<0)
                n = N.vars ;
        else
                n = 1 ;
		
		*ntf = nsf(theElp->gtype) ;
		nbf = theElp->nnds ;
        
        initm(ep.K,ep.n*ep.n) ;
		initm(ep.S,ep.n*ep.n) ;
		initm(ep.J,ep.n*ep.n) ;
		initm(ep.F,ep.n) ;
						
		if(sets.smallH>0)
		{
			if(check_mixed(theElp)==1)
			{
				ni=get_mixedgps(theElp,g);
			}
			else
				ni = get_gspts(theElp->vtype,g) ;
		}
		else
			ni = get_gspts(theElp->vtype,g) ;
		for(j=0;j<theElp->nnds;j++)
			np[j] = theElp->nps[j];
		wt=0.0;
		for(i=0;i<ni;i++)
		{
			wt+=g[i].w;
		}
		for(i=0;i<ni;i++) {
			theElp->gpcode[i] = 1;
		}
		dtmin = 1000000.0;
		for(i=0;i<ni;i++)
		{
			gpIndex = i;
			get_shape(theElp->vtype,&g[i],&bf) ;
			get_shape(theElp->gtype,&g[i],&tf) ;
			get_PGgf(theElp,&gbf,&gtf,&bf,&tf,np,g[i].w) ;
			dt = update_PGKe(&gbf,&gtf,nbf,*ntf,np,nvar,theElp,wt) ;
			if(dt<dtmin)
				dtmin = dt;
		}
		
		if(N.trans > 0)
		{ /* mass matrix lumping */
			if(lumpS == 1)
			{
				for(i=0;i<(*ntf)*n;i++)
				{
					sdiag = 0.0 ;
					for(j=0;j<(nbf)*n;j++)
					{
						sdiag += SE(i,j) ;
						SE(i,j) = 0.0 ;
					}
					SE(i,i) = sdiag ;
				}
			}
			for(i=0;i<(*ntf)*n;i++){
				FE(i) = 0.0;
				for(j=0;j<nbf;j++){
					for(jj=0;jj<n;jj++){
						SE(i,j*n+jj) /= dtmin ;
					}
				}
			}
		}		 		
		return(nbf) ;
}


double  *K_init(nvar,thediagp,code)
int             nvar, code ;
long            *(*thediagp) ;

{
        int     i, j, nff , ind, n, *mpdArray, d;
        double *kp, *dp ;
        unsigned        elsize ;
        long    bigint, *diagp, il ;
        static long nelem;
		struct element *elp;
        
        if(nvar<0)
                n = N.vars ;
        else
                n = 1 ;
        Nukns = N.nodes * n;
        if(code==0) {
                if(*thediagp != NULL)
                        free(*thediagp) ;
                elsize = sizeof(long) ;
                if((*thediagp = (long int *) calloc((int)Nukns,sizeof(long) )) == NULL)
                        return(NULL) ;
                diagp = *thediagp ;
                *diagp = 0 ;
                ind = 0 ;
                gp.iptrs[0]->ui = ind ;
                diagp++ ;
                ind++ ;
                for(i=1;i<n;i++) {
                        nff = i +1;
                        *diagp = * (diagp-1) + nff;
                        diagp++;
                        ind++;
                }

				mpdArray = calloc(N.nodes, sizeof(int));
				elp = gp.El;
				for(i=0;i<N.elms;i++){
					d = elp->nps[2]->n - elp->nps[1]->n;
					if(d > 0){
						if(d > mpdArray[(elp->nps[2]->n) - 1]){
							mpdArray[(elp->nps[2]->n) - 1] = d;
						}
					}
					else{
						d *= -1;
						if(d > mpdArray[(elp->nps[1]->n) - 1]){
							mpdArray[(elp->nps[1]->n) - 1] = d;
						}
					}
					d = elp->nps[1]->n - elp->nps[0]->n;
					if(d > 0){
						if(d > mpdArray[(elp->nps[1]->n) - 1]){
							mpdArray[(elp->nps[1]->n) - 1] = d;
						}
					}
					else{
						d *= -1;
						if(d > mpdArray[(elp->nps[0]->n) - 1]){
							mpdArray[(elp->nps[0]->n) - 1] = d;
						}
					}
					d = elp->nps[2]->n - elp->nps[0]->n;
					if(d > 0){
						if(d > mpdArray[(elp->nps[2]->n) - 1]){
							mpdArray[(elp->nps[2]->n) - 1] = d;
						}
					}
					else{
						d *= -1;
						if(d > mpdArray[(elp->nps[0]->n) - 1]){
							mpdArray[(elp->nps[0]->n) - 1] = d;
						}
					}
					elp = elp->nextelp;
				}
                for(i=1;i<N.nodes;i++) {        
                        nff =  (mpdArray[i]  + 1) * n - (n-1) ;
//                      nff =  (mpd_i(i,nvar)  + 1) * n - (n-1) ;
                        for(j=0;j<n;j++) {
                                gp.iptrs[i]->ui = ind;
                                ind++;
                                *diagp = * (diagp-1) + nff;
                                diagp++;
                                nff +=1;

                        }
                }
                nelem = *(diagp-1) + 1 ;/*calc's # of entries in assembled skyline matrix*/
 				free(mpdArray);
        }
        NelemInK = nelem;
        elsize = sizeof(double) ;
        bigint =  nelem * (long) elsize  ;
        kp =  dp = (double *) malloc(bigint)  ;                       
		if (kp == NULL) {
			printf(" memory allocation failure in K_init, size = %ld\n", bigint);
		} 
		else if (dp == NULL) {
			printf(" memory allocation failure in K_init, size = %ld\n", bigint);
		} else {
                for(il=0;il<nelem;il++)
                        *(dp++) = 0.0 ;
        }
        return(kp) ;
}



double  *L_init(nvar,thediagp,code)
int              nvar,code ;
long            *(*thediagp) ;

{
        int     i;
        double *kp, *dp ;
        unsigned  elsize ;
        long    bigint, nelem, il, *diagp ;
        
        diagp = *thediagp ;
        for(i=0;i<Nukns;i++) /* FAYE: changed to Nukns 88-08-22 */
                        diagp++ ;
        
        nelem = *(diagp-1) + 1 ;
        elsize = sizeof(double) ;
        bigint =  nelem * (long) elsize  ;
        kp =  dp = (double *) malloc(bigint)  ;           
        if(kp == NULL)
                printf(" memory allocation failure in L_Init, size = %ld\n",bigint);
        else {
                for(il=0;il<nelem;il++)
                        *(dp++) = 0.0 ;
				memJacobian = 2.0 * bigint/1024.0/1024.0;
//              printf(" Memory allocated for Jacobian Matrix = %10.3f MBytes\n",memJacobian);
        }
        return(kp) ;
}


double  *F_init(num)
int             num ;

{
        double  *p ;
        unsigned        nelem, elsize ;
        
        nelem = num ;
        elsize = sizeof(double) ;
        p = (double *) calloc(nelem,elsize) ;
        return(p) ;
}

struct fxnodes  *fn_init(nvar)
int             nvar ;

{       
        int             i, j ;
        struct fxnodes  *p, *ap ;
        unsigned        nelem, elsize ;
        
        Nrfixed = 0 ;
        for(i=0;i<N.nodes;i++) {
                if(isvarfixed(i,nvar) == -1) 
                        Nrfixed++ ;
        }
        nelem = Nrfixed ;
        elsize = sizeof(struct fxnodes) ;
        ap = p = (struct fxnodes *) calloc(nelem,elsize) ;
        for(i=0;i<(int)nelem;i++) {
                ap->i = -1 ;
                for(j=0;j<MDOF;j++)
                        ap->j[j] = -1 ;
                ap++ ;
        }
        return(p) ;
}

int             ep_init()

{
        int     i, ns, nsm ;
        struct element  *elp ;
        
        if(ep.K != NULL)
                free(ep.K) ;
        elp = gp.El ;
        nsm = 0 ;
        for(i=0;i<N.elms;i++) {
                nsm = ( (ns = nsf(elp->vtype)) > nsm ) ? ns : nsm ;
                elp = elp->nextelp ;
        }
        nsm *= N.vars ;

/*      printf("Size of Element Matrices = %d\n",nsm) ;*/

        if( (ep.K = (double *) calloc(nsm*nsm, sizeof(double)) ) == NULL )
                return(0) ;
        if(ep.F != NULL)
                free(ep.F) ;
        if( (ep.F = (double *) calloc(nsm, sizeof(double)) ) == NULL )
                return(0) ;
        if(N.trans > 0) {
                if(ep.S != NULL)
                        free(ep.S) ;
                if( (ep.S = (double *) calloc(nsm*nsm, sizeof(double)) ) == NULL )
                        return(0) ;
        }
        if(N.trans > 1) {
                if(ep.M != NULL)
                        free(ep.M) ;
                if( (ep.M = (double *) calloc(nsm*nsm, sizeof(double)) ) == NULL )
                return(0) ;
        }
        if(N.nonlin > 0) {
                if(ep.J != NULL)
                        free(ep.J) ;
                if( (ep.J = (double *) calloc(nsm*nsm, sizeof(double)) ) == NULL )
                return(0) ;
        }
        if(N.nonlin > 0) {      
                if(ep.DM != NULL)
                        free(ep.DM) ;
                if( (ep.DM = (double *) calloc(nsm*nsm*nsm, sizeof(double)) ) == NULL )
                return(0) ;
        
}
        return(nsm) ;
}

int             put_Ke(ns,nt,elmntp,nvar,Kp,Lp,Fp,diagp,code) 
int             ns, code ;
double  *Kp,  *Lp, *Fp ;
struct element *elmntp ;
long    *diagp ;

{
        int     i, nvars, j, ii, jj, frow, fcol, row, col ;
        struct node     *inp, *jnp ;

        if(nvar < 0)
                nvars = N.vars ;
        else
                nvars = 1 ;
        
        for(i=0;i<nt;i++) {
                inp = elmntp->nps[i] ;
                frow = inp->i * nvars ;
                for(j=0;j<ns;j++) {
                        jnp = elmntp->nps[j] ;
                        fcol = jnp->i * nvars;
                        for(ii=0;ii<nvars;ii++) {
                                row = frow + ii ;
                                for(jj=0;jj<nvars;jj++) {
                                        col = fcol + jj ;
/*      printf("\n row =%d, col =%d, SE =%g\n",i*nvars+ii,j*nvars+jj,SE(i*nvars+ii,j*nvars+jj)) ;*/
                                        if( (col >= row) && ( code < 2 ) && ( code > -1 ) ) {
                                                if(N.trans == 0)
                                                        *(Kp + *(diagp+col) - col + row ) += KE(i*nvars+ii,j*nvars+jj) ;
                                                else
                                                        *(Kp + *(diagp+col) - col + row ) += SE(i*nvars+ii,j*nvars+jj) ;
/*printf("row =%d, col =%d, upper =%g\n",row,col,*(Kp + *(diagp+col) - col + row )) ;*/
                                        }
                                        if( (col <= row) && (N.sym != 0 ) ) {
                                                if(N.trans == 0)
                                                        *(Lp + *(diagp+row) - row + col ) += KE(i*nvars+ii,j*nvars+jj) ;
                                                else
                                                        *(Lp + *(diagp+row) - row + col ) += SE(i*nvars+ii,j*nvars+jj) ;
/*printf("row =%d, col =%d, Lower =%lf\n",row,col,*(Lp + *(diagp+row) - row + col )) ;*/
                                        }
                                        if( (col == row) && ( code == -1 ) ) {
                                                if(N.trans == 0)
                                                        *(Kp + row ) += KE(i*nvars+ii,j*nvars+jj) ;
                                                else
                                                        *(Kp + row ) += SE(i*nvars+ii,j*nvars+jj) ;
                                        }
                                        if(  (code == 2 ) || code == -1 ) {
                                                *(Fp  + row ) += - KE(i*nvars+ii,j*nvars+jj) * jnp->u[jj] ;
                                        }
                                }
                        }
                }
                for(ii=0;ii<nvars;ii++) {
                        row = frow + ii ;
                        *(Fp + row) += FE(i*nvars+ii) ;
                }
        }
		return (0);
}
/*
int             calc_BF(nvar,nset)
int             nvar, nset ;

{
        int             i ,assBF();
        
        assBF(nvar,0,nset,0);
        solve(0,3,1) ;
        for(i=0;i<Nrfixed;i++) {
                gp.iptrs[*(BFnodes+i)]->p[5] = *(eqnsets[nset].Fp + i) ;
        }
        return(0) ;
}

int             assBF(varnum,code,nset,elcode)
int             varnum, code, nset, elcode ;

{
        double                          *KBF_init(), *L_init(),*F_init() ;
        struct belement         *belp ;
        int			i,j,ns = 1,ind, prevnode, nsl, itnum ;
        char                                            m1[255], m2[255], m3[255], m4[255] ;
		int bnd_fluxes();
        struct fxnodes *the_fnp ;
               
        bnd_fluxes(0) ;

        itnum = 0 ;
        
        if(code< 2) {
                sprintf(m1," Assembling Boundary Flux Matrix ");
                sprintf(m2," Number of equations to be solved = %d ",Nrfixed) ;
                sprintf(m3," ") ;
                sprintf(m4," ") ;
                StartMes(m1,m2,m3,m4) ;
                
                if(eqnsets[nset].Kp != NULL)
                        free(eqnsets[nset].Kp) ;
                if( (eqnsets[nset].Kp = KBF_init(varnum,&eqnsets[nset].diagp,code)) == NULL )
                        return(-1) ;
                        
                eqnsets[nset].neqns = Nrfixed ;
                if( (BFnodes = ( int *) calloc(Nrfixed+1, sizeof(int)) ) == NULL )
                        return(0) ;

                sprintf(m1," Assembling Boundary Flux Matrix ");
                sprintf(m2," Number of equations to be solved = %d ",Nrfixed) ;
                sprintf(m3," Elements in K matrix - %ld ",*(eqnsets[nset].diagp+Nrfixed-1)+1) ;
                sprintf(m4," ") ;
                UpDateMes(m1,m2,m3,m4) ;
        
                if(eqnsets[nset].Fp != NULL) {
                        free(eqnsets[nset].Fp) ;
                }
                if( ( eqnsets[nset].Fp = F_init(Nrfixed) )  == NULL  ) 
                        return(-1) ;
        }
        belp = gp.B ;
        ind = 0 ;
        nsl = 0 ;
        prevnode = -9 ;
        for(i=0;i<N.belms;i++) {
//                ns = get_bKe(belp,varnum,2,itnum) ;
                if(ns > 0 ) {
                        if(prevnode == (belp->nps[0])->i){ 
                                ind += nsl-1 ;
                                for(j=1;j<belp->nnds;j++)
                                        *(BFnodes+ind+j) = (belp->nps[j])->i;
                        }
                        else {
                                ind += nsl ;
                                for(j=0;j<belp->nnds;j++)
                                        *(BFnodes+ind+j) = (belp->nps[j])->i ;
                        }
                        for(j=0;j<ns;j++) {
                                the_fnp = first_fnp ;
                                while( the_fnp->i != (belp->nps[j])->i)
                                        the_fnp++ ;
                                *(eqnsets[nset].Fp + ind + j) = the_fnp->f ;
                        }
        
                        put_BFKe(ns,ns,belp,varnum,eqnsets[nset].Kp,eqnsets[nset].Fp,eqnsets[nset].diagp,ind) ;
                        nsl = ns ;
                        prevnode = (belp->nps[ns-1])->i ;
                }
                belp = belp->nextbelp ;
        }
        if(code < 2)
                EndMes() ;
        return(0) ;
}

double  *KBF_init(nvar,thediagp,code)
int             nvar, code ;
long            *(*thediagp) ;

{
        int     i, nff;
        double *kp, *dp ;
        unsigned        elsize ;
        long    bigint, nelem,*diagp, il ;
        
        nelem = Nrfixed ;
        if(*thediagp != NULL)
                free(*thediagp) ;
        if((*thediagp = (long int *)calloc((int)nelem,sizeof(long))) == NULL)
                return(NULL) ;
        diagp = *thediagp ;
        *diagp = 0 ;
        diagp++ ;
        for(i=1;i<Nrfixed;i++) {
                if(i>3)
                        nff =  4 ;
                else 
                        nff = i+1 ;
                *diagp = *(diagp-1) + nff ;
                diagp++ ;
        }
        nelem = *(diagp-1) + 1 ;
        elsize = sizeof(double) ;
        bigint =  nelem * (long) elsize  ;
        kp =  dp = (double *) malloc(bigint)  ;             // NON STANDARD FUNCTION CALL 
        for(il=0;il<nelem;il++)
                *(dp++) = 0.0 ;
        return(kp) ;
}

int             put_BFKe(ns,nt,elmntp,nvar,Kp,Fp,diagp,code) 
int             ns, code ;
double  *Kp, *Fp ;
struct element *elmntp ;
long    *diagp ;

{
        int     i, j, row, col ;
        
        for(i=0;i<nt;i++) {
                row = code + i ;
                for(j=0;j<ns;j++) {
                        col = code + j ;
                        if( col >= row) {
                                *(Kp + *(diagp+col) - col + row ) += KE(i,j) ;
                        }
                }
        }
		return (0);
}

int             bnd_fluxes(nvar)
int             nvar ;

{
        int     i, tj ;
        double  q ;
        struct fxnodes *the_fnp ;
        
        the_fnp = first_fnp ;
        for(i=0;i<Nrfixed;i++) {
                q = 0.0 ;
                tj = 0 ;
                while(the_fnp->j[tj] >= 0) {
                        q += the_fnp->k[tj] * gp.iptrs[the_fnp->j[tj]]->u[nvar] ;
                        tj++ ;
                }
                q -= the_fnp->f ;
                the_fnp->f = q ;
/*
                printf(" Flux at node %d = %12.6f \n",the_fnp->i,q) ;
                if(i%10 == 0)
                        scanf_s(" %c",&anychar) ;

                the_fnp++ ;
        }
        return(0) ;
}
*/

int             mpd_i(ind,nvar)
int             ind , nvar;

{
        int     i,j,n,mpd,td,d_i() ;
        struct element  *elmntp ;
        
        mpd = 0 ;
        
        if(N.dims == 1)
                return(3) ;
        elmntp = gp.El ;
        for(i=0;i<N.elms;i++) {
                n = -1 ;
                for(j=0;j<elmntp->nnds;j++) {
                        if(elmntp->nps[j] != NULL) {
                                if(gp.iptrs[ind]->n == (elmntp->nps[j])->n )
                                        n = gp.iptrs[ind]->n ;
                        }
                }
                if( n >= 0 ) {
                        for(j=0;j<elmntp->nnds;j++) {
                                if(elmntp->nps[j] != NULL) {
                                        if((elmntp->nps[j])->n >= 0)
                                                td = d_i(ind,(elmntp->nps[j])->i) ;
                                        mpd = (mpd > td ) ? mpd : td ;
                                }
                        }
                }
                elmntp = elmntp->nextelp ;
        }
        return(mpd) ;
}

int             d_i(n1,n2)
int             n1, n2 ;

{
        return(n1 - n2) ;
}

/*        
int             fix_psi(code)
int             code;

{
        double  psi, psi1, psi2() ;
        struct belement *belp ;
        int     i, j, tovar, fromvar, snode ;
        
        if(code == 0) {
                fromvar = 0 ;
                tovar = 1 ;
        }
        else {
                fromvar = 1 ;
                tovar = 0 ;
        }
        bnd_fluxes(fromvar) ;
        belp = gp.B ;
        psi1 = (belp->nps[0])->u[tovar] ;
        snode = (belp->nps[0])->i ;
        for(i=0;i<N.belms;i++) {
                if(snode < -9999) {
                        psi1 = (belp->nps[0])->u[tovar] ;
                        snode = (belp->nps[0])->i ;
                }
                psi = psi2(belp,psi1,fromvar,tovar);
                if(belp->bcs[tovar] == 1) {
                        for(j=0;j<belp->nnds;j++)
                                (belp->nps[j])->u[tovar] = psi1 + (psi - psi1) * j / (belp->nnds - 1);
                }
                psi1 = psi ;
                if((belp->nps[belp->nnds - 1])->i == snode)
                        snode = -10000 ;
                belp = belp->nextbelp ;
        }
        return(0) ;
}

double  psi2(belp,psi1,fv,tv)
struct belement *belp ;
double  psi1 ;
int             fv, tv ;

{
        struct node             *inp ;
        struct fxnodes  *the_fnp ;
        int             i, j;
        double  psi ;
        
        psi = psi1 ;
        if(belp->bcs[fv] == 0)
                return(psi1) ;
        for(i=0;i<belp->nnds;i++) {
                inp = belp->nps[i] ;
                the_fnp = first_fnp ;
                j = 0 ;
                while( (j<Nrfixed )&& (the_fnp->i != inp->i) ) {
                        the_fnp++ ;
                        j++ ;
                }
                if(j<Nrfixed) {
                        psi += the_fnp->f ;
                        the_fnp->f = 0.0 ;
                }
        }
        return(psi) ;
}
*/

double          get_QUL2(elmntp,nvar)
struct element          *elmntp ;
int             nvar ;

{
        struct element          theEl, *theElp ;
        struct gausspts         g[MDOF] ;
        struct shapefuncs       bf ;
        struct node                     *np[MNSF] ;
        int                     i, j, nbf, ni, ntf, dir ;
        double          J[NDIMS][NDIMS], detJ, ua, x, y, z, err, L2norm, uexact(), get_J() ;
        int ElPick();

        theElp = &theEl ;
        ElPick(elmntp,theElp,&ntf,&dir) ;
        nbf = theElp->nnds ;
        
        if(N.dims == 1) 
                ni = get_gspts(911,g) ;
        if(N.dims == 2) {
                if(theElp->vtype == 229)
                        ni = get_gspts(929,g) ;
                else
                        ni = get_gspts(921,g) ;
        }
        for(j=0;j<theElp->nnds;j++)
                np[j] = theElp->nps[j];
        
        L2norm = 0.0 ;
        for(i=0;i<ni;i++) {
                get_shape(theElp->vtype,&g[i],&bf) ;
                detJ = get_J(np,&bf,J) ;
                ua = x = y = z = 0.0 ;
                for(j=0;j<nbf;j++) {
                        ua += bf.f[j] * np[j]->u[nvar] ; 
                        x += bf.f[j] * np[j]->x[0] ;
                        if(N.dims > 1) 
                                y += bf.f[j] * np[j]->x[1] ;
                        if(N.dims >2) 
                                z += bf.f[j] * np[j]->x[2] ;
                }
                err = ua - uexact(nvar,x,y,z,np[0]->p[0]) ;
                L2norm += err * err * fabs(detJ) * g[i].w ;
        }
        return(L2norm) ;
}

int    get_PGKeNumJ(elmntp,theElp,nvar,code,ntf,itnum)
struct element  *elmntp, *theElp ;
int    nvar, code, *ntf, itnum ;

{
	double	Sm[9][9], Fm[9];
	int		i,j, iNode, iVar, col;
	double	eps=0.00000001;
	int		n, nbf, nm;
	double	*fptr ;
	int		get_PGKe();

	if(nvar<0)
		n = N.vars ;
	else
		n = 1 ;
		
	*ntf = nsf(elmntp->gtype) ;
	nbf = theElp->nnds ;

	if(code > 0 || elmntp->matrices == NULL )
	{
		get_PGKe(elmntp,theElp,nvar,code,ntf,itnum);
		for(i=0;i<9;i++){
			Fm[i] = FE(i);
		}

		for(iNode=0;iNode<3;iNode++){
			for(iVar=0;iVar<3;iVar++){
				elmntp->nps[iNode]->u[iVar] += eps;
				get_PGKe(elmntp,theElp,nvar,code,ntf,itnum);
				col = iNode*3 + iVar;
				for(j=0;j<9;j++){
					Sm[j][col] = -(FE(j) - Fm[j])/eps;
				}
				elmntp->nps[iNode]->u[iVar] -= eps;
			}
		}

		for(i=0;i<9;i++){
				for(j=0;j<9;j++){
					SE(i,j) = Sm[i][j];
				}
				FE(i) = Fm[i];
		}

		if(code != 2)
		{
			if(elmntp->matrices != NULL)
				free(elmntp->matrices) ;
			nm = *ntf  * (nbf*n + 1) *n ;
			if((elmntp->matrices = (double *)calloc(nm ,sizeof(double)) ) == NULL)
				return(-1) ;
			fptr = elmntp->matrices ;
			for(i=0;i<*ntf*n;i++)
				for(j=0;j<nbf*n;j++)
				{
					*fptr = SE(i,j) ;
					fptr++ ;
				}
			for(i=0;i<*ntf*n;i++)
			{
				*fptr = FE(i);
				fptr++ ;
			}			
		}
	}
    else
	{
		fptr = elmntp->matrices ;
		for(i=0;i<*ntf*n;i++)
			for(j=0;j<nbf*n;j++)
			{
				SE(i,j) = *fptr ;
				fptr++ ;
			}
		for(i=0;i<*ntf*n;i++)
		{
			FE(i) = *fptr ;
			fptr++ ;
		}
		
    }
	return (3);
}
int             get_PGKeAnalJ(elmntp,theElp,nvar,code,ntf,itnum)
struct element          *elmntp, *theElp ;
int             nvar, code, *ntf, itnum ;

{
        int		i, j, nbf, nm, n ;
        double	*fptr ; 

		if(nvar<0)
			n = N.vars ;
        else
			n = 1 ;
		
		*ntf = nsf(elmntp->gtype) ;
		nbf = theElp->nnds ;
        
        if(code > 0 || elmntp->matrices == NULL )
		{
			get_PGKe(elmntp,theElp,nvar,code,ntf,itnum);

			if(code != 2)
			{
				if(elmntp->matrices != NULL)
					free(elmntp->matrices) ;
				nm = *ntf  * (nbf*n + 1) *n ;
				if((elmntp->matrices = (double *)calloc(nm ,sizeof(double)) ) == NULL)
					return(-1) ;
				fptr = elmntp->matrices ;
				for(i=0;i<*ntf*n;i++)
					for(j=0;j<nbf*n;j++)
					{
						*fptr = SE(i,j) ;
						fptr++ ;
					}
				for(i=0;i<*ntf*n;i++)
				{
					*fptr = FE(i);
					fptr++ ;
				}
				
			}
        }
        else
		{
			fptr = elmntp->matrices ;
			for(i=0;i<*ntf*n;i++)
				for(j=0;j<nbf*n;j++)
				{
					SE(i,j) = *fptr ;
					fptr++ ;
				}
			for(i=0;i<*ntf*n;i++)
			{
				FE(i) = *fptr ;
				fptr++ ;
			}
			
        }
		
		return(3) ;
}

int             get_PGKe(elmntp,theElp,nvar,code,ntf,itnum)
struct element          *elmntp, *theElp ;
int             nvar, code, *ntf, itnum ;

{
        struct gausspts         g[MGPTS] ;
        struct shapefuncs       bf, tf;
        struct node                     *np[MNSF] ;
        struct govfuncs	gbf, gtf ;
        int		i, j, nbf, ni, n, jj ;
        double	tdt, mtdt, sdiag ,wt;
        double	get_mlt();
		int		initm(), get_PGgf(), check_mixed(), get_mixedgps();

        if(nvar<0)
                n = N.vars ;
        else
                n = 1 ;
		
		*ntf = nsf(elmntp->gtype) ;
		theElp = elmntp ;
		nbf = theElp->nnds ;
        
        initm(ep.K,ep.n*ep.n) ;
	    if(N.trans > 0)
			initm(ep.S,ep.n*ep.n) ;
		if(N.nonlin > 0)
			initm(ep.J,ep.n*ep.n) ;
		initm(ep.F,ep.n) ;
						
		if(sets.smallH>0)
		{
			if(check_mixed(theElp)==1)
			{
				ni=get_mixedgps(theElp,g);
			}
			else
				ni = get_gspts(theElp->vtype,g) ;
		}
		else
			ni = get_gspts(theElp->vtype,g) ;
		for(j=0;j<theElp->nnds;j++)
			np[j] = theElp->nps[j];
		wt=0.0;
		for(i=0;i<ni;i++)
		{
			wt+=g[i].w;
		}
		if(itnum == 0){
			for(i=0;i<ni;i++) {
				theElp->gpcode[i] = 1;
			}
		}

		for(i=0;i<ni;i++)
		{
			gpIndex = i;
			get_shape(theElp->vtype,&g[i],&bf) ;
			get_shape(theElp->gtype,&g[i],&tf) ;
			get_PGgf(theElp,&gbf,&gtf,&bf,&tf,np,g[i].w) ;
			update_PGKe(&gbf,&gtf,nbf,*ntf,np,nvar,theElp,wt) ;
		}
		
		if(N.trans > 0)
		{ /* mass matrix lumping */
			if(lumpS == 1)
			{
				for(i=0;i<(*ntf)*n;i++)
				{
					sdiag = 0.0 ;
					for(j=0;j<(nbf)*n;j++)
					{
						sdiag += SE(i,j) ;
						SE(i,j) = 0.0 ;
					}
					SE(i,i) = sdiag ;
				}
			}
			tdt = tvals.theta; // * tvals.dt ;
			mtdt = (1.0 - tvals.theta) ; // * tvals.dt ;
			for(i=0;i<(*ntf)*n;i++)
			{
				FE(i) *= 1.0; //tvals.dt ;
				for(j=0;j<nbf;j++)
				{
					for(jj=0;jj<n;jj++)
					{
						SE(i,j*n+jj) /= tvals.dt ;
						if(itnum == 0)
							FE(i) += (SE(i,j*n+jj) - mtdt *  KE(i,j*n+jj)) * (theElp->nps[j])->uo[jj] ;
						SE(i,j*n+jj) += tdt * KE(i,j*n+jj) ;
					}
				}
			}
			if(itnum == 0){
				for(i=0;i<(*ntf)*n;i++){
					theElp->FEn[i] = FE(i);
				}
			}
			else{
				for(i=0;i<(*ntf)*n;i++){
					FE(i) = theElp->FEn[i];
				}
			}

			for(i=0;i<(*ntf)*n;i++)
			{
				for(j=0;j<nbf;j++)
				{
					for(jj=0;jj<n;jj++)
					{
						FE(i) -= (SE(i,j*n+jj) * (theElp->nps[j])->u[jj]) ;
						if(jacobianType == 0)
							SE(i,j*n+jj) += tdt * JE(i,j*n+jj) ;
					}
				}
			}
		}
		 		
	return(nbf) ;
}

double  get_mlt(gbfp,np,nbf)
int             nbf;
struct  node    *np[];
struct  govfuncs        *gbfp;

{
        double	depth;
        int		j;
        
        depth = 0.0;
        for(j=0;j<nbf;j++)
                depth += gbfp->f[j]*np[j]->u[0];
        if(depth<tvals.minH) 
                return(0.0);
        return(1.0);
}       
        
double  get_bmlt(gbfp,np,nbf)
int             nbf;
struct  node    *np[];
struct  shapefuncs      *gbfp;

{
        double	depth;
        int		j;
        
        depth = 0.0;
        for(j=0;j<nbf;j++)
                depth += gbfp->f[j]*np[j]->u[0];
        if(depth<tvals.minH) 
                return(0.0);
        return(1.0);
}       

int    get_bKeNumJ(elmntp,theElp,nvar,code,ntf,itnum)
struct belement  *elmntp, *theElp ;
int    nvar, code, *ntf, itnum ;

{
	double	Sm[6][6], Fm[6];
	int		i,j, iNode, iVar, col;
	double	eps=0.00000001;
	int		n, nbf, nm, get_bKe();
	double	*fptr ;
	
	if(nvar < 0)
			n = N.vars ;
		else
			n = 1 ;
		
	*ntf = nsf(elmntp->gtype) ;
	nbf =  nsf(elmntp->vtype) ;

	if(code > 0 || elmntp->matrices == NULL )
	{
		get_bKe(elmntp,theElp,nvar,code,ntf,itnum);
		for(i=0;i<6;i++){
			Fm[i] = FE(i);
		}

		for(iNode=0;iNode<2;iNode++){
			for(iVar=0;iVar<3;iVar++){
				elmntp->nps[iNode]->u[iVar] += eps;
				get_bKe(elmntp,theElp,nvar,code,ntf,itnum);
				col = iNode*3 + iVar;
				for(j=0;j<6;j++){
					Sm[j][col] = -(FE(j) - Fm[j])/eps;
				}
				elmntp->nps[iNode]->u[iVar] -= eps;
			}
		}

		for(i=0;i<6;i++){
			for(j=0;j<6;j++){
				SE(i,j) = Sm[i][j];
			}
			FE(i) = Fm[i];
		}

		if(code != 2)
		{
			if(elmntp->matrices != NULL)
				free(elmntp->matrices) ;
			nm = *ntf  * (nbf*n + 1) *n ;
			if((elmntp->matrices = (double *)calloc(nm ,sizeof(double)) ) == NULL)
				return(-1) ;
			fptr = elmntp->matrices ;
			for(i=0;i<*ntf*n;i++)
				for(j=0;j<nbf*n;j++)
				{
					*fptr = SE(i,j) ;
					fptr++ ;
				}
			for(i=0;i<*ntf*n;i++)
			{
				*fptr = FE(i);
				fptr++ ;
			}			
		}
	}
	else
	{
		fptr = elmntp->matrices ;
		for(i=0;i<*ntf*n;i++)
			for(j=0;j<nbf*n;j++)
			{
				SE(i,j) = *fptr ;
				fptr++ ;
			}
		for(i=0;i<*ntf*n;i++)
		{
			FE(i) = *fptr ;
			fptr++ ;
		}
	}
	return (nbf);
}

int             get_bKeAnalJ(elmntp,theElp,nvar,code,ntf,itnum)
struct belement         *elmntp, *theElp ;
int             nvar, code, *ntf, itnum ;

{
        
        int		i, j, nbf, nm, n ;
        double	*fptr ; 
        
        if(nvar < 0)
			n = N.vars ;
		else
			n = 1 ;
		
        *ntf = nsf(elmntp->gtype) ;  
		nbf =  nsf(elmntp->vtype) ;

        if(code > 0 )
		{
			get_bKe(elmntp,theElp,nvar,code,ntf,itnum);

			if(code != 2 )
			{
				if(elmntp->matrices != NULL)
					free(elmntp->matrices) ;
				nm = *ntf  * (nbf*n + 1) *n ;
				if((elmntp->matrices = (double *)calloc(nm ,sizeof(double)) ) == NULL)
					return(-1) ;
				fptr = elmntp->matrices ;
				for(i=0;i<*ntf*n;i++)
					for(j=0;j<nbf*n;j++)
					{
						*fptr = SE(i,j) ;
						fptr++ ;
					}
				for(i=0;i<*ntf*n;i++)
				{
					*fptr = FE(i);
					fptr++ ;
				}
				
			}
		}
		else
		{
			fptr = elmntp->matrices ;
			for(i=0;i<*ntf*n;i++)
				for(j=0;j<nbf*n;j++)
				{
					SE(i,j) = *fptr ;
					fptr++ ;
				}
			for(i=0;i<*ntf*n;i++)
			{
				FE(i) = *fptr ;
				fptr++ ;
			}
			
        }

		return(nbf) ;
}


int             get_bKe(elmntp,theElp,nvar,code,ntf,itnum)
struct belement         *elmntp, *theElp ;
int             nvar, code, *ntf, itnum ;

{
        
        struct gausspts         g[MGPTS] ;
        struct shapefuncs       fv, ft, *fgp ;
        struct node	*np[MNSF] ;
        int		i, j, ig, nbf, ni, n, jj ;
        double	Jaa, Jab, dS, tdt, mtdt,dx,dy ;
        double	mlt,get_bmlt();
		int		check_bmixed(), get_boundarygps();
        
        if(nvar < 0)
			n = N.vars ;
		else
			n = 1 ;
		nbf =  nsf(elmntp->vtype) ;
        *ntf = nsf(elmntp->gtype) ;       
        theElp = elmntp ;

        initm(ep.K,ep.n*ep.n) ;
		if(N.trans > 0)
			initm(ep.S,ep.n*ep.n) ;
		if(N.nonlin > 0)
			initm(ep.J,ep.n*ep.n) ;
		initm(ep.F,ep.n) ;
		
		if(sets.smallH>0)
		{
			if(check_bmixed(elmntp)==1)
				ni=get_boundarygps(elmntp,g);
			else
				ni = get_gspts(elmntp->vtype,g) ;
		}
		else
			ni = get_gspts(elmntp->vtype,g) ;
		for(j=0;j<elmntp->nnds;j++)
			np[j] = elmntp->nps[j] ;
		mlt=1.0;
		for(ig=0;ig<ni;ig++)
		{
			get_shape(elmntp->vtype,&g[ig],&fv) ;
			get_shape(elmntp->gtype,&g[ig],&ft) ;
			Jaa = Jab = 0.0 ;
			fgp = &fv ;
			for(i=0;i<fgp->dof;i++)
			{
				Jaa += fgp->dfdr[i] * np[i]->x[0] ;
				Jab += fgp->dfdr[i] * np[i]->x[1] ;
			}
			dS = g[ig].w * sqrt(Jaa * Jaa + Jab * Jab) ;
			dx = g[ig].w * Jaa;
			dy = g[ig].w * Jab;
			update_BKe(elmntp,&fv,&ft,nbf,*ntf,np,nvar,dx,dy) ;
		}
		if(N.trans > 0)
		{
			tdt = tvals.theta; // * tvals.dt ;
			mtdt = (1.0 - tvals.theta); // * tvals.dt ;
			for(i=0;i<*(ntf)*n;i++)
			{
				FE(i) *= 1.0; //tvals.dt ;
				for(j=0;j<nbf;j++)
				{
					for(jj=0;jj<n;jj++)
					{
						SE(i,j*n+jj) /= tvals.dt ;
						if(itnum == 0)
							FE(i) += (SE(i,j*n+jj) - mtdt *  KE(i,j*n+jj)) * (theElp->nps[j])->uo[jj] ;
						SE(i,j*n+jj) += tdt * KE(i,j*n+jj) ;
					}
				}
			}
		}

		if(itnum == 0){
			for(i=0;i<(*ntf)*n;i++){
				elmntp->FEn[i] = FE(i);
			}
		}
		else{
			for(i=0;i<(*ntf)*n;i++){
				FE(i) = elmntp->FEn[i];
			}
		}

		for(i=0;i<(*ntf)*n;i++)
		{
			for(j=0;j<nbf;j++)
			{
				for(jj=0;jj<n;jj++)
				{
					FE(i) -= (SE(i,j*n+jj) * (theElp->nps[j])->u[jj]) ;
					if(jacobianType == 0)
						SE(i,j*n+jj) += tdt * JE(i,j*n+jj);
				}
			}
		}

		return(nbf) ;
}

int             ElPick(elmntp,theElp,ntfp,dirp)
struct element  *elmntp, *theElp ;
int             *ntfp, *dirp ;

{
        *dirp = 1 ;
        *ntfp = nsf(elmntp->gtype) ;
        *theElp = *elmntp ;
        return(0) ;
}

int             get_PGgf(elp,gbfp,gtfp,bfp,tfp,np,w)
struct element          *elp ;
struct shapefuncs       *bfp, *tfp ;
struct node                     *np[] ;
struct govfuncs         *gbfp,*gtfp ;
double                          w ;

{
        double  J[NDIMS][NDIMS], get_J() ;
        int             i, j ;
        double  set_BUw();
		int		invertJ();

        gbfp->dof = bfp->dof ;
        gtfp->dof = tfp->dof ;
        
        for(i=0;i<N.params;i++) {
                gbfp->p[i] = elp->p[i] ;
        }
                        i = (gbfp->dof > gtfp->dof) ? gbfp->dof : gtfp->dof ;
                        for(j=0;j<i;j++) {
                                gbfp->f[j] = 0.0  ;
                                gtfp->f[j] = 0.0  ;
                                gbfp->dfdx[j] = 0.0 ;     
                                gtfp->dfdx[j] = 0.0 ;
                                gbfp->dfdy[j] = 0.0 ;
                                gtfp->dfdy[j] = 0.0 ;
                        }
        for(j=0;j<bfp->dof;j++)
                gbfp->f[j] = bfp->f[j]  ;
        for(j=0;j<tfp->dof;j++)
                gtfp->f[j] = tfp->f[j]  ;
        
        gbfp->detJ = get_J(np,bfp,J) ;
                if (gbfp->detJ <= 0.0){
                        printf("Bad element # %d, J = %f\n",elp->n,gbfp->detJ);
                }
        gbfp->w = w ;
        invertJ(J,gbfp->detJ) ;
        
        for(j=0;j<bfp->dof;j++) {
        
                switch (N.dims) {
                
                        case 1 :
                                gbfp->dfdx[j] = J[0][0] * bfp->dfdr[j] ;
                                gbfp->d2fdx[j] = J[0][0] * J[0][0] * bfp->d2fdr[j] ;
                                gbfp->dfdy[j] = 0.0 ;
                                gbfp->dfdz[j] = 0.0 ;
                                break ;
                                
                        case 2 :
                                gbfp->dfdx[j] = J[0][0] * bfp->dfdr[j] + J[0][1] * bfp->dfds[j] ;
                                gbfp->dfdy[j] = J[1][0] * bfp->dfdr[j] + J[1][1] * bfp->dfds[j] ;
                                gbfp->dfdz[j] = 0.0 ;
                                break ;
                                
                        case 3 :
                                gbfp->dfdx[j] = J[0][0] * bfp->dfdr[j] + J[0][1] * bfp->dfds[j] + J[0][2] *bfp->dfdt[j] ;
                                gbfp->dfdy[j] = J[1][0] * bfp->dfdr[j] + J[1][1] * bfp->dfds[j] + J[1][2] *bfp->dfdt[j] ;
                                gbfp->dfdz[j] = J[2][0] * bfp->dfdr[j] + J[2][1] * bfp->dfds[j] + J[2][2] *bfp->dfdt[j] ;
                                break ;
                }
        }
                
        for(j=0;j<tfp->dof;j++) {
        
                switch (N.dims) {
                
                        case 1 :
                                gtfp->dfdx[j] = J[0][0] * tfp->dfdr[j] ;
                                gtfp->dfdy[j] = 0.0 ;
                                gtfp->dfdz[j] = 0.0 ;
                                break ;
                                
                        case 2 :
                                gtfp->dfdx[j] = J[0][0] * tfp->dfdr[j]*set_BUw(elp,j,0) + J[0][1] * tfp->dfds[j] *set_BUw(elp,j,1) ;
                                gtfp->dfdy[j] = J[1][0] * tfp->dfdr[j] *set_BUw(elp,j,0) + J[1][1] * tfp->dfds[j]*set_BUw(elp,j,1)  ;
                                gtfp->dfdz[j] = 0.0 ;
                                break ;
                                
                        case 3 :
                                gtfp->dfdx[j] = J[0][0] * tfp->dfdr[j] + J[0][1] * tfp->dfds[j] + J[0][2] *tfp->dfdt[j] ;
                                gtfp->dfdy[j] = J[1][0] * tfp->dfdr[j] + J[1][1] * tfp->dfds[j] + J[1][2] *tfp->dfdt[j] ;
                                gtfp->dfdz[j] = J[2][0] * tfp->dfdr[j] + J[2][1] * tfp->dfds[j] + J[2][2] *tfp->dfdt[j] ;
                                break ;
                }
        }
        gbfp->detJ *= gbfp->w ;
        return(0) ;
}


int             get_QU1dKe(elmntp,theElp,nvar,code,ntf)
struct element          *elmntp, *theElp ;
int             nvar, code, *ntf ;

{
        
        struct gausspts         g[MGPTS] ;
        struct shapefuncs       bf, tf, *fgp;
        struct node                     *np[MNSF] ;
        struct govfuncs         gbf, gtf, *gbfp, *gtfp ;
        int                                     i, j, ig, nbf, ni, dir ;
        double                          Jaa, Jab, dS, ubar, tdt ;
        double                          BoundK(), BoundF() ;
		int		QU1DPick();
        
        initm(ep.K,ep.n*ep.n) ;
        if(N.trans > 0)
                initm(ep.S,ep.n*ep.n) ;
        initm(ep.F,ep.n) ;
        
        QU1DPick(elmntp,theElp,ntf,&dir) ;
        nbf = theElp->nnds ;
        gbfp = &gbf ;
        gtfp = &gtf ;
        
        if(code > -1 ) {
                ni = get_gspts(theElp->gtype,g) ;
                for(j=0;j<theElp->nnds;j++)
                        np[j] = theElp->nps[j] ;
        
                for(ig=0;ig<ni;ig++) {
                        get_shape(theElp->vtype,&g[ig],&bf) ;
                        get_shape(theElp->gtype,&g[ig],&tf) ;
                        fgp = &bf ;
                        for(i=0;i<N.params;i++) {
                                gbfp->p[i] = elmntp->p[i] ;
                                for(j=0;j<fgp->dof;j++)
                                        gbfp->p[i] += fgp->f[j] * np[j]->p[i] ;
                        }
                        Jaa = Jab = 0.0 ;
                        for(i=0;i<fgp->dof;i++) {       
                                Jaa += fgp->dfdr[i] * np[i]->x[0] ;
                                Jab += fgp->dfdr[i] * np[i]->x[1] ;
                        }
                        dS = sqrt(Jaa * Jaa + Jab * Jab) ;
                        ubar = ((gbfp->p[0] * Jaa) + (gbfp->p[1] * Jab))/dS ;
                        for(i=0;i<*ntf;i++) {
                                for(j=0;j<nbf;j++) {
                                        KE(i,j) += (gbfp->p[2] * tf.dfdr[i] * bf.dfdr[j] / dS
                                                + ubar * tf.f[i] * fgp->dfdr[j]
                                        /*      - gbfp->p[3]*tf.f[i]*bf.f[j]*dS  */  )*g[ig].w ;
                                        SE(i,j) +=      ( tf.f[i] * bf.f[j] * dS) * g[ig].w ;
                                }
                        }
                }
        }
        if(N.trans > 0) {
                tdt = tvals.theta * tvals.dt ;
                for(i=0;i<*ntf;i++) {
                        for(j=0;j<nbf;j++) {
                                SE(i,j) += tdt * KE(i,j) ;
                        }
                }
        }
        return(nbf) ;
}

double dotm(A,B,l,m,n)


int     l,m,n;
double  A[3][3], B[3][3];

{
        register double  t = 0.0;
        int     i;

        for(i=0;i<n;i++)

                t += A[l][i] * B[i][m];

        return(t);
}

double  set_BUw(elp,j,i)
struct  element         *elp;
int             i,j;

{
        double  mult,r1,r2,s1,s2;
        int             nbn,l1,l2;
		int		nbounds(), get_nlcs();
        
        nbn = nbounds(elp->gtype);
        mult = 1.0;
        
        if((elp->nps[j]->fxc)==1) {
                if(j==0&&(elp->nps[1]->fxc == 1)) {
                        l1 = 0;
                        l2 = 1;
                        get_nlcs(elp->gtype,l1,&r1,&s1);
                        get_nlcs(elp->gtype,l2,&r2,&s2);
                        if(r1==r2) {
                                if(i==0) mult = 0.0;
                        }
                        if(s1==s2) {
                                if(i==1) mult = 0.0;
                        }

                }
                
                if(j==0&&(elp->nps[nbn-1]->fxc) == 1) {
                        l1 = 0;
                        l2 = nbn-1;
                        get_nlcs(elp->gtype,l1,&r1,&s1);
                        get_nlcs(elp->gtype,l2,&r2,&s2);
                        if(r1==r2) {
                                if(i==0) mult = 0.0;
                        }
                        if(s1==s2) {
                                if(i==1) mult = 0.0;
                        }

                }
                
                if(j==(nbn-1)&&elp->nps[nbn-2]->fxc==1) {
                        l1 = nbn-1;
                        l2 = nbn-2;
                        get_nlcs(elp->gtype,l1,&r1,&s1);
                        get_nlcs(elp->gtype,l2,&r2,&s2);
                        if(r1==r2) {
                                if(i==0) mult = 0.0;
                        }
                        if(s1==s2) {
                                if(i==1) mult = 0.0;
                        }

                }
                
                if(j==(nbn-1)&&elp->nps[0]->fxc==1) {
                        l1 = nbn-1;
                        l2 = 0;
                        get_nlcs(elp->gtype,l1,&r1,&s1);
                        get_nlcs(elp->gtype,l2,&r2,&s2);
                        if(r1==r2) {
                                if(i==0) mult = 0.0;
                        }
                        if(s1==s2) {
                                if(i==1) mult = 0.0;
                        }

                }

                if(j!=0&&j!=(nbn-1)) {
                        if(elp->nps[j+1]->fxc==1) {
                                l1 = j;
                                l2 = j+1;
                                get_nlcs(elp->gtype,l1,&r1,&s1);
                                get_nlcs(elp->gtype,l2,&r2,&s2);
                                if(r1==r2) {
                                        if(i==0) mult = 0.0;
                                }
                                if(s1==s2) {
                                        if(i==1) mult = 0.0;
                                }

                        }
                        if( elp->nps[j-1]->fxc==1) {
                                l1 = j;
                                l2 = j-1;
                                get_nlcs(elp->gtype,l1,&r1,&s1);
                                get_nlcs(elp->gtype,l2,&r2,&s2);
                                if(r1==r2) {
                                        if(i==0) mult = 0.0;
                                }
                                if(s1==s2) {
                                        if(i==1) mult = 0.0;
                                }

                        }
                }
        }
        return(mult);
}       
        
/* Ashraf       Jan 15, 1994 */
double  update_GWKe(f,v,nbf,ntf,np,vcode,theElp,wt,Sox,Soy ) 
int             nbf, ntf, vcode  ;
struct  element *theElp ;
struct  node    *np[] ;
struct  govfuncs        *f, *v ;
double  Sox,Soy,wt;

{
        int     i, j, n ;
        double  T=1.0;

                        
        if(vcode<0)
                n = N.vars ;
        else
                n = 1 ;

/*      
                v->dfdx[i] = UW * f->dfdx[i] * sqrt(f->detJ/f->w*wt) ;
                v->dfdy[i] = UW * f->dfdy[i] * sqrt(f->detJ/f->w*wt) ;
*/

/*  Compute Parameters Here */
        
                
        for(i=0;i<ntf;i++) {
                for(j=0;j<nbf;j++) {

                        KE(i*n,j*n+1) -= f->dfdx[i]*f->f[j]*f->detJ;
                        KE(i*n,j*n+2) -= f->dfdy[i]*f->f[j]*f->detJ;
                        KE(i*n+1,j*n) -= f->dfdx[i]*f->f[j]*f->detJ;
                        KE(i*n+1,j*n+1) += f->f[i]*f->f[j]/T*f->detJ;
                        KE(i*n+2,j*n) -= f->dfdy[i]*f->f[j]*f->detJ;
                        KE(i*n+2,j*n+2) += f->f[i]*f->f[j]/T*f->detJ;


                }
        }


        for(i=0;i<ntf;i++) {
                for(j=0;j<nbf;j++) {
        
                        SE(i*n,j*n) += f->f[i]*f->f[j]*f->detJ;

                }
        }

        for(i=0;i<ntf;i++) {
        
                                FE(i*n+1) += f->f[i]*Sox*f->detJ;
                                FE(i*n+2) += f->f[i]*Soy*f->detJ;
                
        }

        
        return(0.0) ;
}
        
        

double  update_GWBKe(belp,f,v,nbf,ntf,np,nvar,dx,dy) 
int                     nbf, ntf ;
double                  dx,dy ;
struct node             *np[] ;
struct belement         *belp ;
struct shapefuncs       *f, *v ;

{
        int             i, j, n,l,m;
        double  Qx, Qy, H, Qn, Qt, dS;
        double  AK[3][3],cbcc[3],bcc[3],phi[3];

        if(nvar<0)
                n = N.vars ;
        else
                n = 1 ;

        for(i=0;i<N.vars;i++) {
                cbcc[i] = 0.0;
                bcc[i] = 0.0;
                phi[i] = 0.0;
        }

        for(i=0;i<N.vars;i++) {
                for(j=0;j<N.vars;j++) 
                                AK[i][j] = 0.0;
        }

        
        dS = sqrt( dx*dx + dy*dy );

        H = belp->p[0];
        Qn = belp->p[1];
        Qt = belp->p[2];

        Qx = (-Qn*dy-Qt*dx)/dS;
        Qy = ( Qn*dx-Qt*dy)/dS;


        phi[0]  = H;
        phi[1]  = Qx;
        phi[2]  = Qy;


        switch(belp->bcs[0]) {

        case 0:
                bcc[0] = 0.0;
                bcc[1] = 1.0;
                bcc[2] = 1.0;
                break;

        case 3:
                bcc[0] = 1.;
                bcc[1] = 0.;
                bcc[2] = 0.;
                break;

        }


        for(i=0;i<n;i++)
                cbcc[i] = fabs(bcc[i]-1.);
                


/* AK are the original boundary terms */


        AK[0][0] = 0.0;
        AK[0][1] = dy;
        AK[0][2] = -dx;

        AK[1][0] = dy;          
        AK[1][1] = 0.0;
        AK[1][2] = 0.0;

        AK[2][0] = -dx;
        AK[2][1] = 0.0;
        AK[2][2] = 0.0;
        

        for(i=0;i<ntf;i++) {
                for(j=0;j<nbf;j++) {

                        if (tvals.theta<0.005) {
/*                              SE(i*n,j*n) += elmntp->p[6] * (1.E18 * elmntp->p[0]) ;
                                SE(i*n,j*n+1) += elmntp->p[6] * (1.E18 * elmntp->p[1]) ;
                                SE(i*n+1,j*n) += elmntp->p[6] * (1.E18 * elmntp->p[3]) ;  
                                SE(i*n+1,j*n+1) += elmntp->p[6] * (1.E18 * elmntp->p[4]) ;*/
                        }               
                        else {
                                for(l=0;l<n;l++) {
                                        for(m=0;m<n;m++) {
                                                KE(i*n+l,j*n+m) += AK[l][m] * cbcc[m] * f->f[i] * f->f[j];
                                                FE(i*n+l) -= AK[l][m] * bcc[m] * phi[m]*f->f[i]*f->f[j];    
                                        }
                                }                                       
                        }
                }

        }
        return(0.0);
}       
        
/* Ashraf       Jan 15, 1994 */
        
/* Ashraf Nov 29, 1994 */

void trans_bcs(time)
double  time;
{
        double  depth, omega, period=12.4*3600, depthav=19.05;
        struct belement *belp;
        struct boundaryseg      *bseg;
        int i;
        
        belp=gp.B;
        bseg=gp.BSEG;
        
        for(i=0;i<N.boundarysegs;i++) {
                omega = bseg->p[1];
                depth = bseg->p[0]*cos(omega*time) + depthav;
                while(belp->nps[0]->n!=bseg->nstart)
                        belp=belp->nextbelp;
                while(belp->nps[belp->nnds-1]->n!=bseg->nend){
                        belp->p[0]=depth;
                        belp=belp->nextbelp;
                }
                belp->p[0]=depth;
                bseg=bseg->nextbseg;

        }
}

//	April, 2003  -  P. Steffler
//	Calculates nodal velocities on the basis of a projection of the 
//	discharge intensity and depth distribution as an alternative
//	to velocities based simply on nodal qx, qy, and H.
//	Also performs a projection on depth which has a smoothing effect.

void updateVelocities()
{
	double *M, *FH, *FU, *FV;				// global matrices (vectors actually)
	double ME[3][3], FHE[3], FUE[3], FVE[3];// element matrices
	int iElm, iNode, iTest, iBasis, iRow;	// counters
	int iGP, nGP;							// gauss point counter, number of gauss points
	struct gausspts	g[MGPTS];				// locations and weights for element gauss points
	struct shapefuncs	f;					// element basis functions
	double	qx, qy, H, tice, depth;			// solution variables evaluated at the gauss point
	struct element	*elmP;					// current element pointer
	struct node		*nodeP;					// current node pointer

	// initialize global vectors
	M = calloc(N.nodes,sizeof(double));
	FH = calloc(N.nodes,sizeof(double));
	FU = calloc(N.nodes,sizeof(double));
	FV = calloc(N.nodes,sizeof(double));

	elmP = gp.El;

	// cycle through all the elements
	for(iElm=0;iElm<N.elms;iElm++){

		// initialize element matrices
		for(iTest=0;iTest<3;iTest++){
			for(iBasis=0;iBasis<3;iBasis++){
				ME[iTest][iBasis] = 0.0;
			}
			FUE[iTest] = 0.0;
			FVE[iTest] = 0.0;
			FHE[iTest] = 0.0;
		}

		//build element matrices

		// get the locations and weights of the gauss integration points (fegauss.c)
		if(check_mixed(elmP)==1){
			nGP = get_mixedgps(elmP,g);
		}
		else
			nGP = get_gspts(elmP->vtype,g);

		// cycle through gauss points
		for(iGP=0;iGP<nGP;iGP++){
			//interpolate nodal variable values to gauss point using element shape functions
			get_shape(elmP->vtype,&g[iGP],&f) ;
			H = tice = qx = qy = 0;
	        for(iBasis=0;iBasis<3;iBasis++) {
				nodeP = elmP->nps[iBasis];
                H += f.f[iBasis]*nodeP->u[0] ;
				tice += f.f[iBasis]*nodeP->ice[0] ;
                qx += f.f[iBasis]*nodeP->u[1] ;
                qy += f.f[iBasis]*nodeP->u[2] ;
		    }
			tice = tvals.sg_ice * tice;

			// calculate element matrix contributions for this gauss point
			for(iTest=0;iTest<3;iTest++){
				for(iBasis=0;iBasis<3;iBasis++){
					ME[iTest][iBasis] += f.f[iTest] * f.f[iBasis] * g[iGP].w;
				}
				depth = H - tice;
				if(depth > tvals.gwH){
					FHE[iTest] += f.f[iTest] * depth * g[iGP].w;
					FUE[iTest] += f.f[iTest] * qx/depth * g[iGP].w;
					FVE[iTest] += f.f[iTest] * qy/depth * g[iGP].w;
				}
			}
		}

//		assemble element matrices into global vectors, lumping the ME matrix into a global M vector
		for(iTest=0;iTest<3;iTest++){
			iRow = elmP->nps[iTest]->i;
			for(iBasis=0;iBasis<3;iBasis++){
				*(M + iRow) += ME[iTest][iBasis];
			}
			*(FH + iRow) += FHE[iTest];
			*(FU + iRow) += FUE[iTest];
			*(FV + iRow) += FVE[iTest];
		}

		elmP = elmP->nextelp;
	}

//	update nodal velocities by solving diagonal systems [M]{U} = {FU}
	nodeP = gp.N;
	for(iNode=0;iNode<N.nodes;iNode++){
		nodeP->ud[0] = *(FH + iNode) / *(M + iNode);
		nodeP->ud[1] = *(FU + iNode) / *(M + iNode);
		nodeP->ud[2] = *(FV + iNode) / *(M + iNode);
		nodeP = nodeP->nextnp;
	}

//	clean up memory
	free(M);
	free(FH);
	free(FU);
	free(FV);

	return;
}

void setGWFtoZero(){

	int i;
	struct node	*nodeP;
	nodeP = gp.N;
	for(i=0;i<N.nodes;i++)
	{
		if ((nodeP->u[0]) < 0.0)
		{
			nodeP->u[1]=0.0;
			nodeP->u[2]=0.0;

		}				
		nodeP=nodeP->nextnp;
	}
}

double avgCelerity(){

	int i,noOfWetNodes=0;
	struct node	*nodeP;
	double avgCel=0.0;
	
	nodeP = gp.N;
	for(i=0;i<N.nodes;i++)
	{
		if ((nodeP->u[0]-(nodeP->ice[0]*tvals.sg_ice)) > tvals.gwH)
		{
			noOfWetNodes++;
			avgCel += sqrt(nodeP->u[1]*nodeP->u[1]+nodeP->u[2]*nodeP->u[2])/fabs(nodeP->u[0]-(nodeP->ice[0]*tvals.sg_ice))+ sqrt(GRAV*fabs(nodeP->u[0]-(nodeP->ice[0]*tvals.sg_ice)));

		}	
		nodeP=nodeP->nextnp;
	}
	avgCel=avgCel/noOfWetNodes;
	return avgCel;

}

void findMaxchange(){
	
	double *dp;
	int i,maxVarInaNode;
	double chgInH, chgInqx, chgInqy,maxChgInaNode;
	struct node	*nodeP;

	dp = eqnsets[0].Fp ;

	nodeP = gp.N;

	chgInH=0.0;
	chgInqx=0.0;
	chgInqy=0.0;
	maxChgInaNode=0.0;
	maxChange=0.0;
	maxVarInaNode=0;

	for(i=0;i<N.nodes;i++) 
	{
		chgInH=*dp;
		dp++;
		chgInqx=*dp;
		dp++;
		chgInqy=*dp;
		dp++;

		if((nodeP->u[0]-(nodeP->ice[0]*tvals.sg_ice)) < tvals.gwH || (nodeP->u[0]-(nodeP->ice[0]*tvals.sg_ice)-chgInH) < tvals.gwH)
		{
			if(fabs(chgInH) > 1.0)
				maxChgInaNode = chgInH;
			else
				maxChgInaNode = 0.0; //ignoring all the changes in groundwater portion of the domain

			maxVarInaNode=0;
		}
		else
		{
			if(fabs(chgInH) > fabs(chgInqx) && fabs(chgInH) > fabs(chgInqy)){
				maxChgInaNode = chgInH;
				maxVarInaNode=0;
			}
			else if(fabs(chgInqx)>fabs(chgInqy)){
				maxChgInaNode=chgInqx;
				maxVarInaNode=1;
			}
			else{
				maxChgInaNode=chgInqy;
				maxVarInaNode=2;
			}
		}

		if(fabs(maxChgInaNode)>=fabs(maxChange) )
		{
			maxChange = maxChgInaNode;
			maxnodeP = nodeP;
		}

		nodeP = nodeP->nextnp;
	}	
}

void useHalfdp(nset)
int nset;
{

	double	*dp ;
	int		i,j;

	dp = eqnsets[nset].Fp ;
			
	for(i=0;i<N.nodes;i++) {
		for(j=0;j<N.vars;j++) {
			gp.iptrs[i]->u[j] -= *dp ;
			gp.iptrs[i]->u[j] += 0.5*(*dp);
			dp++;
		}
	}
}
void findMaxchangetemp(RHS)
double	*RHS;
{
	
	double *dp;
	int i,maxVarInaNode;
	double chgInH, chgInqx, chgInqy,maxChgInaNode;
	struct node	*nodeP;

	dp = RHS ;

	nodeP = gp.N;

	chgInH=0.0;
	chgInqx=0.0;
	chgInqy=0.0;
	maxChgInaNode=0.0;
	maxChange=0.0;
	maxVarInaNode=0;

	for(i=0;i<N.nodes;i++) 
	{
		chgInH=*dp;
		dp++;
		chgInqx=*dp;
		dp++;
		chgInqy=*dp;
		dp++;

		if((nodeP->u[0]-(nodeP->ice[0]*tvals.sg_ice)) < tvals.gwH || (nodeP->u[0]-(nodeP->ice[0]*tvals.sg_ice)-chgInH) < tvals.gwH)
		{
			if(fabs(chgInH) > 1.0)
				maxChgInaNode = chgInH;
			else
				maxChgInaNode = 0.0; //ignoring all the changes in groundwater portion of the domain

			maxVarInaNode=0;
		}
		else
		{
			if(fabs(chgInH) > fabs(chgInqx) && fabs(chgInH) > fabs(chgInqy)){
				maxChgInaNode = chgInH;
				maxVarInaNode=0;
			}
			else if(fabs(chgInqx)>fabs(chgInqy)){
				maxChgInaNode=chgInqx;
				maxVarInaNode=1;
			}
			else{
				maxChgInaNode=chgInqy;
				maxVarInaNode=2;
			}
		}

		if(fabs(maxChgInaNode)>=fabs(maxChange) )
		{
			maxChange = maxChgInaNode;
			maxnodeP = nodeP;
		}

		nodeP = nodeP->nextnp;
	}	
}
void useHalfdptemp(RHS)
double *RHS;
{

	double	*dp ;
	int		i,j;

	dp =RHS ;
			
	for(i=0;i<N.nodes;i++) {
		for(j=0;j<N.vars;j++) {
			gp.iptrs[i]->u[j] -= *dp ;
			gp.iptrs[i]->u[j] += 0.5*(*dp);
			dp++;
		}
	}
}



