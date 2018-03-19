//#pragma options cpluscmt	
#include "Fe_PS.h"

#define	KE(A,B) ( *(ep.K + ep.n*(A) + (B)) )
#define	SE(A,B) ( *(ep.S + ep.n*(A) + (B)) )
#define	ME(A,B) ( *(ep.M + ep.n*(A) + (B)) )
#define	JE(A,B) ( *(ep.J + ep.n*(A) + (B)) )
#define	FE(B) ( *(ep.F + B) )

extern struct control  N ;
extern struct pointers gp;
extern struct RegMesh		Mesh ;
extern struct fxnodes	*first_fnp ;
extern int		Nukns ;
extern struct elmpointers ep ;
extern struct transient		tvals ;
extern struct eqnset	eqnsets[4] ; 
extern double			uchange, maxChange ;
extern int		maxNode, maxVar;
extern double	defndp , defelp ;

int		varindex(ndindex, var)
int	ndindex, var ;

{
	int	i, vi ;
	
	vi = gp.iptrs[ndindex]->ui ;
	for(i=0;i<var;i++) {
		if(isvarfixed(ndindex,i) == 0)
			vi++ ;
	}
	return(vi) ;
}

int		isvarfixed(ndindex,var)
int	ndindex, var ;

{
	int	fixcode, rc ;
	
	rc = 0 ;
	fixcode = (1 << var) ;
	if( (fixcode & gp.iptrs[ndindex]->fxc) != 0)
		rc = -1 ;
	return(rc) ;
}

int		one_step(nvar)
int		nvar ;

{
	int		i ;
	
	for(i=0;i<tvals.nsteps;i++) {
		tvals.t += tvals.dt ;
		assemble(nvar,0,0,0) ;
		solve(nvar,1,0) ;
		tvals.dt = tvals.dt * tvals.dtfac ;
	}
	return(0) ;
}
	
int		solve(nvar,code,nset)	//for River2D nvar = -1, code = 2, nset = 0
int		code, nvar, nset ;

{
	int	i, j, ind, solvecode ;
	double	*dp, varsum ;
	
	if(code > 1)
		solvecode = 1 ;
	else
		solvecode = code ;
	if(code < 0) {
		solvecode = 2 ;
		code = - code ;
	}
/*
	printf(" Now solving ...") ;
*/
	if((N.sym ==  0)  ||  (code == 3))
		actcol(eqnsets[nset].Kp,eqnsets[nset].Fp,eqnsets[nset].diagp, eqnsets[nset].neqns, solvecode) ;
	else
		uactcl(eqnsets[nset].Kp,eqnsets[nset].Lp,eqnsets[nset].Fp,eqnsets[nset].diagp, eqnsets[nset].neqns, solvecode) ;
/*
	printf(" All done !\n") ;
*/
	
	dp = eqnsets[nset].Fp ;
	if((code > 0) && (code != 3)) {
		if(nvar >= 0)	//for River2D nvar is always -1 so this if is never called
		{				
			ind = 0 ;
			uchange = 0.0 ;
			for(i=0;i<N.nodes;i++) 
			{
				if(code == 1) {
					uchange += fabs(gp.iptrs[i]->u[nvar] - *dp) ;
					gp.iptrs[i]->u[nvar] = *dp ;
				}
				if(code == 2) {
					*(eqnsets[nset+1].Fp + ind) = *(eqnsets[nset].Fp + gp.iptrs[i]->i) ;
				/*
					gp.iptrs[i]->p[4] = *(eqnsets[nset].Fp + gp.iptrs[i]->i) ;
				*/
				}
				if(code == 4) {
				/*
					printf(" %f  %f\n",gp.iptrs[i]->u[nvar], *dp);
				*/
					gp.iptrs[i]->u[nvar] += *dp ;
				}
				dp++;
				ind++ ;
			}
		} 				//end of if(nvar >= 0)
		else {
			uchange = 0.0 ;
			maxChange = 0.0;
			varsum = 0.0 ;
			for(i=0;i<N.nodes;i++) {
				for(j=0;j<N.vars;j++) 
				{
					if(code == 1)
						gp.iptrs[i]->u[j] = *dp;
					if(code == 2)	//in River2D code is always 2 
					{
						uchange += *dp * *dp;
						if(fabs(*dp) > fabs(maxChange)){
							maxChange = *dp;
							maxNode = i;
							maxVar = j;
						}
						gp.iptrs[i]->u[j] += *dp ;
						varsum += gp.iptrs[i]->u[j] * gp.iptrs[i]->u[j] ;
					}
					dp++;
				}
			}
		if(varsum>0.0)
			uchange = sqrt (uchange/varsum) ;
		else
			uchange = 0.0;
		}
	}
	return(0) ;
}

int		actcol(Kp,Fp,diagp,neqns,code)
double	*Kp, *Fp ;
long	*diagp ;
int		neqns, code ;

{
	int	i, j, nf, id, n, ni, nj ;
	double	dot(), *ip, *idp, t ;
	char	m1[255], m2[255], m3[255], m4[255] ;
	
	for(j=1;j<neqns;j++) {
/*
		if( (j % 10 == 0) && (code < 2)) {
			sprintf(m2,"  %d ",j) ;
			UpDateMes(m1,m2,m3,m4) ;
		}
*/
		nf = *(diagp+j) - *(diagp+j-1) ;
		if(code < 2) {
			ip = Kp + *(diagp+j-1) + 1 ;
			id = j - nf + 1 ;
		
			for(i=1;i<nf;i++) {
				n = ((nj = i-1) < (ni = *(diagp+id) - *(diagp+id-1) - 1)) ? nj : ni ;
				idp = Kp + *(diagp+id++) ;
				*ip++ -= dot(ip-n-1,idp-n,n) ;
			}
			ip = Kp + *(diagp+j-1) + 1 ;
			id = j - nf + 1 ;
			idp = Kp + *(diagp+j) ;
			for(i=1;i<nf;i++) {
				t = *ip ;
				*ip /= *(Kp + *(diagp+id++)) ;
				*idp -= t * *ip++ ;
			}
		}
		ip = Kp + *(diagp+j-1) + 1 ;
		id = j - nf + 1 ;
		if(code > 0)
			*(Fp + j) -= dot(ip,Fp+id,nf-1) ; 
	}
	if(code == 0) {
		EndMes() ;
		return(0) ;
	}
	for(j=0;j<neqns;j++)
		*(Fp + j) /= *(Kp + *(diagp+j)) ;
		
	for(j=neqns-1;j>-1;j--) {
		if(j == 0)
			nf = 1;
		else
			nf = *(diagp+j) - *(diagp+j-1) ;
		for(i=1;i<nf;i++) {
			*(Fp + j - nf + i) -= *(Fp + j) * *(Kp + *(diagp+j-1) + i) ;
		}
	}

	if(code < 2)
		EndMes() ;
	return(0) ;
}

int		uactcl(Kp,Lp,Fp,diagp,neqns,code)
double	*Kp, *Lp, *Fp ;
long	*diagp ;
int		neqns, code ;

{
	int	i, j, nf, id, n, ni, nj ;
	double	dot(), *ip, *iLp, *idp, t ;
	char	m1[255], m2[255], m3[255], m4[255] ;
	for(j=1;j<neqns;j++) {

/*
		if( (j % 10 == 0) && (code < 2)) {
			sprintf(m2,"  %d ",j) ;
			UpDateMes(m1,m2,m3,m4) ;
		}
*/
		nf = *(diagp+j) - *(diagp+j-1) ;
		if(code < 2) {
			iLp = Lp + *(diagp+j-1) + 1 ;
			id = j - nf + 1 ;
		
			for(i=0;i<nf-1;i++) {
				n = ((nj = i) < (ni = *(diagp+id) - *(diagp+id-1) - 1)) ? nj : ni ;
				idp = Kp + *(diagp+id) ;
				id++ ;
				*iLp = (*iLp - dot(iLp-n,idp-n,n))/ *idp ;
//if(*idp<0.01)
//printf("attention:  dividing by %g\n",*idp);
				iLp++ ;
			}
			ip = Kp + *(diagp+j-1) + 1 ;
			id = j - nf + 1 ;
		
			for(i=0;i<nf;i++) {
				n = ((nj = i) < (ni = *(diagp+id) - *(diagp+id-1) - 1)) ? nj : ni ;
				idp = Lp + *(diagp+id) ;
				id++ ;
				*ip -= dot(ip-n,idp-n,n) ;
				ip++ ;
			}
		}
		iLp = Lp + *(diagp+j-1) + 1 ;
		id = j - nf + 1 ;
		if(code > 0)
			*(Fp + j) -= dot(iLp,Fp+id,nf-1) ; 
	}
	if(code == 0) {
		EndMes() ;
		return(0) ;
	}
	for(j=neqns-1;j>-1;j--) {
		*(Fp + j) /= *(Kp + *(diagp+j)) ;
		if(j == 0)
			nf = 1;
		else
			nf = *(diagp+j) - *(diagp+j-1) ;
		for(i=1;i<nf;i++) {
			*(Fp + j - nf + i) -= *(Fp + j) * *(Kp + *(diagp+j-1) + i) ;
		}
	}
	if(code < 2)
		EndMes() ;
	return(0) ;
}


double	dot(p1,p2,n)
double	*p1, *p2 ;
int		n ;

{
	register double	t=0.0 ;
	
	while(n-- > 0) {
			t += *p1++ * *p2++ ;
	}
	return(t) ;
}

int		get_Derivs(elmntp,pt,nvar,ddx,ddy,ddz)
struct element		*elmntp ;
struct gausspts 	*pt ;
int		nvar ;
double	*ddx, *ddy, *ddz ;

{
	struct shapefuncs	fv;
	struct node			*np[MNSF] ;
	struct govfuncs		gf ;
	int					i, j, ns, ni, nm ;
	
	ns =  nsf(elmntp->vtype) ;
	for(j=0;j<elmntp->nnds;j++)
		np[j] = elmntp->nps[j] ;
	get_shape(elmntp->vtype,pt,&fv) ;
	get_gf(elmntp,&gf,&fv,&fv,np,pt->w) ;
	*ddx = *ddy = *ddz = 0.0 ;
	for(j=0;j<elmntp->nnds;j++) {
		*ddx += gf.dfdx[j] * np[j]->u[nvar] ;
		*ddy += gf.dfdy[j] * np[j]->u[nvar] ;
		*ddz += gf.dfdz[j] * np[j]->u[nvar] ;
	}
	*ddx *= gf.p[0] ;
	*ddy *= gf.p[0] ;
	*ddz *= gf.p[0] ;
	return(ns) ;
}

int		get_gf(elp,gfp,fvp,fgp,np,w)
struct element		*elp ;
struct shapefuncs	*fvp, *fgp ;
struct node			*np[MNSF] ;
struct govfuncs		*gfp ;
double				w ;

{
	double	J[NDIMS][NDIMS], get_J() ;
	int		i, j ;
	
	gfp->dof = fvp->dof ;
	
	for(i=0;i<N.params;i++) {
		gfp->p[i] = elp->p[i] ;
		for(j=0;j<fvp->dof;j++)
			gfp->p[i] += fvp->f[j] * np[j]->p[i] ;
	}
	for(j=0;j<fvp->dof;j++)
		gfp->f[j] = fvp->f[j]  ;
	
	gfp->detJ = get_J(np,fgp,J) ;
	gfp->w = w ;
	invertJ(J,gfp->detJ) ;
	
	for(j=0;j<fvp->dof;j++) {
	
		switch (N.dims) {
		
			case 1 :
				gfp->dfdx[j] = J[0][0] * fvp->dfdr[j] ;
				gfp->dfdy[j] = 0.0 ;
				gfp->dfdz[j] = 0.0 ;
				break ;
				
			case 2 :
				gfp->dfdx[j] = J[0][0] * fvp->dfdr[j] + J[0][1] * fvp->dfds[j] ;
				gfp->dfdy[j] = J[1][0] * fvp->dfdr[j] + J[1][1] * fvp->dfds[j] ;
				gfp->dfdz[j] = 0.0 ;
				break ;
				
			case 3 :
				gfp->dfdx[j] = J[0][0] * fvp->dfdr[j] + J[0][1] * fvp->dfds[j] + J[0][2] *fvp->dfdt[j] ;
				gfp->dfdy[j] = J[1][0] * fvp->dfdr[j] + J[1][1] * fvp->dfds[j] + J[1][2] *fvp->dfdt[j] ;
				gfp->dfdz[j] = J[2][0] * fvp->dfdr[j] + J[2][1] * fvp->dfds[j] + J[2][2] *fvp->dfdt[j] ;
				break ;
		}
	}
	gfp->detJ *= gfp->w ;
	return(0) ;
}

double	get_J(np,fgp,J)
struct shapefuncs	*fgp ;
struct node			*np[MNSF] ;
double			J[NDIMS][NDIMS] ;

{
	int	i,j ;
	
	for(i=0;i<N.dims;i++)
		for(j=0;j<N.dims;j++)
			J[i][j] = 0.0 ;
			
	for(i=0;i<fgp->dof;i++) 	
		J[0][0] += fgp->dfdr[i] * np[i]->x[0] ;
	if(N.dims == 1)
		return( J[0][0] ) ;
	
	for(i=0;i<fgp->dof;i++) {
		J[0][1] += fgp->dfdr[i] * np[i]->x[1] ;
		J[1][0] += fgp->dfds[i] * np[i]->x[0] ;
		J[1][1] += fgp->dfds[i] * np[i]->x[1] ;
	}
	if(N.dims == 2) 
		return( J[0][0] * J[1][1] - J[0][1] * J[1][0] ) ;
	
	for(i=0;i<fgp->dof;i++) {
		J[0][2] += fgp->dfdr[i] * np[i]->x[2] ;
		J[1][2] += fgp->dfds[i] * np[i]->x[2] ;
		J[2][2] += fgp->dfdt[i] * np[i]->x[2] ;
		J[2][0] += fgp->dfdt[i] * np[i]->x[0] ;
		J[2][1] += fgp->dfdt[i] * np[i]->x[1] ;
	}
	return( J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1] )
	    -   J[0][1] * (J[1][0] * J[2][2] - J[2][0] * J[1][2] )
		+   J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][1] ) ) ;
}

int		invertJ(J,detJ)
double	J[NDIMS][NDIMS], detJ ;

{
	double	K[NDIMS][NDIMS] ;
	int		i,j ;
	
	for(i=0;i<N.dims;i++)
		for(j=0;j<N.dims;j++)
			K[i][j] = J[i][j] ;
			
	switch (N.dims) {
	
		case 1 :
			J[0][0] = 1.0 / K[0][0] ;
			break ;
			
		case 2 :
			J[0][0] =  K[1][1] / detJ ;
			J[0][1] = -K[0][1] / detJ ;
			J[1][0] = -K[1][0] / detJ ;
			J[1][1] =  K[0][0] / detJ ;
			break ;
			
		case 3 :
			J[0][0] = ( K[1][1] * K[2][2] - K[1][2] * K[2][1] ) / detJ ;
			J[0][1] = ( K[0][2] * K[2][1] - K[0][1] * K[2][2] ) / detJ ;
			J[0][2] = ( K[0][1] * K[1][2] - K[1][1] * K[0][2] ) / detJ ;
			J[1][0] = ( K[1][2] * K[2][0] - K[1][0] * K[2][2] ) / detJ ;
			J[1][1] = ( K[0][0] * K[2][2] - K[0][2] * K[2][0] ) / detJ ;
			J[1][2] = ( K[1][0] * K[0][2] - K[0][0] * K[1][2] ) / detJ ;
			J[2][0] = ( K[1][0] * K[2][1] - K[1][1] * K[2][0] ) / detJ ;
			J[2][1] = ( K[0][1] * K[2][0] - K[0][0] * K[2][1] ) / detJ ;
			J[2][2] = ( K[0][0] * K[1][1] - K[0][1] * K[1][0] ) / detJ ;
			break ;
	}
	return(0) ;
}


int		initm(K,n)
double	*K ;
int		n ;

{

	register int	i ;

	for(i=0;i<n;i++) 
			*(K++) = 0.0 ;
	return(0) ;
}

double		GetElErr(elmntp,nvar)
struct element		*elmntp ;
int		nvar ;

{
/*
	struct gausspts 	g[MGPTS] ;
	struct shapefuncs	fv, fve, fg, *fgp;
	struct node			*np[MNSF] ;
	struct govfuncs		gf ;
	int					i, ii, j, ns, ni ;
	double				J[NDIMS][NDIMS], dfdxe, dfdye,fxiui, fyiui, err, area ;
	
	ns =  nsf(elmntp->vtype) ;
	ni = get_gspts(elmntp->vtype,g) ;
	for(j=0;j<elmntp->nnds;j++)
		np[j] = elmntp->nps[j]] ;
	
	err = area = 0.0 ;
	for(i=0;i<ni;i++) {
		get_shape(elmntp->vtype,&g[i],&fv) ;
		mid_pt(&g[i],&fve) ;
		fgp = &fv ;
		get_gf(elmntp,&gf,&fv,fgp,np,g[i].w) ;
	
		for(ii=0;ii<N.dims;ii++)
			for(j=0;j<N.dims;j++)
				J[ii][j] = 0.0 ;
			
		for(ii=0;ii<fgp->dof;ii++) {	
			J[0][0] += fv.dfdr[ii] * np[ii]->x[0] ;
			J[0][1] += fv.dfdr[ii] * np[ii]->x[1] ;
			J[1][0] += fv.dfds[ii] * np[ii]->x[0] ;
			J[1][1] += fv.dfds[ii] * np[ii]->x[1] ;
		}
		invertJ(J,gf.detJ) ;
		dfdxe = J[0][0] * fve.dfdr[8] + J[0][1] * fve.dfds[8] ;
		dfdye = J[1][0] * fve.dfdr[8] + J[1][1] * fve.dfds[8] ;
		
		fxiui = fyiui = 0.0 ;
		for(j=0;j<fgp->dof;j++) {
			fxiui += gf.dfdx[j] * np[j]->u[nvar] ;
			fyiui += gf.dfdy[j] * np[j]->u[nvar] ;
		}
		err += (fxiui * dfdxe + fyiui * dfdye) * gf.detJ ;
		area += gf.detJ ;
	}
*/
	return(0.0) ;
}

int	mid_pt(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,s ;
	
	r = x->x ;
	s = x->y ;
	
	p->dfdr[8] = -2.0 * r * (1.0 - s * s) ;
	p->dfds[8] = -2.0 * s * (1.0 - r * r) ;
	
	return(0) ;
}

int		NodalFlux(np,nvar,fx,fy,fz)
struct node	*np ;
double		*fx,*fy,*fz ;
int			nvar ;

{
	double	fxe[10], fye[10], fze[10] ;
	int		i, j, in;
	struct element	*elp ;
	struct gausspts	pt ;
	
	elp = gp.El ;
	i = 0 ;
	pt.z = pt.w = 0.0 ;
	while(elp != NULL) {
		in = NdinEl(np->n,&elp) ;
		if(in < 0)
			break ;
		get_nlcs(elp->vtype,in,&(pt.x),&(pt.y)) ;
		get_Derivs(elp,&pt,nvar,&fxe[i],&fye[i],&fze[i]) ;
		i++ ;
		elp = elp->nextelp ;
	}
	*fx = *fy = *fz = 0.0 ;
	for(j=0;j<i;j++) {
		*fx += fxe[j] ;
		*fy += fye[j] ;
		*fz += fze[j] ;
	}
	if(i > 0) {
		*fx /= i ;
		*fy /= i ;
		*fz /= i ;
	}
	return(i) ;
}

int		NdinEl(n,elpp)
int		n ;
struct element	*(*elpp) ;

{
	int	i, rc;
	struct element	*selp ;
	
	rc = 0 ;
	selp = *elpp ;
	while(selp != NULL) {
		for(i=0;i<selp->nnds;i++) {
			if(n == (selp->nps[i])->n) {
				rc = 1 ;
				*elpp = selp ;
				break ;
			}
		}
		if(rc == 1)
			return(i) ;
		selp = selp->nextelp ;
	}
	return(-1) ;
}


/*

double	integrate(elcode,elp,nvar)
int		elcode,  nvar ;
struct element	*elp ;

{
	struct node		*np[MNSF] ;
	struct gausspts	g[MGPTS] ;
	struct shapefuncs	sf ;
	int		i, j, ng ;
	double	integral = 0.0 , t[10] ;
	
	for(i=0;i<elp->nnds;i++)
		np[i] = elp->nps[i] ;
	ng = get_gspts(elp->vtype,g) ;
	for(i=0;i<ng;i++) {
		for(j=0;j<10;j++)
			t[j] = 0.0 ;
		get_shape(elp->vtype,&g[i],&sf) ;
		switch (elcode) {
			case	1 :
				switch (N.dims) {
					case 2 :
						for(j=0;j<sf.dof;j++) {
							t[0] += sf.f[j] * np[j]->u[nvar] ;
							t[1] += sf.dfdr[j] * np[j]->x[0] ;
							t[2] += sf.dfdr[j] * np[j]->x[1] ;
						}
						integral += t[0] * sqrt(t[1]*t[1] + t[2]*t[2]) * g[i].w ;
						break ;
					case 3 :
						break ;
				}
				break ;
			case 0 :
				break ;
		}
	}
	return(integral) ;
}
*/
			
int		get_cline(elp,cval,niv,nvar,tol,xl,yl,csystem)
int		niv, nvar ,csystem;
double	xl[], yl[], tol, cval ;
struct element	*elp ;

{
	int	i, j, nipts, rc ;
	double	x, y, Dx, Dy, dfdr, dfds ;
	struct node	*np[MNSF] ;
	struct shapefuncs		sf ;
	struct gausspts		xg ;
	if(niv == 0)
		tol = -1.0 ;
	nipts = isvalinel(elp,np,cval,nvar,tol,&dfdr,&dfds,xl,yl) ;
	switch (nipts) {
		case	0 :
			break ;
		case 2 :
			if((fabs(xl[0]-xl[1])<.000001)&&(fabs(yl[0]-yl[1])<.00001)){
				nipts = 0;
				break;
			}
			if(csystem==1) {
				Dx = xl[1] - xl[0] ; /*printf("Dx =%g\n",Dx);*/
				Dy =  yl[1] - yl[0] ;
				xl[niv+1] = xl[1] ;
				yl[niv+1] = yl[1] ;
				for(i=1;i<=niv;i++) {
					xl[i] = xl[i-1] + (xl[niv+1] - xl[i-1]) * 1 / (niv - i + 2) ;
					yl[i] = yl[i-1] + (yl[niv+1] - yl[i-1]) * 1 / (niv - i + 2) ;
					find_2Dloc(np,nvar+N.dims,elp->gtype,cval,tol,Dx,Dy,dfdr,dfds,&xl[i],&yl[i]) ;
				}
				xg.z = 0.0 ;
				xg.w = -1.0 ;
				for(i=0;i<=niv+1;i++) {
					xg.x = xl[i] ;
					xg.y = yl[i] ; /*printf(" i=%d\txl=%g\tyl=%g\n",i,xl[i],yl[i]);*/
					rc = get_shape(elp->gtype,&xg,&sf) ;
					xl[i] = yl[i] = 0.0 ;
					for(j=0;j<sf.dof;j++) {
						xl[i] += sf.f[j] * np[j]->x[0] ;
						yl[i] += sf.f[j] * np[j]->x[1] ; /*printf(" i=%d\txl=%g\tyl=%g\n",i,xl[i],yl[i]);*/
					}
				}
			}
			break ;
	}
	
	return(nipts) ;
}

int		get_cline1(elp,cval,ncoord,niv,nvar,tol,xl,yl,zl,bl)
int		ncoord, niv, nvar ;
double		xl[], yl[], zl[], bl[], tol, cval ;
struct element	*elp ;

{
//This function finds the intersection of a coordinate with the mesh
	int	i, j, nipts, rc ;
	double	x, y, Dx, Dy, dfdr, dfds ;
	struct node	*np[MNSF] ;
	struct shapefuncs		sf ;
	struct gausspts		xg ;
	if(niv == 0)
		tol = -1.0 ;
	nipts = isvalinel1(elp,np,cval,nvar,ncoord,tol,&dfdr,&dfds,xl,yl,zl,bl) ;/*printf("nipts =%d\n",nipts);*/
	switch (nipts) {
		case	0 :
			break ;
		case 2 :
			Dx = xl[1] - xl[0] ; /*printf("Dx =%g\n",Dx);*/
			Dy = yl[1] - yl[0] ;
			xl[niv+1] = xl[1] ;
			yl[niv+1] = yl[1] ;
/*for niv!=0			for(i=1;i<=niv;i++) {
				xl[i] = xl[i-1] + (xl[niv+1] - xl[i-1]) * 1 / (niv - i + 2) ;
				find_2Dloc(np,nvar+N.dims,elp->gtype,cval,tol,Dx,Dy,dfdr,dfds,&xl[i],&yl[i]) ;
			}*/
			xg.z = 0.0 ;
			xg.w = -1.0; 
			
/*			if(ncoord==0) {

				for(i=0;i<niv+1;i++) {
					xg.x = xl[i];
					xg.y = yl[i];
					rc = get_shape(elp->gtype,&xg,&sf);
					xl[i] = 0.0;
					yl[i] = 0.0;
					for(j=0;j<sf.dof;j++) { 
						xl[i] += sf.f[j] * np[j]->x[1] ;
						yl[i] += sf.f[j] * np[j]->x[0] ;
					}
				}
			}
			else	*/{
					
				for(i=0;i<=niv+1;i++) {
					xg.x = xl[i] ;
					xg.y = yl[i] ;
					rc = get_shape(elp->gtype,&xg,&sf) ;
					xl[i] = yl[i] = 0.0 ;
					for(j=0;j<sf.dof;j++) {
						xl[i] += sf.f[j] * np[j]->x[0] ; /*printf(" i=%d\txl=%g\n",i,xl[i]);*/
						yl[i] += sf.f[j] * np[j]->x[1] ;
					}
				}
			}
			break ;
	}/*  printf("nipts =%d\n",nipts);*/
	return(nipts);
}
	

int		isvalinel(elp,np,cval,nvar,tol,dfdr,dfds,xl,yl)
int		nvar ;
double		xl[], yl[], tol, cval, *dfdr, *dfds ;
struct node	*np[] ;
struct element	*elp ;

{
	int	j, k, nbps, etyp, ng, next, prev, ncpts, aj[MNSF], ak[MNSF] ;
	double		Dx, Dy, r, s, t, f, fr, fs, ft ;
	struct gausspts	g[MGPTS] ;
	
	nbps = nbounds(elp->gtype) ;
	etyp = elp->gtype ;
	for(j=0;j<elp->nnds;j++) 
		np[j] =  elp->nps[j] ;
	prev = (cval <= np[nbps-1]->u[nvar]) ? 1 : -1 ;
	ncpts = 0 ;
	for(j=0;j<nbps;j++) {
		next = (cval <= np[j]->u[nvar]) ? 1 : -1 ;
		if(next != prev) {
			aj[ncpts] = j ;
			if( j == 0 )
				ak[ncpts] = nbps-1 ;
			else
				ak[ncpts] = j - 1 ;
			ncpts++ ;
		}
		prev = next ;
	}
	if((ncpts > 0) && (tol > 0.0) ) {
		ng = get_gspts(etyp,g) ;
		*dfdr = *dfds = 0.0 ;
		for(j=0;j<ng;j++) {
			r = g[j].x ;
			s = g[j].y ;
			t = g[j].z ;
			get_value(np,nvar,etyp,r,s,t,2,&f,&fr,&fs,&ft) ;
			*dfdr += fr ;
			*dfds += fs ;
		}
		*dfdr /= ng ;
		*dfds /= ng ;
	}
	for(j=0;j<ncpts;j++) { /*printf("j =%d\n",j);*/
			GetIPoint(aj[j],ak[j],np,etyp,cval,nvar,&xl[j],&yl[j],&Dx,&Dy) ;
	/*			printf("xl =%g\tyl =%g\n",xl[j],yl[j]);*/
			if(tol > 0.0)
				find_2Dloc(np,nvar+N.dims,etyp,cval,tol,Dx,Dy,*dfdr,*dfds,&xl[j],&yl[j]) ;
	}
	return(ncpts) ;
}
int		isvalinel1(elp,np,cval,nvar,ncoord,tol,dfdr,dfds,xl,yl,zl,bl)
int		ncoord, nvar ;
double		xl[], yl[], zl[], bl[], tol, cval, *dfdr, *dfds ;
struct node	*np[] ;
struct element	*elp ;

{
	int		j, k, nbps, etyp, ng, next, prev, ncpts, aj[MNSF], ak[MNSF] ;
	double		Dx, Dy, r, s, t, f, fr, fs, ft ;
	struct gausspts	g[MGPTS] ;
	
	nbps = nbounds(elp->gtype) ;
	etyp = elp->gtype ;
	for(j=0;j<elp->nnds;j++) 
		np[j] =  elp->nps[j] ;
	prev = (cval <= np[nbps-1]->x[ncoord]) ? 1 : -1 ;/*printf("xcoord=%g\n",np[nbps-1]->x[ncoord]);*/
	ncpts = 0 ;
	for(j=0;j<nbps;j++) {
		next = (cval <= np[j]->x[ncoord]) ? 1 : -1 ;
		if(next != prev) {
			aj[ncpts] = j ;
			if( j == 0 )
				ak[ncpts] = nbps-1 ;
			else
				ak[ncpts] = j - 1 ;
			ncpts++ ;
		}
		prev = next ;
	}
/*	if((ncpts > 0) && (tol > 0.0) ) {
		ng = get_gspts(etyp,g) ;
		*dfdr = *dfds = 0.0 ;
		for(j=0;j<ng;j++) {
			r = g[j].x ;
			s = g[j].y ;
			t = g[j].z ;
			get_value(np,nvar,etyp,r,s,t,2,&f,&fr,&fs,&ft) ;
			*dfdr += fr ;
			*dfds += fs ;
		}
		*dfdr /= ng ;
		*dfds /= ng ;
	}*/
	for(j=0;j<ncpts;j++) { /* printf(" j =%d\n",j);*/
			GetIPoint1(aj[j],ak[j],np,etyp,cval,ncoord,nvar,&xl[j],&yl[j],&zl[j],&bl[j],&Dx,&Dy) ;
/*				printf("xl =%g\n",xl[j]);*/
/*			if(tol > 0.0)
				find_2Dloc(np,nvar+N.dims,etyp,cval,tol,Dx,Dy,*dfdr,*dfds,&xl[j],&yl[j]) ;*/
	}/*printf("ncpts =%d\n",ncpts);*/
	return(ncpts);
}

int		GetIPoint(j,k,np,eltype,u,n,xp,yp,Dx,Dy) 
double	*xp, *yp, *Dx, *Dy, u ;
int		j, k, n, eltype;
struct node	*np[] ;

{
	double	x1, y1, x2, y2 ;
	
	get_nlcs(eltype,k,&x1,&y1) ;
	get_nlcs(eltype,j,&x2,&y2) ; /*printf("x1 =%g\tx2 =%g\ny1 =%g\ty2 =%g\n",x1,x2,y1,y2);*/
	*xp = x1 + (u - np[k]->u[n])/(np[j]->u[n] - np[k]->u[n])*(x2 - x1) ;
	*yp = y1 + (u - np[k]->u[n])/(np[j]->u[n] - np[k]->u[n])*(y2 - y1) ;
	*Dx = y2 - y1 ;
	*Dy = x1 - x2 ;
	return(0) ;
}
int		GetIPoint1(j,k,np,eltype,cval,ncoord,nvar,xp,yp,zp,bp,Dx,Dy) 
double		*xp, *yp, *zp, *bp, *Dx, *Dy, cval ;
int		j, k, eltype, ncoord, nvar;
struct node	*np[] ;

{
	double	x1, y1, x2, y2 ;
/*	printf(" j=%d\t k=%d\n",j,k);*/
	get_nlcs(eltype,k,&x1,&y1) ;
	get_nlcs(eltype,j,&x2,&y2) ;/* printf("x1 =%g\tx2=%g\ny1 =%g\ty2 =%g\n",x1,x2,y1,y2);*/
/*	if(ncoord==1) {*/
		*xp = x1 + (cval - np[k]->x[ncoord])/(np[j]->x[ncoord] - np[k]->x[ncoord])*(x2-x1);
		*yp = y1 + (cval - np[k]->x[ncoord])/(np[j]->x[ncoord] - np[k]->x[ncoord])*(y2-y1);
		*Dx = y2-y1 ;
		*Dy = x1-x2 ;	
/*		}
		else
			{	
		*xp = y1 + (cval - np[k]->x[ncoord])/(np[j]->x[ncoord] - np[k]->x[ncoord])*(y2 - y1) ;
		*yp = x1 + (cval - np[k]->x[ncoord])/(np[j]->x[ncoord] - np[k]->x[ncoord])*(x2 - x1) ;
		*Dx = x1 - x2;
		*Dy = y2 - y1;
	}*/
	*zp = (np[k]->u[nvar]+np[k]->p[0])+(cval-np[k]->x[ncoord])/(np[j]->x[ncoord]-np[k]->x[ncoord])*((np[j]->u[nvar]+np[j]->p[0])-(np[k]->u[nvar]+np[k]->p[0])) ;
	*bp = np[k]->p[0]+(cval-np[k]->x[ncoord])/(np[j]->x[ncoord]-np[k]->x[ncoord])*(np[j]->p[0]-np[k]->p[0]);

	return(0) ;
}


int		find_2Dloc(np,nvar,eltype,fc,tol,Dr,Ds,dfdr,dfds,r,s)
double	fc, tol, Dr, Ds, dfdr, dfds, *r, *s ;
int		eltype, nvar ;
struct node		*np[] ;

{
	int i, eq, imax = 100 ;
	double	t, f, ddr, dds, dfdt, df, dr, ds, sl ;
	
/*	printf(" in find_2Dloc...") ;*/
	
	t = 0.0 ;
	if(fabs(Dr) >= fabs(Ds)) {
		eq = 2 ;
		sl = -Ds/Dr ;
	}
	else {
		eq = 1 ;
		sl = -Dr/Ds ;
	}
	for(i=0;i<imax;i++) {
		if( get_value(np,nvar,eltype,*r,*s,t,2,&f,&ddr,&dds,&dfdt)  != 0)
			return(-1);
		if(fabs(df = f - fc) < tol) {
/*			printf(" fc %f f %f df %f tol %f\n",fc,f,df,tol);*/
			break ;
		}
		if(eq == 2) {
			ds = df / (ddr*sl + dds) ;
			dr = sl*ds ;
		}
		else {
			dr = df / (ddr + sl*dds) ;
			ds = sl * dr ;
		}
		if(fabs(dr) > 0.5 || fabs(ds) > 0.5) {
			if(eq == 2) {
				ds = df / (dfdr*sl + dfds) ;
				dr = sl*ds ;
			}
			else {
				dr = df / (dfdr + sl*dfds) ;
				ds = sl * dr ;
			}
		}
		if(fabs(dr) > 0.5 || fabs(ds) > 0.5) 
			return(-2) ;
		*r -= dr ;
		*s -= ds ;
/*		printf(" %d ",i) ;*/
	}
/*	printf(" ...out\n") ;*/
	return(i) ;
}

int		get_value(np,nvar,eltype,r,s,t,code,f,dfdr,dfds,dfdt)
double	*f, *dfdr, *dfds, *dfdt, r, s, t ;
int		eltype, nvar, code ;
struct node		*np[] ;

{
	int	i, rc ;
	struct shapefuncs		sf ;
	struct gausspts		x ;
	
	x.x = r ;
	x.y = s ;
	x.z = t ;
	x.w = 0.0  ;
	if(code == 0)
		x.w = -1.0 ;
	rc = get_shape(eltype,&x,&sf) ;
	
	switch (code) {
		case	3 :
			*dfdt = 0.0 ;
			for(i=0;i<sf.dof;i++) 
				*dfdt += sf.dfdt[i] * np[i]->x[nvar] ;
		case 2 :
			*dfds = 0.0 ;
			for(i=0;i<sf.dof;i++) 
				*dfds += sf.dfds[i] * np[i]->x[nvar] ;
		case 1 :
			*dfdr = 0.0 ;
			for(i=0;i<sf.dof;i++) 
				*dfdr += sf.dfdr[i] * np[i]->x[nvar] ;
		case 0 :
			*f = 0.0 ;
			for(i=0;i<sf.dof;i++) 
				*f += sf.f[i] * np[i]->x[nvar] ;
	}
	return(rc) ;
}
