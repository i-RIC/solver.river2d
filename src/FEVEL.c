#include "Fe_PS.h"
#include "ctype.h"
#define	MAXLINE	200

extern struct control  N ;
extern struct pointers gp;
extern struct transient tvals ;

double get_dbl();

int	get_xsecvel(code)
int	code;
{
	struct	element	*elp;
	int		niv,i,j,nnodes,nordered,nvar,get_velocity();
	double	cval,ymin;
	FILE		*outf,*put_fptr();
	
	struct	node	vn[50],ond[50],*fn;
	printf("Input x-coordinate value at which velocities should be printed\n");
	scanf("%lf",&cval);

	if((outf=put_fptr(1))==NULL) 
		return(-1);
	
	nnodes=0;
	nordered=0;
	niv=code;
	nvar=0;
	
	elp=gp.El;
// check nvar
	for(i=0;i<N.elms;i++) {
		nnodes=get_velocity(elp,nvar,niv,cval,vn,nnodes);
//nnodes should be incremented by get_velocity
		elp=elp->nextelp;
	}
printf("nnodes=%d\n",nnodes);
for(i=0;i<0;i++)
printf("i=%d\tx=%g\ty=%g\th=%g\tqx=%g\tuvel=%g\n",i,vn[i].x[0],vn[i].x[1],vn[i].u[0],vn[i].u[1],(vn[i].u[1]/vn[i].u[0]));
		fn=NULL;
	for(i=0;i<nnodes;i++)
		vn[i].fxc=0;
	for(j=0;j<nnodes;j++) {
		ymin=10E8;
		
		for(i=0;i<nnodes;i++) {
			if((vn[i].x[1]<ymin)&&(vn[i].fxc==0)) {
				ymin=vn[i].x[1];
				fn=&vn[i];
			}
		}
//now turn the found node off
		if(fn!=NULL)
			fn->fxc=1;
//now turn off all nodes with the same value of y

		if(fn!=NULL){
			for(i=0;i<nnodes;i++) {
				if((fabs(fn->x[1]-vn[i].x[1])<10E-8)&&(fn!=&vn[i]))
					vn[i].x[1]=1000000.00;
			}
	
			if(fn->x[1]==1000000.00)
				break;
			else {
				ond[nordered]=*fn;
				nordered++;
			}
	
		}
/*	
		ond[nordered]=*fn;
		nordered++;
*/
	}
printf("nordered=%d\n",nordered);
// Now print the results
	fprintf(outf,"y-coord	bl-calc	wsl-calc	depth	u-calc	v-vel\n");
	for(i=0;i<nordered;i++)
	fprintf(outf,"%g\t%g\t%g\t%g\t%g\t%g\n",ond[i].x[1],ond[i].p[0],(ond[i].u[0]+ond[i].p[0]),ond[i].u[0],(ond[i].u[1]/ond[i].u[0]),(ond[i].u[2]/ond[i].u[0]));

	fclose(outf);

	return(nordered);
}
	
int	get_velocity(elp,nvar,niv,cval,vn,nnodes)
struct	element	*elp;
struct	node		vn[];
int		niv,nnodes;
double	cval;

{
	double	xl[MNSF],yl[MNSF],zl[MNSF],bl[MNSF],qxl[MNSF],qyl[MNSF],tol;
	int		i,ncpts,ncoord=0,get_cline2();
	
	tol = 0.02;
	
	if((ncpts=get_cline2(elp,cval,ncoord,niv,nvar,tol,xl,yl,zl,bl,qxl,qyl))==2) {
		for(i=0;i<(ncpts+niv);i++) {
			vn[nnodes].x[0]=xl[i];
			vn[nnodes].x[1]=yl[i];
			vn[nnodes].p[0]=bl[i];
			vn[nnodes].u[0]=zl[i]-bl[i];
			vn[nnodes].u[1]=qxl[i];
			vn[nnodes].u[2]=qyl[i];
			nnodes++;
		}
	}
	return(nnodes);
}

int		get_cline2(elp,cval,ncoord,niv,nvar,tol,xl,yl,zl,bl,qxl,qyl)
int		ncoord, niv, nvar ;
double		xl[], yl[], zl[], bl[], tol, cval,qxl[],qyl[] ;
struct element	*elp ;

{
//This function finds the intersection of a coordinate with the mesh
//It is used to print out velocity vals at a xsec
	int	i, j, nipts, rc, isvalinel2(), get_shape() ;
	double Dx, Dy, dfdr, dfds ;
	struct node	*np[MNSF] ;
	struct shapefuncs		sf ;
	struct gausspts		xg ;
	if(niv == 0)
		tol = -1.0 ;
	nipts = isvalinel2(elp,np,cval,nvar,ncoord,tol,&dfdr,&dfds,xl,yl,zl,bl,qxl,qyl) ;

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
	}
//printf("nipts =%d\n",nipts);
	return(nipts);
}	

int		isvalinel2(elp,np,cval,nvar,ncoord,tol,dfdr,dfds,xl,yl,zl,bl,qxl,qyl)
int		ncoord, nvar ;
double		xl[], yl[], zl[], bl[], tol, cval, *dfdr, *dfds,qxl[],qyl[] ;
struct node	*np[] ;
struct element	*elp ;

{
	int		j, nbps, etyp, next, prev, ncpts, aj[MNSF], ak[MNSF],nbounds(),GetIPoint2() ;
	double		Dx, Dy ;
	
	nbps = nbounds(elp->gtype) ;
	etyp = elp->gtype ;
	for(j=0;j<elp->nnds;j++) 
		np[j] =  elp->nps[j] ;
	prev = (cval <= np[nbps-1]->x[ncoord]) ? 1 : -1;
//printf("xcoord=%g\n",np[nbps-1]->x[ncoord]);
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
			GetIPoint2(aj[j],ak[j],np,etyp,cval,ncoord,nvar,&xl[j],&yl[j],&zl[j],&bl[j],&qxl[j],&qyl[j],&Dx,&Dy) ;
/*				printf("xl =%g\n",xl[j]);*/
/*			if(tol > 0.0)
				find_2Dloc(np,nvar+N.dims,etyp,cval,tol,Dx,Dy,*dfdr,*dfds,&xl[j],&yl[j]) ;*/
	}/*printf("ncpts =%d\n",ncpts);*/
	return(ncpts);
}	

int		GetIPoint2(j,k,np,eltype,cval,ncoord,nvar,xp,yp,zp,bp,qxp,qyp,Dx,Dy) 
double		*xp, *yp, *zp, *bp, *Dx, *Dy, cval,*qxp,*qyp ;
int		j, k, eltype, ncoord, nvar;
struct node	*np[] ;

{
	int		get_nlcs();
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

	*qxp = np[k]->u[1]+(cval-np[k]->x[ncoord])/(np[j]->x[ncoord]-np[k]->x[ncoord])*(np[j]->u[1]-np[k]->u[1]);

	*qyp = np[k]->u[2]+(cval-np[k]->x[ncoord])/(np[j]->x[ncoord]-np[k]->x[ncoord])*(np[j]->u[2]-np[k]->u[2]);


	return(0) ;
}


