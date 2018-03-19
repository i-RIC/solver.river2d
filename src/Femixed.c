#include "Fe_PS.h"

extern struct  transient       tvals;

int     check_mixed(elp)
struct  element *elp;
{
	int     i;
	double  mult;
	
	for(i=0;i<(elp->nnds-1);i++) {
		mult = ((elp->nps[i]->uo[0]-(elp->nps[i]->ice[0]*tvals.sg_ice))-tvals.gwH)*
			   ((elp->nps[i+1]->uo[0]-(elp->nps[i+1]->ice[0]*tvals.sg_ice))-tvals.gwH);
		if(mult<0.0)
			return(1);
	}
	mult=((elp->nps[elp->nnds-1]->uo[0]-(elp->nps[elp->nnds-1]->ice[0]*tvals.sg_ice))-tvals.gwH)*
		 ((elp->nps[0]->uo[0]-(elp->nps[0]->ice[0]*tvals.sg_ice))-tvals.gwH);
	if(mult<0.0)
		return(1);
		
	return(0);
}

int     get_mixedgps(elp,g)
struct  element *elp;
struct  gausspts        g[];

{
	struct  node    left[3],right[3],np[4];
	int     i,nipts,nleft,nright,fac,nnodes,ngauss,n,nvar=0;
	
	double  xl[MNSF],yl[MNSF],dx,dy,m,tol1=10E-7,xtr[3],ytr[3],f[3][3];
	double  areatr(),Area,weight,level;
	
	level=tvals.gwH;


	switch(elp->gtype) {
/* to skip the midside node for h.o. triangles*/
		case    210:
			fac = 1;
			Area = 0.5;
			nnodes=get_3nodes(np);
			break;
		case    211:
			fac = 1;
			Area = 4.0;
			nnodes=get_4nodes(np);
			break;
		case    111:
			fac = 1;
			break;
		case    220:
			fac = 2;
			Area = 0.5;
			nnodes=get_3nodes(np);
			break;
		case    221:
			fac = 2;
			Area = 4.0;
			nnodes=get_4nodes(np);
			break;
		case    121:
			fac = 2;
			break;
	}
/*
	if((nipts=get_cline(elp,tvals.gwH,0,0,0.02,xl,yl,0))!=2) {
		printf("no intersection points found in element # %d\n", elp->n);
		return(0);
	}
*/
	if((nipts=find_intersection(elp,np,xl,yl,nnodes,fac,level,0,&nright,&nleft,right,left))!=2) {
		printf("no intersection points found in element # %d\n", elp->n);
		return(0);
	}





/*
	dx=xl[1]-xl[0];
	dy=yl[1]-yl[0];
	
	nleft=nright=0;
	
	for(i=0;i<nnodes;i++) {
		if(fabs(dx)<tol1) {
			if(dy>0.0) {
				if(np[i].x[0]<xl[0]) {
					left[nleft]=np[i];
					nleft++;
				}
				else {
					right[nright]=np[i];
					nright++;
				}
			}
			else {
				if(np[i].x[0]>xl[0]) {
					left[nleft]=np[i];
					nleft++;
				}
				else {
					right[nright]=np[i];
					nright++;
				}
			}
		}
		else {
			m=dy/dx;
			if(dx>0.0) {
				if((np[i].x[1]-m*np[i].x[0])>(yl[0]-m*xl[0])) {
					left[nleft]=np[i];
					nleft++;
				}
				else {
					right[nright]=np[i];
					nright++;
				}
			}
			else {
				if((np[i].x[1]-m*np[i].x[0])<(yl[0]-m*xl[0])) {
					left[nleft]=np[i];
					nleft++;
				}
				else {
					right[nright]=np[i];
					nright++;
				}
			}
		}
	}
*/
	ngauss = 0;
	switch(nright) {
		case    1:
			xtr[0]=xl[0];
			ytr[0]=yl[0];
			xtr[1]=right[0].x[0];
			ytr[1]=right[0].x[1];
			xtr[2]=xl[1];
			ytr[2]=yl[1];
			ngauss=find_gps(xtr,ytr,g,Area,ngauss);

			break;
		
		case    2:
			xtr[0]=xl[0];
			ytr[0]=yl[0];
			xtr[2]=xl[1];
			ytr[2]=yl[1];

			n=find_bestnode(xl[0],yl[0],right[0].x[0],right[0].x[1],right[1].x[0],right[1].x[1],xl[1],yl[1]);
			xtr[1]=right[n].x[0];
			ytr[1]=right[n].x[1];

			ngauss=find_gps(xtr,ytr,g,Area,ngauss);
/*now find the 2nd triangle*/
			
			xtr[1]=right[0].x[0];
			xtr[2]=right[1].x[0];
			ytr[1]=right[0].x[1];
			ytr[2]=right[1].x[1];

			if(n==0) {
				xtr[0]=xl[1];
				ytr[0]=yl[1];
			}
			else {
				xtr[0]=xl[0];
				ytr[0]=yl[0];
			}
			ngauss=find_gps(xtr,ytr,g,Area,ngauss);
			
			break;
			
		case    3:
/*quality of triangles is not checked
//first triangle:*/
			xtr[0]=xl[1];
			ytr[0]=yl[1];
			xtr[1]=xl[0];
			ytr[1]=yl[0];
			xtr[2]=right[0].x[0];
			ytr[2]=right[0].x[1];
			
			ngauss=find_gps(xtr,ytr,g,Area,ngauss);

/*2nd triangle:*/
			xtr[0]=right[0].x[0];
			xtr[1]=right[1].x[0];
			xtr[2]=xl[1];
			ytr[0]=right[0].x[1];
			ytr[1]=right[1].x[1];
			ytr[2]=yl[1];
			
			ngauss=find_gps(xtr,ytr,g,Area,ngauss);

/*3rd triangle:*/
			xtr[0]=right[1].x[0];
			xtr[1]=right[2].x[0];
			xtr[2]=xl[1];
			ytr[0]=right[1].x[1];
			ytr[1]=right[2].x[1];
			ytr[2]=yl[1];

			ngauss=find_gps(xtr,ytr,g,Area,ngauss);
			
			break;
			

	}
	
	switch(nleft) {
		case    1:
			xtr[0]=xl[0];
			ytr[0]=yl[0];
			xtr[1]=xl[1];
			ytr[1]=yl[1];
			xtr[2]=left[0].x[0];
			ytr[2]=left[0].x[1];
			ngauss=find_gps(xtr,ytr,g,Area,ngauss);

			break;
		
		case    2:
			xtr[0]=xl[0];
			ytr[0]=yl[0];
			xtr[1]=xl[1];
			ytr[1]=yl[1];
			
			n=find_bestnode(xl[0],yl[0],xl[1],yl[1],left[0].x[0],left[0].x[1],left[1].x[0],left[1].x[1]);
			xtr[2]=left[n].x[0];
			ytr[2]=left[n].x[1];
			

			find_gps(xtr,ytr,g,Area,ngauss);
			ngauss+=3;
			
/*now find the 2nd triangle*/
			
			xtr[1]=left[0].x[0];
			xtr[2]=left[1].x[0];
			ytr[1]=left[0].x[1];
			ytr[2]=left[1].x[1];
			
			xtr[0]=xl[n];
			ytr[0]=yl[n];
						

			find_gps(xtr,ytr,g,Area,ngauss);
			ngauss+=3;
			break;

		case    3:
/*quality of triangles is not checked
//first triangle:*/
			xtr[0]=xl[0];
			ytr[0]=yl[0];
			xtr[1]=xl[1];
			ytr[1]=yl[1];
			xtr[2]=left[0].x[0];
			ytr[2]=left[0].x[1];
			
			ngauss=find_gps(xtr,ytr,g,Area,ngauss);

/*2nd triangle:*/
			xtr[0]=left[0].x[0];
			xtr[1]=left[1].x[0];
			xtr[2]=xl[0];
			ytr[0]=left[0].x[1];
			ytr[1]=left[1].x[1];
			ytr[2]=yl[0];
			
			ngauss=find_gps(xtr,ytr,g,Area,ngauss);

/*3rd triangle:*/
			xtr[0]=left[1].x[0];
			xtr[1]=left[2].x[0];
			xtr[2]=xl[0];
			ytr[0]=left[1].x[1];
			ytr[1]=left[2].x[1];
			ytr[2]=yl[0];

			ngauss=find_gps(xtr,ytr,g,Area,ngauss);
			
			break;

			

	}       


	return(ngauss);
}

int     get_3nodes(np)
struct  node    np[];
{
	np[0].x[0] = 1.0;
	np[0].x[1] = 0.0;
	np[1].x[0] = 0.0;
	np[1].x[1] = 1.0;
	np[2].x[0] = 0.0;
	np[2].x[1] = 0.0;
	return(3);
}

int     get_4nodes(np)
struct  node    np[];
{
	np[0].x[0] = 1.0;
	np[0].x[1] = -1.0;
	np[1].x[0] = 1.0;
	np[1].x[1] = 1.0;
	np[2].x[0] = -1.0;
	np[2].x[1] = 1.0;
	np[3].x[0] = -1.0;
	np[3].x[1] = -1.0;
	return(4);
}

double  areatr(x1,y1,x2,y2,x3,y3)
double  x1,y1,x2,y2,x3,y3;

{
	double  ar;
	
	ar = 0.5*((x2-x1)*(y3-y1)
		- (y2-y1)*(x3-x1));
	return(ar);
}

int     find_gps(xtr,ytr,g,Area,ngauss)
double  xtr[],ytr[],Area;
struct  gausspts        g[];
int     ngauss;
{
	int     i,j;
	double  weight,areatr(),f[3][3];
	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
			if(i==j)
				f[i][j]=2.0/3.0;
			else
				f[i][j]=1.0/6.0;
		}
	}
	
	weight = areatr(xtr[0],ytr[0],xtr[1],ytr[1],xtr[2],ytr[2])/3.0;
	for(i=ngauss;i<(ngauss+3);i++) {
		g[i].x=xtr[0]*f[0][i-ngauss]+xtr[1]*f[1][i-ngauss]+xtr[2]*f[2][i-ngauss];
		g[i].y=ytr[0]*f[0][i-ngauss]+ytr[1]*f[1][i-ngauss]+ytr[2]*f[2][i-ngauss];
		g[i].z=0.0;
		g[i].w=weight;
	}
	
	return(ngauss+3);
}

int     find_bestnode(xa,xb,ya,yb,x0,y0,x1,y1)
double  xa,xb,ya,yb,x0,y0,x1,y1;
{
	double  alpha0,alpha1,dist();
	
	alpha0=areatr(xa,xb,ya,yb,x0,y0)/(dist(xa,ya,xb,yb)+dist(xb,yb,x0,y0)+dist(x0,y0,xa,ya));
	alpha1=areatr(xa,xb,ya,yb,x1,y1)/(dist(xa,ya,xb,yb)+dist(xb,yb,x1,y1)+dist(x1,y1,xa,ya));
	
	if(alpha0>alpha1)
		return(0);
	else
		return(1);
}

double  dist(x1,y1,x2,y2)
double  x1,y1,x2,y2;
{
	double  d;
	d=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	return(d);
}


int     find_intersection(elp,np,xl,yl,nnodes,fac,level,nvar,nright,nleft,right,left)
struct  element *elp;
struct  node    np[],right[],left[];
double  xl[],yl[],level;
int             nnodes,fac,nvar,*nright,*nleft;
{
	int     i,ncut,nlast,side,j,nr,nl;
	double  mult,factor;
	ncut=0;
	nr=0;
	side=0;
	left[0]=np[0];
	nl=1;
	
	for(i=0;i<(nnodes*fac);i++) {
		if(i==(nnodes-1)*fac)
			nlast=0;
		else
			nlast=i+fac;
	
		mult=((elp->nps[i]->uo[nvar]-(elp->nps[i]->ice[0]*tvals.sg_ice))-level)*((elp->nps[nlast]->uo[nvar]-(elp->nps[nlast]->ice[0]*tvals.sg_ice))-level);

		if(mult>0.0){
/*no intersection*/
			if(side==0) {
/*we are on the left side*/
				if(nr==0){
/*this is our first time on the left side*/
					left[nl]=np[nlast/fac];
					nl++;
				}
				else{
/*we are returning to the left side -> insert*/
					if(nlast!=0) {
						for(j=nl;j>0;j--)
							left[j]=left[j-1];
						left[0]=np[nlast/fac];
						nl++;
					}
				}
			}
			else{
/* we are on the right side*/
				right[nr]=np[nlast/fac];
				nr++;
			}
		}       
				
				
		if(mult==0.0) {
			if((elp->nps[i]->uo[nvar]-(elp->nps[i]->ice[nvar]*tvals.sg_ice))==level) {
				xl[ncut]=np[i/fac].x[0];
				yl[ncut]=np[i/fac].x[1];
				ncut++;

				side=abs(1-side);
				if(i==0)
					nl--;
				if(nlast!=0){
					if(side==0) {
						if(nr==0){
							left[nl]=np[nlast/fac];
							nl++;
						}
						else {
							if(nlast!=0) {
								for(j=nl;j>0;j--)
									left[j]=left[j-1];
								left[0]=np[nlast/fac];
								nl++;
							}
						}
					}
					else {
/*we are on the right side*/
						right[nr]=np[nlast/fac];
						nr++;
					}
				}
			}
		}
		if(mult<0.0) {
			factor=fabs((elp->nps[i]->uo[nvar]-(elp->nps[i]->ice[nvar]*tvals.sg_ice))-level)/
				(fabs((elp->nps[i]->uo[nvar]-(elp->nps[i]->ice[nvar]*tvals.sg_ice))-level)+
				fabs((elp->nps[nlast]->uo[nvar]-(elp->nps[nlast]->ice[0]*tvals.sg_ice))-level));
			xl[ncut]=np[i/fac].x[0]+(np[nlast/fac].x[0]-np[i/fac].x[0])*factor;
			yl[ncut]=np[i/fac].x[1]+(np[nlast/fac].x[1]-np[i/fac].x[1])*factor;
			
			ncut++;
			side=abs(side-1);
			if(side==0) {
/* we must be returning to left -> insert*/
				if(nlast!=0) {
					for(j=nl;j>0;j--)
						left[j]=left[j-1];
					left[0]=np[nlast/fac];
					nl++;
				}
			}
			else {
/* we have switched from left to right*/
				right[nr]=np[nlast/fac];
				nr++;
			}
		}
	}       
	
	*nright=nr;
	*nleft=nl;

	return(ncut);
}


int     check_bmixed(belp)
struct  belement        *belp;
{
	int     i;
	double  mult;
	
	for(i=0;i<(belp->nnds-1);i++) {
		mult = ((belp->nps[i]->uo[0]-(belp->nps[i]->ice[0]*tvals.sg_ice))-tvals.gwH)*
				((belp->nps[i+1]->uo[0]-(belp->nps[i+1]->ice[0]*tvals.sg_ice))-tvals.gwH);
		if(mult<0.0)
			return(1);
	}
		
	return(0);
}

int     get_boundarygps(belp,g)
struct  belement        *belp;
struct  gausspts        g[];

{
	int     ngauss=0;
	double  r[2],ri;
	
	ri=2.0*fabs((belp->nps[0]->uo[0]-(belp->nps[0]->ice[0]*tvals.sg_ice))-tvals.gwH)/
				(fabs((belp->nps[0]->uo[0]-(belp->nps[0]->ice[0]*tvals.sg_ice))-tvals.gwH)+
				fabs((belp->nps[(belp->nnds-1)]->uo[0]-(belp->nps[belp->nnds-1]->ice[0]*tvals.sg_ice))-tvals.gwH));
	
	r[0]=-1.0;
	r[1]=ri-1.0;
	ngauss=find_bgps(r,g,ngauss);
	r[0]=ri-1.0;
	r[1]=1.0;
	ngauss=find_bgps(r,g,ngauss);
	
	return(ngauss);
}

int     find_bgps(r,g,ngauss)
double  r[];
struct  gausspts        g[];
int             ngauss;
{
	double  f[2];
	
	f[0]=0.5*(1.0+1.0/sqrt(3.0));
	f[1]=0.5*(1.0-1.0/sqrt(3.0));
	
	g[ngauss].x=r[0]*f[0]+r[1]*f[1];
	g[ngauss].y=0.0;
	g[ngauss].z=0.0;
	g[ngauss].w=(r[1]-r[0])/2.0;
	
	ngauss++;
	
	g[ngauss].x=r[0]*f[1]+r[1]*f[0];
	g[ngauss].y=0.0;
	g[ngauss].z=0.0;
	g[ngauss].w=(r[1]-r[0])/2.0;
	
	ngauss++;
	
	return(ngauss);
}
	
