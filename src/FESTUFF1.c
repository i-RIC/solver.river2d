//#pragma options cpluscmt	
#include "Fe_PS.h"


extern struct control  N ;


int	shqu_2l1(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,s,mr,ms ;
	
	r = x->x ;
	s = x->y ;
	if( (r<-1.00001) || (r>1.00001) || (s<-1.00001) || (s>1.00001) ) return(4) ;
	
	mr = 1.0 - r ;
	ms = 1.0 - s ;
	
	p->dof = 4 ;
	
	p->f[0] = r * ms ;
	p->f[1] = r * s ;
	p->f[2] = mr * s ;
	p->f[3] = mr * ms ;
	
	if(x->w < 0.0)
		return(0) ;
	
	p->dfdr[0] =  ms ;
	p->dfdr[1] =  s ;
	p->dfdr[2] = -s ;
	p->dfdr[3] = -ms ;
	
	p->dfds[0] = -r ;
	p->dfds[1] =  r ;
	p->dfds[2] =  mr ;
	p->dfds[3] = -mr ;
	
	return(0) ;
}


int	shqu9_2l2(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,s,mr,ms,pr,ps,t,tr,ts,pmr,pms ;
	
	r = x->x ;
	s = x->y ;
	if( (r<-1.00001) || (r>1.00001) || (s<-1.00001) || (s>1.00001) ) return(9) ;
	
	t = 0.25 ;
	tr = 2.0 * r ;
	ts = 2.0 * s ;
	mr = 1.0 - r ;
	ms = 1.0 - s ;
	pr = 1.0 + r ;
	ps = 1.0 + s ;
	pmr = pr * mr ;
	pms = ps * ms ;
	
	p->dof = 9 ;
	
	p->f[8] = -t * pr * ms * r * s ;
	p->f[1] = t * pr * ps * r * s ;
	p->f[4] = -t * r * mr * s * ps ;
	p->f[6] = t * r * mr * s * ms ;
	p->f[0] = 0.5 * r * pr * pms ;
	p->f[2] = 0.5 * pmr * s * ps ;
	p->f[5] = -0.5 * r * mr * pms ;
	p->f[7] = -0.5 * pmr * s * ms ;
	p->f[3] = pmr * pms ;
	
	if(x->w < 0.0)
		return(0) ;
	
	p->dfdr[8] = -t * s * ms * (tr + 1.0) ;
	p->dfdr[1] = t * s * ps * (1.0 + tr) ;
	p->dfdr[4] = -t * s * ps * (1.0 - tr) ;
	p->dfdr[6] = t * s * ms * (1.0 - tr) ;
	p->dfdr[0] =  0.5 * pms * (1.0 + tr) ;
	p->dfdr[2] = -r * s * ps ;
	p->dfdr[5] = -0.5 * (1.0 - tr) * pms ;
	p->dfdr[7] =  r * s * ms ;
	p->dfdr[3] = -2.0 * r * pms ;
	
	p->dfds[8] = -t * r * pr * (1.0 - ts) ;
	p->dfds[1] = t * r * pr * (1.0 + ts) ;
	p->dfds[4] = -t * r * mr * (1.0 + ts) ;
	p->dfds[6] = t * r * mr * (1.0 - ts) ;
	p->dfds[0] = -s * r * pr ;
	p->dfds[2] = 0.5 * pmr * (1.0 + ts) ;
	p->dfds[5] = s * r * mr ;
	p->dfds[7] = -0.5 * pmr * (1.0 - ts) ;
	p->dfds[3] = -2.0 * s * pmr ;
	
	return(0) ;
}


int	shqu8_2l2(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,s,mr,ms,pr,ps,t,tr,ts,pmr,pms ;
	
	r = x->x ;
	s = x->y ;
	if( (r<-1.00001) || (r>1.00001) || (s<-1.00001) || (s>1.00001) ) return(8) ;
	
	t = 0.25 ;
	tr = 2.0 * r ;
	ts = 2.0 * s ;
	mr = 1.0 - r ;
	ms = 1.0 - s ;
	pr = 1.0 + r ;
	ps = 1.0 + s ;
	pmr = pr * mr ;
	pms = ps * ms ;
	
	p->dof = 8 ;
	
	p->f[7] = -0.5 * ms * r * s ;
	p->f[1] = t * pr * ps * r * s ;
	p->f[4] = -0.5 * r * mr * s ;
	p->f[0] = 0.5 * r * pr * pms ;
	p->f[2] = 0.5 * pmr * s * ps ;
	p->f[5] = -0.5 * r * mr * ms ;
	p->f[6] = -0.5 * mr * s * ms ;
	p->f[3] = pmr * pms ;
	
	if(x->w < 0.0)
		return(0) ;
	
	p->dfdr[7] = -0.5 * ms * s  ;
	p->dfdr[1] = t * s * ps * (1.0 + tr) ;
	p->dfdr[4] = -0.5 * s * (1.0 - tr) ;
	p->dfdr[0] =  0.5 * pms * (1.0 + tr) ;
	p->dfdr[2] = -r * s * ps ;
	p->dfdr[5] = -0.5 * (1.0 - tr) * ms ;
	p->dfdr[6] =  0.5 * s * ms ;
	p->dfdr[3] = -2.0 * r * pms ;
	
	p->dfds[7] = 0.5 * r * (1.0 - ts) ;
	p->dfds[1] = t * r * pr * (1.0 + ts) ;
	p->dfds[4] = -0.5 * r * mr ;
	p->dfds[0] = -s * r * pr ;
	p->dfds[2] = 0.5 * pmr * (1.0 + ts) ;
	p->dfds[5] = 0.5 * r * mr ;
	p->dfds[6] = -0.5 * mr * (1.0 - ts) ;
	p->dfds[3] = -2.0 * s * pmr ;
	
	return(0) ;
}

int	shqu6_2l2(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,s,mr,ms,pr,ps,t,tr,ts,pmr,pms ;
	
	r = x->x ;
	s = x->y ;
	if( (r<-1.00001) || (r>1.00001) || (s<-1.00001) || (s>1.00001) ) return(6) ;
	
	t = 0.25 ;
	tr = 2.0 * r ;
	ts = 2.0 * s ;
	mr = 1.0 - r ;
	ms = 1.0 - s ;
	pr = 1.0 + r ;
	ps = 1.0 + s ;
	pmr = pr * mr ;
	pms = ps * ms ;
	
	p->dof = 6 ;
	
	p->f[1] = 0.5 * pr * r * s ;
	p->f[4] = -0.5 * r * mr * s ;
	p->f[0] = 0.5 * r * pr * ms ;
	p->f[2] = pmr * s  ;
	p->f[5] = -0.5 * r * mr * ms ;
	p->f[3] = pmr * ms ;
	
	if(x->w < 0.0)
		return(0) ;
	
	p->dfdr[1] =  0.5 * s * (1.0 + tr) ;
	p->dfdr[4] = -0.5 * s * (1.0 - tr) ;
	p->dfdr[0] =  0.5 * ms * (1.0 + tr) ;
	p->dfdr[2] = -2.0 * r * s ;
	p->dfdr[5] = -0.5 * (1.0 - tr) * ms ;
	p->dfdr[3] = -2.0 * r * ms ;
	
	p->dfds[1] = 0.5 * pr * r ;
	p->dfds[4] = -0.5 * r * mr ;
	p->dfds[0] = -0.5 * r * pr ;
	p->dfds[2] = pmr ;
	p->dfds[5] = 0.5 * r * mr ;
	p->dfds[3] = -pmr ;
	
	return(0) ;
}



