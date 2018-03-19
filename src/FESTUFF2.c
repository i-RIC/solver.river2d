//#pragma options cpluscmt	
/*			the shape functions for the cubic upwinding functions   239		*/
#include "Fe_PS.h"


extern struct control  N ;



int shcu2_2x1(x,p)
struct gausspts *x;
struct shapefuncs	*p;

{
	double r, s, t, rm, sm, rp, sp, rp3;

	r = x->x ;
	s = x->y ;
	if ( (r<-3.00001) || (r>1.00001) || (s<-1.00001) || (s>1.00001) ) return(6) ;

	t = 0.125 ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	sp = 1.0 + s ;
	sm = 1.0 - s ;
	rp3 = 3.0 + r ;

	p->dof = 6 ;

	p->f[3] = t * rp3 * sm * rm ;
	p->f[0] = t / 2.0 * rp3 * rp * sm ;
	p->f[1] = t / 2.0 * rp * sp * rp3 ;
	p->f[2] = t * rm * sp * rp3 ;
	p->f[4] = t / -2.0 * rp * sp * rm ;
	p->f[5] = t / -2.0 * rp * sm * rm ;

	if (x->w < 0.0)
		return(0) ;

	p->dfdr[3] = -2.0 * t * sm * rp ;
	p->dfdr[0] = t * sm * ( r + 2.0 ) ;
	p->dfdr[1] = t * sp * ( r + 2.0 ) ;
	p->dfdr[2] = -2.0 * t * sp * rp ;
	p->dfdr[4] = t * sp * r ;
	p->dfdr[5] = t * sm * r ;

	p->dfds[3] = -t * rp3 * rm ;
	p->dfds[0] = -t / 2.0 * rp3 * rp ;
	p->dfds[1] = t / 2.0 * rp3 * rp ;
	p->dfds[2] = t * rp3 * rm ;
	p->dfds[4] = -t / 2.0 * rp * rm ;
	p->dfds[5] = t / 2.0 * rp * rm ;

	return (0);
}


int shcu3_2x2(x,p)
struct gausspts *x;
struct shapefuncs	*p;

{
	double r, s, t, rm, sm, rp, sp, rp3, sp3;

	r = x->x ;
	s = x->y ;
	if ( (r<-3.00001) || (r>1.00001) || (s<-3.00001) || (s>1.00001) ) return(8) ;

	t = 0.0625 ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	sp = 1.0 + s ;
	sm = 1.0 - s ;
	rp3 = 3.0 + r ;
	sp3 = 3.0 + s ;

	p->dof = 8 ;

	p->f[3] = t * rp3 * sm * rm * sp3 ;
	p->f[0] = t / 2.0 * rp3 * rp * sm * sp3 ;
	p->f[1] = t / 4.0 * rp * sp * rp3 * sp3 ;
	p->f[2] = t / 2.0 * rm * sp * rp3 * sp3 ;
	p->f[4] = -t * rp * sp * rm ;
	p->f[5] = -t * rp * sm * rm ;
	p->f[6] = -t * sm * rm * sp ;
	p->f[7] = -t * rp * sm * sp ;


	if (x->w < 0.0)
		return(0) ;

	p->dfdr[3] = -2.0 * t * sm * rp * sp3;
	p->dfdr[0] = t * sm * ( r + 2.0 ) * sp3;
	p->dfdr[1] = t / 2.0 * sp * ( r + 2.0 ) * sp3;
	p->dfdr[2] = -t * sp * rp * sp3;
	p->dfdr[4] = 2. * t * sp * r ;
	p->dfdr[5] = 2. * t * sm * r ;
	p->dfdr[6] = t * sm * sp;
	p->dfdr[7] = -t * sm * sp;


	p->dfds[3] = -2.0 * t * sp * rp3 * rm ;
	p->dfds[0] = -t * sp * rp3 * rp ;
	p->dfds[1] = t / 2.0 * ( s + 2.0 ) * rp3 * rp ;
	p->dfds[2] = t * ( s + 2.0 ) * rp3 * rm ;
	p->dfds[4] = -t * rp * rm ;
	p->dfds[5] = t * rp * rm ;
	p->dfds[6] = t * 2.0 * s * rm ;
	p->dfds[7] = t * 2.0 * s * rp ;

	return (0);
}

int shcu4_2x2(x,p)
struct gausspts *x;
struct shapefuncs	*p;

{
	double r, s, t, rm, sm, rp, sp, rp3, sp3;

	r = x->x ;
	s = x->y ;
	if ( (r<-3.00001) || (r>1.00001) || (s<-3.00001) || (s>1.00001) ) return(9) ;

	t = 0.0625 ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	sp = 1.0 + s ;
	sm = 1.0 - s ;
	rp3 = 3.0 + r ;
	sp3 = 3.0 + s ;

	p->dof = 9 ;

	p->f[3] = t * rp3 * sm * rm * sp3 ;
	p->f[0] = t / 2.0 * rp3 * rp * sm * sp3 ;
	p->f[1] = t / 4.0 * rp * sp * rp3 * sp3 ;
	p->f[2] = t / 2.0 * rm * sp * rp3 * sp3 ;
	p->f[4] = t / -4.0 * rp * sp * rm * sp3 ;
	p->f[5] = t / -2.0 * rp * sm * rm * sp3 ;
	p->f[7] = t / -2.0 * rp3 * sm * rm * sp ;
	p->f[8] = t / -4.0 * rp * sm * rp3 * sp ;
	p->f[6] = t / 4.0 * rp * sm * rm * sp ;


	if (x->w < 0.0)
		return(0) ;

	p->dfdr[3] = -2.0 * t * sm * rp * sp3;
	p->dfdr[0] = t * sm * ( r + 2.0 ) * sp3;
	p->dfdr[1] = t / 2.0 * sp * ( r + 2.0 ) * sp3;
	p->dfdr[2] = -t * sp * rp * sp3;
	p->dfdr[4] = t / 2.0 * sp * r * sp3;
	p->dfdr[5] = t * sm * r * sp3;
	p->dfdr[7] = t * sm * rp * sp;
	p->dfdr[8] = -t / 2.0 * sm * ( r + 2.0 ) * sp;
	p->dfdr[6] = -t / 2.0 * sm * r * sp;


	p->dfds[3] = -2.0 * t * sp * rp3 * rm ;
	p->dfds[0] = -t * sp * rp3 * rp ;
	p->dfds[1] = t / 2.0 * ( s + 2.0 ) * rp3 * rp ;
	p->dfds[2] = t * ( s + 2.0 ) * rp3 * rm ;
	p->dfds[4] = -t / 2.0 * ( s + 2.0 ) * rp * rm ;
	p->dfds[5] = t * sp * rp * rm ;
	p->dfds[7] = t * s * rp3 * rm ;
	p->dfds[8] = t / 2.0 * s * rp * rp3 ;
	p->dfds[6] = -t / 2.0 * s * rp * rm ;

	return (0);
}

int shcu3_3x1(x,p)
struct gausspts *x;
struct shapefuncs	*p;

{
	double r, s, t, rm, sm, rp, sp, rp3, rp5, r2p8 ;

	r = x->x ;
	s = x->y ;
	if ( (r<-5.00001) || (r>1.00001) || (s<-1.00001) || (s>1.00001) ) return(8) ;

	t = 1.0 / 32. ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	sp = 1.0 + s ;
	sm = 1.0 - s ;
	rp3 = 3.0 + r ;
	rp5 = 5.0 + r ;

	p->dof = 8 ;

	p->f[3] = t * rp3 * sm * rm * rp5 ;
	p->f[0] = t / 3.0 * rp3 * rp * sm * rp5 ;
	p->f[1] = t / 3.0 * rp * sp * rp3 * rp5 ;
	p->f[2] = t * rm * sp * rp3 * rp5 ;
	p->f[4] = -t * rp * sp * rm * rp5 ;
	p->f[5] = t / 3.0 * rp * sp * rm * rp3 ;
	p->f[6] = t / 3.0 * rp * sm * rm * rp3 ;
	p->f[7] = -t * rp * sm * rm * rp5 ;

	if (x->w < 0.0)
		return(0) ;

	r2p8 =	2.0 * r + 8.0 ;

	p->dfdr[3] = t * sm * ( (-rp5 * rp3) + rm * r2p8 ) ;
	p->dfdr[0] = t / 3.0 * sm * ( ( rp5 * rp3) + rp * r2p8 ) ;
	p->dfdr[1] = t / 3.0 * sp * ( ( rp5 * rp3) + rp * r2p8 ) ;
	p->dfdr[2] = t * sp * ( (-rp5 * rp3) + rm * r2p8 ) ;
	p->dfdr[4] = -t * sp * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[5] = t / 3.0 * sp * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[6] = t / 3.0 * sm * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[7] = -t * sm * ( -2.0 * r * rp5 + rp * rm ) ;

	p->dfds[3] = -t * rp3 * rm * rp5 ;
	p->dfds[0] = -t / 3.0 * rp3 * rp * rp5 ;
	p->dfds[1] = t / 3.0 * rp * rp3 * rp5 ;
	p->dfds[2] = t * rm * rp3 * rp5 ;
	p->dfds[4] = -t * rp * rm * rp5 ;
	p->dfds[5] = t / 3.0 * rp * rm * rp3 ;
	p->dfds[6] = -t / 3.0 * rp * rm * rp3 ;
	p->dfds[7] = t * rp * rm * rp5 ;

	return (0);
}

int shcu4_3x2(x,p)
struct gausspts *x;
struct shapefuncs	*p;

{
	double r, s, t, rm, sm, rp, sp, rp3, rp5, r2p8, sp3 ;

	r = x->x ;
	s = x->y ;
	if ( (r<-5.00001) || (r>1.00001) || (s<-3.00001) || (s>1.00001) ) return(10) ;

	t = 1.0 / 64. ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	sp = 1.0 + s ;
	sm = 1.0 - s ;
	rp3 = 3.0 + r ;
	rp5 = 5.0 + r ;
	sp3 = 3.0 + s ;

	p->dof = 10 ;

	p->f[3] = t * rp3 * sm * rm * rp5 * sp3 ;
	p->f[0] = t / 3.0 * rp3 * rp * sm * rp5 * sp3 ;
	p->f[1] = t / 6.0 * rp * sp * rp3 * rp5 * sp3 ;
	p->f[2] = t / 2.0 * rm * sp * rp3 * rp5 * sp3 ;
	p->f[4] = -2.0 * t * rp * sp * rm * rp5 ;
	p->f[5] = 2.0 * t / 3.0 * rp * sp * rm * rp3 ;
	p->f[6] = 2.0 * t / 3.0 * rp * sm * rm * rp3 ;
	p->f[7] = -2.0 * t * rp * sm * rm * rp5 ;
	p->f[8] = -4.0 * t * sp * rm * sm ;
	p->f[9] = -4.0 * t * sp * rp * sm ;

	if (x->w < 0.0)
		return(0) ;

	r2p8 =	2.0 * r + 8.0 ;

	p->dfdr[3] = t * sm * sp3 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfdr[0] = t / 3.0 * sm * sp3 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfdr[1] = t / 6.0 * sp * sp3 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfdr[2] = t / 2.0 * sp * sp3 * ( (-rp5 * rp3) + rm * r2p8 ) ;
	p->dfdr[4] = -2.0 * t * sp * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[5] = 2.0 * t / 3.0 * sp * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[6] = 2.0 * t / 3.0 * sm * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[7] = -2.0 * t * sm * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[8] = 4.0 * t * sm * sp ;
	p->dfdr[9] = -4.0 * t * sm * sp ;

	p->dfds[3] = -t * 2.0 * rp3 * sp * rm * rp5 ;
	p->dfds[0] = -t * 2.0 / 3.0 * rp3 * rp * sp * rp5 ;
	p->dfds[1] = t / 3.0 * rp * rp3 * rp5 * ( s + 2.0 ) ;
	p->dfds[2] = t * rm * rp3 * rp5 * ( s + 2.0 ) ;
	p->dfds[4] = -2.0 * t * rp * rm * rp5 ;
	p->dfds[5] = 2.0 * t / 3.0 * rp * rm * rp3 ;
	p->dfds[6] = -2.0 * t / 3.0 * rp * rm * rp3 ;
	p->dfds[7] = 2.0 * t * rp * rm * rp5 ;
	p->dfds[8] = 8.0 * t * s * rm ;
	p->dfds[9] = 8.0 * t * s * rp ;

	return (0);
}

int sh4_3x2b(x,p)
struct gausspts *x;
struct shapefuncs	*p;

{
	double r, s, t, rm, sm, rp, sp, rp3, rp5, r2p8, sp3 ;

	r = x->y ;
	s = x->x ;
	if ( (r<-5.00001) || (r>1.00001) || (s<-3.00001) || (s>1.00001) ) return(10) ;

	t = 1.0 / 64. ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	sp = 1.0 + s ;
	sm = 1.0 - s ;
	rp3 = 3.0 + r ;
	rp5 = 5.0 + r ;
	sp3 = 3.0 + s ;

	p->dof = 10 ;

	p->f[3] = t * rp3 * sm * rm * rp5 * sp3 ;
	p->f[2] = t / 3.0 * rp3 * rp * sm * rp5 * sp3 ;
	p->f[1] = t / 6.0 * rp * sp * rp3 * rp5 * sp3 ;
	p->f[0] = t / 2.0 * rm * sp * rp3 * rp5 * sp3 ;
	p->f[4] = -2.0 * t * rp * sp * rm * rp5 ;
	p->f[5] = 2.0 * t / 3.0 * rp * sp * rm * rp3 ;
	p->f[6] = 2.0 * t / 3.0 * rp * sm * rm * rp3 ;
	p->f[7] = -2.0 * t * rp * sm * rm * rp5 ;
	p->f[8] = -4.0 * t * sp * rm * sm ;
	p->f[9] = -4.0 * t * sp * rp * sm ;

	if (x->w < 0.0)
		return(0) ;

	r2p8 =	2.0 * r + 8.0 ;

	p->dfds[3] = t * sm * sp3 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfds[2] = t / 3.0 * sm * sp3 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfds[1] = t / 6.0 * sp * sp3 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfds[0] = t / 2.0 * sp * sp3 * ( (-rp5 * rp3) + rm * r2p8 ) ;
	p->dfds[4] = -2.0 * t * sp * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfds[5] = 2.0 * t / 3.0 * sp * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfds[6] = 2.0 * t / 3.0 * sm * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfds[7] = -2.0 * t * sm * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfds[8] = 4.0 * t * sm * sp ;
	p->dfds[9] = -4.0 * t * sm * sp ;

	p->dfdr[3] = -t * 2.0 * rp3 * sp * rm * rp5 ;
	p->dfdr[2] = -t * 2.0 / 3.0 * rp3 * rp * sp * rp5 ;
	p->dfdr[1] = t / 3.0 * rp * rp3 * rp5 * ( s + 2.0 ) ;
	p->dfdr[0] = t * rm * rp3 * rp5 * ( s + 2.0 ) ;
	p->dfdr[4] = -2.0 * t * rp * rm * rp5 ;
	p->dfdr[5] = 2.0 * t / 3.0 * rp * rm * rp3 ;
	p->dfdr[6] = -2.0 * t / 3.0 * rp * rm * rp3 ;
	p->dfdr[7] = 2.0 * t * rp * rm * rp5 ;
	p->dfdr[8] = 8.0 * t * s * rm ;
	p->dfdr[9] = 8.0 * t * s * rp ;

	return (0);
}

int shcu5_3x2(x,p)
struct gausspts *x;
struct shapefuncs	*p;

{
	double r, s, t, rm, sm, rp, sp, rp3, rp5, r2p8, sp3 ;

	r = x->x ;
	s = x->y ;
	if ( (r<-5.00001) || (r>1.00001) || (s<-3.00001) || (s>1.00001) ) return(11) ;

	t = 1.0 / 64. ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	sp = 1.0 + s ;
	sm = 1.0 - s ;
	rp3 = 3.0 + r ;
	rp5 = 5.0 + r ;
	sp3 = 3.0 + s ;

	p->dof = 11 ;

	p->f[3] = t * rp3 * sm * rm * rp5 * sp3 ;
	p->f[0] = t / 3.0 * rp3 * rp * sm * rp5 * sp3 ;
	p->f[1] = t / 6.0 * rp * sp * rp3 * rp5 * sp3 ;
	p->f[2] = t / 2.0 * rm * sp * rp3 * rp5 * sp3 ;
	p->f[4] = -t / 2.0 * rp * sp * rm * rp5 * sp3 ;
	p->f[5] = 2.0 * t / 3.0 * rp * sp * rm * rp3 ;
	p->f[6] = 2.0 * t / 3.0 * rp * sm * rm * rp3 ;
	p->f[7] = -t * rp * sm * rm * rp5 * sp3 ;
	p->f[8] = t * rp * sp * rm * sm ;
	p->f[9] = -2.0 * t * rp3 * sp * rm * sm ;
	p->f[10] = -t * rp3 * sp * rp * sm ;

	if (x->w < 0.0)
		return(0) ;

	r2p8 =	2.0 * r + 8.0 ;

	p->dfdr[3] = t * sm * sp3 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfdr[0] = t / 3.0 * sm * sp3 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfdr[1] = t / 6.0 * sp * sp3 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfdr[2] = t / 2.0 * sp * sp3 * ( (-rp5 * rp3) + rm * r2p8 ) ;
	p->dfdr[4] = -t / 2.0 * sp * sp3 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[5] = 2.0 * t / 3.0 * sp * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[6] = 2.0 * t / 3.0 * sm * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[7] = -t * sm * sp3 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[8] = t * sp * sm * ( -2.0 * r ) ;
	p->dfdr[9] = -t * 2.0 * sm * sp * ( -2.0 * rp ) ;
	p->dfdr[10] = -t * sp * sm * ( 2.0 * r + 4.0 ) ;

	p->dfds[3] = -t * 2.0 * rp3 * sp * rm * rp5 ;
	p->dfds[0] = -t * 2.0 / 3.0 * rp3 * rp * sp * rp5 ;
	p->dfds[1] = t / 3.0 * rp * rp3 * rp5 * ( s + 2.0 ) ;
	p->dfds[2] = t * rm * rp3 * rp5 * ( s + 2.0 ) ;
	p->dfds[4] = -t * rp * rm * rp5 * ( s + 2.0 ) ;
	p->dfds[5] = 2.0 * t / 3.0 * rp * rm * rp3 ;
	p->dfds[6] = -2.0 * t / 3.0 * rp * rm * rp3 ;
	p->dfds[7] = 2.0 * t * rp * rm * rp5 * sp ;
	p->dfds[8] = -2.0 * t * rp * s * rm ;
	p->dfds[9] = 4.0 * t * rp3 * s * rm ;
	p->dfds[10] = 2.0 * t * rp3 * s * rp ;

	return (0);
}

int sh5_3x2b(x,p)
struct gausspts *x;
struct shapefuncs	*p;

{
	double r, s, t, rm, sm, rp, sp, rp3, rp5, r2p8, sp3 ;

	r = x->y ;
	s = x->x ;
	if ( (r<-5.00001) || (r>1.00001) || (s<-3.00001) || (s>1.00001) ) return(11) ;

	t = 1.0 / 64. ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	sp = 1.0 + s ;
	sm = 1.0 - s ;
	rp3 = 3.0 + r ;
	rp5 = 5.0 + r ;
	sp3 = 3.0 + s ;

	p->dof = 11 ;

	p->f[3] = t * rp3 * sm * rm * rp5 * sp3 ;
	p->f[2] = t / 3.0 * rp3 * rp * sm * rp5 * sp3 ;
	p->f[1] = t / 6.0 * rp * sp * rp3 * rp5 * sp3 ;
	p->f[0] = t / 2.0 * rm * sp * rp3 * rp5 * sp3 ;
	p->f[4] = -t / 2.0 * rp * sp * rm * rp5 * sp3 ;
	p->f[5] = 2.0 * t / 3.0 * rp * sp * rm * rp3 ;
	p->f[6] = 2.0 * t / 3.0 * rp * sm * rm * rp3 ;
	p->f[7] = -t * rp * sm * rm * rp5 * sp3 ;
	p->f[8] = t * rp * sp * rm * sm ;
	p->f[9] = -2.0 * t * rp3 * sp * rm * sm ;
	p->f[10] = -t * rp3 * sp * rp * sm ;

	if (x->w < 0.0)
		return(0) ;

	r2p8 =	2.0 * r + 8.0 ;

	p->dfds[3] = t * sm * sp3 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfds[2] = t / 3.0 * sm * sp3 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfds[1] = t / 6.0 * sp * sp3 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfds[0] = t / 2.0 * sp * sp3 * ( (-rp5 * rp3) + rm * r2p8 ) ;
	p->dfds[4] = -t / 2.0 * sp * sp3 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfds[5] = 2.0 * t / 3.0 * sp * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfds[6] = 2.0 * t / 3.0 * sm * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfds[7] = -t * sm * sp3 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfds[8] = t * sp * sm * ( -2.0 * r ) ;
	p->dfds[9] = -t * 2.0 * sm * sp * ( -2.0 * rp ) ;
	p->dfds[10] = -t * sp * sm * ( 2.0 * r + 4.0 ) ;

	p->dfdr[3] = -t * 2.0 * rp3 * sp * rm * rp5 ;
	p->dfdr[2] = -t * 2.0 / 3.0 * rp3 * rp * sp * rp5 ;
	p->dfdr[1] = t / 3.0 * rp * rp3 * rp5 * ( s + 2.0 ) ;
	p->dfdr[0] = t * rm * rp3 * rp5 * ( s + 2.0 ) ;
	p->dfdr[4] = -t * rp * rm * rp5 * ( s + 2.0 ) ;
	p->dfdr[5] = 2.0 * t / 3.0 * rp * rm * rp3 ;
	p->dfdr[6] = -2.0 * t / 3.0 * rp * rm * rp3 ;
	p->dfdr[7] = 2.0 * t * rp * rm * rp5 * sp ;
	p->dfdr[8] = -2.0 * t * rp * s * rm ;
	p->dfdr[9] = 4.0 * t * rp3 * s * rm ;
	p->dfdr[10] = 2.0 * t * rp3 * s * rp ;

	return (0);
}

