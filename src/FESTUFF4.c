/*			the shape functions for the cubic upwinding functions   239		*/
#include "Fe_PS.h"


extern struct control  N ;

int shcu7_3x3(x,p)
struct gausspts *x;
struct shapefuncs	*p;

{
	double r, s, t, rm, sm, rp, sp, rp3, rp5, r2p8, sp3, sp5, s2p8 ;

	r = x->x ;
	s = x->y ;
	if ( (r<-5.00001) || (r>1.00001) || (s<-5.00001) || (s>1.00001) ) return(14) ;

	t = 1.0 / 256. ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	sp = 1.0 + s ;
	sm = 1.0 - s ;
	rp3 = 3.0 + r ;
	rp5 = 5.0 + r ;
	sp3 = 3.0 + s ;
	sp5 = 5.0 + s ;

	p->dof = 14 ;

	p->f[3] = t * rp3 * sm * rm * rp5 * sp3 * sp5 ;
	p->f[0] = t / 3.0 * rp3 * rp * sm * rp5 * sp3 * sp5 ;
	p->f[1] = t / 9.0 * rp * sp * rp3 * rp5 * sp3 * sp5 ;
	p->f[2] = t / 3.0 * rm * sp * rp3 * rp5 * sp3 * sp5 ;
	p->f[4] = -2.0 * t * rp * sp * rm * rp5 * sp3 ;
	p->f[5] = 2.0 * t / 3.0 * rp * sp * sp3 * rm * rp3 ;
	p->f[6] = 4.0 * t / 3.0 * rp * sm * sp3 * rm * rp3 ;
	p->f[7] = -4.0 * t * rp * sm * rm * rp5 * sp3 ;
	p->f[8] = -2.0 / 3.0 * t * rp * sp * rm * sm * rp3 ;
	p->f[9] = 2.0 * t * rp * sp * rm * sm * rp5 ;
	p->f[10] = -t * rp3 * sp * rm * rp5 * sm * sp5 ;
	p->f[11] = 8.0 * t / 3.0 * sp * rm * sm * sp3 ;
	p->f[12] = 8.0 * t / 3.0 * sp * rp * sm * sp3 ;
	p->f[13] = -t / 3.0 * rp3 * rp5 * sp * rp * sm * sp5 ;

	if (x->w < 0.0)
		return(0) ;

	r2p8 =	2.0 * r + 8.0 ;
	s2p8 =	2.0 * s + 8.0 ;

	p->dfdr[3] = t * sm * sp3 * sp5 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfdr[0] = t / 3.0 * sm * sp3 * sp5 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfdr[1] = t / 9.0 * sp * sp3 * sp5 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfdr[2] = t / 3.0 * sp * sp3 * sp5 * ( (-rp5 * rp3) + rm * r2p8 ) ;
	p->dfdr[4] = -2.0 * t * sp * sp3 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[5] = 2.0 * t / 3.0 * sp * sp3 * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[6] = 4.0 * t / 3.0 * sm * sp3 * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[7] = -4.0 * t * sm * sp3 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[8] = -2.0 / 3.0 * t * sp * sm * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[9] = 2.0 * t * sm * sp * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[10] = -t * sm * sp * sp5 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfdr[11] = -8.0 * t / 3.0 * sm * sp * sp3 ;
	p->dfdr[12] = 8.0 * t / 3.0 * sp * sm * sp3  ;
	p->dfdr[13] = -t / 3.0 * sp * sm * sp5 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;

	p->dfds[3] = t * rp3 * rm * rp5 * ( -( sp5 * sp3 ) + sm * s2p8 ) ;
	p->dfds[0] = t / 3.0 * rp3 * rp * rp5 * ( -( sp5 * sp3 ) + sm * s2p8 ) ;
	p->dfds[1] = t / 9.0 * rp * rp3 * rp5 * ( ( sp5 * sp3 ) + sp * s2p8 ) ;
	p->dfds[2] = t / 3.0 * rm * rp3 * rp5 * ( ( sp5 * sp3 ) + sp * s2p8 ) ;
	p->dfds[4] = -2.0 * t * rp * rm * rp5 * ( 2.0 * s + 4.0 ) ;
	p->dfds[5] = 2.0 * t / 3.0 * rp * rm * rp3 * ( 2.0 * s + 4.0 ) ;
	p->dfds[6] = 4.0 * t / 3.0 * rp * rm * rp3 * (-2.0 * sp );
	p->dfds[7] = -4.0 * t * rp * rm * rp5 * ( -2.0 * sp ) ;
	p->dfds[8] = -2.0 / 3.0 * t * rp * rm * rp3 * ( -2.0 * s) ;
	p->dfds[9] = 2.0 * t * rp * rm * rp5 * ( -2.0 * s ) ;
	p->dfds[10] = -t * rp3 * rm * rp5 * ( -2.0 * s * sp5 + sp * sm ) ;
	p->dfds[11] = 8.0 * t / 3.0 * rm * ( -2.0 * s * sp3 + sp * sm ) ;
	p->dfds[12] = 8.0 * t / 3.0 * rp * ( -2.0 * s * sp3 + sp * sm ) ;
	p->dfds[13] = -t / 3.0 * rp3 * rp * rp5 * ( -2.0 * s * sp5 + sp * sm ) ;

	return (0);
}

int sh7_3x3b(x,p)
struct gausspts *x;
struct shapefuncs	*p;

{
	double r, s, t, rm, sm, rp, sp, rp3, rp5, r2p8, sp3, sp5, s2p8 ;

	r = x->y ;
	s = x->x ;
	if ( (r<-5.00001) || (r>1.00001) || (s<-5.00001) || (s>1.00001) ) return(14) ;

	t = 1.0 / 256. ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	sp = 1.0 + s ;
	sm = 1.0 - s ;
	rp3 = 3.0 + r ;
	rp5 = 5.0 + r ;
	sp3 = 3.0 + s ;
	sp5 = 5.0 + s ;

	p->dof = 14 ;

	p->f[3] = t * rp3 * sm * rm * rp5 * sp3 * sp5 ;
	p->f[2] = t / 3.0 * rp3 * rp * sm * rp5 * sp3 * sp5 ;
	p->f[1] = t / 9.0 * rp * sp * rp3 * rp5 * sp3 * sp5 ;
	p->f[0] = t / 3.0 * rm * sp * rp3 * rp5 * sp3 * sp5 ;
	p->f[4] = -2.0 * t * rp * sp * rm * rp5 * sp3 ;
	p->f[5] = 2.0 * t / 3.0 * rp * sp * sp3 * rm * rp3 ;
	p->f[6] = 4.0 * t / 3.0 * rp * sm * sp3 * rm * rp3 ;
	p->f[7] = -4.0 * t * rp * sm * rm * rp5 * sp3 ;
	p->f[8] = -2.0 / 3.0 * t * rp * sp * rm * sm * rp3 ;
	p->f[9] = 2.0 * t * rp * sp * rm * sm * rp5 ;
	p->f[10] = -t * rp3 * sp * rm * rp5 * sm * sp5 ;
	p->f[11] = 8.0 * t / 3.0 * sp * rm * sm * sp3 ;
	p->f[12] = 8.0 * t / 3.0 * sp * rp * sm * sp3 ;
	p->f[13] = -t / 3.0 * rp3 * rp5 * sp * rp * sm * sp5 ;

	if (x->w < 0.0)
		return(0) ;

	r2p8 =	2.0 * r + 8.0 ;
	s2p8 =	2.0 * s + 8.0 ;

	p->dfds[3] = t * sm * sp3 * sp5 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfds[2] = t / 3.0 * sm * sp3 * sp5 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfds[1] = t / 9.0 * sp * sp3 * sp5 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfds[0] = t / 3.0 * sp * sp3 * sp5 * ( (-rp5 * rp3) + rm * r2p8 ) ;
	p->dfds[4] = -2.0 * t * sp * sp3 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfds[5] = 2.0 * t / 3.0 * sp * sp3 * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfds[6] = 4.0 * t / 3.0 * sm * sp3 * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfds[7] = -4.0 * t * sm * sp3 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfds[8] = -2.0 / 3.0 * t * sp * sm * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfds[9] = 2.0 * t * sm * sp * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfds[10] = -t * sm * sp * sp5 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfds[11] = -8.0 * t / 3.0 * sm * sp * sp3 ;
	p->dfds[12] = 8.0 * t / 3.0 * sp * sm * sp3  ;
	p->dfds[13] = -t / 3.0 * sp * sm * sp5 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;

	p->dfdr[3] = t * rp3 * rm * rp5 * ( -( sp5 * sp3 ) + sm * s2p8 ) ;
	p->dfdr[2] = t / 3.0 * rp3 * rp * rp5 * ( -( sp5 * sp3 ) + sm * s2p8 ) ;
	p->dfdr[1] = t / 9.0 * rp * rp3 * rp5 * ( ( sp5 * sp3 ) + sp * s2p8 ) ;
	p->dfdr[0] = t / 3.0 * rm * rp3 * rp5 * ( ( sp5 * sp3 ) + sp * s2p8 ) ;
	p->dfdr[4] = -2.0 * t * rp * rm * rp5 * ( 2.0 * s + 4.0 ) ;
	p->dfdr[5] = 2.0 * t / 3.0 * rp * rm * rp3 * ( 2.0 * s + 4.0 ) ;
	p->dfdr[6] = 4.0 * t / 3.0 * rp * rm * rp3 * (-2.0 * sp );
	p->dfdr[7] = -4.0 * t * rp * rm * rp5 * ( -2.0 * sp ) ;
	p->dfdr[8] = -2.0 / 3.0 * t * rp * rm * rp3 * ( -2.0 * s) ;
	p->dfdr[9] = 2.0 * t * rp * rm * rp5 * ( -2.0 * s ) ;
	p->dfdr[10] = -t * rp3 * rm * rp5 * ( -2.0 * s * sp5 + sp * sm ) ;
	p->dfdr[11] = 8.0 * t / 3.0 * rm * ( -2.0 * s * sp3 + sp * sm ) ;
	p->dfdr[12] = 8.0 * t / 3.0 * rp * ( -2.0 * s * sp3 + sp * sm ) ;
	p->dfdr[13] = -t / 3.0 * rp3 * rp * rp5 * ( -2.0 * s * sp5 + sp * sm ) ;

	return (0);
}

int shcu8_3x3(x,p)
struct gausspts *x;
struct shapefuncs	*p;

{
	double r, s, t, rm, sm, rp, sp, rp3, rp5, r2p8, sp3, sp5, s2p8 ;

	r = x->x ;
	s = x->y ;
	if ( (r<-5.00001) || (r>1.00001) || (s<-5.00001) || (s>1.00001) ) return(15) ;

	t = 1.0 / 256. ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	sp = 1.0 + s ;
	sm = 1.0 - s ;
	rp3 = 3.0 + r ;
	rp5 = 5.0 + r ;
	sp3 = 3.0 + s ;
	sp5 = 5.0 + s ;

	p->dof = 15 ;



	p->f[3] = t * rp3 * sm * rm * rp5 * sp3 * sp5 ;
	p->f[0] = t / 3.0 * rp3 * rp * sm * rp5 * sp3 * sp5 ;
	p->f[1] = t / 9.0 * rp * sp * rp3 * rp5 * sp3 * sp5 ;
	p->f[2] = t / 3.0 * rm * sp * rp3 * rp5 * sp3 * sp5 ;
	p->f[4] = -t / 3.0 * rp * sp * rm * rp5 * sp3 * sp5 ;
	p->f[5] = 2.0 * t / 3.0 * rp * sp * sp3 * rm * rp3 ;
	p->f[6] = 4.0 * t / 3.0 * rp * sm * sp3 * rm * rp3 ;
	p->f[7] = -t * rp * sm * rm * rp5 * sp3 * sp5 ;
	p->f[8] = -2.0 / 3.0 * t * rp * sp * rm * sm * rp3 ;
	p->f[9] = t * rp * sp * rm * sm * sp5 * rp5 ;
	p->f[10] = -t * rp3 * sp * rm * rp5 * sm * sp5 ;
	p->f[11] = -2.0 / 3.0 * t * sp * rm * rp * sm * sp3 ;
	p->f[12] = 4.0 * t / 3.0 * sp * rp3 * rm * sm * sp3 ;
	p->f[13] = 2.0 * t / 3.0 * sp * rp3 * rp * sm * sp3 ;
	p->f[14] = -t / 3.0 * rp3 * rp5 * sp * rp * sm * sp5 ;

	if (x->w < 0.0)
		return(0) ;
	r2p8 =	2.0 * r + 8.0 ;
	s2p8 =	2.0 * s + 8.0 ;

	p->dfdr[3] = t * sm * sp3 * sp5 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfdr[0] = t / 3.0 * sm * sp3 * sp5 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfdr[1] = t / 9.0 * sp * sp3 * sp5 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfdr[2] = t / 3.0 * sp * sp3 * sp5 * ( (-rp5 * rp3) + rm * r2p8 ) ;
	p->dfdr[4] = -t / 3.0 * sp5 * sp * sp3 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[5] = 2.0 * t / 3.0 * sp3 * sp * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[6] = 4.0 * t / 3.0 * sp3 * sm * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[7] = -t * sm * sp3 * sp5 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[8] = -2.0 / 3.0 * t * sp * sm * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[9] = t * sm * sp * sp5 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[10] = -t * sm * sp * sp5 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfdr[11] = -2.0 / 3.0 * t * sm * sp * sp3 * ( -2.0 * r ) ;
	p->dfdr[12] = -2.0 * t / 3.0 * sm * sp * sp3 * rp ;
	p->dfdr[13] = 2.0 * t / 3.0 * sp * sm * sp3 * ( 2.0 * r + 4.0 ) ;
	p->dfdr[14] = -t / 3.0 * sp * sm * sp5 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;

	p->dfds[3] = t * rp3 * rm * rp5 * ( -( sp5 * sp3 ) + sm * s2p8 ) ;
	p->dfds[0] = t / 3.0 * rp3 * rp * rp5 * ( -( sp5 * sp3 ) + sm * s2p8 ) ;
	p->dfds[1] = t / 9.0 * rp * rp3 * rp5 * ( ( sp5 * sp3 ) + sp * s2p8 ) ;
	p->dfds[2] = t / 3.0 * rm * rp3 * rp5 * ( ( sp5 * sp3 ) + sp * s2p8 ) ;
	p->dfds[4] = -t / 3.0 * rp * rm * rp5 * ( ( sp5 * sp3 ) + sp * s2p8 ) ;
	p->dfds[5] = 2.0 * t / 3.0 * rp * rm * rp3 * ( 2.0 * s + 4.0 ) ;
	p->dfds[6] = 4.0 * t / 3.0 * rp * rm * rp3 * ( -2.0 * sp ) ;
	p->dfds[7] = -t * rp * rm * rp5 * ( -( sp5 * sp3 ) + sm * s2p8 ) ;
	p->dfds[8] = -2.0 / 3.0 * t * rp * rm * rp3 * ( -2.0 * s) ;
	p->dfds[9] = t * rp * rm * rp5 * ( -2.0 * s * sp5 + sp * sm ) ;
	p->dfds[10] = -t * rp3 * rm * rp5 * ( -2.0 * s * sp5 + sp * sm ) ;
	p->dfds[11] = -2.0 / 3.0 * t * rm * rp * ( -2.0 * s * sp3 + sp * sm ) ;
	p->dfds[12] = 4.0 * t / 3.0 * rm * rp3 * ( -2.0 * s * sp3 + sp * sm ) ;
	p->dfds[13] = 2.0 * t / 3.0 * rp * rp3 * ( -2.0 * s * sp3 + sp * sm ) ;
	p->dfds[14] = -t / 3.0 * rp3 * rp * rp5 * ( -2.0 * s * sp5 + sp * sm ) ;

	return (0);
}

