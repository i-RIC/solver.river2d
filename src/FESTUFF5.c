/*			the shape functions for the cubic upwinding functions   239		*/
#include "Fe_PS.h"


extern struct control  N ;

int shcu9_3x3(x,p)
struct gausspts *x;
struct shapefuncs	*p;

{
	double r, s, t, rm, sm, rp, sp, rp3, rp5, r2p8, sp3, sp5, s2p8 ;

	r = x->x ;
	s = x->y ;
	if ( (r<-5.00001) || (r>1.00001) || (s<-5.00001) || (s>1.00001) ) return(16) ;

	t = 1.0 / 256. ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	sp = 1.0 + s ;
	sm = 1.0 - s ;
	rp3 = 3.0 + r ;
	rp5 = 5.0 + r ;
	sp3 = 3.0 + s ;
	sp5 = 5.0 + s ;

	p->dof = 16 ;

	p->f[3] = t * rp3 * sm * rm * rp5 * sp3 * sp5 ;
	p->f[0] = t / 3.0 * rp3 * rp * sm * rp5 * sp3 * sp5 ;
	p->f[1] = t / 9.0 * rp * sp * rp3 * rp5 * sp3 * sp5 ;
	p->f[2] = t / 3.0 * rm * sp * rp3 * rp5 * sp3 * sp5 ;
	p->f[4] = -t / 3.0 * rp * sp * rm * rp5 * sp3 * sp5 ;
	p->f[5] = t / 9.0 * rp * sp * rm * rp3 * sp3 * sp5 ;
	p->f[6] = t / 3.0 * rp * sm * rm * rp3 * sp3 * sp5 ;
	p->f[7] = -t * rp * sm * rm * rp5 * sp3 * sp5 ;
	p->f[8] = -t / 3.0 * rp * sm * rm * rp3 * sp * sp5 ;
	p->f[9] = t * rp * sp * rm * rp5 * sm * sp5 ;
	p->f[10] = -t * rp3 * sp * rm * rp5 * sm * sp5 ;
	p->f[11] = t / 9.0 * rp3 * sp * rp * rm * sm * sp3 ;
	p->f[12] = -t / 3.0 * rp * sp * rm * rp5 * sm * sp3 ;
	p->f[13] = t / 3.0 * rp3 * sp * rm * rp5 * sm * sp3 ;
	p->f[14] = t / 9.0 * rp3 * sp * rp * rp5 * sm * sp3 ;
	p->f[15] = -t / 3.0 * rp3 * sp * rp * rp5 * sm * sp5 ;

	if (x->w < 0.0)
		return(0) ;

	r2p8 =	2.0 * r + 8.0 ;
	s2p8 =	2.0 * s + 8.0 ;


	p->dfdr[3] = t * sm * sp3 * sp5 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfdr[0] = t / 3.0 * sm * sp3 * sp5 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfdr[1] = t / 9.0 * sp * sp3 * sp5 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfdr[2] = t / 3.0 * sp * sp3 * sp5 * ( (-rp5 * rp3) + rm * r2p8 ) ;
	p->dfdr[4] = -t / 3.0 * sp * sp3 * sp5 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[5] = t / 9.0 * sp * sp3 * sp5 * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[6] = t / 3.0 * sm * sp3 * sp5 * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[7] = -t * sm * sp3 * sp5 * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[8] = -t / 3.0 * sm * sp * sp5 * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[9] = t * sm * sp5 * sp * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[10] = -t * sm * sp * sp5 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfdr[11] = t / 9.0 * sm * sp * sp3 * ( -2.0 * r * rp3 + rp * rm ) ;
	p->dfdr[12] = -t / 3.0 * sm * sp3 * sp * ( -2.0 * r * rp5 + rp * rm ) ;
	p->dfdr[13] = t / 3.0 * sm * sp * sp3 * ( (-rp5 * rp3 ) + rm * r2p8 ) ;
	p->dfdr[14] = t / 9.0 * sp * sm * sp3 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;
	p->dfdr[15] = -t / 3.0 * sp * sm * sp5 * ( ( rp5 * rp3 ) + rp * r2p8 ) ;

	p->dfds[3] = t * rp3 * rm * rp5 * ( -( sp5 * sp3 ) + sm * s2p8 ) ;
	p->dfds[0] = t / 3.0 * rp3 * rp * rp5 * ( -( sp5 * sp3 ) + sm * s2p8 ) ;
	p->dfds[1] = t / 9.0 * rp * rp3 * rp5 * ( ( sp5 * sp3 ) + sp * s2p8 ) ;
	p->dfds[2] = t / 3.0 * rm * rp3 * rp5 * ( ( sp5 * sp3 ) + sp * s2p8 ) ;
	p->dfds[4] = -t / 3.0 * rp * rm * rp5 * ( ( sp5 * sp3 ) + sp * s2p8 ) ;
	p->dfds[5] = t / 9.0 * rp * rm * rp3  * ( ( sp5 * sp3 ) + sp * s2p8 ) ;
	p->dfds[6] = t / 3.0 * rp * rm * rp3  * ( -( sp5 * sp3 ) + sm * s2p8 ) ;
	p->dfds[7] = -t * rp * rm * rp5  * ( -( sp5 * sp3 ) + sm * s2p8 ) ;
	p->dfds[8] = -t / 3.0 * rp * rm * rp3 * ( -2.0 * s * sp5 + sp * sm ) ;
	p->dfds[9] = t * rp * rm * rp5 * ( -2.0 * s * sp5 + sp * sm ) ;
	p->dfds[10] = -t * rp3 * rm * rp5 * ( -2.0 * s * sp5 + sp * sm ) ;
	p->dfds[11] = t / 9.0 * rp3 * rp * rm * ( -2.0 * s * sp3 + sp * sm ) ;
	p->dfds[12] = -t / 3.0 * rp * rm * rp5 * ( -2.0 * s * sp3 + sp * sm ) ;
	p->dfds[13] = t / 3.0 * rp3 * rm * rp5 * ( -2.0 * s * sp3 + sp * sm ) ;
	p->dfds[14] = t / 9.0 * rp3 * rp * rp5 * ( -2.0 * s * sp3 + sp * sm ) ;
	p->dfds[15] = -t / 3.0 * rp3 * rp * rp5 * ( -2.0 * s * sp5 + sp * sm ) ;

	return (0);

}
