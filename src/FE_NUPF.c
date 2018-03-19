//#pragma options cpluscmt	
#include "Fe_PS.h"


extern struct control  	N ;
extern double UW ;

int		SU1DPick(elmntp,theElp,ntfp,dirp)
struct element	*elmntp, *theElp ;
int		*ntfp, *dirp ;

{
	int	j;
	double	q01 ;
	struct node			*np[MNSF] ;
	
	theElp->vtype = 111 ;
	*ntfp = 2;
	for(j=0;j<elmntp->nnds;j++)
		np[j] = elmntp->nps[j] ;
		
	for(j=0;j<N.params;j++)
		theElp->p[j] = elmntp->p[j] ;
		
	q01 = (np[1]->u[1]-np[0]->u[1])
			/(np[1]->u[0]-np[0]->u[0]) ;

/*	this part used to send gi only */
	q01 = 1.0 ;	 

/*	this part used to send vi+gi only: so we need to know which way to upwind */

	if (q01 >= 0.0) {
			theElp->gtype = 119 ;
			theElp->nnds = 2 ;
			theElp->nps[0] = elmntp->nps[0] ;
			theElp->nps[1] = elmntp->nps[1] ;
			return(0) ;
	}
	else {
			theElp->gtype = 118 ;
			theElp->nnds = 2 ;
			theElp->nps[0] = elmntp->nps[0] ;
			theElp->nps[1] = elmntp->nps[1] ;
			return(0) ;
	}


}

int		QU1DPick(elmntp,theElp,ntfp,dirp)
struct element	*elmntp, *theElp ;
int		*ntfp, *dirp ;

{
	int	j;
	double	q01 ;
	struct node			*np[MNSF] ;
	
	theElp->gtype = 111 ;
	*ntfp = 3;
	for(j=0;j<elmntp->nnds;j++)
		np[j] = elmntp->nps[j] ;
		
	for(j=0;j<N.params;j++)
		theElp->p[j] = elmntp->p[j] ;
		
	q01 = np[0]->p[0] * (np[1]->x[0]-np[0]->x[0])
			+ np[0]->p[1] *(np[1]->x[1]-np[0]->x[1]) ;
	q01 = 1.0;

	if (q01 >= 0.0) {
		if(elmntp->nps[2] == NULL ){
			theElp->vtype = theElp->gtype = 111 ;
			*ntfp = 2;
			theElp->nnds = 2 ;
			theElp->nps[0] = elmntp->nps[0] ;
			theElp->nps[1] = elmntp->nps[1] ;
			return(0) ;
		}
		else {
			theElp->vtype = 129 ;
			theElp->gtype = 129 ;
			theElp->nnds = 3 ;
			theElp->nps[0] = elmntp->nps[0] ;
			theElp->nps[1] = elmntp->nps[1] ;
			theElp->nps[2] = elmntp->nps[2] ;
			return(1) ;
		}
	}
	else {
		if(elmntp->nps[3] == NULL ){
			theElp->vtype = theElp->gtype = 111 ;
			theElp->nnds = 2 ;
			theElp->nps[0] = elmntp->nps[0] ;
			theElp->nps[1] = elmntp->nps[1] ;
			return(0) ;
		}
		else  {
			theElp->vtype = 128 ;
			theElp->gtype = 111 ;
			theElp->nnds = 3 ;
			theElp->nps[0] = elmntp->nps[0] ;
			theElp->nps[1] = elmntp->nps[1] ;
			theElp->nps[2] = elmntp->nps[3] ;
			return(1) ;
		}
	}
}
		

int		CU1DPick(elmntp,theElp,ntfp,dirp)
struct element	*elmntp, *theElp ;
int		*ntfp, *dirp ;

{
	int	j;
	double	q01 ;
	struct node			*np[MNSF] ;
	
	theElp->gtype = 111 ;
	*ntfp = 2;
	for(j=0;j<elmntp->nnds;j++)
		np[j] = elmntp->nps[j] ;
		
	for(j=0;j<N.params;j++) 
		theElp->p[j] = elmntp->p[j] ;
		
	q01 = (np[0]->p[0]+np[1]->p[0])/2 ;

	if (q01 >= 0.0) {
		if((elmntp->nps[2] == NULL) || (elmntp->nps[3] == NULL) ){
			theElp->vtype = theElp->gtype = 111 ;
			theElp->nnds = 2 ;
			theElp->nps[0] = elmntp->nps[0] ;
			theElp->nps[1] = elmntp->nps[1] ;
			return(0) ;
		}
		else {
			theElp->vtype = 139 ;
			theElp->gtype = 111 ;
			*ntfp = 2;
			theElp->nnds = 4 ;
			theElp->nps[0] = elmntp->nps[0] ;
			theElp->nps[1] = elmntp->nps[1] ;
			theElp->nps[2] = elmntp->nps[2] ;
			theElp->nps[3] = elmntp->nps[3] ;
			return(1) ;
		}
	}
	else {
		if((elmntp->nps[4] == NULL) || (elmntp->nps[5] == NULL) ){
			theElp->vtype = theElp->gtype = 111 ;
			theElp->nnds = 2 ;
			theElp->nps[0] = elmntp->nps[0] ;
			theElp->nps[1] = elmntp->nps[1] ;
			return(0) ;
		}
		else  {
			theElp->vtype = 138 ;
			theElp->gtype = 111 ;
			*ntfp = 2;
			theElp->nnds = 4 ;
			theElp->nps[0] = elmntp->nps[0] ;
			theElp->nps[1] = elmntp->nps[1] ;
			theElp->nps[2] = elmntp->nps[4] ;
			theElp->nps[3] = elmntp->nps[5] ;
			return(1) ;
		}
	}
/*	return(-1) ;  */
}


int		QUelPick(elmntp,theElp,ntfp,dirp)
struct element	*elmntp, *theElp ;
int		*ntfp, *dirp ;

{
	int	j;
	double	q01, q12, q23, q30 ;
	struct node			*np[MNSF] ;
	
	theElp->gtype = 219 ;
	*ntfp = 4;
	for(j=0;j<elmntp->nnds;j++)
		np[j] = elmntp->nps[j] ;
		
	for(j=0;j<N.params;j++) 
		theElp->p[j] = elmntp->p[j] ;

	q01 = (np[0]->p[0]+np[1]->p[0])/2*(np[0]->x[1]-np[1]->x[1])
			- (np[0]->p[1]+np[1]->p[1])/2*(np[0]->x[0]-np[1]->x[0]) ;
	q12 = (np[1]->p[0]+np[2]->p[0])/2*(np[1]->x[1]-np[2]->x[1])
			- (np[1]->p[1]+np[2]->p[1])/2*(np[1]->x[0]-np[2]->x[0]) ;
	q23 = (np[2]->p[0]+np[3]->p[0])/2*(np[2]->x[1]-np[3]->x[1])
			- (np[2]->p[1]+np[3]->p[1])/2*(np[2]->x[0]-np[3]->x[0]) ;
	q30 = (np[3]->p[0]+np[0]->p[0])/2*(np[3]->x[1]-np[0]->x[1])
			- (np[3]->p[1]+np[0]->p[1])/2*(np[3]->x[0]-np[0]->x[0]) ;
/*
	q01 = -np[1]->u[0] + np[0]->u[0] ;
	q12 = -np[2]->u[0] + np[1]->u[0] ;
	q23 = -np[3]->u[0] + np[2]->u[0] ;
	q30 = -np[0]->u[0] + np[3]->u[0] ;
*/	

	if((q23>=q01) && (q30>=q12)) {
		if((elmntp->nps[14] == NULL) || (elmntp->nps[15] == NULL) ){
			if((elmntp->nps[11] == NULL) || (elmntp->nps[12] == NULL) ){
				theElp->vtype = theElp->gtype = 211 ;
				theElp->nnds = 4 ;
				theElp->nps[0] = elmntp->nps[0] ;
				theElp->nps[1] = elmntp->nps[1] ;
				theElp->nps[2] = elmntp->nps[2] ;
				theElp->nps[3] = elmntp->nps[3] ;
				return(0) ;
			}
			else {
				theElp->vtype = 228 ;
				theElp->nnds = 6 ;
				theElp->nps[0] = elmntp->nps[0] ;
				theElp->nps[1] = elmntp->nps[1] ;
				theElp->nps[2] = elmntp->nps[2] ;
				theElp->nps[3] = elmntp->nps[3] ;
				theElp->nps[4] = elmntp->nps[11] ;
				theElp->nps[5] = elmntp->nps[12] ;
				return(1) ;
			}
		}
		else {
			if((elmntp->nps[11] == NULL) || (elmntp->nps[12] == NULL) ){
				theElp->vtype = 228 ;
				theElp->nnds = 6 ;
				theElp->nps[0] = elmntp->nps[1] ;
				theElp->nps[1] = elmntp->nps[2] ;
				theElp->nps[2] = elmntp->nps[3] ;
				theElp->nps[3] = elmntp->nps[0] ;
				theElp->nps[4] = elmntp->nps[14] ;
				theElp->nps[5] = elmntp->nps[15] ;
				return(1) ;
			}
			else if (elmntp->nps[13] == NULL) {
				theElp->vtype = 224 ;
				theElp->nnds = 8 ;
				theElp->nps[0] = elmntp->nps[0] ;
				theElp->nps[1] = elmntp->nps[1] ;
				theElp->nps[2] = elmntp->nps[2] ;
				theElp->nps[3] = elmntp->nps[3] ;
				theElp->nps[4] = elmntp->nps[11] ;
				theElp->nps[5] = elmntp->nps[12] ;
				theElp->nps[7] = elmntp->nps[14] ;
				theElp->nps[8] = elmntp->nps[15] ;
				return(2) ;
			}
			else {
				theElp->vtype = 229 ;
				theElp->nnds = 9 ;
				theElp->nps[0] = elmntp->nps[0] ;
				theElp->nps[1] = elmntp->nps[1] ;
				theElp->nps[2] = elmntp->nps[2] ;
				theElp->nps[3] = elmntp->nps[3] ;
				theElp->nps[4] = elmntp->nps[11] ;
				theElp->nps[5] = elmntp->nps[12] ;
				theElp->nps[6] = elmntp->nps[13] ;
				theElp->nps[7] = elmntp->nps[14] ;
				theElp->nps[8] = elmntp->nps[15] ;
				return(3) ;
			}
		}
	}
	if((q23>=q01) && (q30<q12)) {
		if((elmntp->nps[11] == NULL) || (elmntp->nps[12] == NULL) ){
			if((elmntp->nps[8] == NULL) || (elmntp->nps[9] == NULL) ){
				theElp->vtype = 211 ;
				theElp->nnds = 4 ;
				theElp->nps[0] = elmntp->nps[0] ;
				theElp->nps[1] = elmntp->nps[1] ;
				theElp->nps[2] = elmntp->nps[2] ;
				theElp->nps[3] = elmntp->nps[3] ;
				return(0) ;
			}
			else {
				theElp->vtype = 228 ;
				theElp->nnds = 6 ;
				theElp->nps[0] = elmntp->nps[3] ;
				theElp->nps[1] = elmntp->nps[0] ;
				theElp->nps[2] = elmntp->nps[1] ;
				theElp->nps[3] = elmntp->nps[2] ;
				theElp->nps[4] = elmntp->nps[8] ;
				theElp->nps[5] = elmntp->nps[9] ;
				return(1) ;
			}
		}
		else {
			if((elmntp->nps[8] == NULL) || (elmntp->nps[9] == NULL) ){
				theElp->vtype = 228 ;
				theElp->nnds = 6 ;
				theElp->nps[0] = elmntp->nps[0] ;
				theElp->nps[1] = elmntp->nps[1] ;
				theElp->nps[2] = elmntp->nps[2] ;
				theElp->nps[3] = elmntp->nps[3] ;
				theElp->nps[4] = elmntp->nps[11] ;
				theElp->nps[5] = elmntp->nps[12] ;
				return(1) ;
			}
			else if (elmntp->nps[10] == NULL) {
				theElp->vtype = 224 ;
				theElp->nnds = 8 ;
				theElp->nps[0] = elmntp->nps[3] ;
				theElp->nps[1] = elmntp->nps[0] ;
				theElp->nps[2] = elmntp->nps[1] ;
				theElp->nps[3] = elmntp->nps[2] ;
				theElp->nps[4] = elmntp->nps[8] ;
				theElp->nps[5] = elmntp->nps[9] ;
				theElp->nps[7] = elmntp->nps[11] ;
				theElp->nps[8] = elmntp->nps[12] ;
				return(2) ;
			}
			else {
				theElp->vtype = 229 ;
				theElp->nnds = 9 ;
				theElp->nps[0] = elmntp->nps[3] ;
				theElp->nps[1] = elmntp->nps[0] ;
				theElp->nps[2] = elmntp->nps[1] ;
				theElp->nps[3] = elmntp->nps[2] ;
				theElp->nps[4] = elmntp->nps[8] ;
				theElp->nps[5] = elmntp->nps[9] ;
				theElp->nps[6] = elmntp->nps[10] ;
				theElp->nps[7] = elmntp->nps[11] ;
				theElp->nps[8] = elmntp->nps[12] ;
				return(3) ;
			}
		}
	}
	if((q23<q01) && (q30<q12)) {
		if((elmntp->nps[8] == NULL) || (elmntp->nps[9] == NULL) ){
			if((elmntp->nps[5] == NULL) || (elmntp->nps[6] == NULL) ){
				theElp->vtype = 211 ;
				theElp->nnds = 4 ;
				theElp->nps[0] = elmntp->nps[0] ;
				theElp->nps[1] = elmntp->nps[1] ;
				theElp->nps[2] = elmntp->nps[2] ;
				theElp->nps[3] = elmntp->nps[3] ;
				return(0) ;
			}
			else {
				theElp->vtype = 228 ;
				theElp->nnds = 6 ;
				theElp->nps[0] = elmntp->nps[2] ;
				theElp->nps[1] = elmntp->nps[3] ;
				theElp->nps[2] = elmntp->nps[0] ;
				theElp->nps[3] = elmntp->nps[1] ;
				theElp->nps[4] = elmntp->nps[5] ;
				theElp->nps[5] = elmntp->nps[6] ;
				return(1) ;
			}
		}
		else {
			if((elmntp->nps[5] == NULL) || (elmntp->nps[6] == NULL) ){
				theElp->vtype = 228 ;
				theElp->nnds = 6 ;
				theElp->nps[0] = elmntp->nps[3] ;
				theElp->nps[1] = elmntp->nps[0] ;
				theElp->nps[2] = elmntp->nps[1] ;
				theElp->nps[3] = elmntp->nps[2] ;
				theElp->nps[4] = elmntp->nps[8] ;
				theElp->nps[5] = elmntp->nps[9] ;
				return(1) ;
			}
			else if (elmntp->nps[7] == NULL) {
				theElp->vtype = 224 ;
				theElp->nnds = 8 ;
				theElp->nps[0] = elmntp->nps[2] ;
				theElp->nps[1] = elmntp->nps[3] ;
				theElp->nps[2] = elmntp->nps[0] ;
				theElp->nps[3] = elmntp->nps[1] ;
				theElp->nps[4] = elmntp->nps[5] ;
				theElp->nps[5] = elmntp->nps[6] ;
				theElp->nps[6] = elmntp->nps[8] ;
				theElp->nps[7] = elmntp->nps[9] ;
				return(2) ;
			}
			else {
				theElp->vtype = 229 ;
				theElp->nnds = 9 ;
				theElp->nps[0] = elmntp->nps[2] ;
				theElp->nps[1] = elmntp->nps[3] ;
				theElp->nps[2] = elmntp->nps[0] ;
				theElp->nps[3] = elmntp->nps[1] ;
				theElp->nps[4] = elmntp->nps[5] ;
				theElp->nps[5] = elmntp->nps[6] ;
				theElp->nps[6] = elmntp->nps[7] ;
				theElp->nps[7] = elmntp->nps[8] ;
				theElp->nps[8] = elmntp->nps[9] ;
				return(3) ;
			}
		}
	}
	if((q23<q01) && (q30>=q12)) {
		if((elmntp->nps[5] == NULL) || (elmntp->nps[6] == NULL) ){
			if((elmntp->nps[14] == NULL) || (elmntp->nps[15] == NULL) ){
				theElp->vtype = 211 ;
				theElp->nnds = 4 ;
				theElp->nps[0] = elmntp->nps[0] ;
				theElp->nps[1] = elmntp->nps[1] ;
				theElp->nps[2] = elmntp->nps[2] ;
				theElp->nps[3] = elmntp->nps[3] ;
				return(0) ;
			}
			else {
				theElp->vtype = 228 ;
				theElp->nnds = 6 ;
				theElp->nps[0] = elmntp->nps[1] ;
				theElp->nps[1] = elmntp->nps[2] ;
				theElp->nps[2] = elmntp->nps[3] ;
				theElp->nps[3] = elmntp->nps[0] ;
				theElp->nps[4] = elmntp->nps[14] ;
				theElp->nps[5] = elmntp->nps[15] ;
				return(1) ;
			}
		}
		else {
			if((elmntp->nps[14] == NULL) || (elmntp->nps[15] == NULL) ){
				theElp->vtype = 228 ;
				theElp->nnds = 6 ;
				theElp->nps[0] = elmntp->nps[2] ;
				theElp->nps[1] = elmntp->nps[3] ;
				theElp->nps[2] = elmntp->nps[0] ;
				theElp->nps[3] = elmntp->nps[1] ;
				theElp->nps[4] = elmntp->nps[5] ;
				theElp->nps[5] = elmntp->nps[6] ;
				return(1) ;
			}
			else if (elmntp->nps[4] == NULL) {
				theElp->vtype = 224 ;
				theElp->nnds = 8 ;
				theElp->nps[0] = elmntp->nps[1] ;
				theElp->nps[1] = elmntp->nps[2] ;
				theElp->nps[2] = elmntp->nps[3] ;
				theElp->nps[3] = elmntp->nps[0] ;
				theElp->nps[4] = elmntp->nps[14] ;
				theElp->nps[5] = elmntp->nps[15] ;
				theElp->nps[6] = elmntp->nps[5] ;
				theElp->nps[7] = elmntp->nps[6] ;
				return(2) ;
			}
			else {
				theElp->vtype = 229 ;
				theElp->nnds = 9 ;
				theElp->nps[0] = elmntp->nps[1] ;
				theElp->nps[1] = elmntp->nps[2] ;
				theElp->nps[2] = elmntp->nps[3] ;
				theElp->nps[3] = elmntp->nps[0] ;
				theElp->nps[4] = elmntp->nps[14] ;
				theElp->nps[5] = elmntp->nps[15] ;
				theElp->nps[6] = elmntp->nps[4] ;
				theElp->nps[7] = elmntp->nps[5] ;
				theElp->nps[8] = elmntp->nps[6] ;
				return(3) ;
			}
		}
	}
	return(-1) ;
}

//int		CU2DPick(elmntp,theElp,ntfp,dirp)
//struct element	*elmntp, *theElp ;
//int		*ntfp, *dirp ;
//
//{
//	int	j;
//	double	q01, q12, q23, q30 ;
//	struct node	*np[MNSF] ;
//	int node_hands[MNSF] ;
//
//	theElp->gtype = 211 ;
//	*ntfp = 4;
//	for(j=0;j<elmntp->nnds;j++)
//		np[j] = elmntp->nps[j] ;
//
//	for(j=0;j<N.params;j++)
//		theElp->p[j] = elmntp->p[j] ;
///*
//	q01 = (np[0]->p[0]+np[1]->p[0])/2*(np[0]->x[1]-np[1]->x[1])
//			- (np[0]->p[1]+np[1]->p[1])/2*(np[0]->x[0]-np[1]->x[0]) ;
//	q12 = (np[1]->p[0]+np[2]->p[0])/2*(np[1]->x[1]-np[2]->x[1])
//			- (np[1]->p[1]+np[2]->p[1])/2*(np[1]->x[0]-np[2]->x[0]) ;
//	q23 = (np[2]->p[0]+np[3]->p[0])/2*(np[2]->x[1]-np[3]->x[1])
//			- (np[2]->p[1]+np[3]->p[1])/2*(np[2]->x[0]-np[3]->x[0]) ;
//	q30 = (np[3]->p[0]+np[0]->p[0])/2*(np[3]->x[1]-np[0]->x[1])
//			- (np[3]->p[1]+np[0]->p[1])/2*(np[3]->x[0]-np[0]->x[0]) ;
//*/
//	q01 = -np[1]->u[0] + np[0]->u[0] ;
//	q12 = -np[2]->u[0] + np[1]->u[0] ;
//	q23 = -np[3]->u[0] + np[2]->u[0] ;
//	q30 = -np[0]->u[0] + np[3]->u[0] ;
//	
//
//	if((q23>=q01) && (q30>=q12)) {
//		node_hands[0] = 0 ;
//		node_hands[1] = 1 ;
//		node_hands[2] = 2 ;
//		node_hands[3] = 3 ;
//		node_hands[4] = 11 ;
//		node_hands[5] = 28 ;
//		node_hands[6] = 29 ;
//		node_hands[7] = 12 ;
//		node_hands[8] = 30 ;
//		node_hands[9] = 13 ;
//		node_hands[10] = 14 ;
//		node_hands[11] = 15 ;
//		node_hands[12] = 31 ;
//		node_hands[13] = 32 ;
//		node_hands[14] = 33 ;
//		node_hands[15] = 34 ;
//	}
//	if((q23>=q01) && (q30<q12)) {
//		node_hands[0] = 3 ;
//		node_hands[1] = 0 ;
//		node_hands[2] = 1 ;
//		node_hands[3] = 2 ;
//		node_hands[4] = 8 ;
//		node_hands[5] = 23 ;
//		node_hands[6] = 24 ;
//		node_hands[7] = 9 ;
//		node_hands[8] = 25 ;
//		node_hands[9] = 10 ;
//		node_hands[10] = 11 ;
//		node_hands[11] = 12 ;
//		node_hands[12] = 26 ;
//		node_hands[13] = 27 ;
//		node_hands[14] = 28 ;
//		node_hands[15] = 29 ;
//	}
//	if((q23<q01) && (q30<q12)) {
//		node_hands[0] = 2 ;
//		node_hands[1] = 3 ;
//		node_hands[2] = 0 ;
//		node_hands[3] = 1 ;
//		node_hands[4] = 5 ;
//		node_hands[5] = 18 ;
//		node_hands[6] = 19 ;
//		node_hands[7] = 6 ;
//		node_hands[8] = 20 ;
//		node_hands[9] = 7 ;
//		node_hands[10] = 8 ;
//		node_hands[11] = 9 ;
//		node_hands[12] = 21 ;
//		node_hands[13] = 22 ;
//		node_hands[14] = 23 ;
//		node_hands[15] = 24 ;
//	}
//	if((q23<q01) && (q30>=q12)) {
//		node_hands[0] = 1 ;
//		node_hands[1] = 2 ;
//		node_hands[2] = 3 ;
//		node_hands[3] = 0 ;
//		node_hands[4] = 14 ;
//		node_hands[5] = 33 ;
//		node_hands[6] = 34 ;
//		node_hands[7] = 15 ;
//		node_hands[8] = 35 ;
//		node_hands[9] = 4 ;
//		node_hands[10] = 5 ;
//		node_hands[11] = 6 ;
//		node_hands[12] = 16 ;
//		node_hands[13] = 17 ;
//		node_hands[14] = 18 ;
//		node_hands[15] = 19 ;
//	}
//	reall_nodes_CU2D( elmntp, theElp, node_hands) ;
//
//	return(-1) ;
//}


int reall_nodes_CU2D(elmntp, theElp, n_h)
struct element	*elmntp, *theElp ;
int	n_h[MNSF] ;
{

		if((elmntp->nps[n_h[12]] == NULL)){
			theElp->vtype = 238 ;
			theElp->nnds = 15 ;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[4]] ;
			theElp->nps[5] = elmntp->nps[n_h[5]] ;
			theElp->nps[6] = elmntp->nps[n_h[6]] ;
			theElp->nps[7] = elmntp->nps[n_h[7]] ;
			theElp->nps[8] = elmntp->nps[n_h[8]] ;
			theElp->nps[9] = elmntp->nps[n_h[9]] ;
			theElp->nps[10] = elmntp->nps[n_h[10]] ;
			theElp->nps[11] = elmntp->nps[n_h[13]] ;
			theElp->nps[12] = elmntp->nps[n_h[14]] ;
			theElp->nps[13] = elmntp->nps[n_h[15]] ;
			theElp->nps[14] = elmntp->nps[n_h[11]] ;
		
		}
		else {
			theElp->vtype = 239 ;
			theElp->nnds = 16 ;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[4]] ;
			theElp->nps[5] = elmntp->nps[n_h[5]] ;
			theElp->nps[6] = elmntp->nps[n_h[6]] ;
			theElp->nps[7] = elmntp->nps[n_h[7]] ;
			theElp->nps[8] = elmntp->nps[n_h[8]] ;
			theElp->nps[9] = elmntp->nps[n_h[9]] ;
			theElp->nps[10] = elmntp->nps[n_h[10]] ;
			theElp->nps[11] = elmntp->nps[n_h[12]] ;
			theElp->nps[12] = elmntp->nps[n_h[13]] ;
			theElp->nps[13] = elmntp->nps[n_h[14]] ;
			theElp->nps[14] = elmntp->nps[n_h[15]] ;
			theElp->nps[15] = elmntp->nps[n_h[11]] ;
		
		}
		if ((elmntp->nps[n_h[13]] == NULL)){
			theElp->vtype = 2371 ;
			theElp->nnds = 14 ;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
  			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[4]] ;
			theElp->nps[5] = elmntp->nps[n_h[5]] ;
			theElp->nps[6] = elmntp->nps[n_h[6]] ;
			theElp->nps[7] = elmntp->nps[n_h[7]] ;
			theElp->nps[8] = elmntp->nps[n_h[8]] ;
			theElp->nps[9] = elmntp->nps[n_h[9]] ;
			theElp->nps[10] = elmntp->nps[n_h[10]] ;
			theElp->nps[11] = elmntp->nps[n_h[14]] ;
			theElp->nps[12] = elmntp->nps[n_h[15]] ;
			theElp->nps[13] = elmntp->nps[n_h[11]] ;
		
		}
		if ((elmntp->nps[n_h[8]] == NULL)){
			theElp->vtype = 2372 ;
			theElp->nnds = 14 ;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
  			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[11]] ;
			theElp->nps[5] = elmntp->nps[n_h[15]] ;
			theElp->nps[6] = elmntp->nps[n_h[14]] ;
			theElp->nps[7] = elmntp->nps[n_h[10]] ;
			theElp->nps[8] = elmntp->nps[n_h[13]] ;
			theElp->nps[9] = elmntp->nps[n_h[9]] ;
			theElp->nps[10] = elmntp->nps[n_h[7]] ;
			theElp->nps[11] = elmntp->nps[n_h[6]] ;
			theElp->nps[12] = elmntp->nps[n_h[5]] ;
			theElp->nps[13] = elmntp->nps[n_h[4]] ;

		}
		if ((elmntp->nps[n_h[14]] == NULL) || (elmntp->nps[n_h[15]] == NULL)){
			theElp->vtype = 2361;
			theElp->nnds = 12;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
  			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[4]] ;
			theElp->nps[5] = elmntp->nps[n_h[5]] ;
			theElp->nps[6] = elmntp->nps[n_h[6]] ;
			theElp->nps[7] = elmntp->nps[n_h[7]] ;
			theElp->nps[8] = elmntp->nps[n_h[8]] ;
			theElp->nps[9] = elmntp->nps[n_h[9]] ;
			theElp->nps[10] = elmntp->nps[n_h[10]] ;
			theElp->nps[11] = elmntp->nps[n_h[11]] ;
	
		}
		if ((elmntp->nps[n_h[5]] == NULL) || (elmntp->nps[n_h[6]] == NULL)){
			theElp->vtype = 2363;
			theElp->nnds = 12;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
  			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[11]] ;
			theElp->nps[5] = elmntp->nps[n_h[15]] ;
			theElp->nps[6] = elmntp->nps[n_h[14]] ;
			theElp->nps[7] = elmntp->nps[n_h[10]] ;
			theElp->nps[8] = elmntp->nps[n_h[13]] ;
			theElp->nps[9] = elmntp->nps[n_h[9]] ;
			theElp->nps[10] = elmntp->nps[n_h[7]] ;
			theElp->nps[11] = elmntp->nps[n_h[4]] ;
	
		}
		if ((elmntp->nps[n_h[8]] == NULL) && ((elmntp->nps[n_h[14]] == NULL)
												|| (elmntp->nps[n_h[15]] == NULL))){
			theElp->vtype = 2351;
			theElp->nnds = 11 ;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
  			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[4]] ;
			theElp->nps[5] = elmntp->nps[n_h[5]] ;
			theElp->nps[6] = elmntp->nps[n_h[6]] ;
			theElp->nps[7] = elmntp->nps[n_h[7]] ;
			theElp->nps[8] = elmntp->nps[n_h[9]] ;
			theElp->nps[9] = elmntp->nps[n_h[10]] ;
			theElp->nps[10] = elmntp->nps[n_h[11]] ;
	
		}
		if ((elmntp->nps[n_h[13]] == NULL) && ((elmntp->nps[n_h[5]] == NULL)
				|| (elmntp->nps[n_h[6]] == NULL))){
			theElp->vtype = 2353;
			theElp->nnds = 11 ;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
  			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[11]] ;
			theElp->nps[5] = elmntp->nps[n_h[15]] ;
			theElp->nps[6] = elmntp->nps[n_h[14]] ;
			theElp->nps[7] = elmntp->nps[n_h[10]] ;
			theElp->nps[8] = elmntp->nps[n_h[9]] ;
			theElp->nps[9] = elmntp->nps[n_h[7]] ;
			theElp->nps[10] = elmntp->nps[n_h[4]] ;
		}
		if((elmntp->nps[n_h[8]] == NULL) && (elmntp->nps[n_h[13]] == NULL)){
			theElp->vtype = 2362 ;
			theElp->nnds = 13 ;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
  			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[4]] ;
			theElp->nps[5] = elmntp->nps[n_h[5]] ;
			theElp->nps[6] = elmntp->nps[n_h[6]] ;
			theElp->nps[7] = elmntp->nps[n_h[7]] ;
			theElp->nps[8] = elmntp->nps[n_h[9]] ;
			theElp->nps[9] = elmntp->nps[n_h[10]] ;
			theElp->nps[10] = elmntp->nps[n_h[14]] ;
			theElp->nps[11] = elmntp->nps[n_h[15]] ;
			theElp->nps[12] = elmntp->nps[n_h[11]] ;

		}
		if((elmntp->nps[n_h[9]] == NULL)){
		 	theElp->vtype = 2352 ;
			theElp->nnds = 12 ;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[4]] ;
			theElp->nps[5] = elmntp->nps[n_h[5]] ;
			theElp->nps[6] = elmntp->nps[n_h[6]] ;
			theElp->nps[7] = elmntp->nps[n_h[7]] ;
			theElp->nps[8] = elmntp->nps[n_h[10]] ;
			theElp->nps[9] = elmntp->nps[n_h[14]] ;
			theElp->nps[10] = elmntp->nps[n_h[15]] ;
			theElp->nps[11] = elmntp->nps[n_h[11]] ;
		}
	
		else if(((elmntp->nps[n_h[6]] == NULL)	|| (elmntp->nps[n_h[5]] == NULL))
				&& ((elmntp->nps[n_h[14]] == NULL)	|| (elmntp->nps[n_h[15]] == NULL))){
			theElp->vtype = 232 ;
			theElp->nnds = 9 ;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[4]] ;
			theElp->nps[5] = elmntp->nps[n_h[7]] ;
			theElp->nps[6] = elmntp->nps[n_h[9]] ;
			theElp->nps[7] = elmntp->nps[n_h[10]] ;
			theElp->nps[8] = elmntp->nps[n_h[11]] ;
					
		}
		if(((elmntp->nps[n_h[14]] == NULL)	|| (elmntp->nps[n_h[15]] == NULL))
		 && (elmntp->nps[n_h[9]] == NULL)){
			theElp->vtype = 2341 ;
			theElp->nnds = 10 ;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[4]] ;
			theElp->nps[5] = elmntp->nps[n_h[5]] ;
			theElp->nps[6] = elmntp->nps[n_h[6]] ;
			theElp->nps[7] = elmntp->nps[n_h[7]] ;
			theElp->nps[8] = elmntp->nps[n_h[10]] ;
			theElp->nps[9] = elmntp->nps[n_h[11]] ;
		
			if((elmntp->nps[n_h[6]] == NULL)	|| (elmntp->nps[n_h[5]] == NULL)){
				theElp->vtype = 2311 ;
				theElp->nnds = 8 ;
				theElp->nps[0] = elmntp->nps[n_h[0]] ;
				theElp->nps[1] = elmntp->nps[n_h[1]] ;
				theElp->nps[2] = elmntp->nps[n_h[2]] ;
				theElp->nps[3] = elmntp->nps[n_h[3]] ;
				theElp->nps[4] = elmntp->nps[n_h[4]] ;
				theElp->nps[5] = elmntp->nps[n_h[7]] ;
				theElp->nps[6] = elmntp->nps[n_h[10]] ;
				theElp->nps[7] = elmntp->nps[n_h[11]] ;
			}
		}
		if(((elmntp->nps[n_h[6]] == NULL)	|| (elmntp->nps[n_h[5]] == NULL))
		 && (elmntp->nps[n_h[9]] == NULL)){
			theElp->vtype = 2342 ;
			theElp->nnds = 10 ;
			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[11]] ;
			theElp->nps[5] = elmntp->nps[n_h[15]] ;
			theElp->nps[6] = elmntp->nps[n_h[14]] ;
			theElp->nps[7] = elmntp->nps[n_h[10]] ;
			theElp->nps[8] = elmntp->nps[n_h[7]] ;
			theElp->nps[9] = elmntp->nps[n_h[4]] ;
		}
		if((elmntp->nps[n_h[10]] == NULL) || (elmntp->nps[n_h[11]] == NULL)){
			theElp->vtype = 233 ;
			theElp->nnds = 8 ;
  			theElp->nps[0] = elmntp->nps[n_h[0]] ;
			theElp->nps[1] = elmntp->nps[n_h[1]] ;
			theElp->nps[2] = elmntp->nps[n_h[2]] ;
			theElp->nps[3] = elmntp->nps[n_h[3]] ;
			theElp->nps[4] = elmntp->nps[n_h[4]] ;
			theElp->nps[5] = elmntp->nps[n_h[5]] ;
			theElp->nps[6] = elmntp->nps[n_h[6]] ;
			theElp->nps[7] = elmntp->nps[n_h[7]] ;
		
			if ((elmntp->nps[n_h[6]] == NULL)	|| (elmntp->nps[n_h[5]] == NULL)){
				theElp->vtype = 230;
				theElp->nnds = 6 ;
				theElp->nps[0] = elmntp->nps[n_h[0]] ;
				theElp->nps[1] = elmntp->nps[n_h[1]] ;
				theElp->nps[2] = elmntp->nps[n_h[2]] ;
				theElp->nps[3] = elmntp->nps[n_h[3]] ;
				theElp->nps[4] = elmntp->nps[n_h[4]] ;
				theElp->nps[5] = elmntp->nps[n_h[7]] ;
			}
		}
	
		if((elmntp->nps[n_h[4]] == NULL) || (elmntp->nps[n_h[7]] == NULL)){
			theElp->vtype = 233 ;
			theElp->nnds = 8 ;
  		 	theElp->nps[0] = elmntp->nps[n_h[1]] ;
		 	theElp->nps[1] = elmntp->nps[n_h[2]] ;
		 	theElp->nps[2] = elmntp->nps[n_h[3]] ;
		 	theElp->nps[3] = elmntp->nps[n_h[0]] ;
		 	theElp->nps[4] = elmntp->nps[n_h[10]] ;
		 	theElp->nps[5] = elmntp->nps[n_h[14]] ;
		 	theElp->nps[6] = elmntp->nps[n_h[15]] ;
		 	theElp->nps[7] = elmntp->nps[n_h[11]] ;
		
			if ((elmntp->nps[n_h[14]] == NULL)	|| (elmntp->nps[n_h[15]] == NULL)){
				theElp->vtype = 230;
				theElp->nnds = 6 ;
				theElp->nps[0] = elmntp->nps[n_h[1]] ;
				theElp->nps[1] = elmntp->nps[n_h[2]] ;
				theElp->nps[2] = elmntp->nps[n_h[3]] ;
				theElp->nps[3] = elmntp->nps[n_h[0]] ;
				theElp->nps[4] = elmntp->nps[n_h[10]] ;
				theElp->nps[5] = elmntp->nps[n_h[11]] ;
			}
		}
		if(((elmntp->nps[n_h[10]] == NULL)	|| (elmntp->nps[n_h[11]] == NULL))
		 	&& ((elmntp->nps[n_h[7]] == NULL)	|| (elmntp->nps[n_h[4]] == NULL))){
		 	
		 	theElp->vtype = 211;
		 	theElp->nnds = 4 ;
		 	theElp->nps[0] = elmntp->nps[n_h[0]] ;
		 	theElp->nps[1] = elmntp->nps[n_h[1]] ;
		 	theElp->nps[2] = elmntp->nps[n_h[2]] ;
		 	theElp->nps[3] = elmntp->nps[n_h[3]] ;
		}			

	return(-1) ;
}
