//#pragma options cpluscmt	
#include "Fe_PS.h"
#include "ctype.h"

extern struct control  	N ;
extern struct pointers 	gp;
extern struct RegMesh	Mesh ;
extern int	Nndso, Nelmso, Nbelmso ;

int		n1L, n2L ;
double	x0, x1, yzero, yone, L, ar, defndp, defelp ;

int		MakeMesh(theMesh)
struct	RegMesh	*theMesh ;

{
	int		err,i, j, size ;
	long int	bigint ;
	struct node	*anodep, *np ;
		
	x0 =-1.0 ;
	yzero = -1.0 ;
	x1 = -1.0 + theMesh->nbx * 2.0 ;
	yone = -1.0 + theMesh->nby * 2.0 ;
	n1L = theMesh->ny ;
	n2L = theMesh->nx ;
	L = (yone - yzero)/n1L ;
	ar = (x1 - x0)/(3.0 * n2L) ;
	if(theMesh->eltype == 216 || theMesh->eltype == 226 || theMesh->eltype == 227)
		n2L += 1 ;
	if(gp.N != NULL) {
		free(gp.N) ;
		free(gp.iptrs) ;
		clear_data() ;
	}
	
	j = 0 ;
	for(i=0;i<1;i++) {
		j++ ;
		switch (theMesh->eltype) {
			case 111 :
				N.nodes = n2L+1 ;
				N.dims = 1 ;
				break ;
			case 119 :
				N.nodes = n2L+1 ;
				N.dims = 1 ;
				break ;
			case 121 :
				N.nodes = 2*n2L+1 ;
				N.dims = 1 ;
				break ;
			case 129 :
				N.nodes = n2L + 1;
				N.dims = 1 ;
				break ;
			case 139 :
				N.nodes = n2L + 1;
				N.dims = 1 ;
				break ;
			case 131 :
				N.nodes = 3*n2L+1 ;
				N.dims = 1 ;
				break ;
			case 210 :
				N.nodes = (n1L+1) * (n2L+1) ;
				break ;
			case 211 :
				N.nodes = (n1L+1) * (n2L+1) ;
				break ;
			case 216 :
				N.nodes = (n1L+1) * (n2L+1) ;
				break ;
			case 221 :
				N.nodes = 3*n1L*n2L + 2*(n1L + n2L) + 1 ;
				break ;
			case 226 :
				N.nodes = 3*n1L*n2L + 2*(n1L + n2L) + 1 ;
				break ;
			case 220 :
				N.nodes = (2*n1L + 1) * (2*n2L + 1) ;
				break ;
			case 222 :
				N.nodes = (2*n1L + 1) * (2*n2L + 1) ;
				break ;
			case 227 :
				N.nodes = (2*n1L + 1) * (2*n2L + 1) ;
				break ;
			case 229 :
				N.nodes = (n1L+1) * (n2L+1) ;
				break ;
			case 231 :
				N.nodes = 5*n1L*n2L + 3*(n1L + n2L) + 1 ;
				break ;
			default :
				return(-1) ;
		}
		bigint = (N.nodes*5)/4 * sizeof(struct node *) ; 
		if(N.nodes < 2 || (gp.iptrs = (struct node **) malloc( bigint )) == NULL) break ;
		bigint = (long int) N.nodes * (long int) sizeof(struct node) ;
		size = sizeof(struct node) ;
		printf("memory requested = %ld, %d\n",bigint,size) ;
		if(N.nodes < 2 || (gp.N = (struct node *) malloc( bigint )) == NULL) break ;
		j++ ;
		N.elms = n1L * n2L ;
		if((theMesh->eltype == 210) || (theMesh->eltype == 220) )
			N.elms *= 2 ;
		bigint = (long int) N.elms * (long int) sizeof(struct element) ;
		size = sizeof(struct element) ;
		printf("memory requested = %ld, %d\n",bigint,size) ;
		if(gp.El != NULL)
			free(gp.El) ;
		if(N.elms < 1 || (gp.El = (struct element *) malloc( bigint) ) == NULL) break ;
		j++ ;
		N.belms = 2 * (n1L + n2L) ;
		if((theMesh->eltype == 111)||(theMesh->eltype == 121)||(theMesh->eltype == 119)||
			(theMesh->eltype == 129)||(theMesh->eltype == 131)||
			(theMesh->eltype == 139))
				N.belms = 2 ;
		bigint = (long int) N.belms * (long int) sizeof(struct belement) ;
		if(gp.B != NULL)
			free(gp.B) ;
		if(N.belms < 1 || (gp.B = (struct belement *) malloc( bigint )) == NULL) break ;
		Nndso = N.nodes ;
		Nelmso = N.elms ;
		Nbelmso = N.belms ;
 		put_control(stdout) ;
		
		for(i=1;i<=n2L;i++)
			for(j=1;j<=n1L;j++)
				make_elmnt(i,j,theMesh->eltype) ;

		set_bvalues() ;/*in feio.c*/ 
		j = 0 ;
	}



	err = j ;
	if(err != 0)
		printf(" Control Variable Error %d\n",err) ;
	anodep = gp.N ;
/*
	for(i=0;i<N.nodes;i++) {
		printf(" index %d name %d\n",anodep->i,anodep->n) ;
		printf("       %d      %d\n",gp.iptrs[i]->i,gp.iptrs[i]->n) ;
		anodep = anodep->nextnp ;
	}
*/

	return(err);
}

int		make_elmnt(i,j,eltype)
int		i, j, eltype ;

{
	int	 k, index, teltype  ;
	double 	x, y ;
	struct node  *nodep, *alloc_node() ;
	struct element	*elmntp, *alloc_elm() ;
	struct belement	*belmntp, *alloc_belm() ;
	
	teltype = eltype ;
	if(eltype == 210)
		teltype = 211 ;
	if(eltype == 220)
		teltype = 222 ;
	if(eltype == 119)
		teltype = 111 ;
	switch (teltype) {
		case 111 :
			index = i - 1 ;
			x = x0 + 3 * (i-1) * ar ;
			y = yzero ;
			nodep = alloc_node(index,x,y) ;
			if(i == 1) {
				index = 0 ;
				belmntp = alloc_belm(index,1) ;
				belmntp->nps[0] = gp.N ; 
				leftb(belmntp) ;
			}
			if(i == n2L) {
				index = i ;
				x = x1 ;
				y = yzero ;
				nodep = alloc_node(index,x,y) ;
				index = 1 ;
				belmntp = alloc_belm(index,1) ;
				belmntp->nps[0] = gp.N + i ; 
				rightb(belmntp) ;
				belmntp->nextbelp = NULL ;
			}
			index = i - 1 ;
			elmntp = alloc_elm(index,eltype) ;
			if(eltype == 119){
				elmntp->vtype = 111;
				elmntp->gtype = 119;
			}
			elmntp->nps[0] = gp.N + i-1 ; 
			elmntp->nps[1] = gp.N + i ; 
			if(index == N.elms-1)
				elmntp->nextelp = NULL ;
			break ;
			
		case 121 :
			index = 2*(i - 1) ;
			x = x0 + 3.0 * (i-1) * ar ;
			y = yzero ;
			nodep = alloc_node(index,x,y) ;
			index ++ ;
			x +=  1.5  * ar ;
			y = yzero ;
			nodep = alloc_node(index,x,y) ;
			if(i == 1) {
				index = 0 ;
				belmntp = alloc_belm(index,1) ;
				belmntp->nps[0] = gp.N ; 
				leftb(belmntp) ;
			}
			if(i == n2L) {
				index = 2*i ;
				x = x1 ;
				y = yzero ;
				nodep = alloc_node(index,x,y) ;
				index = 1 ;
				belmntp = alloc_belm(index,1) ;
				belmntp->nps[0] = gp.N + 2*i ; 
				rightb(belmntp) ;
				belmntp->nextbelp = NULL ;
			}
			index = i - 1 ;
			elmntp = alloc_elm(index,eltype) ;
			elmntp->nps[0] = gp.N + 2*(i-1) ; 
			elmntp->nps[1] = gp.N + 2*i - 1 ; 
			elmntp->nps[2] = gp.N + 2*i ; 
			if(index == N.elms-1)
				elmntp->nextelp = NULL ;
			break ;
			
		case 129 :
			index = i - 1 ;
			x = x0 + 3.0 * (i-1) * ar ;
			y = yzero ;
			nodep = alloc_node(index,x,y) ;
			if(i == 1) {
				index = 0 ;
				belmntp = alloc_belm(index,1) ;
				belmntp->nps[0] = gp.N ; 
				leftb(belmntp) ;
			}
			if(i == n2L) {
				index = i ;
				x = x1 ;
				y = yzero ;
				nodep = alloc_node(index,x,y) ;
				index = 1 ;
				belmntp = alloc_belm(index,1) ;
				belmntp->nps[0] = gp.N + i ; 
				rightb(belmntp) ;
				belmntp->nextbelp = NULL ;
			}
			index = i - 1 ;
			elmntp = alloc_elm(index,eltype) ;
			elmntp->nps[0] = gp.N + i-1 ; 
			elmntp->nps[1] = gp.N + i ;
			if(i != n2L) 
				elmntp->nps[3] = gp.N + i + 1 ;
			else 
				elmntp->nps[3] = NULL ;
			if(i != 1) 
				elmntp->nps[2] = gp.N + i - 2 ; 
			else 
				elmntp->nps[2] = NULL ;
			if(index == N.elms-1)
				elmntp->nextelp = NULL ;
			break ;
			
		case 139 :
			index = i - 1 ;
			x = x0 + 3 * (i-1) * ar ;
			y = yzero ;
			nodep = alloc_node(index,x,y) ;
			if(i == 1) {
				index = 0 ;
				belmntp = alloc_belm(index,1) ;
				belmntp->nps[0] = gp.N ; 
				leftb(belmntp) ;
			}
			if(i == n2L) {
				index = i ;
				x = x1 ;
				y = yzero ;
				nodep = alloc_node(index,x,y) ;
				index = 1 ;
				belmntp = alloc_belm(index,1) ;
				belmntp->nps[0] = gp.N + i ; 
				rightb(belmntp) ;
				belmntp->nextbelp = NULL ;
			}
			index = i - 1 ;
			elmntp = alloc_elm(index,eltype) ;
			elmntp->nps[0] = gp.N + i-1 ; 
			elmntp->nps[1] = gp.N + i ;
			if(i != n2L) 
				elmntp->nps[3] = gp.N + i + 1 ;
			else 
				elmntp->nps[3] = NULL ;
			if(i == 1) {
				elmntp->nps[2] = NULL ;
				elmntp->nps[3] = NULL ;
			}
			else if(i == 2) {
				elmntp->nps[2] = gp.N + i - 2 ; 
				elmntp->nps[3] = NULL ;
			}
			else  {
				elmntp->nps[2] = gp.N + i - 2 ; 
				elmntp->nps[3] = gp.N + i - 3  ;
			}
			if(i == n2L) {
				elmntp->nps[4] = NULL ;
				elmntp->nps[5] = NULL ;
			}
			else if(i == n2L-1) {
				elmntp->nps[4] = gp.N + i + 1 ; 
				elmntp->nps[5] = NULL ;
			}
			else  {
				elmntp->nps[4] = gp.N + i + 1 ; 
				elmntp->nps[5] = gp.N + i + 2  ;
			}
			
			if(index == N.elms-1)
				elmntp->nextelp = NULL ;
			break ;
			
		case 211 :
			index = (i-1) * (n1L+1) + j - 1 ;
			if(i == 1)
				x = x0 ;
			else 
				x = x0 + 3.0 * (i-1) * ar ;
			y = yzero + (j-1) * L *1.0;
//printf("i=%d\tj=%d\tx=%g\ty=%g\n",i,j,x,y);
			nodep = alloc_node(index,x,y) ;
			if(j == 1) {
				index = i - 1 ;
				belmntp = alloc_belm(index,111) ;
				belmntp->nps[0] = gp.N + (i-1) * (n1L + 1)  ; 
				belmntp->nps[1] = gp.N + i * (n1L + 1)  ; 
				bottomb(belmntp) ;
			}
			if(i == 1) {
				index = 2 * n2L + 2 * n1L - j ;
				belmntp = alloc_belm(index,111) ;
				belmntp->nps[0] = gp.N + j ; 
				belmntp->nps[1] = gp.N + j - 1 ; 
				leftb(belmntp) ;
				if(j == 1)
					belmntp->nextbelp = NULL ;
			}
			if(j == n1L) {
				index = (i-1) * (n1L+1) + j ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
				index = 2 * n2L + n1L - i ;
				belmntp = alloc_belm(index,111) ;
				belmntp->nps[0] = gp.N + (1 + i ) * (n1L + 1) -1 ; 
				belmntp->nps[1] =  gp.N + i  * (n1L + 1) -1 ; 
				topb(belmntp) ;
			}
			if(i == n2L) {
				index = i * (n1L+1) + j - 1 ;
				x = x1 ;
				y = yzero + (j-1) * L *1.0;
				nodep = alloc_node(index,x,y) ;
				index =  n2L + j - 1 ;
				belmntp = alloc_belm(index,111) ;
				belmntp->nps[0] = gp.N + n2L * (n1L + 1) + j -1 ; 
				belmntp->nps[1] = gp.N + n2L * (n1L + 1) + j ; 
				rightb(belmntp) ;
			}
			if(i == n2L && j == n1L) {
				index = i * (n1L+1) + j  ;
				x = x1 ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
			}
			if(eltype == 211) {
				index = (i - 1) * n1L + j - 1;
				elmntp = alloc_elm(index,eltype) ;
				elmntp->nps[3] = gp.N + (i-1) * (n1L + 1) + j -1 ; 
				elmntp->nps[0] = gp.N + i * (n1L + 1) + j-1 ; 
				elmntp->nps[1] = gp.N + i * (n1L + 1) + j ; 
				elmntp->nps[2] = gp.N + (i-1) * (n1L + 1) + j ;
			}
			else {
				index = 2 * ((i - 1) * n1L + j - 1) ;
				elmntp = alloc_elm(index,eltype) ;
				elmntp->nps[0] = gp.N + (i-1) * (n1L + 1) + j - 1 ; 
				elmntp->nps[1] = gp.N + i * (n1L + 1) + j - 1 ; 
				elmntp->nps[2] = gp.N + (i-1) * (n1L + 1) + j ;
				index++ ;
				elmntp = alloc_elm(index,eltype) ;
				elmntp->nps[2] = gp.N + i * (n1L + 1) + j - 1 ; 
				elmntp->nps[0] = gp.N + i * (n1L + 1) + j ; 
				elmntp->nps[1] = gp.N + (i-1) * (n1L + 1) + j ;
			}
			if(index == N.elms-1)
				elmntp->nextelp = NULL ;
			break ;
			
		case 229 :
			index = (i-1) * (n1L+1) + j - 1 ;
			if(i == 1)
				x = x0 ;
			else 
				x = x0 + 3.0 * (i-1) * ar ;
			y = yzero + (j-1) * L ;
			nodep = alloc_node(index,x,y) ;
			if(j == 1) {
				index = i - 1 ;
				belmntp = alloc_belm(index,111) ;
				belmntp->nps[0] = gp.N + (i-1) * (n1L + 1) ; 
				belmntp->nps[1] = gp.N + i * (n1L + 1) ; 
				bottomb(belmntp) ;
			}
			if(i == 1) {
				index = 2 * n2L + 2 * n1L - j ;
				belmntp = alloc_belm(index,111) ;
				belmntp->nps[0] = gp.N + j; 
				belmntp->nps[1] = gp.N + j - 1 ; 
				leftb(belmntp) ;
				if(j == 1)
					belmntp->nextbelp = NULL ;
			}
			if(j == n1L) {
				index = (i-1) * (n1L+1) + j ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
				index = 2 * n2L + n1L - i ;
				belmntp = alloc_belm(index,111) ;
				belmntp->nps[0] = gp.N + (1 + i ) * (n1L + 1) - 1  ; 
				belmntp->nps[1] =  gp.N + i  * (n1L + 1) - 1 ; 
				topb(belmntp) ;
			}
			if(i == n2L) {
				index = i * (n1L+1) + j - 1 ;
				x = x1 ;
				y = yzero + (j-1) * L ;
				nodep = alloc_node(index,x,y) ;
				index =  n2L + j - 1 ;
				belmntp = alloc_belm(index,111) ;
				belmntp->nps[0] = gp.N + n2L * (n1L + 1) + j -1 ; 
				belmntp->nps[1] = gp.N + n2L * (n1L + 1) + j ; 
				rightb(belmntp) ;
			}
			if(i == n2L && j == n1L) {
				index = i * (n1L+1) + j  ;
				x = x1 ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
			}
			index = (i - 1) * n1L + j - 1;
			
			elmntp = alloc_elm(index,eltype) ;
				
			elmntp->nps[3] = gp.N + (i-1) * (n1L + 1) + j - 1; 
			elmntp->nps[0] = gp.N + i * (n1L + 1) + j - 1; 
			elmntp->nps[1] = gp.N + i * (n1L + 1) + j ; 
			elmntp->nps[2] = gp.N + (i-1) * (n1L + 1) + j ;
			if(i != n2L && j != 1)
				elmntp->nps[4] = gp.N + (i+1) * (n1L + 1) + j-2 ;
			else
				elmntp->nps[4] = NULL ;
			if(i != n2L)
				elmntp->nps[5] = gp.N + (i+1) * (n1L + 1) + j-1 ;
			else
				elmntp->nps[5] = NULL ;
			if(i != n2L)
				elmntp->nps[6] = gp.N + (i+1) * (n1L + 1) + j ;
			else
				elmntp->nps[6] = NULL ;
			if(i != n2L && j != n1L)
				elmntp->nps[7] = gp.N + (i+1) * (n1L + 1) + j+1 ;
			else
				elmntp->nps[7] = NULL ;
			if(j != n1L)
				elmntp->nps[8] = gp.N + i * (n1L + 1) + j + 1 ;
			else
				elmntp->nps[8] = NULL ;
			if(j != n1L)
				elmntp->nps[9] = gp.N + (i-1) * (n1L + 1) + j + 1 ;
			else
				elmntp->nps[9] = NULL ;
			if(i != 1 && j != n1L)
				elmntp->nps[10] = gp.N + (i-2) * (n1L + 1) + j + 1 ;
			else
				elmntp->nps[10] = NULL ;
			if(i != 1 )
				elmntp->nps[11] = gp.N + (i-2) * (n1L + 1) + j ;
			else
				elmntp->nps[11] = NULL ;
			if(i != 1 )
				elmntp->nps[12] = gp.N + (i-2) * (n1L + 1) + j - 1 ;
			else
				elmntp->nps[12] = NULL ;
			if(i != 1 && j != 1)
				elmntp->nps[13] = gp.N + (i-2) * (n1L + 1) + j - 2 ;
			else
				elmntp->nps[13] = NULL ;
			if(j != 1)
				elmntp->nps[14] = gp.N + (i-1) * (n1L + 1) + j - 2 ;
			else
				elmntp->nps[14] = NULL ;
			if(j != 1)
				elmntp->nps[15] = gp.N + i * (n1L + 1) + j - 2 ;
			else
				elmntp->nps[15] = NULL ;
			
			if(index == N.elms-1)
				elmntp->nextelp = NULL ;
			break ;
			
		case 221 :
			index = (i-1) * (3*n1L+2) + 2*j - 2 ;
			x = x0 + 3 * (i-1) * ar ;
			y = yzero + (j-1) * L ;
			nodep = alloc_node(index,x,y) ;
			index++ ;
			y = yzero + (j-1) *L + 0.5 *L ;
			nodep = alloc_node(index,x,y) ;
			index = (i-1) * (3*n1L+2) + 2*(n1L+1) + j - 2 ;
			x = x0 + 3 * (i-1) * ar + 1.5 * ar ;
			y = yzero + (j-1) * L ;
			nodep = alloc_node(index,x,y) ;
			if(j == 1) {
				index = i - 1 ;
				belmntp = alloc_belm(index,121) ;
				belmntp->nps[0] = gp.N + (i-1) * (3*n1L + 2)  ; 
				belmntp->nps[1] = gp.N + (i-1) * (3*n1L+2) + 2*(n1L+1) - 1 ; 
				belmntp->nps[2] = gp.N + i * (3*n1L + 2) ; 
				bottomb(belmntp) ;
			}
			if(i == 1) {
				index = 2 * n2L + 2 * n1L - j ;
				belmntp = alloc_belm(index,121) ;
				belmntp->nps[0] = gp.N + 2*j ; 
				belmntp->nps[1] = gp.N + 2*j -1; 
				belmntp->nps[2] = gp.N + 2*j - 2 ; 
				leftb(belmntp) ;
				if(j == 1)
					belmntp->nextbelp = NULL ;
			}
			if(j == n1L) {
				index = (i-1) * (3*n1L+2) + 2*n1L ;
				x = x0 + 3 * (i-1) * ar ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
				nodep->u[1] = 1.0 ;
				index = i * (3*n1L+2) - 1  ;
				x = x0 + 3 * (i-1) * ar + 1.5*ar ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
				nodep->u[1] = 1.0 ;
				index = 2 * n2L + n1L - i ;
				belmntp = alloc_belm(index,121) ;
				belmntp->nps[0] = gp.N + i * (3*n1L+2) + 2*n1L  ; 
				belmntp->nps[1] = gp.N + i * (3*n1L+2) - 1 ; 
				belmntp->nps[2] = gp.N + (i-1) * (3*n1L+2) + 2*n1L ; 
				topb(belmntp) ;
			}
			if(i == n2L) {
				index = i * (3*n1L+2) + 2 * j - 2 ;
				x = x1 ;
				y = yzero + (j-1) * L ;
				nodep = alloc_node(index,x,y) ;
				index = i * (3*n1L+2) + 2*j - 1 ;
				x = x1 ;
				y = yzero + (j-1) * L + 0.5 * L ;
				nodep = alloc_node(index,x,y) ;
				index =  n2L + j - 1 ;
				belmntp = alloc_belm(index,121) ;
				belmntp->nps[0] = gp.N + n2L * (3*n1L + 2) + 2*j - 2 ; 
				belmntp->nps[1] = gp.N + n2L * (3*n1L + 2) + 2*j - 1; 
				belmntp->nps[2] = gp.N + n2L * (3*n1L + 2) + 2*j ;
				rightb(belmntp) ;
			}
			if(i == n2L && j == n1L) {
				index = N.nodes - 1  ;
				x = x1 ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
			}
			index = (i - 1) * n1L + j - 1;
			elmntp = alloc_elm(index,eltype) ;
			if(index == N.elms-1)
				elmntp->nextelp = NULL ;
			elmntp->nps[0] = gp.N + i * (3*n1L + 2) + 2*j - 2 ; 
			elmntp->nps[1] = gp.N + i * (3*n1L + 2) + 2*j -1; 
			elmntp->nps[2] = gp.N + i * (3*n1L + 2) + 2*j ; 
			elmntp->nps[3] = gp.N + (i-1) * (3*n1L + 2) + 2*n1L + j + 1 ; 
			elmntp->nps[4] = gp.N + (i-1) * (3*n1L + 2) + 2*j ; 
			elmntp->nps[5] = gp.N + (i-1) * (3*n1L + 2) + 2*j - 1 ; 
			elmntp->nps[6] = gp.N + (i-1) * (3*n1L + 2) + 2*j - 2 ; 
			elmntp->nps[7] = gp.N + (i-1) * (3*n1L + 2) + 2*n1L + j ; 
			break ;
			
			
		case 222 :
			index = (2 * i - 2) * (2*n1L+1) + 2*j - 2 ;
			x = x0 + 3 * (i-1) * ar ;
			y = yzero + (j-1) * L ;
			nodep = alloc_node(index,x,y) ;
			index++ ;
			y = yzero + (j-1) *L + 0.5 *L ;
			nodep = alloc_node(index,x,y) ;
			index = (2 * i - 1) * (2*n1L+1) + 2*j - 2 ;
			x = x0 + 3 * (i-1) * ar + 1.5 * ar ;
			y = yzero + (j-1) * L ;
			nodep = alloc_node(index,x,y) ;
			index++ ;
			y = yzero + (j-1) *L + 0.5 *L ;
			nodep = alloc_node(index,x,y) ;
			if(j == 1) {
				index = i - 1 ;
				belmntp = alloc_belm(index,121) ;
				belmntp->nps[0] = gp.N + (2 * i - 2) * (2*n1L+1) ; 
				belmntp->nps[1] = gp.N + (2 * i - 1) * (2*n1L+1) ; 
				belmntp->nps[2] = gp.N + (2 * i) * (2*n1L+1) ; 
				bottomb(belmntp) ;
			}
			if(i == 1) {
				index = 2 * n2L + 2 * n1L - j ;
				belmntp = alloc_belm(index,121) ;
				belmntp->nps[0] = gp.N + 2*j ; 
				belmntp->nps[1] = gp.N + 2*j - 1; 
				belmntp->nps[2] = gp.N + 2*j - 2 ; 
				leftb(belmntp) ;
				if(j == 1)
					belmntp->nextbelp = NULL ;
			}
			if(j == n1L) {
				index = (2 * i - 1) * (2*n1L+1) - 1 ;
				x = x0 + 3 * (i-1) * ar ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
				nodep->u[1] = 1.0 ;
				index = (2 * i) * (2*n1L+1) - 1  ;
				x = x0 + 3 * (i-1) * ar + 1.5*ar ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
				nodep->u[1] = 1.0 ;
				index = 2 * n2L + n1L - i ;
				belmntp = alloc_belm(index,121) ;
				belmntp->nps[0] = gp.N + (2 * i + 1) * (2*n1L+1) - 1 ; 
				belmntp->nps[1] = gp.N + (2 * i) * (2*n1L+1) - 1 ; 
				belmntp->nps[2] = gp.N + (2 * i - 1) * (2*n1L+1) - 1 ; 
				topb(belmntp) ;
			}
			if(i == n2L) {
				index = (2 * i) * (2*n1L+1) + 2*j - 2 ;
				x = x1 ;
				y = yzero + (j-1) * L ;
				nodep = alloc_node(index,x,y) ;
				index = (2 * i) * (2*n1L+1) + 2*j - 1 ;
				x = x1 ;
				y = yzero + (j-1) * L + 0.5 * L ;
				nodep = alloc_node(index,x,y) ;
				index =  n2L + j - 1 ;
				belmntp = alloc_belm(index,121) ;
				belmntp->nps[0] = gp.N + (2 * i) * (2*n1L+1) + 2*j - 2 ; 
				belmntp->nps[1] = gp.N + (2 * i) * (2*n1L+1) + 2*j - 1 ; 
				belmntp->nps[2] = gp.N + (2 * i) * (2*n1L+1) + 2*j ;
				rightb(belmntp) ;
			}
			if(i == n2L && j == n1L) {
				index = N.nodes - 1  ;
				x = x1 ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
			}
			if(eltype == 222) {
				index = (i - 1) * n1L + j - 1;
				elmntp = alloc_elm(index,eltype) ;
				elmntp->nps[0] = gp.N + (2 * i) * (2*n1L+1) + 2*j - 2 ; 
				elmntp->nps[1] = gp.N + (2 * i) * (2*n1L+1) + 2*j - 1 ; 
				elmntp->nps[2] = gp.N + (2 * i) * (2*n1L+1) + 2*j ; 
				elmntp->nps[3] = gp.N + (2 * i - 1) * (2*n1L+1) + 2*j ; 
				elmntp->nps[4] = gp.N + (2 * i - 2) * (2*n1L+1) + 2*j ; 
				elmntp->nps[5] = gp.N + (2 * i - 2) * (2*n1L+1) + 2*j - 1 ; 
				elmntp->nps[6] = gp.N + (2 * i - 2) * (2*n1L+1) + 2*j - 2 ; 
				elmntp->nps[7] = gp.N + (2 * i - 1) * (2*n1L+1) + 2*j - 2 ; 
				elmntp->nps[8] = gp.N + (2 * i - 1) * (2*n1L+1) + 2*j - 1 ; 
			}
			else	{
				index = 2 * ((i - 1) * n1L + j - 1) ;
				elmntp = alloc_elm(index,eltype) ;
				elmntp->nps[2] = gp.N + (2 * i) * (2*n1L+1) + 2*j - 2 ; 
				elmntp->nps[4] = gp.N + (2 * i - 2) * (2*n1L+1) + 2*j ; 
				elmntp->nps[5] = gp.N + (2 * i - 2) * (2*n1L+1) + 2*j - 1 ; 
				elmntp->nps[0] = gp.N + (2 * i - 2) * (2*n1L+1) + 2*j - 2 ; 
				elmntp->nps[1] = gp.N + (2 * i - 1) * (2*n1L+1) + 2*j - 2 ; 
				elmntp->nps[3] = gp.N + (2 * i - 1) * (2*n1L+1) + 2*j - 1 ; 
				index++ ;
				elmntp = alloc_elm(index,eltype) ;
				elmntp->nps[0] = gp.N + (2 * i) * (2*n1L+1) + 2*j - 2 ; 
				elmntp->nps[1] = gp.N + (2 * i) * (2*n1L+1) + 2*j - 1 ; 
				elmntp->nps[2] = gp.N + (2 * i) * (2*n1L+1) + 2*j ; 
				elmntp->nps[3] = gp.N + (2 * i - 1) * (2*n1L+1) + 2*j ; 
				elmntp->nps[4] = gp.N + (2 * i - 2) * (2*n1L+1) + 2*j ; 
				elmntp->nps[5] = gp.N + (2 * i - 1) * (2*n1L+1) + 2*j - 1; 
			}	
			if(index == N.elms-1)
				elmntp->nextelp = NULL ;
			break ;
			
			
		case 231 :
			index = (i-1) * (5*n1L+3) + 3*j - 3 ;
			x = x0 + 3 * (i-1) * ar ;
			y = yzero + (j-1) * L ;
			nodep = alloc_node(index,x,y) ;
			index++ ;
			y = yzero + (j-1) *L + 0.33333333 *L ;
			nodep = alloc_node(index,x,y) ;
			index++ ;
			y = yzero + (j-1) *L + 0.66666667 *L ;
			nodep = alloc_node(index,x,y) ;
			
			index = (i-1) * (5*n1L+3) + 3*n1L + j ;
			x = x0 + 3 * (i-1) * ar + 1.0 * ar ;
			y = yzero + (j-1) * L ;
			nodep = alloc_node(index,x,y) ;
			index = (i-1) * (5*n1L+3) + 4*n1L + j + 1 ;
			x = x0 + 3 * (i-1) * ar + 2.0 * ar ;
			y = yzero + (j-1) * L ;
			nodep = alloc_node(index,x,y) ;
			if(j == 1) {
				index = i - 1 ;
				belmntp = alloc_belm(index,131) ;
				belmntp->nps[0] = gp.N + (i-1) * (5*n1L + 3) ; 
				belmntp->nps[1] = gp.N + (i-1) * (5*n1L+3) + 3*n1L + 1 ; 
				belmntp->nps[2] = gp.N + (i-1) * (5*n1L+3) + 4*n1L + 2 ; 
				belmntp->nps[3] = gp.N + i * (5*n1L + 3) ; 
				bottomb(belmntp) ;
			}
			if(i == 1) {
				index = 2 * n2L + 2 * n1L - j ;
				belmntp = alloc_belm(index,131) ;
				belmntp->nps[0] = gp.N + 3*j ; 
				belmntp->nps[1] = gp.N + 3*j - 1; 
				belmntp->nps[2] = gp.N + 3*j - 2 ; 
				belmntp->nps[3] = gp.N + 3*j - 3 ; 
				leftb(belmntp) ;
				if(j == 1)
					belmntp->nextbelp = NULL ;
			}
			if(j == n1L) {
				index = (i-1) * (5*n1L+3) + 3*n1L ;
				x = x0 + 3 * (i-1) * ar ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
				nodep->u[1] = 1.0 ;
				index = (i-1) * (5*n1L+3) + 4*n1L + 1  ;
				x = x0 + 3 * (i-1) * ar + 1.0*ar ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
				nodep->u[1] = 1.0 ;
				index = (i-1) * (5*n1L+3) + 5*n1L + 2  ;
				x = x0 + 3 * (i-1) * ar + 2.0*ar ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
				nodep->u[1] = 1.0 ;
				index = 2 * n2L + n1L - i ;
				belmntp = alloc_belm(index,131) ;
				belmntp->nps[0] = gp.N + i * (5*n1L+3) + 3*n1L ; 
				belmntp->nps[1] = gp.N + i * (5*n1L+3) - 1 ; 
				belmntp->nps[2] = gp.N + (i-1) * (5*n1L+3) + 4*n1L + 1 ; 
				belmntp->nps[3] = gp.N + (i-1) * (5*n1L+3) + 3*n1L ; 
				topb(belmntp) ;
			}
			if(i == n2L) {
				index = i * (5*n1L+3) + 3 * j - 3 ;
				x = x1 ;
				y = yzero + (j-1) * L ;
				nodep = alloc_node(index,x,y) ;
				index++ ;
				x = x1 ;
				y = yzero + (j-1) * L + 0.33333333 * L ;
				nodep = alloc_node(index,x,y) ;
				index++ ;
				x = x1 ;
				y = yzero + (j-1) * L + 0.666666667 * L ;
				nodep = alloc_node(index,x,y) ;
				index =  n2L + j - 1 ;
				belmntp = alloc_belm(index,131) ;
				belmntp->nps[0] = gp.N + n2L * (5*n1L + 3) + 3*j - 3 ; 
				belmntp->nps[1] = gp.N + n2L * (5*n1L + 3) + 3*j - 2 ; 
				belmntp->nps[2] = gp.N + n2L * (5*n1L + 3) + 3*j - 1;
				belmntp->nps[3] = gp.N + n2L * (5*n1L + 3) + 3*j ; 
				rightb(belmntp) ;
			}
			if(i == n2L && j == n1L) {
				index = N.nodes - 1  ;
				x = x1 ;
				y = yone ;
				nodep = alloc_node(index,x,y) ;
			}
			index = (i - 1) * n1L + j - 1;
			elmntp = alloc_elm(index,eltype) ;
			if(index == N.elms-1)
				elmntp->nextelp = NULL ;
			elmntp->nps[0] = gp.N + i * (5*n1L + 3) + 3*j - 3 ; 
			elmntp->nps[1] = gp.N + i * (5*n1L + 3) + 3*j - 2 ; 
			elmntp->nps[2] = gp.N + i * (5*n1L + 3) + 3*j - 1; 
			elmntp->nps[3] = gp.N + i * (5*n1L + 3) + 3*j ; 
			elmntp->nps[4] = gp.N + (i-1) * (5*n1L + 3) + 4*n1L + j + 2 ; 
			elmntp->nps[5] = gp.N + (i-1) * (5*n1L + 3) + 3*n1L + j + 1 ; 
			elmntp->nps[6] = gp.N + (i-1) * (5*n1L + 3) + 3*j ; 
			elmntp->nps[7] = gp.N + (i-1) * (5*n1L + 3) + 3*j - 1; 
			elmntp->nps[8] = gp.N + (i-1) * (5*n1L + 3) + 3*j - 2 ; 
			elmntp->nps[9] = gp.N + (i-1) * (5*n1L + 3) + 3*j - 3 ; 
			elmntp->nps[10] = gp.N + (i-1) * (5*n1L + 3) + 3*n1L + j ; 
			elmntp->nps[11] = gp.N + (i-1) * (5*n1L + 3) + 4*n1L + j + 1 ; 
			break ;
						
	}

	return(0) ;
}

struct node *alloc_node(index,x,y)
int		index ;
double	x, y ;

{
	struct node	*nodep ;
	int	k ;
	
	nodep = gp.N + index ;
	nodep->n = index + 1  ;
	nodep->i = index  ;
	nodep->fxc = 0 ;
	gp.iptrs[index] = nodep ;
	
	nodep->x[0] = x ;
	nodep->x[1] = y ;
	for(k=0;k<N.params;k++)
		nodep->p[k] = defndp ;
	for(k=0;k<N.vars;k++)
		nodep->u[k] = 0.00 ;
	nodep->nextnp = nodep + 1 ;
/*	nodevalues(nodep) ;

	printf("Node %d, i = %d,  n = %d,\n",index,nodep->i,nodep->n);
	printf("Node %d, i = %d,  n = %d,\n",index,(gp.iptrs[index])->i,(gp.iptrs[index])->n);
*/
	return(nodep) ;
}

struct element		*alloc_elm(index, eltype)
int		index, eltype ;

{
	int   k ;
	struct element  *elmntp ;
	
	elmntp = gp.El + index ;
	elmntp->n = index + 1 ;
	elmntp->vtype = eltype;
	elmntp->gtype = eltype;
	elmntp->nnds =  nsf(elmntp->vtype) ;
	if(eltype == 229) {
		elmntp->gtype = 211;
		elmntp->nnds = 16 ;
	}
	if(eltype == 129) {
		elmntp->gtype = 111;
		elmntp->nnds = 4 ;
	}
	if(eltype == 139) {
		elmntp->gtype = 111;
		elmntp->nnds = 6 ;
	}
	for(k=0;k<N.params;k++)
		elmntp->p[k] = defelp ;
	elmntp->matrices = NULL ;
	elmntp->nextelp = elmntp + 1 ;
	elvalues(elmntp) ;

	return(elmntp) ;
}

struct belement		*alloc_belm(index, eltype)
int		index, eltype ;

{
	int		j, k ;
	struct belement  *belmntp ;
	
	belmntp = gp.B + index ;
	belmntp->n = index + 1 ;
	belmntp->vtype = eltype;
	belmntp->gtype = eltype;
	belmntp->nnds =  nsf(belmntp->vtype) ;
	for(k=0;k<N.bparams;k++)
		belmntp->p[k] = defelp ;
	for(j=0;j<N.vars;j++) {
		belmntp->bcs[j] = 0 ;
	}
	belmntp->matrices = NULL ;
	belmntp->nextbelp = belmntp + 1 ;

	return(belmntp) ;
}

int		map_mesh(theblock,maptype,rownum,colnum)
int		maptype, rownum, colnum ;
struct blockmap	*theblock ;

{
	struct node	*np ;
	struct gausspts	g ;
	struct shapefuncs	sf ;
	int		i, j ;
//rownum is the block number in the y-direction
//colnum is the block number in the x-direction
	np = gp.N ;
	g.z = 0.0 ;
	g.w = -1.0 ;
	for(i=0;i<N.nodes;i++) {
		g.x = np->x[0] - 2.0 * colnum ;
		g.y = np->x[1] - 2.0 * rownum ;
		if(get_shape(maptype,&g,&sf) == 0) { 
			np->x[0] = 0.0 ;
			np->x[1] = 0.0 ;
			for(j=0;j<sf.dof;j++) {
				np->x[0] += sf.f[j] * theblock->xnodes[j] ;
				np->x[1] += sf.f[j] * theblock->ynodes[j] ;
			}
		}
//printf("node=%d\tx=%g\ty=%g\n",np->n,np->x[0],np->x[1]);
		np = np->nextnp ;
	}
	return(0) ;
}

