//#pragma options cpluscmt        
#include "Fe_PS.h"
#include "ctype.h"
#include "iriclib.h"
#include "cgnslib.h"
#include <stdio.h>
//#include <direct.h>

#define MAXLINE 200

extern struct control  N ;
extern struct pointers gp;
extern struct transient tvals ;
extern struct settings  sets;
extern double           UW,DIF, epsilon1, epsilon3;

double get_dbl();

int     Nndso, Nelmso, Nbelmso, Nbsego;

int             input(fp)
FILE*             fp ;

{
        int             err,i ;
        FILE    *fptr ;
        
        if((fptr = fp ) != NULL) {
                if( (err = get_control(fptr)) == 0) {
                        put_control(stdout) ;
			
			for(i=0;i<N.nodes;i++) 
                                get_node(fptr,i) ;


                        for(i=0;i<N.elms;i++)
                                get_elm(fptr,i) ;

                        for(i=0;i<N.belms;i++)
                                get_belm(fptr,i) ;

                        find_neighbour();

                        for(i=0;i<N.boundarysegs;i++)
                                get_bseg(fptr,i) ;
                        if(fptr != stdin) fclose(fptr) ;
                                set_bvalues() ;
                }

                else
		  printf(" Control Variable Error %d\n",err) ;
		//if(fptr != stdin) fclose(fptr) ; modified by C.Frias June 15 2010
        }
        else {
                printf(" Unable to open that file\n") ;
                err = -1 ;
                }

        return(err);
}

int             get_control(f)
FILE    *f ;

{
        int     i, j, k, in, ii ;
        long            bigint ;
        
        j = 0 ;
        for(i=0;i<1;i++) {
                j++ ;
                N.trans = get_int(f) ;
                N.meshtype = get_int(f);
                if(N.trans < 0 || N.trans > 1) break ;
                j++ ;
                if(N.trans == 1)
                        get_trans(f) ;
                N.dims = get_int(f) ;
                if(N.dims < 1 || N.dims > NDIMS) break ;
                j++ ;
                N.vars = get_int(f) ;
                if(N.vars < 1 || N.vars > NVARS) break ;
                j++ ;
                for(ii=0;ii<((N.vars+1)*N.vars);ii++)
                        N.Keqns[ii] = get_int(f) ;
                j++ ;
                N.params = get_int(f) ;
                if(N.params < 0 || N.params > NPARAMS) break ;
                j++ ;
                N.bparams = get_int(f) ;
                if(N.bparams < 0 || N.bparams > NBPARAMS) break ;
                j++ ;
                if(gp.N != NULL) {
                        clear_data() ;
                        free(gp.N) ;
                        free(gp.iptrs) ;
                }
                N.nodes = get_int(f) ;
                Nndso = N.nodes ;
                bigint = (5 * N.nodes)/4 * sizeof(struct node *) ;
                if(N.nodes < 2 || (gp.iptrs = (struct node **) malloc(bigint)) == NULL) break ;
                bigint = (long int) N.nodes * (long int) sizeof(struct node) ;
                if(N.nodes < 2 || (gp.N = (struct node *) malloc(bigint)) == NULL) break ;
                j++ ;
                N.elms = get_int(f) ;
                Nelmso = N.elms ;
                bigint = (long int) N.elms * (long int) sizeof(struct element) ;
                if(gp.El != NULL)
                        free(gp.El) ;
                if(N.elms < 1 || (gp.El = (struct element *) malloc(bigint)) == NULL) break ;
                j++ ;
                N.belms = get_int(f) ;
                Nbelmso = N.belms ;
                bigint = (long int) N.belms * (long int) sizeof(struct belement) ;
                if(gp.B != NULL)
                        free(gp.B) ;
                if(N.belms < 1 || (gp.B = (struct belement *) malloc( bigint)) == NULL) break ;
                
                j++ ;
                N.boundarysegs = get_int(f) ;
                Nbsego = N.boundarysegs ;
                if(N.boundarysegs>0) {
                        bigint = (long int) N.boundarysegs * (long int) sizeof(struct boundaryseg) ;
                        if(gp.BSEG != NULL)
                                free(gp.BSEG) ;
                        if(N.boundarysegs < 1 || (gp.BSEG = (struct boundaryseg *) malloc( bigint)) == NULL) break ;
                }
                j = 0 ;
        }
        return(j) ;
}

int             get_trans(f)
FILE *f ;
{
        double get_dbl() ;

        tvals.nsteps = get_int(f) ;
        tvals.dtfac = get_dbl(f) ;
        tvals.t = get_dbl(f) ;
        tvals.dt = get_dbl(f) ;
        tvals.theta = get_dbl(f) ;
        UW = get_dbl(f) ;
        DIF     = get_dbl(f);
        sets.latitude=get_dbl(f);
        tvals.S=get_dbl(f);
        sets.diffusivewave=get_int(f);
        sets.uwJ=get_dbl(f);
        sets.plotcode=get_int(f);
        sets.transbcs=get_int(f);
        sets.maxitnum=get_int(f);
        sets.smallH=get_int(f);
        sets.JE=get_int(f);
//      tvals.minH=get_dbl(f);
		epsilon1 = get_dbl(f);
        tvals.gwH=get_dbl(f);
//      tvals.GWD=get_dbl(f);
		epsilon3 = get_dbl(f);
        tvals.T=get_dbl(f);
		tvals.sg_ice = 0.92;
		tvals.itnum = 0;
        return(0) ;
}       


int             clear_data() 

{
        int             i ;
        struct node     *np1, *np2 ;
        struct element  *elp1, *elp2 ;
        struct belement *belp1, *belp2 ;
        
        np1 = gp.N ;
        for(i=0;i<N.nodes;i++) {
                np2 = np1->nextnp ;
                if(i >= Nndso)
                        free(np1) ;
                np1 = np2 ;
        }
        belp1 = gp.B ;
        for(i=0;i<N.belms;i++) {
                belp2 = belp1->nextbelp ;
                if(i >= Nbelmso)
                        free(belp1) ;
                belp1 = belp2 ;
        }
        elp1 = gp.El ;
        for(i=0;i<N.elms;i++) {
                elp2 = elp1->nextelp ;
                if(elp1->matrices != NULL)
                        free(elp1->matrices) ;
                if(i >= Nelmso)
                        free(elp1) ;
                elp1 = elp2 ;
        }
        return(0) ;
}
        
int             get_node(f,i)
FILE    *f ;
int             i ;

{
        int      j  ;
        double  get_dbl() ;
        struct node  *nodep ;
        
        nodep = gp.N + i ;
        nodep->n = get_int(f)  ;
        nodep->i = i ;
        nodep->fxc = 0 ;
        gp.iptrs[i] = nodep ;
        for(j=0;j<N.dims;j++) 
                nodep->x[j] = get_dbl(f) ;
        for(j=N.dims;j<NDIMS;j++) 
                nodep->x[j] = 0.0 ;
        for(j=0;j<N.params;j++)
                nodep->p[j] = get_dbl(f) ;
        for(j=0;j<N.vars;j++)
                nodep->u[j] = get_dbl(f) ;
		for(j=0;j<N.params;j++)
				nodep->ice[j] = 0.0;
        nodep->nextnp = nodep + 1 ;
        return(0) ;
}

int             get_elm(f,i)
FILE    *f ;

{
        int  j, jj, nnds, n, rc ;
        struct element  *elmntp ;
        struct node             *np ;
        
        elmntp = gp.El + i ;
        rc = 0 ;
        elmntp->n = get_int(f) ;
        elmntp->vtype = get_int(f);
        elmntp->gtype = get_int(f);
        elmntp->nnds =  nsf(elmntp->vtype) ;
        nnds = nsf(elmntp->gtype);
        if(nnds > elmntp->nnds) elmntp->nnds = nnds ;
        if(elmntp->vtype == 229)
                elmntp->nnds = 16 ;
        if(elmntp->vtype == 129)
                elmntp->nnds = 4 ;
        if(elmntp->vtype == 139)
                elmntp->nnds = 6 ;
        for(j=0;j<elmntp->nnds;j++) {
                if((n = get_int(f)) >= 0) { 
 /*                      np = gp.N ;
                        for(jj=0;jj<N.nodes;jj++) {
                                if(np->n == n)
                                        break ;
                                np = np->nextnp ;
                        } 
                        if(np == NULL) {
                                printf(" Non-existent node %d referred to in element %d\n",n,elmntp->n) ;
                                rc = -1 ;
                        }
                        elmntp->nps[j] = np ;
*/
						elmntp->nps[j] = gp.N + (n-1);
                }
                else
                        elmntp->nps[j] = NULL ;
        }
        for(j=0;j<N.params;j++)
                elmntp->p[j] = get_dbl(f) ;
        elmntp->matrices = NULL ;
        if(i < N.elms - 1)
                elmntp->nextelp = elmntp + 1 ;
        else
                elmntp->nextelp = NULL ;
        return(rc) ;
}

int             get_belm(f,i)
FILE    *f ;

{
        int  j, jj, nnds, n, rc ;
        struct belement  *belmntp ;
        struct node             *np;
        
        belmntp = gp.B + i ;
        rc = 0 ;
        belmntp->n = get_int(f) ;
        belmntp->vtype = get_int(f);
        belmntp->gtype = get_int(f);
        belmntp->nnds =  nsf(belmntp->vtype) ;
        nnds = nsf(belmntp->gtype);
        if(nnds > belmntp->nnds) belmntp->nnds = nnds ;
        for(j=0;j<belmntp->nnds;j++) {
                n = get_int(f) ;
 /*               np = gp.N ;
                for(jj=0;jj<N.nodes;jj++) {
                        if(np->n == n)
                                break ;
                        np = np->nextnp ;
                } 
                if(np == NULL) {
                        printf(" Non-existent node %d referred to in element %d\n",n,belmntp->n) ;
                        rc = -1 ;
                }
                belmntp->nps[j] = np ;
*/
 				belmntp->nps[j] = gp.N + (n-1);
       }
        /*
        np->fxc = 1;
         aaa */

        for(j=0;j<N.bparams;j++)
                belmntp->p[j] = get_dbl(f) ;
        for(j=0;j<N.vars;j++) {
                belmntp->bcs[j] = get_int(f) ;
        }
        belmntp->matrices = NULL ;
        belmntp->elp=NULL;
        belmntp->bseg=NULL;
        if(i < N.belms - 1)
                belmntp->nextbelp = belmntp + 1;
        else
                belmntp->nextbelp = gp.B ;
        return(rc) ;
}

int             get_bseg(f,i)
FILE    *f ;

{
        int  j, jj, nnds, n;
        struct boundaryseg  *bseg ;
        struct node             *np;
        double  get_dbl() ;
		char c;

        bseg = gp.BSEG + i ;
        bseg->n = get_int(f) ;
        
        bseg->bcs[0] = get_int(f);
        bseg->p[0] = get_dbl(f);
        bseg->p[1] = get_dbl(f);
        bseg->nstart = get_int(f);
        bseg->nend = get_int(f);

		c = getc(f);
		if(isdigit(c))
			ungetc(c,f);
		else
		{
		while (c != '\n')
			c = getc(f);
		}
		    
        if(i < N.boundarysegs - 1)
                bseg->nextbseg = bseg + 1;
        else
                bseg->nextbseg = NULL ;
        return(0) ;
}


double  get_dbl(fptr)
FILE    *fptr ;

{
        char    s[25] ;
        int             n ;
        
        n  = getnumstr(fptr,s) ;
        return( atof(s) ) ;
}

int             get_int(fptr)
FILE    *fptr ;

{
        char    s[25] ;
        int             n ;
        double   d;
        
        n = getnumstr(fptr,s) ;
        d = atof(s) ;
        return( (int)  d ) ;
}

int     getnumstr(fptr,s)
char    *s ;
FILE    *fptr ;

{
        register int  c ;
        register char *cs ;
        
        cs = s ;
        while( (c = getc(fptr)) != EOF) {
                if((isdigit(c) != 0) || (c == '.') || (c == '-') || (c == '+')) {
                        *cs++ = c ;
                        while( (c = getc(fptr)) != EOF) {
                                if((isdigit(c) != 0) || (c == '.') || (c == 'e') || (c == 'E') || (c == '-') || (c == '+') )
                                        *cs++ = c ;
                                else
                                        break ;
                        }
                        *cs = '\0' ;
                        break ;
                }
        }
        return(cs - s) ;
}

int             set_bvalues()

{
        int             i, j ;
        struct node             *nptr ;
        struct  belement        *belmntp;
        
        j = 0 ;
        nptr = gp.N ;
        for(i=0;i<N.nodes;i++) {
                nptr->i = i ;
                gp.iptrs[i] = nptr ;
                nptr = nptr->nextnp ;
        }
        
        belmntp = gp.B;
        for(i=0;i<N.belms;i++) {
                for(j=0;j<belmntp->nnds;j++)
        /*              if(belmntp->bcs[0]==0)*/
                                belmntp->nps[j]->fxc = 1;
                                
                belmntp = belmntp->nextbelp;
        }
        if(N.boundarysegs>0)
                update_BCs();
        N.ukns = N.nodes * N.vars ;
        printf(" Number of unknowns = %d\n",N.ukns) ;
        return(j) ;
}

int             output(nin) //Original output routine
int             nin ;

{
        int             err,i ;
        FILE    *fptr, *put_fptr() ;
        struct node     *np ;
        struct element  *elp ;
        struct belement *belp ;
        struct  boundaryseg     *bseg;

        if((fptr = put_fptr(nin) ) != NULL) {
                if( (err = put_control(fptr)) == 0) {
                        fprintf(fptr,"\n Node Information \n");
                        fprintf(fptr,"\n Node #, Coordinates, Parameters, Variables\n\n");
                        np = gp.N ;
                        for(i=0;i<N.nodes;i++) {
//if(np->u[0]>0.05)
                                put_node(fptr,i,np) ;
                                np = np->nextnp ;
                        }
                        fprintf(fptr,"\n Element Information \n");
                        fprintf(fptr,"\n Element #, vtype, gtype, nodes\n\n");
                        elp = gp.El ;
                        for(i=0;i<N.elms;i++) {
                                put_elm(fptr,i,elp) ;
                                elp = elp->nextelp ;
                        }
			fprintf(fptr,"\n Boundary Element #, vtype, gtype, nodes, boundary condition codes\n\n");
                        belp = gp.B ;
                        for(i=0;i<N.belms;i++) {
                                put_belm(fptr,i,belp) ;
                                belp = belp->nextbelp ;
                        }
                        fprintf(fptr,"\n Boundary Seg #,Boundary type,stage,QT,start node #,end node #\n\n");
                        bseg = gp.BSEG ;
                        for(i=0;i<N.boundarysegs;i++) {
                                put_bseg(fptr,i,bseg) ;
                                bseg = bseg->nextbseg ;
                        }
                        if(fptr != stdin) fclose(fptr) ;
                }
                else {
                        printf(" Control Variable Error %d\n",err) ;
                        fclose(fptr) ;
                }
        }
        else
                printf(" Unable to open that file\n") ;
        return(err);
}




int             perimeter(nin)
int             nin ;

{
        int             err,i ;
        FILE    *fptr, *put_fptr() ;
        struct belement *belp ;

        if((fptr = put_fptr(nin) ) != NULL) {
                fprintf(fptr,"\n Boundary Node Information \n");
                fprintf(fptr,"\n Node #, Coordinates, Parameters, Variables\n\n");
                belp = gp.B ;
                for(i=0;i<N.belms;i++) {
                        put_node(fptr,i,belp->nps[0]) ;
                        belp = belp->nextbelp ;
                }
                if(fptr != stdin) fclose(fptr) ;
        }
        else
                printf(" Unable to open that file\n") ;
        return(err);
}

int             put_control(f) //ORIGINAL put_control routine
FILE    *f ;

{
        int     i, ii ;
        
        fprintf(f," Transient analysis = %d\n",N.trans) ;
        fprintf(f," Mesh type = %d\n",N.meshtype) ;
        if(N.trans == 1){
                fprintf(f," Number of Time Steps = %d\n",tvals.nsteps) ;
                fprintf(f," Delta t Acceleration Factor = %5.3f\n",tvals.dtfac) ;
                fprintf(f," Time = %5.3f\n",tvals.t) ;
                fprintf(f," Delta t = %5.3f\n", tvals.dt) ;
                fprintf(f," Theta = %5.3f\n", tvals.theta) ;
                fprintf(f," UW = %5.3f\n", UW) ;
                fprintf(f," DIF = %5.3f\n", DIF) ;
                fprintf(f," Latitude \t\t = %5.3f\t degrees\n", sets.latitude) ;
                fprintf(f," Old equal New option\t = %d\t zero for old not equal new\n", sets.oldeqnew) ;
                fprintf(f," Diffusive wave Solution = %d\t zero for fully dynamic only\n", sets.diffusivewave) ;
                fprintf(f," UW Jacobian terms included = %5.3f\t zero for not included\n", sets.uwJ) ;
                fprintf(f," Plot Code \t\t\t\t\t= %d\t zero for xsec one for contour two for velocity and three for threeD\n", sets.plotcode) ;
                fprintf(f," Transient Boundary Condition = %d\t zero for Steady BCs\n", sets.transbcs) ;
                fprintf(f," Maximum Number of Iterations = %d\n", sets.maxitnum) ;
                fprintf(f," Small Depthes Occur \t= %d\t zero for no small depth calculations\n", sets.smallH) ;
                fprintf(f," Jacobian Terms included = %d\t zero for not included\n", sets.JE) ;
                fprintf(f," Minimum Depth for Diffusive Wave Calculation= %5.3f\n", tvals.minH) ;
                fprintf(f," Minimum Depth for Groundwater Flow Calculation= %5.3f\n", tvals.gwH) ;
                fprintf(f," Groundwater Flow Artificial Diffusion = %5.3f\n", tvals.GWD) ;
                fprintf(f," Transmissivity of Aquifer\t\t= %5.3f\n", tvals.T);


        }
        fprintf(f," Dimensions = %d\n",N.dims) ;
        fprintf(f," Number of Variables = %d\n",N.vars) ;
        fprintf(f," [K] governing equation numbers \n") ;
                for(i=0;i<N.vars;i++) {
                        for(ii=0;ii<N.vars;ii++) {
                                fprintf(f,"\t %d ",N.Keqns[i*(N.vars+1) + ii])  ;
                        }
                        fprintf(f,"\t\t %d \n",N.Keqns[i*(N.vars+1) + N.vars]) ;
                }
        fprintf(f," Number of Parameters = %d\n",N.params) ;
        fprintf(f," Number of Boundary Parameters = %d\n",N.bparams) ;
        fprintf(f," Number of Nodes = %d\n",N.nodes) ;
        fprintf(f," Number of Elements = %d\n",N.elms) ;
        fprintf(f," Number of  Boundary Elements = %d\n",N.belms) ;
        fprintf(f," Number of  Boundary Segments = %d\n",N.boundarysegs) ;
        fprintf(f, " \n") ;
        
        return(0) ;
}




int             put_node(f,i,nodep) //ORIGINAL put_node routine
FILE    *f ;
int             i ;
struct node  *nodep ;

{
        int     j ;
        
        put_int(f,nodep->n)  ;
//      put_int(f,nodep->fxc)  ;

        for(j=0;j<N.dims;j++) 
                put_dbl(f,(double) nodep->x[j]) ;
        for(j=0;j<N.params;j++)
                put_dbl(f,(double) nodep->p[j]) ;
        for(j=0;j<N.vars;j++)
                put_dbl(f,(double) nodep->u[j]) ;
        fprintf(f," \n");
        return(0) ;
}

int             put_elm(f,i,elmntp) //ORIGINAL pu_elm routine
int             i ;
FILE    *f ;
struct element *elmntp ;

{
        int             j;
        
        put_int(f,elmntp->n) ;
        put_int(f,elmntp->vtype);
        put_int(f,elmntp->gtype);
        for(j=0;j<elmntp->nnds;j++) {
                if(elmntp->nps[j] != NULL)
                        put_int(f,(elmntp->nps[j])->n) ;
                else
                        put_int(f,-1) ;
        }
        for(j=0;j<N.params;j++)
                put_dbl(f,(double) elmntp->p[j]) ;
        fprintf(f," \n");
        return(0) ;
}

int             put_belm(f,i,belmntp)
int             i ;
FILE    *f ;
struct belement *belmntp ;

{
        int             j;

        put_int(f,belmntp->n) ;
        put_int(f,belmntp->vtype);
        put_int(f,belmntp->gtype);
        for(j=0;j<belmntp->nnds;j++)
                put_int(f,(belmntp->nps[j])->n) ;
        for(j=0;j<N.bparams;j++)
                put_dbl(f,(double) belmntp->p[j]) ;
        for(j=0;j<N.vars;j++)
                put_int(f,belmntp->bcs[j]) ;
        fprintf(f," \n");
        return(0) ;
}

int             put_bseg(f,i,bseg)
int             i ;
FILE    *f ;
struct boundaryseg *bseg ;

{
        int             j;

        put_int(f,bseg->n) ;
        put_int(f,bseg->bcs[0]);
        for(j=0;j<2;j++)
                put_dbl(f,(double) bseg->p[j]) ;
        put_int(f,bseg->nstart) ;
        put_int(f,bseg->nend) ;
        fprintf(f," \n");
        return(0) ;
}

int     put_dbl(fptr,d)
FILE    *fptr ;
double  d;

{
        fprintf(fptr,"\t%f",d) ;
        return(0) ;
}

int             put_int(fptr,i)
FILE    *fptr ;
int             i ;

{
        fprintf(fptr,"\t%d",i) ;
        return(0) ;
}

int     get_IFG4dat(nin)
int     nin;

{
        FILE    *inf,*outf,*get_fptr(),*put_fptr(),*meshf;
        char    line[MAXLINE];
        int     c,i,n,ioc[20],xsecnum,xs,iartus,nwus,ndus;
        double  f2m,cf2cm,location,weight,sozf,slope,breakloc,dnum;
        double  get_limiteddouble(),ltotal,ymax,sumb,sumH;
        struct  xsec    *xsp,*newxsp,*alloc_xsec();
        struct  datapoint       *alloc_dp(),*dp;
        double  dwsl;
        int             iartds,nwds,ndds;
        
        N.xsecs=0;
        N.datapts=0;
        ltotal=0.0;
        ymax=-1000.0;

        printf(" should an upstream artificial inflow sec be generated? Enter 1 for yes\n");
        scanf("%d",&iartus);
        if(iartus==1){
                printf("At how many depthes upstream should the artificial sec be located?\n");
                scanf("%d",&ndus);
                printf("At how many widths upstream should the artificial sec be located?\n");
                scanf("%d",&nwus);
        }
        printf(" should a downstream artificial inflow sec be generated? Enter 1 for yes\n");
        scanf("%d",&iartds);
        if(iartds==1){
                printf("At how many depthes upstream should the artificial sec be located?\n");
                scanf("%d",&ndds);
                printf("At how many widths upstream should the artificial sec be located?\n");
                scanf("%d",&nwds);
        }
        printf("Enter Minimum allowable depth in meters\n");
        scanf("%lf",&dwsl) ;
        if((inf=get_fptr(1))==NULL)
                return(-1);
        
        if((outf=put_fptr(1))==NULL)
                return(-1);
// read and print the two title lines

        for(i=0;i<2;i++) {
                n=user_getline(line,inf);
//              putline(line,outf,n);
        }
// read the IOC line (22 IOC values)
        if((c=getc(inf))=='I')
                getlinename(9,inf);
        n=1;
        for(i=1;i<23;i++)
                if((c=getc(inf))!='\n') {
                        ioc[i]=c;
                        n++;
                }
                else
                        break;
//      for(i=1;i<n;i++)
//              putc(ioc[i],outf);
        if(c!='\n') {
                user_getline(line,inf);
//              putc('\n',outf);
        }
//skip the BMAX,NMAX and NSLP lines
        if((ioc[14]==3)||(ioc[14]==5))
                if((c=getc(inf))=='B')
                        user_getline(line,inf);
        if((ioc[15]==1)||(ioc[15]==2))
                if((c=getc(inf))=='N')
                        user_getline(line,inf);
        if(ioc[16]==1)
                if((c=getc(inf))=='N')
                        user_getline(line,inf);

//              check for conversion    after going to end of the above line
        
        f2m=1.0;
        cf2cm=1.0;
        if(c=getc(inf)=='M') 
                user_getline(line,inf);
        else {
                f2m=0.3048;
                cf2cm=0.028317;
                user_getline(line,inf);
        }
//now read the QARD lines (note that the first character might have been read)
//  for now skip all the QARD lins

//      putline(line,outf,n);
        while((c=getc(inf))=='Q')
                user_getline(line,inf);
                
//now skip the DEPTH lines
        if(c=='D'){
                user_getline(line,inf);
                while((c=getc(inf))=='D')
                        user_getline(line,inf);
        }
//now read the XSEC line
        xsp=alloc_xsec();
        gp.X=xsp;
        
        while(c!='E') {
                if(c=='X'){
                        getlinename(3,inf);
                        xsecnum = get_xsecnum(6,inf);
                }
                location        = get_limiteddouble(10,inf)*f2m;        
                weight=get_limiteddouble(4,inf);
                sozf=get_limiteddouble(10,inf)*f2m;
                user_getline(line,inf);
        

                xsp->n=xsecnum;
                xsp->loc=location;
                xsp->weight=weight;
                xsp->sozf=sozf;
                ltotal+=location;
printf("xsecnum=%d\n",xsp->n);

printf("location=%g\n",ltotal); 
printf("sozf=%g\n",sozf);
printf("weight=%g\n",xsp->weight);      
        
        
//now read the xsec geometry

                while((c=getc(inf))!='N') {
                        if((xs=get_xsecnum(9,inf))==xsecnum) {
                                while((dnum=get_limiteddouble(nin,inf))!=-1) {
                                        if(dnum==-2) {
                                                user_getline(line,inf);
                                                break;
                                        }
                                        xsp->dp[xsp->ndatapts]=alloc_dp();
                                        xsp->dp[xsp->ndatapts]->b=get_limiteddouble(nin,inf)*f2m;
                                        xsp->dp[xsp->ndatapts]->y=dnum*f2m;
                                        xsp->ndatapts++;
                                        }
                        }
                }
//skip all the NS lines for now
                if(c=='N') {
                        user_getline(line,inf);
                        while((c=getc(inf))=='N')
                                user_getline(line,inf);
                }
//skip the WSL line for now
                if(c=='W') {
                        user_getline(line,inf);
                        while((c=getc(inf))=='W')
                                user_getline(line,inf);
                }
                if(c=='C') {
                        getlinename(4,inf);
                        xsecnum = get_xsecnum(5,inf);
                        xsp->wsl=get_limiteddouble(10,inf)*f2m;
                        xsp->Q=get_limiteddouble(10,inf)*cf2cm;
//fprintf(outf,"wsl=%g\tQ=%g\n",xsp->wsl,xsp->Q);
                }
//now skip all the other CAL lines and the VEL lines
                while((c=getc(inf))!='X') {
                        if(c=='E')
                                break;
                        user_getline(line,inf);
                }
//now calculate ymax    
                if(xsp->dp[xsp->ndatapts-1]->y>ymax)
                        ymax=xsp->dp[xsp->ndatapts-1]->y;

printf("number of data pts=%d\n",xsp->ndatapts);
//for(i=0;i<xsp->ndatapts;i++) {
//fprintf(outf,"%d\t%g\t%g\n",i,xsp->dp[i]->y,xsp->dp[i]->b);

//}

                if(c=='X') {
                        newxsp=alloc_xsec();
                        xsp->nextxsp=newxsp;
                        xsp=newxsp;
                }

        }
        
        if(iartus==1) {
                
                        
                newxsp=alloc_xsec();
                xsp->nextxsp=newxsp;
                newxsp->ndatapts=xsp->ndatapts;
                
//now calculate parameters of last xsec
                for(i=0;i<xsp->ndatapts;i++) {
                        dp=xsp->dp[i];
                        dp->h=xsp->wsl-dp->b;
                        if(dp->h>dwsl) {
                                if(xsp->nwet==0)
                                        xsp->ifirst=i;
                                xsp->nwet++;
                                xsp->ilast=i;
                        }
                }
//              xsp->ilast=xsp->ifirst+xsp->nwet-1;
                xsp->nwet=xsp->ilast-xsp->ifirst+1;
                
                sumb=0.0;
                sumH=0.0;
                for(i=xsp->ifirst;i<(xsp->ilast+1);i++) {
                        sumb+=xsp->dp[i]->b;
                        sumH+=(xsp->wsl-xsp->dp[i]->b);
                }
                sumb/=xsp->nwet;
                sumH/=xsp->nwet;
                
                
                
                
                
                xsp->nwet=0;


                for(i=0;i<xsp->ndatapts;i++) {
                        newxsp->dp[i]=alloc_dp();
                        newxsp->dp[i]->b=sumb;
                        newxsp->dp[i]->y=xsp->dp[i]->y;
                }

//set the bed level of the data points before ifirst and after ilast
// to a very high value in order to exclude these points
                for(i=0;i<xsp->ifirst;i++) 
                        newxsp->dp[i]->b=10000.0;
                for(i=(xsp->ilast+1);i<xsp->ndatapts;i++)
                        newxsp->dp[i]->b=10000.0;
printf("sumH=%g\tsumb=%g\n",sumH,sumb);         
                        
                newxsp->wsl=xsp->wsl;
                newxsp->Q=xsp->Q;
                newxsp->wsl=xsp->wsl;
                newxsp->n=xsp->n+1;
                
                if((ndus*sumH)>(nwus*(xsp->dp[xsp->ndatapts-1]->y-xsp->dp[0]->y)))
                        newxsp->loc=(ndus*sumH);
                else
                        newxsp->loc=(nwus*(xsp->dp[xsp->ndatapts-1]->y-xsp->dp[0]->y));
                ltotal+=newxsp->loc;
        }
        

        if(iartds==1) {
                
                        
                newxsp=alloc_xsec();
                xsp=gp.X;
                gp.X=newxsp;
                newxsp->nextxsp=xsp;
                newxsp->ndatapts=xsp->ndatapts;
                
//now calculate parameters of last xsec
                for(i=0;i<xsp->ndatapts;i++) {
                        dp=xsp->dp[i];
                        dp->h=xsp->wsl-dp->b;
                        if(dp->h>dwsl) {
                                if(xsp->nwet==0)
                                        xsp->ifirst=i;
                                xsp->nwet++;
                        }
                }
                xsp->ilast=xsp->ifirst+xsp->nwet-1;
                
                sumb=0.0;
                sumH=0.0;
                for(i=xsp->ifirst;i<(xsp->ilast+1);i++) {
                        sumb+=xsp->dp[i]->b;
                        sumH+=(xsp->wsl-xsp->dp[i]->b);
                }
                sumb/=xsp->nwet;
                sumH/=xsp->nwet;
                
                
                
                
                
                xsp->nwet=0;


                for(i=0;i<xsp->ndatapts;i++) {
                        newxsp->dp[i]=alloc_dp();
                        newxsp->dp[i]->b=sumb;
                        newxsp->dp[i]->y=xsp->dp[i]->y;
                }

//set the bed level of the data points before ifirst and after ilast
// to a very high value in order to exclude these points
                for(i=0;i<xsp->ifirst;i++) 
                        newxsp->dp[i]->b=10000.0;
                for(i=(xsp->ilast+1);i<xsp->ndatapts;i++)
                        newxsp->dp[i]->b=10000.0;
printf("sumH=%g\tsumb=%g\n",sumH,sumb);         
                        
                newxsp->wsl=xsp->wsl;
                newxsp->Q=xsp->Q;
//              newxsp->n=xsp->n+1;
                
                
                
                if((ndds*sumH)>(nwds*(xsp->dp[xsp->ndatapts-1]->y-xsp->dp[0]->y)))
                        xsp->loc=(ndds*sumH);
                else
                        xsp->loc=(nwds*(xsp->dp[xsp->ndatapts-1]->y-xsp->dp[0]->y));
                ltotal+=newxsp->loc;
                
//      now renumber all the xsecs
                xsp=gp.X;
                xsp->n=1;
                xsp=xsp->nextxsp;
                for(i=0;i<(N.xsecs-1);i++){
                        xsp->n++;
                        xsp=xsp->nextxsp;
                }
                xsp=gp.X;
printf("Here are the new xsec numbers total of %d\n",N.xsecs);
for(i=0;i<N.xsecs;i++){
printf("i=%d\txsp->n=%d\n",i,xsp->n);
xsp=xsp->nextxsp;
}
        }
        
        
printf("ltotal=%g\tymax=%g\n",ltotal,ymax);                     
        complete_data(ltotal,ymax,dwsl);
                                
        if((meshf=put_fptr(1))==NULL)
                return(-1);
        put_interpolationdata(outf);
        put_meshdata(meshf);
        fclose(outf);
        fclose(meshf);
        return(0);
}
int     get_IFG4xsecs(nin)
int     nin;

{
        FILE    *inf,*outf,*get_fptr(),*put_fptr(),*meshf;
        char    line[MAXLINE];
        int     c,i,n,ioc[20],xsecnum,xs,iartus,nwus,ndus,j;
        double  f2m,cf2cm,location,weight,sozf,slope,breakloc,dnum;
        double  get_limiteddouble(),ltotal,ymax,sumb,sumH;
        struct  xsec    *xsp,*newxsp,*alloc_xsec(),*temp;
        struct  datapoint       *alloc_dp(),*dp;
        double  xcoord;
        
        N.xsecs=0;
        N.datapts=0;
        ltotal=0.0;
        ymax=-1000.0;

//      printf("Enter Minimum allowable depth in meters\n");
        if((inf=get_fptr(1))==NULL)
                return(-1);
        
        if((outf=put_fptr(1))==NULL)
                return(-1);
// read and print the two title lines

        for(i=0;i<2;i++) {
                n=user_getline(line,inf);
//              putline(line,outf,n);
        }
// read the IOC line (22 IOC values)
        if((c=getc(inf))=='I')
                getlinename(9,inf);
        n=1;
        for(i=1;i<23;i++)
                if((c=getc(inf))!='\n') {
                        ioc[i]=c;
                        n++;
                }
                else
                        break;
//      for(i=1;i<n;i++)
//              putc(ioc[i],outf);
        if(c!='\n') {
                user_getline(line,inf);
//              putc('\n',outf);
        }
//skip the BMAX,NMAX and NSLP lines
        if((ioc[14]==3)||(ioc[14]==5))
                if((c=getc(inf))=='B')
                        user_getline(line,inf);
        if((ioc[15]==1)||(ioc[15]==2))
                if((c=getc(inf))=='N')
                        user_getline(line,inf);
        if(ioc[16]==1)
                if((c=getc(inf))=='N')
                        user_getline(line,inf);

//              check for conversion    after going to end of the above line
        
        f2m=1.0;
        cf2cm=1.0;
        if(c=getc(inf)=='M') 
                user_getline(line,inf);
        else {
                f2m=0.3048;
                cf2cm=0.028317;
                user_getline(line,inf);
        }
//now read the QARD lines (note that the first character might have been read)
//  for now skip all the QARD lins

//      putline(line,outf,n);
        while((c=getc(inf))=='Q')
                user_getline(line,inf);
                
//now skip the DEPTH lines
        if(c=='D'){
                user_getline(line,inf);
                while((c=getc(inf))=='D')
                        user_getline(line,inf);
        }
//now read the XSEC line
        xsp=alloc_xsec();
        gp.X=xsp;
        
        while(c!='E') {
                if(c=='X'){
                        getlinename(3,inf);
                        xsecnum = get_xsecnum(6,inf);
                }
                location        = get_limiteddouble(10,inf)*f2m;        
                weight=get_limiteddouble(4,inf);
                sozf=get_limiteddouble(10,inf)*f2m;
                user_getline(line,inf);
        

                xsp->n=xsecnum;
                xsp->loc=location;
                xsp->weight=weight;
                xsp->sozf=sozf;
                ltotal+=location;
printf("xsecnum=%d\n",xsp->n);

printf("location=%g\n",ltotal); 
printf("sozf=%g\n",sozf);
printf("weight=%g\n",xsp->weight);      
        
        
//now read the xsec geometry

                while((c=getc(inf))!='N') {
                        if((xs=get_xsecnum(9,inf))==xsecnum) {
                                while((dnum=get_limiteddouble(nin,inf))!=-1) {
                                        if(dnum==-2) {
                                                user_getline(line,inf);
                                                break;
                                        }
                                        xsp->dp[xsp->ndatapts]=alloc_dp();
                                        xsp->dp[xsp->ndatapts]->b=get_limiteddouble(nin,inf)*f2m;
                                        xsp->dp[xsp->ndatapts]->y=dnum*f2m;
                                        xsp->ndatapts++;
                                        }
                        }
                }
//skip all the NS lines for now
                if(c=='N') {
                        user_getline(line,inf);
                        while((c=getc(inf))=='N')
                                user_getline(line,inf);
                }
//skip the WSL line for now
                if(c=='W') {
                        user_getline(line,inf);
                        while((c=getc(inf))=='W')
                                user_getline(line,inf);
                }
                if(c=='C') {
                        getlinename(4,inf);
                        xsecnum = get_xsecnum(5,inf);
                        xsp->wsl=get_limiteddouble(10,inf)*f2m;
                        xsp->Q=get_limiteddouble(10,inf)*cf2cm;
//fprintf(outf,"wsl=%g\tQ=%g\n",xsp->wsl,xsp->Q);
                }
//now skip all the other CAL lines and the VEL lines
        
                j=0;
                while((c=getc(inf))!='X') {
                        if(c=='E')
                                break;
                        if(c=='V') {
                                getlinename(3,inf);
                                xsecnum=get_xsecnum(6,inf);
                                for(i=0;i<12;i++) {
                                        xsp->dp[j]->u=get_limiteddouble(nin,inf)*f2m;
                                        printf(" u =%g\t i=%d\tj=%d\tsec #=%d\n",xsp->dp[j]->u,i,j,xsp->n);
                                        j++;
                                        if(j>=(xsp->ndatapts-1))
                                                break;
                                }
                        }
                        
                        
                        
                        
                        
                        
                        
                        else
                                user_getline(line,inf);
                }
//now calculate ymax    

printf("number of data pts=%d\n",xsp->ndatapts);
//for(i=0;i<xsp->ndatapts;i++) {
//fprintf(outf,"%d\t%g\t%g\n",i,xsp->dp[i]->y,xsp->dp[i]->b);

//}

                if(c=='X') {
                        newxsp=alloc_xsec();
                        xsp->nextxsp=newxsp;
                        xsp=newxsp;
                }

        }
//Now calculate the proper x-coords of the cross sec (starting at us)

printf("out of while loop\n");
        xsp=gp.X;
        temp=gp.X;
        for(i=0;i<(N.xsecs-1);i++) 
                temp=temp->nextxsp;
        temp->nextxsp=gp.X;
        
        
        
        
        xcoord=ltotal;
        for(i=0;i<N.xsecs;i++) {
                xcoord -= xsp->loc;
                xsp->xcoord = xcoord;
                temp=xsp;
                xsp = xsp->nextxsp;
                xsp->previousxsp=temp;

        }
        
        
        
        
printf("ltotal=%g\tymax=%g\n",ltotal,ymax);                     
                                
        put_xsecdata(outf);
        fclose(outf);
        return(0);
}

int     put_xsecdata(f)
FILE    *f;
{
        int     i,j;
        struct  xsec    *xsp;
        
        
        fprintf(f,"node x       y       bl      ks                      u[three]        bcp(three)      bcs\n");

//put info for outflow boundary
        xsp=gp.X;
        xsp=xsp->previousxsp;
        for(i=0;i<N.xsecs;i++) {
                fprintf(f,"xross section number = %d\t\tlocation=%g\n",xsp->n,xsp->xcoord);
                fprintf(f,"y    bl-meas wsl-meas        depth   u-meas\n");     
        
                for(j=0;j<xsp->ndatapts;j++) {
                        if(xsp->dp[j]->u<0.0)
                                xsp->dp[j]->u=0.0;
                        put_veldata(f,xsp,xsp->dp[j],j+1);
                }
                xsp=xsp->previousxsp;
        }
        return(0);
}


int     put_veldata(f,xsp,dp,j)
FILE    *f;
struct  xsec    *xsp;
struct  datapoint       *dp;
int     j;
{
        fprintf(f,"%g\t%g\t%g\t%g\t%g\n",
                        dp->y,dp->b,xsp->wsl,(xsp->wsl-dp->b),dp->u);
                
        return(0);
}







int     put_meshdata(f)
FILE    *f;
{
        int     i,j;
        int     bcin=1,bcout=3;
        struct  xsec    *xsp;
        double shift=0.001;
        
        xsp=gp.X;
        xsp=xsp->previousxsp;
        
        fprintf(f,"number of nodes per element = 3\n"); 
        fprintf(f,"number of loops = 1\n");
        fprintf(f,"number of ordering segments = 1\n");
        fprintf(f,"order segment pairs = \t");
        fprintf(f,"%d\t%d\t",(N.xsecs+1),(N.xsecs+2));
        fprintf(f,"\n");
        fprintf(f,"node x       y       bl      ks                      u[three]        bcp(three)      bcs\n");

//put info for first corner node on outflow boundary

        xsp=gp.X;
        
        put_nodedata(f,xsp->dp[xsp->ifirst],-shift,shift);
        put_bdata(f,xsp->dp[xsp->ifirst],bcout);

//put information for last node on edge (no flow bc)
        put_nodedata(f,xsp->dp[(xsp->ilast)],-shift,-shift);
        put_bdata(f,xsp->dp[(xsp->ilast)],0);
        fprintf(f,"\n");
//put info of upper no flow boundary
        xsp=xsp->nextxsp;
        for(j=0;j<(N.xsecs-2);j++) {
                put_nodedata(f,xsp->dp[xsp->ilast],0.0,-shift);
                put_bdata(f,xsp->dp[xsp->ilast],0);
                xsp=xsp->nextxsp;
        }
        fprintf(f,"\n");
//put info of inflow boundary

        xsp=gp.X;
        xsp=xsp->previousxsp;
        
        put_nodedata(f,xsp->dp[xsp->ilast],shift,-shift);
        put_bdata(f,xsp->dp[xsp->ilast],bcin);
        
//put information for last node on edge (no flow bc)
        put_nodedata(f,xsp->dp[(xsp->ifirst)],shift,shift);
        put_bdata(f,xsp->dp[(xsp->ifirst)],0);
        fprintf(f,"\n");
//put info on lower no flow boundary
        xsp=xsp->previousxsp;
        for(j=0;j<(N.xsecs-2);j++) {
                put_nodedata(f,xsp->dp[xsp->ifirst],0.0,shift);
                put_bdata(f,xsp->dp[xsp->ifirst],0);
printf("xsp->n=%d\txcoord=%g\txsp->loc=%g\n",xsp->n,xsp->xcoord,xsp->loc);
                xsp=xsp->previousxsp;

        }
        
        fprintf(f,"-1\n");
        fprintf(f,"-1\n");
        
        return(0);
}

int     put_interpolationdata(f)
FILE    *f;
{
        int     i,j;
        int     bcin=1,bcout=3;
        struct  xsec    *xsp;
        
        xsp=gp.X;
        xsp=xsp->previousxsp;
        
        fprintf(f,"number of nodes per element = 3\n"); 
        fprintf(f,"number of loops = 1\n");
        fprintf(f,"node x       y       bl      ks                      u[three]        bcp(three)      bcs\n");

//put info for outflow boundary
        xsp=gp.X;
        for(i=xsp->ifirst;i<(xsp->ilast);i++) {
//              xsp->dp[i]->hb=(xsp->dp[i]->h+xsp->dp[i+1]->h)/2.0;
                put_nodedata(f,xsp->dp[i],0.0,0.0);
                put_bdata(f,xsp->dp[i],bcout);
        }
//put information for last node on edge (no flow bc)
        put_nodedata(f,xsp->dp[(xsp->ilast)],0.0,0.0);
//        put_bdata(f,xsp->dp[(xsp->ilast)],0.0,0.0);
        fprintf(f,"\n");
//put info of upper no flow boundary
        xsp=xsp->nextxsp;
        for(j=0;j<(N.xsecs-2);j++) {
                put_nodedata(f,xsp->dp[xsp->ilast],0.0,0.0);
                put_bdata(f,xsp->dp[xsp->ilast],0);
                xsp=xsp->nextxsp;
        }
        fprintf(f,"\n");
//put info of inflow boundary
        for(i=xsp->ilast;i>(xsp->ifirst);i--) {
//              xsp->dp[i]->hb=(xsp->dp[i]->h+xsp->dp[i-1]->h)/2.0;
                put_nodedata(f,xsp->dp[i],0.0,0.0);
                put_bdata(f,xsp->dp[i],bcin);
        }
//put information for last node on edge (no flow bc)
        put_nodedata(f,xsp->dp[(xsp->ifirst)],0.0,0.0);
        put_bdata(f,xsp->dp[(xsp->ifirst)],0);
        fprintf(f,"\n");
//put info on lower no flow boundary
        xsp=xsp->previousxsp;
        for(j=0;j<(N.xsecs-2);j++) {
                put_nodedata(f,xsp->dp[xsp->ifirst],0.0,0.0);
                put_bdata(f,xsp->dp[xsp->ifirst],0);
printf("xsp->n=%d\txcoord=%g\txsp->loc=%g\n",xsp->n,xsp->xcoord,xsp->loc);
                xsp=xsp->previousxsp;

        }
        fprintf(f,"-1\n");
//put info for interior nodes at each x-sec
        xsp=gp.X;
        xsp=xsp->nextxsp;
        for(j=0;j<(N.xsecs-2);j++) {
                for(i=(xsp->ifirst+1);i<(xsp->ilast-1);i++) {
                        put_nodedata(f,xsp->dp[i],0.0,0.0);
                        fprintf(f,"\n");
                }
                xsp=xsp->nextxsp;
                fprintf(f,"\n");
        }
        fprintf(f,"-1\n");
        
        return(0);
}


int     put_bdata(f,dp,bc)
FILE    *f;
struct  datapoint       *dp;
int     bc;
{
        fprintf(f,"%g\t%g\t%g\t%d\n",dp->hb,dp->qxb,dp->qy,bc);
        return(0);
}


int     put_nodedata(f,dp,dx,dy)
FILE    *f;
struct  datapoint       *dp;
double  dx,dy;
{
        fprintf(f,"%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t",dp->n,(dp->xcoord+dx)
                        ,(dp->ycoord+dy),dp->b,dp->man,dp->h,dp->qx,dp->qy);
                
        return(0);
}

int     complete_data(ltotal,ymax,dwsl)
double  ltotal;
double  ymax,dwsl;
{
        int     i,j;
        double  xcoord,dum,factor;
        struct  xsec    *xsp,*temp;
        struct  datapoint       *dp,*alloc_dp();
        double  h1,h2,dy;

        xsp=gp.X;
        temp=gp.X;
        for(i=0;i<(N.xsecs-1);i++) 
                temp=temp->nextxsp;
        temp->nextxsp=gp.X;
        
        
        
        
        xcoord=ltotal;
        for(i=0;i<N.xsecs;i++) {
                xcoord -= xsp->loc;
                xsp->xcoord = xcoord;
                temp=xsp;
                xsp = xsp->nextxsp;
                xsp->previousxsp=temp;

        }
        xsp=gp.X;
        for(j=0;j<N.xsecs;j++) {
                for(i=0;i<xsp->ndatapts;i++) {
                        dp=xsp->dp[i];
                        dp->xcoord=xsp->xcoord;
                        dp->ycoord=ymax/2.0+dp->y-(xsp->dp[xsp->ndatapts-1]->y)/2.0;
                        dp->h=xsp->wsl-dp->b;
                        if(dp->h>dwsl) {
                                if(xsp->nwet==0)
                                        xsp->ifirst=i;
                                xsp->nwet++;
                                xsp->ilast=i;
//                              xsp->sigmahp+=pow(dp->h,5.0/3.0);
                        }
                }
                if((xsp->ifirst+xsp->nwet-1)!=xsp->ilast)
                        printf("WARNING:Dry points exist in x-sec %d\n",xsp->n);
                xsp=xsp->nextxsp;
        }
//Now insert two interpolated data points which have a depth =dwsl if required
printf("1\n");
        xsp=gp.X;

        for(j=0;j<N.xsecs;j++) {
        if(xsp->ifirst!=0){
printf("1.02\n");
printf("xsp->ifirst=%d\n",xsp->ifirst);
                h1=-(xsp->dp[xsp->ifirst-1]->h-dwsl);
                h2=xsp->dp[xsp->ifirst]->h-dwsl;
printf("h1=%g\th2=%g\n",h1,h2);
                dy=(xsp->dp[xsp->ifirst]->ycoord-xsp->dp[xsp->ifirst-1]->ycoord)*h2/(h1+h2);
//printf("xsp=%g\th1=%g\th2=%g\tdy=%g\n",xsp->xcoord,h1,h2,dy);
printf("1.1\n");
                if(dy>0.1*(xsp->dp[xsp->ifirst]->ycoord-xsp->dp[xsp->ifirst-1]->ycoord)) {
                        xsp->ndatapts++;
                        xsp->nwet++;
                        for(i=(xsp->ndatapts-1);i>xsp->ifirst;i--) 
                                xsp->dp[i]=xsp->dp[i-1];
printf("1.2\n");
        
                        xsp->dp[xsp->ifirst]=alloc_dp();
printf("1.3\n");

                        xsp->dp[xsp->ifirst]->h=dwsl;
                        xsp->dp[xsp->ifirst]->b=xsp->wsl-dwsl;
                        xsp->dp[xsp->ifirst]->ycoord=xsp->dp[xsp->ifirst+1]->ycoord-dy;
                        xsp->dp[xsp->ifirst]->xcoord=xsp->xcoord;
                        xsp->ilast++;
        
                }
//printf("ifirst=%d\tndatapts=%d\tilast=%d\tnwet=%d\n",xsp->ifirst,xsp->ndatapts,xsp->ilast,xsp->nwet);
//for(i=0;i<xsp->ndatapts;i++)
//printf("i=%d\tdatapoint=%d\th=%g\n",i,xsp->dp[i],xsp->dp[i]->h);
//Now for ilast
printf("2.5\n");
                
                h1=xsp->dp[xsp->ilast]->h-dwsl;
                h2=-(xsp->dp[xsp->ilast+1]->h-dwsl);
                dy=(xsp->dp[xsp->ilast+1]->ycoord-xsp->dp[xsp->ilast]->ycoord)*h1/(h1+h2);
                if(dy>(0.1*(xsp->dp[xsp->ilast+1]->ycoord-xsp->dp[xsp->ilast]->ycoord))) {
                        xsp->ndatapts++;
                        xsp->nwet++;
                        for(i=(xsp->ndatapts-1);i>(xsp->ilast+1);i--) 
                                xsp->dp[i]=xsp->dp[i-1];
                        xsp->ilast++;
                        xsp->dp[xsp->ilast]=alloc_dp();
                        xsp->dp[xsp->ilast]->h=dwsl;
                        xsp->dp[xsp->ilast]->b=xsp->wsl-dwsl;
                        xsp->dp[xsp->ilast]->ycoord=xsp->dp[xsp->ilast-1]->ycoord+dy;
                        xsp->dp[xsp->ilast]->xcoord=xsp->xcoord;
                }
                
                xsp=xsp->nextxsp;
        }
        }
        
printf("2\n");
        
/*      
        xsp=gp.X;
        for(j=0;j<N.xsecs;j++) {
                for(i=xsp->ifirst;i<(xsp->ilast+1);i++) {
                        dp=xsp->dp[i];
                        dp->qx=xsp->Q/(xsp->dp[xsp->ilast]->y-xsp->dp[xsp->ifirst]->y);
                        dp->qxb=dp->qx;
                        dp->man=0.17;
                }
                xsp=xsp->nextxsp;
        }
*/
//      calculate qxb and hb for inflow boundary
        xsp=gp.X;
        xsp=xsp->previousxsp;
        
        for(i=xsp->ilast;i>xsp->ifirst;i--) {
                xsp->dp[i]->hb=(xsp->dp[i]->h+xsp->dp[i-1]->h)/2.0;
                if(xsp->dp[i]->hb>0.0)
                        xsp->sigmaHp+=(pow(xsp->dp[i]->hb,5.0/3.0)*(xsp->dp[i]->y
                                                                                -xsp->dp[i-1]->y));
        }
        factor=xsp->Q/xsp->sigmaHp;
        for(i=xsp->ilast;i>xsp->ifirst;i--) {
                if(xsp->dp[i]->hb>0.0)
                        xsp->dp[i]->qxb=pow(xsp->dp[i]->hb,5.0/3.0)*factor;
                else
                        xsp->dp[i]->qxb=0.0;
        }
//calculate hb for outflow boundary

        xsp=gp.X;
        for(i=xsp->ifirst;i<xsp->ilast;i++)
                xsp->dp[i]->hb=(xsp->dp[i]->h+xsp->dp[i+1]->h)/2.0;
                
//      calculate estimate for qx at all data points
// note that this calculation is not very exact,
// because hb is used for calculation of sigmahp, while
// h is used to calculate qx
        xsp=gp.X;
//printf("Q=%g\tsigmahp=%g\n",xsp->Q,xsp->sigmahp);
        for(j=0;j<N.xsecs;j++) {
                for(i=(xsp->ifirst+1);i<(xsp->ilast+1);i++) {
                        xsp->dp[i]->hb=(xsp->dp[i]->h+xsp->dp[i-1]->h)/2.0;
                        if(xsp->dp[i]->hb>0.0)
                                xsp->sigmahp+=(pow(xsp->dp[i]->hb,5.0/3.0)*(xsp->dp[i]->y
                                                                                -xsp->dp[i-1]->y));
//printf("i=%d\thb=%g\thb5thirds=%g\twidth=%g\tsigmahp=%g\th=%g\n",i,xsp->dp[i]->hb,pow(xsp->dp[i]->hb,5.0/3.0),(xsp->dp[i]->y-xsp->dp[i-1]->y),xsp->sigmahp,xsp->dp[i]->h);
                        
                }
                factor=xsp->Q/xsp->sigmahp;
//printf("x-coord=%g\tfactor=%g\tsigmahptimesw=%g\n",xsp->xcoord,factor,xsp->sigmahp);
                for(i=xsp->ifirst;i<(xsp->ilast+1);i++) {
                        if(xsp->dp[i]->h>0.0)
                                xsp->dp[i]->qx=pow(xsp->dp[i]->h,5.0/3.0)*factor;
                        else
                                xsp->dp[i]->qx=0.0;
                        xsp->dp[i]->man=0.31;
                }
                
                xsp=xsp->nextxsp;
        }
printf("3\n");
        
        return(0);
}






struct  xsec    *alloc_xsec()

{
        struct xsec     *xsp;
        int     size;
        
        size = sizeof(struct xsec);
        xsp = (struct xsec *) malloc(size);
        xsp->nextxsp=NULL;
        xsp->ndatapts=0;
        xsp->nwet=0;
        xsp->sigmahp=0.0;
        xsp->sigmaHp=0.0;
        xsp->qx=0.0;
        N.xsecs++;
        return(xsp);
}

struct  datapoint       *alloc_dp()

{
        struct datapoint        *dp;
        int     size;
        
        size = sizeof(struct datapoint);
        dp = (struct datapoint *) malloc(size);
        dp->n=N.datapts;
        dp->nwet=0;
        dp->b=0.0;
        dp->y=0.0;
        dp->man=0.0;
        dp->xcoord=0.0;
        dp->ycoord=0.0;
        dp->h=-1.0;
        dp->u=0.0;
        dp->v=0.0;
        dp->weight=0.0;
        dp->hb=0.0;
        dp->qxb=0.0;
        dp->qx=0.0;
        dp->qy=0.0;
        N.datapts++;
        return(dp);
}

double  get_limiteddouble(n,f)
int     n;
FILE    *f;
{
        char    s[25],*cs;
        int     i,c,k;
        double  d;
        
        for(i=0;i<25;i++)
                s[i]=' ';
        cs=s;
        k=0;
        for(i=0;i<n;i++) {
                if((c=getc(f))!='\n') {
                        if((isdigit(c)!=0) || (c=='.')||(c=='-')||(c=='+')) 
                                *cs++=c;
                }
                else {
                        k=1;
                        break;
                }
        }
        if((k==1)||(cs==s)) {
                if(k==1)
                        d=-1;
                if((cs==s)&&(k==0))
                        d=-2;
        }
        else
                d=atof(s);
        return(d);

}
        







int     get_xsecnum(n,f)
int     n;
FILE    *f;

{       
        char    s[25],*cs;
        int     i,c;
        double  d;
        
        for(i=0;i<25;i++)
                s[i]=' ';
        cs=s;
        for(i=0;i<n;i++) {
                c=getc(f);
                if((isdigit(c)!=0) || (c=='.')||(c=='-')||(c=='+')) 
                        *cs++=c;
        }
        d=atof(s);
        return((int) d);

}
        


int     getlinename(n,f)
int     n;
FILE    *f;

{
        int     i;
        
        for(i=0;i<n;i++)
                getc(f);
        return(0);
}
                
                
                
                
                
                
int     user_getline(line,f)
char    *line;
FILE    *f;
{
        int     c;
        char    *cs;
        
        cs=line;
        
        while((c=getc(f))!=EOF) {
                *cs++=c;
                if(c=='\n') 
                        break;  
        }
        return(cs-line);
}

int     putline(line,f,n)
char    *line;
FILE    *f;
int     n;
{
        int     i;
        for(i=0;i<n;i++) 
//              c=line[i]
                putc(line[i],f);
        return(0);
}

int             outxsec(nin)
int             nin ;

{
        int             err,i ;
        FILE    *fptr, *put_fptr() ;
        struct node     *np ;
        double  coord;
        
        printf("Enter coordinate Value\n");
        scanf(" %lf",&coord);
        if((fptr = put_fptr(nin) ) != NULL) {
                fprintf(fptr,"\n Distance (m), Depth (m), qx, qy ,BL , Stage, U, V\n");
                np = gp.N ;
                for(i=0;i<N.nodes;i++) {
                        if(np->x[nin]==coord)
                                putxsec_node(fptr,i,np,nin) ;
                        np = np->nextnp ;
                }
                if(fptr != stdin) fclose(fptr) ;
        }
        else
                printf(" Unable to open that file\n") ;
        return(err);
}

int             putxsec_node(f,i,nodep,nin)
FILE    *f ;
int             i,nin ;
struct node  *nodep ;

{
        int     j ;
        
        put_dbl(f,(double) nodep->x[abs(nin-1)]) ;
        for(j=0;j<N.vars;j++)
                put_dbl(f,(double) nodep->u[j]) ;
        put_dbl(f,(double) nodep->p[0]) ;
        put_dbl(f,(double) (nodep->p[0]+nodep->u[0]));
        put_dbl(f,(double) (nodep->u[1]/nodep->u[0]));
        put_dbl(f,(double) (nodep->u[2]/nodep->u[0]));
        fprintf(f," \n");
        return(0) ;
}




