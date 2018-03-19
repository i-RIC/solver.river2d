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

int             tecplot_output(nin, s) //modified by C.Frias Jun 16 2010 for tecplot output
int             nin ; 
char            s[63];//added june 16


{
        int             err,i ;
        FILE    *fptr, *put_fptr() ;
        struct node     *np ;
        struct element  *elp ;
        struct belement *belp ;
        struct  boundaryseg     *bseg;

        if((fptr = fopen(s, "a") ) != NULL) { //modified june 16
	  if( (err = tecplot_put_control(fptr)) == 0) {
		  np=gp.N;

                        for(i=0;i<N.nodes;i++) {

                                tecplot_put_node(fptr,i,np) ;
                                np = np->nextnp ;
                        }
                 
                        fprintf(fptr,"\n");

                        elp = gp.El ;
                        for(i=0;i<N.elms;i++) {
                                tecplot_put_elm(fptr,i,elp) ;
                                elp = elp->nextelp ;
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

int             tecplot_put_control(f) //put control routine modified by C.Frias June 16 2010
FILE    *f ;


{

        fprintf(f,"TITLE = \"RIVER 2D OUTPUT FOR TECPLOT\" \n");
        fprintf(f,"VARIABLES = \"X\", \"Y\", \"Z\", \"k\", \"code\", \"DEPTH\", \"U\", \"V\" , \"tau_bx\", \"WSelev\" , \"U2\" , \"V2\" , \"tau_bx2\" \n");
        fprintf(f, "ZONE T=\"time=%.2f\", DATAPACKING=POINT, NODES= %d, ELEMENTS= %d, ZONETYPE=FETRIANGLE \n", tvals.t, N.nodes, N.elms);      
        return(0) ;
}

int             tecplot_put_node(f,i,nodep) //put_node routine modified by C.Frias Jun 16 2010
FILE    *f ;
int             i ;
struct node  *nodep ;

{
          int     j ;
	  double Cs;
	  double rho=1000;
	  double g=9.81;
	  double ks;
	  double tau_bx;
	  double U,V,H;
	  double water_elevation,bed_elevation;

	  tau_bx = 0.0;

        
        for(j=0;j<N.dims;j++) 
                put_dbl(f,(double) nodep->x[j]) ;
        for(j=0;j<N.params;j++)
                put_dbl(f,(double) nodep->p[j]) ;
        for(j=0;j<N.vars;j++)
                put_dbl(f,(double) nodep->u[j]) ;

	ks=nodep->p[1];
	U=nodep->u[1];
	V=nodep->u[2];
	H=nodep->u[0];
	bed_elevation=nodep->p[0];

	Cs=5.75*log (12* H / ks);
	if(H>0)
	  {
	tau_bx=sqrt(U*U +V*V) * U / (Cs *Cs) * rho ; //calculates the bed shear stress;
	put_dbl(f,(double) tau_bx);
	  }
	else
	put_dbl(f,0);

	water_elevation=H+bed_elevation;//calculates the water surface elevation
	put_dbl(f,(double) water_elevation);

	put_dbl(f,(double) U / H); // outputto see if u[] is actually a vector of discharges
	put_dbl(f, (double) V /H); //output to see if u[] is actually a vector of discharges
	put_dbl(f, (double) tau_bx / (H*H)); //output to see if u[] is actually a vector of discharges

        fprintf(f," \n");
        return(0) ;
}

int             tecplot_put_elm(f,i,elmntp) //put_elm routine modified by C.Frias June 16 2010
int             i ;
FILE    *f ;
struct element *elmntp ;

{
        int             j;
        
        for(j=0;j<elmntp->nnds;j++) {
                if(elmntp->nps[j] != NULL)
                        put_int(f,(elmntp->nps[j])->n) ;
                else
                        put_int(f,-1) ;
        }
        fprintf(f," \n");
        return(0) ;
}