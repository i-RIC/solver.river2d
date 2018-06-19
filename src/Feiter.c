#include "Fe_PS.h"
#include <time.h>

#define KE(A,B) ( *(ep.K + ep.n*(A) + (B)) )
#define SE(A,B) ( *(ep.S + ep.n*(A) + (B)) )
#define ME(A,B) ( *(ep.M + ep.n*(A) + (B)) )
#define JE(A,B) ( *(ep.J + ep.n*(A) + (B)) )
#define FE(B)   ( *(ep.F + B) )
#define Vp(A,B) ( *(eqnsets[0].Vp + eqnsets[0].m*(A) + (B)) )
#define Hp(A,B) ( *(eqnsets[0].Hp + eqnsets[0].m*(A) + (B)) )
#define GRAV    9.806
#define E	2.718

extern struct control          N ;
extern struct pointers         gp;
extern struct transient        tvals ;
extern int                     Nukns ;
extern struct elmpointers      ep ;
extern struct eqnset           eqnsets[4] ; 
extern struct settings sets;
extern int						lumpS;
extern double                  uchange, maxChange;
double						   memGMRES;
extern int   				   outputts, outputtsnum;
extern double				   goalt;
extern double	resL2, resLinf;
extern int		resmaxi;
extern int		jacobianType;
int				preconditionerType;
extern void	updateVelocities();
extern int		PCILU_gmres();
extern int		PCJAC_gmres();
extern void	findMaxchange();//added by mostafa Nov 01,09
extern void	setGWFtoZero();// added by mostafa Nov 01,09
extern double avgCelerity(); // added by mostafa Nov 02,09
extern void useHalfdp(nset);
extern struct  node   *maxnodeP;//added by mostafa Nov 01,09

/*** This is the good version of steadynew DO NO EDIT  ***/
int     steadynew(nvar,elcode,tol,dtmax,fp,i,m,k,gmrestol )
int     nvar,elcode,i ;
int		m,k ; //m = No. of Krylov Vectors or No. of steps before restart k = No. of gmres iterations
double  tol,dtmax,gmrestol ;
FILE* fp;

{
        int     itnum, err ;
        FILE    *fptr;
        double  fac, allowChg=0.0, avgCel;
       	int		update_BCs(), reset_u(), set_ptou(), gmres_init()/*, PC_gmres()*/;
        double	avgCelerity(), test_outflow();
		void	setGWFtoZero(), findMaxChange();

		N.trans = 1 ;
        lumpS = 1;
		
		if((fptr = fp ) != NULL)
		{
			if(i == 1)
				fprintf(fptr, "m (number of steps before restart) = %d\nk (number of gmres iterations) = %d\ngmres tolerance = %f\n", m, k, gmrestol);
			tvals.t += tvals.dt ;
			update_BCs();
			itnum = 0 ;

			if( (err = gmres_init(0, m)) != 0) return(err);
			eqnsets[0].assembleFlag = 1;
			if(preconditionerType == 0)
				err = PCJAC_gmres(0, itnum, nvar, elcode, m, k, gmrestol/*, fptr*/);
			else
				err = PCILU_gmres(0, itnum, nvar, elcode, m, k, gmrestol/*, fptr*/);

			setGWFtoZero();
			findMaxchange();
			avgCel = avgCelerity();
			uchange=fabs(maxChange)/avgCel;
			allowChg=tvals.uc*avgCel;

/*PS*/      if((fabs(maxChange) > 1.5 * allowChg) || (err>=k))  //if maxChg is too big reject iteration and write info to logfile											
			{
				fprintf(fptr,"\n\tmaxChange > 1.5*allowChange Time Step Rejected\n");
				fprintf(fptr,"# %d  t= %f  dt= %f allowablechange=%E maxchange=%E solution change= %E\n",
								i,tvals.t,tvals.dt,allowChg, maxChange,uchange);
				fprintf(fptr,"\tmax change occur at node %d \tx = %f \ty = %f\n\t\t from h = %f qx = %f qy = %f \n\t\t to h = %f qx = %f qy = %f\n\n",
                      maxnodeP->n,maxnodeP->x[0],maxnodeP->x[1],maxnodeP->uo[0],
					  maxnodeP->uo[1],maxnodeP->uo[2],
                      maxnodeP->u[0],maxnodeP->u[1],maxnodeP->u[2]);

				tvals.t -= tvals.dt;
				fac = 0.5 * allowChg/fabs(maxChange);  //lower the time step
/*PS*/			if(err>=k) fac = 0.5;
				reset_u(0); //reset the values at the nodes to the values 
							//at the previous time step
			}
			else if(fabs(maxChange) > allowChg) // if true, then reduce the time step but accept iteration and write info to logfile
			{
				fprintf(fptr,"\n\tmaxChg > allowChg Time Step Accepted but deltat reduced\n");
				fprintf(fptr,"# %d  t= %f  dt= %f  allowablechange= %E  maxchange= %E solution change= %E\n",
								i,tvals.t,tvals.dt,allowChg,maxChange,uchange);
				fprintf(fptr,"\tmax change occur at node %d \tx = %f \ty = %f\n\t\t from h = %f qx = %f qy = %f \n\t\t to h = %f qx = %f qy = %f\n\n",
                      maxnodeP->n,maxnodeP->x[0],maxnodeP->x[1],maxnodeP->uo[0],
					  maxnodeP->uo[1],maxnodeP->uo[2],
                      maxnodeP->u[0],maxnodeP->u[1],maxnodeP->u[2]);

				
				set_ptou(0); //copies the values at the nodes this time step to the space 
							 //allocated for the those at the previous time step
				updateVelocities();
				fac = 0.5*allowChg/fabs(maxChange);
			}
			else //amount of change is acceptable; write results to logfile and increase time step
			{
				fprintf(fptr,"\n\tmaxChg < allowChg Time Step Accepted\n");
				fprintf(fptr,"# %d  t= %f  dt= %f  allowablechange= %E  maxchange= %E uchange= %E\n",
								i,tvals.t,tvals.dt,allowChg,maxChange,uchange);
				set_ptou(0); //copies the values at the nodes this time step to the space 
							 //allocated for the those at the previous time step
				updateVelocities();
				fac = allowChg/fabs(maxChange); //factor = goalchange (default 0.05) / uchange 
				if(fac >1.5) fac = 1.5;
/*PS*/			if((err> k/2) && (fac > 1.0)) fac = 1.0;
			}
			tvals.dt *=fac;
			if(tvals.dt > dtmax)
				tvals.dt = dtmax;
		}
        return(0);
}

int     transientnew(nvar,elcode,tol,goaldt,fp,i,m,k,gmrestol)
int     nvar,elcode,i ;
int		m,k ; //m = No. of Krylov Vectors or No. of steps before restart k = No. of gmres iterations
double  tol,gmrestol, goaldt ;
FILE* fp;

{
        int     j, itnum, err, gmrestotalitnum;
        FILE    *fptr, *put_fptr() ;
        double  fac,test_outflow();
		int		update_BCs(), set_ptou(), gmres_init(),PCJAC_gmres();
		double	olddt, oldt;
		void	setGWFtoZero(), findMaxchange(), useHalfdp();
		double	maxChangeOld;
        
        N.trans = 1 ;
        lumpS = 0;  //should be 0 for transient problems but could try 1 to see if any difference

		outputts = 0;
        
        if((fptr = fp ) != NULL)
		{ 
			if(i == 1)
				fprintf(fptr, "m (number of steps before restart) = %d\nk (number of gmres iterations) = %d\ngmres tolerance = %f\n", m, k, gmrestol);
			olddt = 0.0;
			oldt  = tvals.t;
			tvals.t +=tvals.dt;
			if (tvals.t > goalt)
			{
				olddt = tvals.dt;
				tvals.dt = goalt - oldt;
				tvals.t = goalt;
			}
			update_BCs();
			maxChange = 1.0;
			tvals.itnum = 0;
			itnum = 0;
			gmrestotalitnum = 0;
			for(j=0;j<sets.maxitnum;j++) 
			{  
				if(fabs(maxChange)<tol)
					break ;
				if( (err = gmres_init(0, m)) != 0) return(err);
				eqnsets[0].assembleFlag = 1;
				if(preconditionerType == 0)
					err = PCJAC_gmres(0, itnum, nvar, elcode, m, k, gmrestol/*, fptr*/);
				else
					err = PCILU_gmres(0, itnum, nvar, elcode, m, k, gmrestol/*, fptr*/);
				gmrestotalitnum += err;
				tvals.itnum += 1 ;
				itnum += 1;
				setGWFtoZero();
				maxChangeOld = maxChange;
				findMaxchange();
				fprintf(fptr,"itnum = %d\t maxChange = %f\t maxNode = %d\n", tvals.itnum, maxChange,maxnodeP->n);
				
					if(maxChange == 0.0){
					maxChange = 1.0;				
					break;
				}
				if(fabs(maxChange) > tvals.uc*avgCelerity())
					break;

				if(fabs(maxChange) > fabs(0.9*maxChangeOld)){
					useHalfdp(0);
					findMaxchange();
				}
			
			}

			if((fabs(maxChange)>tol) || (err>=k))//reject the time step
			{					
				fprintf(fptr,"\n\tdt = %f didn't work\n",tvals.dt);
				fprintf(fptr,"\tMaximum change at node %d \tx = %f \ty = %f\n\t\t from h = %f qx = %f qy = %f \n\t\t to   h = %f qx = %f qy = %f\n",
                      maxnodeP->n,maxnodeP->x[0],maxnodeP->x[1],maxnodeP->uo[0],
					  maxnodeP->uo[1],maxnodeP->uo[2],
                      maxnodeP->u[0],maxnodeP->u[1],maxnodeP->u[2]);

				tvals.t -= tvals.dt;
				fac = 0.5;
/*PS*/			if(err>=k) fac = 0.5;
				reset_u(0); //reset the values at the nodes to the values 
							//at the previous time step
				tvals.dt *= fac;
			}
			else
			{            
				fprintf(fptr,"time step # %d\tt = %fs\tdt = %fs\ttotal inflow = %g\ttotal outflow = %g\t# of N.R. iterations = %d\t average # of gmres iterations = %d Maxchange%g\n",i,tvals.t,tvals.dt,test_outflow(1),test_outflow(2),tvals.itnum,gmrestotalitnum/itnum,maxChange);               
				set_ptou(0) ;
				updateVelocities();
				fac = 1.5;
				if((err> k/2) && (fac > 1.0)) fac = 1.0;
				tvals.dt *= fac;
				if (olddt != 0.0)
					tvals.dt = olddt;
				if (tvals.t == goalt)
				{
					goalt += goaldt;
					outputts = 1;
					outputtsnum += 1;
				}
			}
			if(tvals.dt > goaldt)
				tvals.dt = goaldt;
        }
        return(0) ;
}

int gmres_init(nset, m)
int nset, m;
{
		double	*F_init(), initm();	
		int n, i;
		int ep_init();
		struct  node    *nodep;
		double res;

		n = N.vars ;
		Nukns = N.nodes * n;
		eqnsets[nset].neqns = Nukns ;
		eqnsets[nset].m = m;

		if (eqnsets[nset].assembleFlag == 0)
		{
			if((eqnsets[nset].Kp != NULL) && (eqnsets[nset].assembleFlag == 0) )	//freeing up memory if direct 
                        free(eqnsets[nset].Kp) ;					//solver was called before gmres

			if((eqnsets[nset].Lp != NULL) && (eqnsets[nset].assembleFlag == 0) )	//freeing up memory if direct 
					   free(eqnsets[nset].Lp) ;						//solver was called before gmres

			if((eqnsets[nset].Fp != NULL) && (eqnsets[nset].assembleFlag == 0))
			    free(eqnsets[nset].Fp) ;
       
	        if( ( eqnsets[nset].Fp = F_init(Nukns) )  == NULL  ) 
				return(-1) ;

			if((eqnsets[nset].Bp != NULL) && (eqnsets[nset].assembleFlag == 0))
		        free(eqnsets[nset].Bp) ;
        
	        if( ( eqnsets[nset].Bp = F_init(Nukns) )  == NULL  ) 
				return(-1) ;

			if((eqnsets[nset].Qp != NULL) && (eqnsets[nset].assembleFlag == 0))
		        free(eqnsets[nset].Qp) ;
        
	        if( ( eqnsets[nset].Qp = F_init(Nukns) )  == NULL  ) 
				return(-1) ;

			if((eqnsets[nset].Rp != NULL) && (eqnsets[nset].assembleFlag == 0))
		        free(eqnsets[nset].Rp) ;
	        
	        if( ( eqnsets[nset].Rp = F_init(m+1) )  == NULL  ) 
				return(-1) ;

			if((eqnsets[nset].Dp != NULL) && (eqnsets[nset].assembleFlag == 0))
		        free(eqnsets[nset].Dp) ;

			if( ( eqnsets[nset].Dp = F_init(Nukns) )  == NULL  ) 
				return(-1) ;

			if((eqnsets[nset].Hp != NULL) && (eqnsets[nset].assembleFlag == 0))
			    free(eqnsets[nset].Hp) ;

			if( ( eqnsets[nset].Hp = (double *) calloc((m+1)*m, sizeof(double)) ) == NULL )
				return(-1) ;

			if((eqnsets[nset].Vp != NULL) && (eqnsets[nset].assembleFlag == 0))
			    free(eqnsets[nset].Vp) ;

			if( ( eqnsets[nset].Vp = (double *) calloc(Nukns*m, sizeof(double)) ) == NULL )
				return(-1) ;

			if((eqnsets[nset].Mp != NULL) && (eqnsets[nset].assembleFlag == 0))
		        free(eqnsets[nset].Mp) ;
	        
	        if( ( eqnsets[nset].Mp = F_init(Nukns) )  == NULL  ) 
				return(-1) ;

			if((eqnsets[nset].Yp != NULL) && (eqnsets[nset].assembleFlag == 0))
		        free(eqnsets[nset].Mp) ;
	        
	        if( ( eqnsets[nset].Yp = F_init(Nukns) )  == NULL  ) 
				return(-1) ;

			if( (ep.n = ep_init() ) == 0 )
				 return(-1) ;

			//allocate memory for the block diagonals
			nodep = gp.N;
			for(i=0;i<N.nodes;i++)
			{
				if( ( nodep->blockDiagKp = (double *) calloc(n*n, sizeof(double)) ) == NULL )
				return(-1) ;
				nodep = nodep->nextnp;
			}

			memGMRES = ((720+sizeof(struct element))*N.elms + sizeof(struct node)*N.nodes + (336+sizeof(struct belement))*N.belms)/1024.0/1024.0;
		}

	//	if (eqnsets[nset].assembleFlag == 0)  						
	//	{
			for(i=0;i<eqnsets[nset].neqns;i++)		// these are the delta phis 
					*(eqnsets[nset].Fp + i) = 0.0 ;	// only want to set to zero at 
	//	}											// first iteration

		for(i=0;i<eqnsets[nset].neqns;i++)
 				*(eqnsets[nset].Bp + i) = 0.0 ;

		for(i=0;i<eqnsets[nset].neqns;i++)
				*(eqnsets[nset].Qp + i) = 0.0 ;

		for(i=0;i<(m+1);i++)
				*(eqnsets[nset].Rp + i) = 0.0 ;

		for(i=0;i<eqnsets[nset].neqns;i++)
				*(eqnsets[nset].Dp + i) = 0.0 ;

		for(i=0;i<((m+1)*m);i++)
				*(eqnsets[nset].Hp + i) = 0.0 ;

		for(i=0;i<(Nukns*m);i++)
				*(eqnsets[nset].Vp + i) = 0.0 ;

		for(i=0;i<eqnsets[nset].neqns;i++)
				*(eqnsets[nset].Mp + i) = 0.0 ;

		for(i=0;i<eqnsets[nset].neqns;i++)
				*(eqnsets[nset].Yp + i) = 0.0 ;

		//reset the block diagonals to zero
		nodep = gp.N;
		for(i=0;i<N.nodes;i++)
		{
			initm(nodep->blockDiagKp,n*n);
			nodep = nodep->nextnp;
		}

		resL2 = 0.0;
		resLinf = 00;
		resmaxi = 0;
        for(i=0;i<eqnsets[nset].neqns;i++)
		{
                res = *(eqnsets[nset].Fp + i);
				resL2 += res*res;
				if(fabs(res) > fabs(resLinf)){
					resLinf = res;
					resmaxi = i;
				}
		}

		resL2 = sqrt(resL2);

		return (0);
}
int             test_elK_file(nelm,nvar,fp) // used to degug
int             nelm, nvar ;				// prints the element
FILE* fp;									// matrices to file


{
        int             i,j, ii, jj,err, ntf, n, itnum=0 ;
        struct element  *elp, tel ;
		int				ep_init(), get_PGKe();

        fprintf(fp,"\n Function test_elK \n");
        set_ptou(0) ;  
		updateVelocities();
        if(nvar < 0)
                n = N.vars ;
        else
                n = 1 ;
        if( (ep.n = ep_init() ) == 0 )
                return(-1) ;
        elp = gp.El ;
        for(i=0;i<N.elms;i++) {
                if(elp->n == nelm)
                        break ;
                elp = elp->nextelp ;
        }fprintf(fp," before getPGKe\n");
        /* get_PGKe returns nbf: # of basis functions */
        err = get_PGKe(elp,&tel,-1,2,&ntf,0) ;
        fprintf(fp,"\n  K matrix and F vector for element %d\n",nelm) ;
        for(j=0;j<err;j++) {
                        fprintf(fp," %d",tel.nps[j]->i) ;
        }
        for(i=0;i<ntf;i++) {
                for(ii=0;ii<n;ii++) { 
                        fprintf(fp,"\n %d ",tel.nps[i]->i) ;
                        for(j=0;j<err;j++) {
                                for(jj=0;jj<n;jj++)
                                        fprintf(fp," %12.4g",KE(i*n+ii,j*n+jj)) ;
                        }

                }
        }
        if(N.trans > 0) {
                fprintf(fp,"\n  S matrix  for element %d\n",nelm) ;
                for(i=0;i<ntf;i++) { 
                        for(ii=0;ii<n;ii++) { 
                                fprintf(fp,"\n ") ;
                                for(j=0;j<err;j++) {
                                        for(jj=0;jj<n;jj++)
                                                fprintf(fp," %11.3g",SE(i*n+ii,j*n+jj)) ;
                                }
                        }
                }
        }
        printf("\n\n") ;
        return(err) ;
}

int JacobiPC(Min, itnum, varnum, elcode)
double *Min;
int itnum, varnum, elcode; 
{
	int i, k, ii, outVectrow;
	struct element          *elp;
    struct belement         *belp;
	int						ntf, ns;
	struct node				*outVectinp;
	int						get_PGKeAnalJ(), get_bKeAnalJ();
	int						get_PGKeNumJ(), get_bKeNumJ();


	//note: does not reinitialize Min before assembling
	//but this should not be a problem unless Min is being
	//updated within a call to the gmres solver - Min will
	//be reinitialized at the beginning of each time step

	elp = gp.El ;
        
	for(k=0;k<N.elms;k++)
	{
		if(jacobianType == 0)
			ns = get_PGKeAnalJ(elp,elp,varnum,elcode,&ntf,itnum) ;
		else
 			ns = get_PGKeNumJ(elp,elp,varnum,elcode,&ntf,itnum) ;
		for(i = 0; i < ntf; i ++)
		{
			outVectinp = elp->nps[i];
			outVectrow = outVectinp->i * N.vars;
			for (ii = 0; ii < N.vars; ii++)
			{					
				Min[outVectrow + ii] += SE(i*N.vars+ii,i*N.vars+ii);
			}
		}
        elp = elp->nextelp ;
    }

    belp = gp.B ;

    for(k=0;k<N.belms;k++)
	{ 
		if(jacobianType == 0)
			ns = get_bKeAnalJ(belp,belp,varnum,elcode,&ntf,itnum) ;
		else
			ns = get_bKeNumJ(belp,belp,varnum,elcode,&ntf,itnum) ;
		for(i = 0; i < ntf; i ++)
		{
			outVectinp = belp->nps[i];
			outVectrow = outVectinp->i * N.vars;
			for (ii = 0; ii < N.vars; ii++)
			{
				Min[outVectrow + ii] += SE(i*N.vars+ii,i*N.vars+ii);
			}
		}
        belp = belp->nextbelp ;
    }

	for( k = 0; k < Nukns; k++) Min[k] = 1.0 / Min[k];

	return(0);
}



