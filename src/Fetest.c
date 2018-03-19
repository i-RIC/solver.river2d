//#pragma options cpluscmt        
#include "Fe_PS.h"

#define KE(A,B) ( *(ep.K + ep.n*(A) + (B)) )
#define SE(A,B) ( *(ep.S + ep.n*(A) + (B)) )
#define ME(A,B) ( *(ep.M + ep.n*(A) + (B)) )
#define FE(B) ( *(ep.F + B) )
#define TRULY -1

extern  struct  control N ;
extern  struct  pointers        gp;
extern  struct  RegMesh Mesh ;
extern  struct  fxnodes *first_fnp ;
extern  int             Nukns ;
extern  struct  elmpointers     ep ;
extern  struct  transient       tvals ;
extern  struct  settings        sets;
extern  struct  eqnset  eqnsets[4];
void	updateVelocities();


int             get_L2norm(nvar)
int             nvar ;

{
        int     i ;
        struct element  *anElp ;
        struct node     *anp ;
        double  L2 = 0.0, get_QUL2(), maxerr, err, uexact() ;
        
        anElp = gp.El ;
        for(i=0;i<N.elms;i++) {
                L2 += get_QUL2(anElp,nvar) ;
                anElp = anElp->nextelp ;
        }
        printf(" L2 norm is %20.14g\n",sqrt(L2)) ;
        anp = gp.N ;
        maxerr = 0.0 ;
        for(i=0;i<N.nodes;i++) {
                err = fabs(uexact(nvar,anp->x[0],0.0,0.0,anp->p[0]) - anp->u[nvar]) ;
                if(err > maxerr)
                        maxerr = err ;
                anp = anp->nextnp ;
        }
        printf(" inf norm is %20.14g\n",maxerr) ;
        return(0) ;
}

int             list_vars(nvar)
int             nvar ;

{
        int     i;
        double  uexact() ;
        
        printf(" Node\t Value \n\n");
        for(i=0;i<N.nodes;i++) {
                printf(" %d\t %20.14g  %20.14g\n",gp.iptrs[i]->n,gp.iptrs[i]->u[nvar],
                  uexact(nvar,gp.iptrs[i]->x[0],0.0,0.0,gp.iptrs[i]->p[2]));
        }
        return ;
}               

int             set_param(nparam)
int             nparam ;

{
        int     i;
        double  p ;
        
        printf(" Input value ");
        scanf(" %lf",&p) ;
        for(i=0;i<N.nodes;i++) {
                gp.iptrs[i]->p[nparam] = p;
        }
        printf(" \n") ;
        return ;
}
                
int             set_ptou(nparam)
int             nparam ;

{
        /* function to set old values (uo, uoo) */
        int     i;
        double  p ;
        
        for(i=0;i<N.nodes;i++) {
                update_p(gp.iptrs[i]);
        }
        return ;
}
                
int             reset_u(nparam)
int             nparam ;

{
        /* function to reset new values to old values (uo, uoo) */
        int     i,j;
        double  p ;
        
        for(i=0;i<N.nodes;i++) {
                for(j=0;j<N.vars;j++) {
                        gp.iptrs[i]->u[j] = gp.iptrs[i]->uo[j];
                }
        }
        return ;
}
                
int             test_elK(nelm,nvar)
int             nelm ;

{
        int             i,j, ii, jj,err, ntf, n, itnum=0 ;
        struct element  *elp, tel ;
                printf("\n Function test_elK \n");
        
        set_ptou(0) ;   /* Faye 89-02-07 */
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
        }printf(" before getPGKe\n");
        /* get_PGKe returns nbf: # of basis functions */
        err = get_PGKe(elp,&tel,-1,2,&ntf,itnum) ;
        printf("\n  K matrix and F vector for element %d\n",nelm) ;
        for(j=0;j<err;j++) {
                        printf(" %d",tel.nps[j]->i) ;
        }
        for(i=0;i<ntf;i++) {
                for(ii=0;ii<n;ii++) { 
                        printf("\n %d ",tel.nps[i]->i) ;
                        for(j=0;j<err;j++) {
                                for(jj=0;jj<n;jj++)
                                        printf(" %12.4g",KE(i*n+ii,j*n+jj)) ;
                        }
/*                      printf(" %11.3g",FE(i*n+ii)) ;*/
                }
        }
        if(N.trans > 0) {
                printf("\n  S matrix  for element %d\n",nelm) ;
                for(i=0;i<ntf;i++) { 
                        for(ii=0;ii<n;ii++) { 
                                printf("\n ") ;
                                for(j=0;j<err;j++) {
                                        for(jj=0;jj<n;jj++)
                                                printf(" %11.3g",SE(i*n+ii,j*n+jj)) ;
                                }
                        }
                }
        }
        printf("\n\n") ;
        return(err) ;
}

int             test_globalK(numit,nvar)
int             numit,nvar;

{
        int             i,j, ii, jj,err, ntf, n ;
        int             varnum=-1,code=0,nset=0,elcode=2, itnum=0;      struct element  *elp, tel ;
        FILE    *fptr, *put_fptr() ;

        tvals.dt = 1.0 ;        
        tvals.theta = 0.5 ;

        printf("\n Function test_globalK \n");
        printf("theta=%g\tdt=%g\n",tvals.theta,tvals.dt);
        set_ptou(0) ;   /* Faye 89-02-07 */
 		updateVelocities();
       
        assemble(varnum,code,nset,elcode,itnum);
        if(numit !=0) {
                for (i=0;i<numit;i++)
                solve(varnum,2,0);
                assemble(varnum,code,nset,elcode,i+1);
        }
        if((fptr = put_fptr(1) ) != NULL){ 
/*              for(i=0; i<*(eqnsets[nset].diagp + Nukns-1) + 1; i++)
                        fprintf(fptr,"\t%d\t\t%g\n",i,*(eqnsets[nset].Kp+i));
*/
                fprintf(fptr,"Upper Matrix\n\n");
                fprintf(fptr,"\t%g\n",*(eqnsets[nset].Kp));
                fprintf(fptr,"\n");
                fprintf(fptr,"-1\n");
                for(j=1;j<Nukns;j++){
                        ii = *(eqnsets[nset].diagp+j) - *(eqnsets[nset].diagp+j-1);
                        for(jj=0;jj<ii;jj++)                                                                    fprintf(fptr,"\t%g\n",*(eqnsets[nset].Kp+*(eqnsets[nset].diagp+j-1) +1+jj));
                        fprintf(fptr,"\n");
                fprintf(fptr,"-1\n");

                }
                fprintf(fptr,"Lower Matrix\n\n");
                fprintf(fptr,"\t%g\n",*(eqnsets[nset].Lp));
                fprintf(fptr,"\n");
//              fprintf(fptr,"-1\n");

                for(j=0;j<Nukns;j++){
                        ii = *(eqnsets[nset].diagp+j) - *(eqnsets[nset].diagp+j-1);
                        for(jj=0;jj<ii;jj++)                                                                    fprintf(fptr,"\t%g\t",*(eqnsets[nset].Lp+*(eqnsets[nset].diagp+j-1) +1+jj));
                        fprintf(fptr,"\n");
                fprintf(fptr,"-1\n");

                }
                fprintf(fptr,"F vector\n\n");
                for(j=0;j<Nukns;j++)                                                    fprintf(fptr,"\t%g\n",*(eqnsets[nset].Fp+j));
                


        }
        fclose(fptr);           

        return(0) ;

}

int             test_globalK1(numit,nvar)
int             numit,nvar;
/* this function prints the global matrix in a format, suitable for the function switch_2D */
{
        int             i,j, ii, jj,err, ntf, n ;
        int             varnum=-1,code=0,nset=0,elcode=2, itnum=0;      struct element  *elp, tel ;
        FILE    *fptr, *put_fptr() ;

        tvals.dt = 1.0 ;        
        tvals.theta = 0.5 ;

        printf("\n Function test_globalK1 \n");
        printf("theta=%g\tdt=%g\n",tvals.theta,tvals.dt);
        set_ptou(0) ;   /* Faye 89-02-07 */
 		updateVelocities();
       
        assemble(varnum,code,nset,elcode,itnum);
        if(numit !=0) {
                for (i=0;i<numit;i++)
                solve(varnum,2,0);
                assemble(varnum,code,nset,elcode,i+1);
        }
        if((fptr = put_fptr(1) ) != NULL){ 
                fprintf(fptr,"Upper Matrix\n\n");
                fprintf(fptr,"\t%g\n",*(eqnsets[nset].Kp));
                fprintf(fptr,"9999.9999\n");
                for(j=1;j<Nukns;j++){
                        ii = *(eqnsets[nset].diagp+j) - *(eqnsets[nset].diagp+j-1);
                        for(jj=0;jj<ii;jj++)                                                                    fprintf(fptr,"\t%g\n",*(eqnsets[nset].Kp+*(eqnsets[nset].diagp+j-1) +1+jj));
                        fprintf(fptr,"9999.9999\n");
                }
                fprintf(fptr,"Lower Matrix\n\n");
                fprintf(fptr,"\t%g\n",*(eqnsets[nset].Lp));
                fprintf(fptr,"9999.9999\n");
                for(j=1;j<Nukns;j++){
                        ii = *(eqnsets[nset].diagp+j) - *(eqnsets[nset].diagp+j-1);
                        for(jj=0;jj<ii;jj++)                                                                    fprintf(fptr,"\t%g\n",*(eqnsets[nset].Lp+*(eqnsets[nset].diagp+j-1) +1+jj));
                        fprintf(fptr,"9999.9999\n");
                }
                fprintf(fptr,"F vector\n\n");
                for(j=0;j<Nukns;j++)                                                    fprintf(fptr,"\t%g\n",*(eqnsets[nset].Fp+j));
                


        }
        fclose(fptr);           

        return(0) ;

}

/*
int             switch_2D(opt,nnode)
int             opt, nnode;

{
        int             i,j,k,l,m;
        int             ndiag,icount,Nukns1d,NVARS1d,Nnodes1d,istart;
        double  Col[Nukns], K2D[Nukns][Nukns];
        double  K2[N.nodes][N.nodes][NVARS][NVARS], K1[N.nodes/2][N.nodes/2][2][2];
        double  K1D[N.nodes/2*(NVARS-1)][N.nodes/2*(NVARS-1)];
        double  F2[Nukns],F1D[N.nodes];
        FILE    *ft,*fi,*fo,*get_fptr(),*put_fptr();
        double  get_dbl();
        
//  initialize all matrices
        NVARS1d = 2;
        Nnodes1d =N.nodes/2;
        Nukns1d = Nnodes1d * 2;
        
        for(i=0;i<Nukns;i++){
                for(j=0;j<Nukns;j++) 
                        K2D[i][j] = 0.0;
        }
        for(i=0;i<N.nodes;i++){
                for(j=0;j<N.nodes;j++){
                        for(m=0;m<NVARS;m++) {
                                for(l=0;l<NVARS;l++)
                                        K2[i][j][m][l] = 0.0;
                        }
                }
        }
        for(i=0;i<Nukns1d;i++){
                for(j=0;j<Nukns1d;j++) 
                        K1D[i][j] = 0.0;
        }
        for(i=0;i<Nnodes1d;i++){
                for(j=0;j<Nnodes1d;j++){
                        for(m=0;m<NVARS1d;m++) {
                                for(l=0;l<NVARS1d;l++)
                                        K1[i][j][m][l] = 0.0;
                        }
                }
        }
        for(i=0;i<Nukns;i++)
                F2[i] = 0.0;
        for(i=0;i<Nukns1d;i++)
                F1D[i] = 0.0;
        printf("\n Function to switch a 2D global matrix to a 1D global matrix \n");
        set_ptou(0);
		updateVelocities();
// ft=put_fptr(1);
        if((fi = get_fptr(1)) != NULL) {
                ndiag =  Nukns;
//              printf("ndiag=%d\n",ndiag);
                for(j=0;j<ndiag;j++){
                        for(i=0;i<ndiag;i++)
                                Col[i] = 0.0;
                        i=icount=0;
                        while ((Col[i] = get_dbl(fi)) != 9999.9999) {
                                icount++;
                                i++;
                        }
                        for(i=0;i<(icount);i++){
                                istart=j-icount+1;
                                K2D[i+istart][j] = Col[i];
                        }
//fprintf(ft,"%g\t",Col[i]);
                                
                }

                for(j=0;j<ndiag;j++){
                        for(i=0;i<ndiag;i++)
                                Col[i] = 0.0;
                        i=icount=0;
                        while((Col[i] = get_dbl(fi)) != 9999.9999) {
                                 i++;
                                icount++;
                        }
                        for(i=0;i<(icount);i++)
                                K2D[j][j-(icount-i)+1] = Col[i];
                }
                for(j=0;j<ndiag;j++){
                        F2[j] = get_dbl(fi);
//                      printf("%g\t",F2[j]);
                }
                for(i=0;i<ndiag;i++){
                        for(j=0;j<ndiag;j++){
                                K2[i/3][j/3][i%3][j%3] = K2D[i][j];
//fprintf(ft,"%g\t",K2[i/3][j/3][i%3][j%3]);
                        }
//              fprintf(ft,"\n");
                }
//o.k.
//  Build 1d matrix from K2 (assume a one-element wide 2d domain)
                for(i=0;i<Nnodes1d;i++){
                        for(j=0;j<Nnodes1d;j++){
                                for(k=0;k<NVARS1d;k++){
                                        for(m=0;m<NVARS1d;m++){
                                                K1[i][j][k][m] = K2[2*i][2*j][k][m] + K2[2*i+1][2*j][k][m] + K2[2*i][2*j+1][k][m] + K2[2*i+1][2*j+1][k][m];
                                            }
                                }
                        }
                }
        }
        for(i=0;i<Nukns1d;i++) {
                for(j=0;j<Nukns1d;j++){
                        K1D[i][j] = K1[i/2][j/2][i%2][j%2];
//  fprintf(ft,"%g\t",K1D[i][j]);
                }
// fprintf(ft,"\n");
        }
        for(i=0;i<Nukns1d;i++){
                F1D[i] = F2[3*i-2*(i%2)] + F2[3*(i+1)-2*(i%2)];
//              printf("%g\n",F1D[i]);
        }
        if((fo=put_fptr(1)) !=NULL) {
                fprintf(fo,"Upper Matrix\n\n");
                for(i=0;i<Nukns1d;i++) {
                        for(j=0;j<=i;j++)
                                fprintf(fo,"\t%g\n",K1D[j][i]);
                        fprintf(fo,"\n");                       
                }
                fprintf(fo,"Lower Matrix\n\n");
                for(i=0;i<Nukns1d;i++) {
                        for(j=0;j<=i;j++)
                                fprintf(fo,"\t%g\t",K1D[i][j]);
                        fprintf(fo,"\n");
                }
                fprintf(fo,"F-vector\n\n");
                for(i=0;i<Nukns1d;i++)
                        fprintf(fo,"\t%g\n",F1D[i]);
        }

        fclose(fo);
        return(0);
}
*/

int     matrixform(unknowns)
int     unknowns;
{
        FILE    *f, *get_fptr();
        FILE    *fptr, *put_fptr() ;

        double  matrix[12][12],temp[12],get_dbl();
        int             i,j,height,fill;
        double  num;
        
        f=get_fptr(1);
        for(i=0;i<unknowns;i++) {
                temp[i]=0.0;
                for(j=0;j<unknowns;j++)
                        matrix[i][j]=0.0;
        }
        for(i=0;i<unknowns;i++) {
                for(j=0;j<unknowns;j++)
                        temp[j]=0.0;
                j=0;
                height=0;
                while((num=get_dbl(f))!=-1) {
                        temp[j]=num;
                        j++;
                        height++;
                }
                fill=i+1-height;
                for(j=0;j<(i+1-fill);j++)
                        matrix[j+fill][i]=temp[j];
        }

        for(i=0;i<unknowns;i++) {
//printf("i=%d\n",i);
                for(j=0;j<unknowns;j++)
                        temp[j]=0.0;
                j=0;
                height=0;
                while((num=get_dbl(f))!=-1) {
                        temp[j]=num;
                        j++;
                        height++;
//printf("j=%d\n",j);
                }
                fill=i+1-height;
                for(j=0;j<(i+1-fill);j++)
                        matrix[i][j+fill]=temp[j];
        }

        
        if((fptr = put_fptr(1) ) != NULL){ 
                fprintf(fptr,"{");
                for(i=0;i<unknowns;i++) {
                        fprintf(fptr,"{");
                        for(j=0;j<unknowns;j++){
                                if(fabs(matrix[i][j])<0.00001)
                                        matrix[i][j]=0.0;
                                if(j!=(unknowns-1))
                                        fprintf(fptr,"%g,",matrix[i][j]);
                                else
                                        fprintf(fptr,"%g",matrix[i][j]);
                        }
                        if(i!=(unknowns-1))
                                fprintf(fptr,"},\n");
                        else
                                fprintf(fptr,"}\n");
                }
                fprintf(fptr,"}");
                

        }
        fclose(fptr);           
        return(0);

}
                        
                


        
int             test_belK(nelm,nvar)
int             nelm ;

{
        int             i,j, ii, jj,err, ntf, n, itnum = 0 ;
        struct belement *elp, tel ;
        
        if(nvar < 0)
                n = N.vars ;
        else
                n = 1 ;
        if( (ep.n = ep_init() ) == 0 )
                return(-1) ;
        elp = gp.B ;
        for(i=0;i<N.belms;i++) {
                if(elp->n == nelm)
                        break ;
                elp = elp->nextbelp ;
        }
        err = get_bKe(elp,&tel,nvar,2,&ntf,itnum) ;
        printf("\n  K matrix and F vector for element %d type %d\n",nelm,tel.vtype) ;
        for(j=0;j<err;j++) {
                        printf("\t %d",tel.nps[j]->i) ;
        }
        for(i=0;i<ntf;i++) {
                for(ii=0;ii<n;ii++) { 
                        printf("\n %d ",tel.nps[i]->i) ;
                        for(j=0;j<err;j++) {
                                for(jj=0;jj<n;jj++)
                                        printf(" %9.4g",KE(i*n+ii,j*n+jj)) ;
                        }
                        printf(" %12.4g",FE(i*n+ii)) ;
                }
        }
        printf("\n\n") ;
        printf("\n  S matrix for element %d type %d\n",nelm,tel.vtype) ;
        for(j=0;j<err;j++) {
                        printf("\t %d",tel.nps[j]->i) ;
        }
        for(i=0;i<ntf;i++) {
                for(ii=0;ii<n;ii++) { 
                        printf("\n %d ",tel.nps[i]->i) ;
                        for(j=0;j<err;j++) {
                                for(jj=0;jj<n;jj++)
                                        printf(" %9.4g",SE(i*n+ii,j*n+jj)) ;
                        }
                }
        }
        printf("\n\n") ;
        return(err) ;
}

int             test_contour(code)
int             code ;

{
        int     nelm,niv,nvar, ni, i, ntf, dir ;
        double  tol, xl[10], yl[10], cval ;
        struct element  *elp, tel ;
        
//      printf(" Input el#, #iv, #var, value, tol \n") ;
//      scanf(" %d %d %d %lf %lf",&nelm,&niv,&nvar,&cval,&tol) ;
printf("input xl[0]\n");
scanf("%lf",xl);
printf("xl[0]=%g\n",xl[0]);
printf("input xl[1]\n");
scanf("%lf",xl);
printf("xl[1]=%g\n",xl[1]);

return(0);
        elp = gp.El ;
        for(i=0;i<N.elms;i++) {
                if(elp->n == nelm)
                        break ;
                elp = elp->nextelp ;
        }
        ElPick(elp,&tel,&ntf,&dir) ;
        ni = get_cline(&tel,cval,niv,nvar,tol,xl,yl,1) ;
        printf(" number of intersections = %d\n",ni) ;
        if(ni == 2) {
                printf(" contour line\n") ;
                for(i=0;i<niv+2;i++)
                        printf(" %f  %f \n",xl[i],yl[i]) ;
        }
        return(0) ;
}
                
int             test_2Delm(el_type)
int             el_type ;

{
        int     i ;
        struct gausspts         x ;
        struct shapefuncs       sf ;    
        
        printf(" Input local coordinates r, s \n") ;
        scanf(" %lf %lf",&x.x,&x.y) ;
        x.z = 0.0 ;
        x.w = 1.0 ;
        printf("\n   r = %10.5f  s = %10.5f\n",x.x,x.y) ;
        if( get_shape(el_type,&x,&sf) == 0 ){
          printf("\n Function number       f          dfdr          dfds\n\n");
          for(i=0;i<sf.dof;i++){
            printf(" %20d    %10.5f    %10.5f     %10.5f\n",i+1,sf.f[i],sf.dfdr[i],sf.dfdr[i]);
          }
        }
        return(0) ;
}

int             TestElErr(elnum)
int             elnum ;

{
        int     i ;
        struct element          *elp ; 
        double  GetElErr(), err ;
        
        elp = gp.El ;
        for(i=0;i<N.elms;i++) {
                if(elp->n == elnum)
                        break ;
                elp = elp->nextelp ;
        }
        err = GetElErr(elp,0) ;
        printf("  error = %f  \n",err) ;
        return(0) ;
}
                
int             test_trans(nvar)
int             nvar ;

{
        int     i, nt ;
        double  tol ;
        
        printf(" Time = %lf \n", tvals.t) ;
        printf(" Input theta, delta t and magnification factor \n") ;
        scanf(" %lf %lf %lf",&tvals.theta,&tvals.dt,&tvals.dtfac) ;
        printf(" Input # of time steps, tolerance and # of steps between printouts\n") ;
        scanf(" %d %lf %d",&tvals.nsteps,&tol,&tvals.iter) ;
//      
//      
        switch(sets.plotcode){          
                case 0 :
                        printf(" Input coordinate # & coordinate value \n") ;
                        scanf(" %d %lf",&tvals.crdnum,&tvals.crdval) ;
                        break ;
                case 1:
                        printf("Input variable number, min and max coordinate values and number of contours\n");
                        scanf("%d %lf %lf %d",&sets.plotn[0],&sets.plotv[0],&sets.plotv[1],&sets.plotn[1]);
                        break;
                case 2:
//                      printf("Input minimum Depth for Velocity Calculation\n");
//                      scanf("%lf",&sets.plotv[0]);
                        sets.plotv[0]=0.1;
                        break;
                case 3 :
                        printf("Input xo and yo, gama and beta, and ndx and ndy\n");
                        scanf("%lf %lf %lf %lf %d %d",&sets.plotv[0],&sets.plotv[1],&sets.plotv[2],&sets.plotv[3],&sets.plotn[0],&sets.plotn[1]);
                        break ;
        }

        transient(nvar,2,tol) ;
        return(0) ;
}
                
int             test_steady(nvar)
int             nvar ;

{
        int     i, nt ;
        double  tol, dtmax ;
        
        printf(" Time = %f \n", tvals.t) ;
        tvals.theta = 1.0;
        tvals.dtfac = 1.0;
        tol = 0.0001;
        tvals.iter = 1;
        if(nvar == 1){
                printf(" Input first delta t, max delta t, goal uchange, and model final time \n") ;
                scanf(" %lf %lf %lf %lf",&tvals.dt,&dtmax,&tvals.uc,&tvals.tfinal) ;
        }
        else{
                printf(" Input max delta t and model final time \n") ;
                scanf(" %lf %lf",&dtmax,&tvals.tfinal) ;
        }
 
        
//      
//      
        switch(sets.plotcode){          
                case 0 :
                        printf(" Input coordinate # & coordinate value \n") ;
                        scanf(" %d %lf",&tvals.crdnum,&tvals.crdval) ;
                        break ;
                case 1:
                        printf("Input variable number, min and max coordinate values and number of contours\n");
                        scanf("%d %lf %lf %d",&sets.plotn[0],&sets.plotv[0],&sets.plotv[1],&sets.plotn[1]);
                        break;
                case 2:
//                      printf("Input minimum Depth for Velocity Calculation\n");
//                      scanf("%lf",&sets.plotv[0]);
                        sets.plotv[0]=0.1;
                        break;
                case 3 :
                        printf("Input xo and yo, gama and beta, and ndx and ndy\n");
                        scanf("%lf %lf %lf %lf %d %d",&sets.plotv[0],&sets.plotv[1],&sets.plotv[2],&sets.plotv[3],&sets.plotn[0],&sets.plotn[1]);
                        break ;
        }

        steady(-1,2,tol,dtmax) ;
        return(0) ;
}

int             test_time(nvar)
int             nvar ;

{
        
        tvals.t = nvar ;
        printf(" Time = %lf \n", tvals.t) ;
        return(0) ;
}

int             test_smooth(nvar)
int             nvar ;

{
        struct node             *np0, *npp, *npm ;
        int i, j;

        for(j=0;j<3;j++){
                npm = gp.N ;
                np0 = npm -> nextnp ;
                npp = np0 -> nextnp ;
        
                for (i=1;i<N.nodes-1;i++){ 
                        np0->x[nvar] = np0->x[nvar] + (npm->x[nvar] - 2*np0->x[nvar] + npp->x[nvar])/2 ;
                        npm = np0 ;
                        np0 = npp ;
                        npp = npp -> nextnp ;
                }                       
        }
        return(0) ;
}

double          test_consL2(elmntp,ident,etype)
int     ident, etype;
struct  element *elmntp;
{
        struct element theEl,*theElp;
        struct gausspts g[MDOF];
        struct shapefuncs bf;
        struct node *np[MNSF];

        int i,j,nbf,ni,nm,ntf,dir;
        double J[NDIMS][NDIMS],detJ,ua,x,quantity,get_J(),Qx,Qy;
        theElp = &theEl;
        ElPick(elmntp,theElp,&ntf,&dir);
        nbf = theElp->nnds;

        ni = get_gspts(etype,g);

        for(j = 0;j < theElp->nnds;j++)
                np[j] = theElp->nps[j];

        quantity = 0.0;

        for(i=0;i<ni;i++) {
                get_shape(theElp->vtype,&g[i],&bf) ;
                detJ = get_J(np,&bf,J) ;
                ua = 0.0;

                for(j=0;j<nbf;j++) 
                        ua += bf.f[j]*np[j]->u[ident];

                quantity += ua * fabs(detJ) * g[i].w ;
        }
        return(quantity);
}

double          test_cons(etype)
int             etype ;

{
        struct node             *np0, *npp ;
        struct element          *elp ;
        double vol , momx, momy, test_consL2();
        int i, j;

        momx = momy = 0.0 ;
        vol = 0.0 ;
        elp = gp.El ;
        for(i=0;i<(N.elms);i++){
                vol += test_consL2(elp,0,etype);  
                momx += test_consL2(elp,1,etype);
                momy += test_consL2(elp,2,etype);

                elp = elp -> nextelp ;
        }
        printf(" Total volume      = %g \n Total x-momentum = %g \n Total y-momentum = %g\n",vol,momx,momy) ;
        return(0) ;
}

int             test_mesh(eltype)
int     eltype ;

{
        int     i, ni ;
        struct element  *elp ;
        
        Mesh.eltype = eltype ;
        printf(" Input nx, nbx, ny, nby \n") ;
        scanf(" %d %d %d %d",&Mesh.nx,&Mesh.nbx,&Mesh.ny,&Mesh.nby) ;
        ni = MkLapMesh(&Mesh) ;/*in fecore.c*/

        printf(" return code = %d\n",ni) ;
        return(0) ;
}


int             test_map(maptype)
int             maptype ;

{
        int     i, n, k, l, numrows, numcols ;
        struct  blockmap        *theblock ;
        struct  node            *np;    
        Mesh.maptype = maptype ;
        if(Mesh.firstmap != NULL)
                free(Mesh.firstmap) ;
        theblock = Mesh.firstmap ;
        for(k=0;k<Mesh.nby;k++) {
                for(l=0;l<Mesh.nbx;l++) {
                        printf(" Input macro element nodes x,y coord. pairs\n") ;
                        n = nsf(maptype) ;
                        for(i=0;i<n;i++) {
                                scanf(" %f %f",&(theblock->xnodes[i]),&(theblock->ynodes[i])) ;
                        }
                        map_mesh(theblock,Mesh.maptype,k,l) ;
                        theblock++ ;
                }
        }
        np = gp.N;
        for(i=0;i<N.nodes;i++) {
                nodevalues(np);
                np = np->nextnp;
        }




        return(0) ;
}

double	test_outflow(code)	// netflow (code =0) inflow (code = 1), outflow (code = 2)
int		code;
	
{
        struct belement 	*belp ;
        struct gausspts  	g[MGPTS] ;
        struct shapefuncs 	fv;
        struct node       	*np[MNSF] ;
        double  discharge = 0.0, dx, dy, qx, qy;
        int i, j, nbf, ni, ig, ib;
 
		switch(code){
		case 0 :	//netflow

			belp = gp.B ;
			for(ib=0;ib<N.belms;ib++){
				if( belp->bcs[0] > 2){ //(belp->bcs[0] == 3) || (belp->bcs[0] == 5)) {
           			nbf =  nsf(belp->vtype);
					for(j=0;j<nbf;j++)
                		np[j] = belp->nps[j] ;
           			ni = get_gspts(belp->vtype,g) ;
					for(ig=0;ig<ni;ig++) {
						get_shape(belp->vtype,&g[ig],&fv) ;
						dx = dy = qx = qy = 0.0 ;
						for(i=0;i<fv.dof;i++) {       
							dx += fv.dfdr[i] * np[i]->x[0] ;
							dy += fv.dfdr[i] * np[i]->x[1] ;
							qx += fv.f[i] * np[i]->u[1];
							qy += fv.f[i] * np[i]->u[2];
						}
						discharge += (qx*dy - qy*dx) * g[ig].w;
					}
				}
				else{ //belp->bcs[0] == 1
					discharge -= belp->p[1] * sqrt((belp->nps[1]->x[0]-belp->nps[0]->x[0])*
												   (belp->nps[1]->x[0]-belp->nps[0]->x[0])
												  +(belp->nps[1]->x[1]-belp->nps[0]->x[1])*
												   (belp->nps[1]->x[1]-belp->nps[0]->x[1]));
					}
				belp = belp->nextbelp ;
			}
			return discharge;
		
		case 1:	//inflow

			belp = gp.B ;
			for(ib=0;ib<N.belms;ib++){
				if( belp->bcs[0] == 1){
           			discharge += belp->p[1] * sqrt((belp->nps[1]->x[0]-belp->nps[0]->x[0])*
												   (belp->nps[1]->x[0]-belp->nps[0]->x[0])
												  +(belp->nps[1]->x[1]-belp->nps[0]->x[1])*
												   (belp->nps[1]->x[1]-belp->nps[0]->x[1]));
				}
				belp = belp->nextbelp ;
			}
			return discharge;
		
		case 2:	//outflow

			belp = gp.B ;
			for(ib=0;ib<N.belms;ib++){
				if( belp->bcs[0] > 2){ //(belp->bcs[0] == 3) || (belp->bcs[0] == 5)) {
           			nbf =  nsf(belp->vtype);
					for(j=0;j<nbf;j++)
                		np[j] = belp->nps[j] ;
           			ni = get_gspts(belp->vtype,g) ;
					for(ig=0;ig<ni;ig++) {
						get_shape(belp->vtype,&g[ig],&fv) ;
						dx = dy = qx = qy = 0.0 ;
						for(i=0;i<fv.dof;i++) {       
							dx += fv.dfdr[i] * np[i]->x[0] ;
							dy += fv.dfdr[i] * np[i]->x[1] ;
							qx += fv.f[i] * np[i]->u[1];
							qy += fv.f[i] * np[i]->u[2];
						}
						discharge += (qx*dy - qy*dx) * g[ig].w;
					}
				}
				belp = belp->nextbelp ;
			}
			return discharge;
		}
		return 0.0;
}


double	test_volume()
	
{
        struct shapefuncs 	fv;
        struct node       	*np[MNSF] ;
		struct element		*elp;
        double  areatr(), volume = 0.0, area, depth;
        int i, j, nbf, ni, ig, ib;
 
		elp = gp.El;
		while(elp != NULL){

			area = areatr(elp->nps[0]->x[0],elp->nps[0]->x[1],
				elp->nps[1]->x[0],elp->nps[1]->x[1],
				elp->nps[2]->x[0],elp->nps[2]->x[1]);
			depth = (elp->nps[0]->u[0] + elp->nps[1]->u[0] + elp->nps[1]->u[0])/3.0;
			if(depth > 0.0)
				volume += area*depth;
			elp = elp->nextelp;
		}
		return volume;
}

int             test_bvalues(code)
int             code;

{
        struct belement *belp ;
        double  discharge, elevation=0.0, dx, dy, qx, qy;
        int i, bcode=0, nsegs=0, segnum = 0;

        belp = gp.B ;
        for(i=0;i<N.belms;i++) {
                if( (belp->bcs[0] == 3) || (belp->bcs[0] == 5)) {
                        if( (bcode != 3) && (bcode != 5) ){
                                elevation = 0.0;
                                nsegs = 0;
                                discharge = 0.0;
                                segnum += 1;
                        }
                        dx = belp->nps[1]->x[0] - belp->nps[0]->x[0]; 
                        dy = belp->nps[1]->x[1] - belp->nps[0]->x[1];
                        qx = (belp->nps[0]->u[1] + belp->nps[1]->u[1]) / 2.0; 
                        qy = (belp->nps[0]->u[2] + belp->nps[1]->u[2]) / 2.0;
                        discharge += qx*dy - qy*dx;
                        elevation += (belp->nps[0]->u[0] + belp->nps[1]->u[0]) / 2.0; 
                        elevation += (belp->nps[0]->p[0] + belp->nps[1]->p[0]) / 2.0; 
                        nsegs += 1;
                }
                else {
                        if( (bcode == 3) || (bcode == 5)){
                                printf("\n Discharge on outflow boundary segment %d = %f\n",segnum, discharge);
                                elevation /= nsegs;
                                printf("\n Average elevation on outflow segment %d = %f\n",segnum,elevation);
                        }
                }
                if( belp->bcs[0] == 1) {
                        if( bcode != 1){
                                elevation = 0.0;
                                nsegs = 0;
                                segnum += 1;
                        }
                        elevation += (belp->nps[0]->u[0] + belp->nps[1]->u[0]) / 2.0; 
                        elevation += (belp->nps[0]->p[0] + belp->nps[1]->p[0]) / 2.0; 
                        nsegs += 1;
                }
                else {
                        if( bcode == 1){
                                elevation /= nsegs;
                                printf("\n Average elevation on inflow segemnt %d = %f\n",segnum,elevation);
                        }
                }
                bcode = belp->bcs[0]; 
                belp = belp->nextbelp ;
        }
        return(0);
}

int             set_bc(code)
int             code ;

{
        int     bestart, beend,i, j, bctype ;
        double  bvalue ;
        struct belement *belp ;
        
        while (TRULY) {
                printf(" Input first, last boundary elements, code and value\n") ;
                scanf(" %d %d %d %lf",&bestart,&beend,&bctype,&bvalue) ;
                if(bestart <= 0)
                        break ;
                belp = gp.B ;
                for(i=0;i<N.belms;i++) {
                        if(belp->n == bestart)
                                break ;
                        belp = belp->nextbelp ;
                }
                do {
                        if(bctype == 0) {
                                belp->bcs[0] = 0 ;
                                belp->bcs[1] = 1 ;
                                for(j=0;j<belp->nnds;j++) {
                                        (belp->nps[j])->u[1] = bvalue ;
                                }
                        }
                        else {
                                belp->bcs[0] = 1 ;
                                belp->bcs[1] = 0 ;
                                for(j=0;j<belp->nnds;j++) {
                                        (belp->nps[j])->u[0] = bvalue ;
                                }
                        }
                        if(belp->n == beend) 
                                break;
                        belp = belp->nextbelp ;
                }
                while(belp->n <= N.belms) ;
        }
        return(0) ;
}

int             test_UMesh(ndisc)
int             ndisc ;

{
        int     ni;

        ni = MkLapUMesh(ndisc) ;
        return(0) ;
}

int             test_mixedgps(nelement)
int             nelement ;

{
        int     i,ni;
        struct  element *elp;
        struct  gausspts        g[MGPTS];
        elp=gp.El;
        for(i=0;i<N.elms;i++) {
                if(elp->n==nelement) 
                        break;
                else
                        elp=elp->nextelp;
        }
        ni = get_mixedgps(elp,g) ;
for(i=0;i<ni;i++) 
printf("g.x=%g\tg.y=%g\tg.w=%g\n",g[i].x,g[i].y,g[i].w);
        return(1) ;
}
int             test_boundgps(nelement)
int             nelement ;

{
        int     i,ni;
        struct  belement        *elp;
        struct  gausspts        g[MGPTS];
        elp=gp.B;
        for(i=0;i<N.belms;i++) {
                if(elp->n==nelement) 
                        break;
                else
                        elp=elp->nextbelp;
        }
        if(check_bmixed(elp)==1)
                ni = get_boundarygps(elp,g) ;
        else {
                printf("no intersection\n");
                return(0);
        }
for(i=0;i<ni;i++) 
printf("g.x=%g\tg.y=%g\tg.w=%g\n",g[i].x,g[i].y,g[i].w);
        return(1) ;
}

