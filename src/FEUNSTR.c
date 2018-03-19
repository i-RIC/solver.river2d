//#pragma options cpluscmt        
#include                "Fe_PS.h"

extern  struct  control N;
extern  struct  pointers                gp;
struct  UnstrMesh       *theMesh;
extern  int             Nndso,Nelmso,Nbelmso;
struct  unode   *lastnodep;
struct  bsegment        *lastbsp;
extern double   defndp,defelp;
double  min,secondmin;
double  xc[20];


int     MkLapUMesh(num)
int             num;
{
        int     size;
        int     ndiscr,nmax;
        size = sizeof(struct UnstrMesh);
        if(theMesh!=NULL)
                free(theMesh);
        theMesh = (struct UnstrMesh *) malloc(size);
        if(num<0)       {
                printf(" Input number of discretizations in the y-direction\n");
                scanf("%d",&ndiscr);
        }
        else    
                ndiscr = num;
        printf(" Input number of closest neighbours to be considered\n");
        scanf("%d",&nmax);      
        N.trans = theMesh->con.trans = 1;
        N.dims  = theMesh->con.dims     = 2;
        N.vars  = theMesh->con.vars     = 3;
        N.Keqns[0] = theMesh->con.Keqns[0] = 1;
        N.Keqns[1] = theMesh->con.Keqns[1] = 1;
        N.params = theMesh->con.params = 3 ;
        N.bparams = theMesh->con.bparams = 7 ;
        theMesh->nd = ndiscr;
//change this in order to read the number of nodes per element 
// from the input file and set eltype (210 for 3 and 220 for 6 nodes)
        theMesh->eltype = 210;
        theMesh->nmax=nmax;
        defndp = 0.00 ;
        defelp = 0.00 ;
        if(num<0) 
                MakeUNodes(num);
        else
                MakeUMesh() ;
        return(0) ;
}


MakeUMesh()

{
        double  maxdim,mindim,find_min(),find_max();
        void    place_intnodes(),alloc_regnode(),make_elements(),set_previousunp();
        void    make_frontsegs(),make_tempbelements();
        void    smooth(),ordermesh();
        double  x,y,c;

        struct unode    *nodep;
        int i,j,ncuts,err,nsegs;
struct  orderseg        *osp;
struct  node    *np;
struct  element *elp;   
        
 //       input_bn();

        maxdim = find_max(1,-1);
        mindim = find_min(1,-1);
        theMesh->ymax = maxdim;
        theMesh->ymin = mindim;

        theMesh->dy = (maxdim-mindim)/theMesh->nd;

        make_tempbelements(theMesh->nd);
        make_frontsegs();
        N.bnodes=N.unodes;
        if(theMesh->nd!=0)
//                place_intnodes();
        nodep=gp.UN;
//        set_previousunp();
        make_elements(1);
        smooth(7);
        ordermesh();
//osp=gp.OR;
//for(i=0;i<N.ordersegs;i++){
//printf("order segment nodes are%d\t%d\n",osp->nps[0]->n,osp->nps[1]->n);
//osp=osp->nextosp;
//}

//np=gp.N;
//for(i=0;i<N.nodes;i++){
//printf("node %d\tx = %g\ty =%g\n",np->n,np->x[0],np->x[1]);
//np=np->nextnp;
//}     

elp=gp.El;
for(i=0;i<N.elms;i++) {
printf("Elm# =%d\tNodes =%d\t%d\t%d\n",elp->n,elp->nps[0]->n,elp->nps[1]->n,elp->nps[2]->n);
elp=elp->nextelp;
}
}

MakeUNodes(num)
int     num;

{
        double  maxdim,mindim,find_min(),find_max();
        void    place_intnodes(),set_previousunp();
        void    make_frontsegs(),make_tempbelements();

struct unode    *nodep;
        int i,j,ncuts,err,nsegs;
struct  tempbelement    *belp;
struct  frontseg                        *fsp;
//printf("0.5\n");
        
        input_bn(num);
//printf("0.75\n");

        maxdim = find_max(1,num);
        mindim = find_min(1,num);
        theMesh->ymax = maxdim;
        theMesh->ymin = mindim;
        
        if(theMesh->nd!=0)
                theMesh->dy = (maxdim-mindim)/theMesh->nd;
//printf("0.8\n");
        make_tempbelements(num);
//printf("0.9\n");
//belp=gp.TB;
//for(i=0;i<N.belms;i++) {
//printf("i=%d\tbelp =%d\tnextbelp =%d\n",i,belp,belp->nexttbp);
//belp=belp->nexttbp;
//}
//printf("1\n");
        make_frontsegs();
//printf("2\n");
        
        if(num==-1)
                N.bnodes=N.unodes;
        else
                N.ibnodes=N.uinodes;
/*
nodep=gp.UN;
for(i=0;i<N.unodes;i++) {
printf("i=%d\tnode=%d\tbelp1=%d\tbelp2=%d\n",i,nodep->n,nodep->belp[0],nodep->belp[1]);
nodep=nodep->nextunp;

}
*/
        if(theMesh->nd!=0)
                place_intnodes(num);
//printf("3\n");

        set_previousunp(num);
//printf("4\n");

//fsp=gp.FS;
//for(i=0;i<N.frontsegs;i++) {
//printf("i=%d\tfsp =%d\tnextfsp =%d\n",i,fsp,fsp->nextfsp);
//fsp=fsp->nextfsp;
//}


}

int     input_bn(num)
int     num;

{
        FILE    *fptr, *get_fptr();
        int     i,j,err;
        int     get_int();
        void    get_datapoints();
        struct  unode   *zeronp;
struct  unode   *np;    
        if((fptr = get_fptr(1) ) != NULL) {
                theMesh->nnodes = get_int(fptr);
                theMesh->nloops = get_int(fptr);
                if(num==-1) 
                        zeronp=gp.UN;
                else
                        zeronp=gp.UIN;
                if(num==-1) {
                        theMesh->nordersegs = get_int(fptr);
                        for(i=0;i<(theMesh->nordersegs*2);i++)
                                theMesh->ordernodes[i] = get_int(fptr);
                }
                if(zeronp != NULL)      {
                        free(zeronp);
                        clear_data();
                }
                for(i=0;i<theMesh->nloops;i++) 
                        get_bnodes(fptr,i,num);
                get_datapoints(fptr,num);
        }
        else {
                printf(" Unable to open that file\n") ;
                err = -1 ;
        }

//np=gp.UN;
//for(i=0;i<N.unodes;i++) {
//printf("n=%d\tu[0]=%g\tp[0]=%g\n",np->n,np->u[0],np->p[0]);
//np=np->nextunp;
//}
        return(err);

}


int     get_bnodes(f,i,num)
FILE    *f;
int     i,num;

{       
        int     j,k,fxc,icount,inter;
        int     size,gpsize;    
        double  get_dbl(),x,y;
        int     get_int();
        void    get_nodevals(),get_bvals(),free_bsegs();
        struct  unode   *nodep, *alloc_tempnode(), *firstnodep;
        struct  bsegment        *alloc_bseg(),*bsp;
        static  struct  unode   *zeronp;
        
        if((gp.BS!=NULL)&&(i==0))
                free_bsegs();
        fxc = 1;
        icount = 0;
        if(num==-1) {
                zeronp=gp.UN;
                inter=1;
        }
        else {
                zeronp=gp.UIN;
                inter=0;
        }
//printf("20\n");
        if((zeronp==NULL)&&(j=get_int(f)!=-1)) {
                x = get_dbl(f);
                y = get_dbl(f);
                zeronp = alloc_tempnode(fxc,x,y,inter);
                gp.BS = alloc_bseg();
                lastnodep = zeronp;
                get_nodevals(f,lastnodep);      
                lastbsp = gp.BS;
                get_bvals(f,lastbsp);
                firstnodep = zeronp;            
                icount++;
        }
        
        while((j = get_int(f)) != -1) {
                x = get_dbl(f);
                y = get_dbl(f);
                
                nodep = lastnodep;
//printf("21\n");
                
                lastnodep = alloc_tempnode(fxc,x,y,inter);
                get_nodevals(f,lastnodep);
                if(icount==0) {
                        firstnodep=lastnodep;
                        bsp=lastbsp;
                        lastbsp=alloc_bseg();
                        get_bvals(f,lastbsp);
                        bsp->nextbsp=lastbsp;

                }
                else {
                        bsp= lastbsp;
//                        lastbsp=alloc_bseg(N.bsegs);
                        get_bvals(f,lastbsp);
                        bsp->nps[0] = nodep;
                        bsp->nps[1] = lastnodep;
                        bsp->nextbsp = lastbsp;
                }
                nodep->nextunp=lastnodep;
                icount++;
        }
        bsp = lastbsp;
        bsp->nps[0] = lastnodep;
        bsp->nps[1] = firstnodep;
//printf("number of nodes=%d\n",N.unodes);      
        if(num==-1)
                gp.UN=zeronp;
        else
                gp.UIN=zeronp;
        return(0);
}

        

struct node *alloc_unode(bc,x,y)
int             bc ;
double  x, y ;

{
        struct node     *nodep ;
        int     k,size,gpsize ;
        
        size=sizeof(struct node);
        gpsize=sizeof(struct node *);

        nodep = (struct node *) malloc(size);
        N.nodes++;
        
        nodep->n = N.nodes  ;
        nodep->fxc =  bc;

//check 
//      gp.iptrs[index] = nodep ;
        
        nodep->x[0] = x ;
        nodep->x[1] = y ;
//check 
//      nodevalues(nodep) ;
        return(nodep) ;
}

struct unode *alloc_tempnode(bc,x,y,inter)
int             bc ;
double  x, y ;
int             inter;
{
        struct unode    *nodep ;
        int     k,size,gpsize ;
        
        size=sizeof(struct unode);
        gpsize=sizeof(struct unode *);

        nodep = (struct unode *) malloc(size);
        if(inter==0) {
                N.uinodes++;
                nodep->n=N.uinodes;
        }
        else {
                N.unodes++;
                nodep->n=N.unodes;
        }

        nodep->fxc =  bc;

        nodep->loc = bc;
        
        nodep->x[0] = x ;
        nodep->x[1] = y ;
        nodep->nbelms=0;
        nodep->rnp = NULL;

        return(nodep) ;
}


double find_min(ndim,num)
        int     ndim,num;

{       
        int     i,nnodes;
        double ymin;
        struct unode *nodep,*zeronp;
        ymin=10E8;
        if(num==-1){
                zeronp=gp.UN;
                nnodes=N.unodes;
        }
        else {
                zeronp=gp.UIN;
                nnodes=N.uinodes;
        }
        nodep = zeronp;
        for(i=0;i<nnodes;i++) {
                if(nodep->x[ndim] < ymin)
                        ymin = nodep->x[ndim];
                nodep=nodep->nextunp;
        }
        return(ymin);
}

double find_max(ndim,num)
        int     ndim,num;

{       
        int     i,nnodes;
        double ymax;
        struct unode *nodep,*zeronp;
        ymax=-10E8;
        if(num==-1){
                zeronp=gp.UN;
                nnodes=N.unodes;
        }
        else {
                zeronp=gp.UIN;
                nnodes=N.uinodes;
        }

        nodep = zeronp;
        for(i=0;i<nnodes;i++) {
                if(nodep->x[ndim] > ymax)
                        ymax = nodep->x[ndim];
                nodep=nodep->nextunp;
        }
        return(ymax);
}

struct  bsegment        *alloc_bseg()

{
        struct  bsegment        *bs;
        int     size,k;
        
        size = sizeof(struct bsegment);
        bs = (struct bsegment *) malloc(size);
        
        N.bsegs++;
        bs->n = N.bsegs;
        
        for(k=0;k<3;k++)
                bs->p[k] = defelp;
                
                
        return(bs);
}

struct tempbelement             *alloc_tbelement()

{
        int             j ;
        struct tempbelement  *belmntp ;
        long    int     size;
        
        size = sizeof(struct tempbelement);
        
        belmntp = (struct tempbelement *) malloc(size) ;
        N.belms++;
        belmntp->n = N.belms ;  
        
        return(belmntp) ;
}

struct belement         *alloc_ubelement()

{
        int             j,k ;
        struct belement  *belmntp ;
        long    int     size;
        
        size = sizeof(struct belement);
        
        belmntp = (struct belement *) malloc(size) ;    
        for(k=0;k<N.bparams;k++)
                belmntp->p[k] = defelp ;
        for(j=0;j<N.vars;j++) {
                belmntp->bcs[j] = 0 ;
        }
        belmntp->matrices = NULL ;
        
        return(belmntp) ;
}

struct closenode                *alloc_cnp(np,dist)
struct  unode   *np;
double  dist;

{
        struct closenode  *cnp ;
        long    int     size;
        
        size = sizeof(struct closenode);        
        cnp = (struct closenode *) malloc(size) ;
        
        cnp->node=np;
        cnp->node->incnp=1;
        cnp->dist=dist;
        cnp->nlamda=0;
        N.closest++;
        
        return(cnp) ;
}

void    make_tempbelements(num)
int     num;

{
        int             i,j,nbelms,nnodes,numnodes,inter;
        double  l,dl,x1,x2,y1,y2,dx,dy,x,y,ratio,ri,rn;
        struct  bsegment *bsp;
        struct  unode   *alloc_tempnode(),*np,*zeronp;
        struct  tempbelement    *alloc_tbelement(),*belp,*lastbelp,*temp;
        void            alloc_nodalbelp(),free_tbelements(),alloc_bcvals();
        void            find_tempnodevals();
        
        free_tbelements();
printf("0.84\n");       
printf("number of nodes at top of make_tempb=%d\n",N.unodes);
        if(num==-1) {
                zeronp=gp.UN;
                inter=1;
                numnodes=N.unodes;
        }
        else {
                zeronp=gp.UIN;
                inter=0;
                numnodes=N.uinodes;
        }
        np=zeronp;
        for(i=0;i<numnodes;i++) {
                np->nbelms=0;
                np=np->nextunp;
        }
printf("0.85\n");
        bsp = gp.BS;
        for(j=0;j<N.bsegs;j++) {
                x1= bsp->nps[0]->x[0];
                x2= bsp->nps[1]->x[0];
                y1= bsp->nps[0]->x[1];
                y2= bsp->nps[1]->x[1];
        
                if(theMesh->nd==0) 
                        nbelms=1;
                else {  
                        l=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

                        nbelms = l/theMesh->dy;

                        if ((l/theMesh->dy-nbelms)> 0.6)
                                nbelms++;
                        if(nbelms==0)
                                nbelms=1;
                }
                nnodes=nbelms-1;
                
                if(gp.TB==NULL) {
                        gp.TB = alloc_tbelement();
                        lastbelp = gp.TB;
                }

                else    {
                        lastbelp = alloc_tbelement();
                        belp->nexttbp = lastbelp;
                }

                if(nbelms <= 1) {
                        belp = lastbelp ;
                        belp->nps[0] = bsp->nps[0];
                        belp->nps[1] = bsp->nps[1];
                        alloc_nodalbelp(belp);
                        alloc_bcvals(bsp,belp);
                }
                else {
                        dl = l/nbelms;
                        dx=dl*(x2-x1)/l;
                        dy=dl*(y2-y1)/l;
                
                        np=zeronp;
                        if(num==-1) {
                                numnodes=N.unodes;
                        }
                        else {
                                numnodes=N.uinodes;
                        }

                        for(i=0;i<(numnodes-1);i++) 
                                np=np->nextunp;
//                      np=np->previousunp;
                
                        for(i=0;i<nnodes;i++) {
                                x = x1+(i+1)*dx;
                                y = y1+(i+1)*dy;
                                                
                                lastnodep=alloc_tempnode(1,x,y,inter);
                                ri=i;
                                rn = nbelms;
                                ratio = (ri+1)/rn;
                                find_tempnodevals(ratio,bsp,lastnodep);
                                
                                belp=lastbelp;
                                lastbelp = alloc_tbelement();
                        

                                belp->nexttbp = lastbelp;

                                if(i==0)
                                        belp->nps[0]=bsp->nps[0];
                                else
                                        belp->nps[0]=np;

                                belp->nps[1]=lastnodep;
                                alloc_nodalbelp(belp);
                                alloc_bcvals(bsp,belp);

                                np->nextunp = lastnodep;

                                np = lastnodep;
                        }
                        belp = lastbelp;
                        belp->nps[0] = np;
                        belp->nps[1] = bsp->nps[1];
                        alloc_nodalbelp(belp);
                        alloc_bcvals(bsp,belp);
                }

                bsp = bsp->nextbsp;

        }
        belp->nexttbp = gp.TB;
printf("number of nodes in belm=%d\n",N.unodes);        
}


int     get_cuts(h)

double  h;
{
        double  x1,x2,y1,y2,one;
        int     j,ncuts;
        struct  bsegment        *bsp;
        
        ncuts = 0;
        bsp = gp.BS;
        
        for(j=0;j<N.bsegs;j++) {
        
                x1 = bsp->nps[0]->x[0];
                y1 = bsp->nps[0]->x[1];
                x2 = bsp->nps[1]->x[0];
                y2 = bsp->nps[1]->x[1];

                one=(y1-h)*(y2-h);
                if (one<0||((one==0)&&((h>y1)||(h>y2)))) {

                        xc[ncuts] = x1+(h-y1)*(x2-x1)/(y2-y1);
                        ncuts++;
                }
                bsp = bsp->nextbsp;
        }
        return(ncuts);
}

void    sortdouble(double v[],int left,int right)
{
        int i,last;
        void    swap(double v[],int i,int j);
        
        if(left >= right)
                return;
        swap(v,left,(left+right)/2);
        last = left;
        for (i= left+1;i<=right;i++)
                if(v[i] < v[left])
                        swap(v,++last,i);
        swap(v,left,last);
        sortdouble(v,left,last-1);
        sortdouble(v,last+1,right);
}

void swap(double v[],int i,int j)
{
        double  temp;
        
        temp = v[i];
        v[i] = v[j];
        v[j] = temp;
}

void    make_elements(num)
int     num;

{
        int             i,j,k,check_intersect(),nnodes,nmax,ntime,nfound,found;
        struct  unode   *np,*ap,*bp,*cp,*zeronp;
        struct  frontseg        *fsp,*alloc_frontseg(),*lastfsp,*temp;
        struct  tempelement     *alloc_tempelement(),*lasttelp,*elmntp,*zerop;
        struct  lamda   *lm;
        double  xa,ya,xb,yb,dx,dy,m,get_lamda();
        double  al,almax,alpha();
        void            update_front(),free_telements(),free_closenodes();
        struct  closenode       *cnp,*scnp,*fnp;
        
        nmax=theMesh->nmax;

        if(num==0) {
                zeronp=gp.UIN;
                nnodes=N.uinodes;
                zerop=gp.IE;
                free_telements(zerop,num);
                gp.IE=NULL;
        }
        if(num==1) {
                zeronp=gp.UN;
                nnodes=N.unodes;
                zerop=gp.TE;
                free_telements(zerop,num);
                gp.TE=NULL;
        }
        zerop=NULL;

//printf("2\n");
        if(N.frontsegs<3)
                return;
        np=zeronp;
        for(i=0;i<nnodes;i++){
                np->nneighbours=0;
                np=np->nextunp;
        }

        while(gp.FS != NULL) {
//printf("1\n");

                fsp=gp.FS;
                fsp=fsp->previousfsp;
                cp=NULL;
                found=0;
//printf("2\n");
                if(gp.C!=NULL) 
                        free_closenodes();
                np=zeronp;
                for(i=0;i<nnodes;i++) {
                        np->reject=0;
                        np->incnp=0;
                        np=np->nextunp;
                }
                
                
                if(N.frontsegs==3) {
//printf("2.5\n");
                        fsp=gp.FS;
                        fsp=fsp->previousfsp;
                        ap=fsp->nps[0];
                        bp=fsp->nps[1];
                        cp=fsp->nextfsp->nps[1];
                        
                        if(zerop == NULL) {
                                zerop = alloc_tempelement(ap,bp,cp,num);
                                lasttelp = zerop;
                        }
                        else {
                                elmntp = lasttelp;
                                lasttelp = alloc_tempelement(ap,bp,cp,num);
                                elmntp->nexttelp = lasttelp;
                                lasttelp->nexttelp=zerop;
                        }
                        update_front(ap,bp,cp);
                        return;
                }

                else {  
//printf("3\n");
        
                        ap = fsp->nps[0];
                        bp = fsp->nps[1];
                        
                        xa = ap->x[0];
                        ya = ap->x[1];
                        xb = bp->x[0];
                        yb = bp->x[1];
        
                        dx = xb - xa;
                        dy = yb - ya;
        
                        ntime = 0;
                        
                        while(found==0) {
                                ntime+=1;
                                np = zeronp;
//printf("ntime=%d\n",ntime);                   
//printf("fsp 0=%d\tfsp 1=%d\n",fsp->nps[0],fsp->nps[1]);
//printf("dx=%g\n",dx);
//printf("dy=%g\n",dy);

                                for(i=0;i<nnodes;i++) {
                        
//printf("np=%d\treject=%d\n",np->n,np->reject);
//printf("node=%d\n",np->n);

//add to this if statement a check of np->reject
                                        if((np!=ap)&&(np!=bp)&&(np->reject!=1)&&(np->incnp!=1)) {
                                        
                                                if(fabs(dx)<10E-7) {
                                                        if(dy>0) {
                                                                if(np->x[0]<xa) {
                                                                        if((ntime==1)&&(num==0))
                                                                                findclosestedgenodes(fsp,np,ntime,nmax);
                                                                        else
                                                                                findclosestneighbours(fsp,np,ntime,nmax);
                                                                }
                                                        }
                                                        else {
                                                                if(np->x[0]>xa) {
                                                                        if((ntime==1)&&(num==0))
                                                                                findclosestedgenodes(fsp,np,ntime,nmax);
                                                                        else
                                                                                findclosestneighbours(fsp,np,ntime,nmax);
                                                                }
                                                        }
                                                }
                                                else {
                                                        m = dy/dx;
                                                        if(dx>0) {
//                                                              if(fabs((np->x[1]-m*np->x[0])-(ya-m*xa))>10E-7) 
                                                                if((np->x[1]-m*np->x[0])>(ya-m*xa)) {
                                                                        if((ntime==1)&&(num==0))
                                                                                findclosestedgenodes(fsp,np,ntime,nmax);
                                                                        else
                                                                                findclosestneighbours(fsp,np,ntime,nmax);
                                                                }
                                                        }
                                                        else {
//printf("before\n");

//                                                              if(fabs((np->x[1]-m*np->x[0])-(ya-m*xa))>10E-7) 
                                                                if((np->x[1]-m*np->x[0])<(ya-m*xa)) {
                                                                        if((ntime==1)&&(num==0))
                                                                                findclosestedgenodes(fsp,np,ntime,nmax);
                                                                        else
                                                                                findclosestneighbours(fsp,np,ntime,nmax);
                                                                }
//printf("after\n");
                                                        }
                                                }
                                        }
                                        np = np->nextunp;
/*
cnp=gp.C;
printf("a=%d\tb=%d\t",ap->n,bp->n);
for(i=0;i<N.closest;i++) {
printf("node= %d\treject=%d\t",cnp->node->n,cnp->node->reject);
cnp=cnp->nextcnp;

}
*/
                                }
//so far we should have pointers to all the closest nodes       
//printf("4\n");        

//printf("N    closest=%d\n",N.closest);


//cnp=gp.C;
//printf("a=%d\tb=%d\t",ap->n,bp->n);
//for(i=0;i<N.closest;i++) {
//printf("node= %d\treject=%d\t",cnp->node->n,cnp->node->reject);
//cnp=cnp->nextcnp;
//}
//printf("\n");
//printf("after cnp=cnp-nextcnp\n");
//if(ap->n==96){
//printf("its 96\n");
//cnp=gp.C;
//for(i=0;i<N.closest;i++) {
//lm=cnp->firstl;
//for(j=0;j<cnp->nlamda;j++){
//printf("node=%d\t,lamda=%g\n",cnp->node->n,lm->l);
//lm=lm->nextl;
//}
//cnp=cnp->nextcnp;
//}
//}

//printf("10\n");

                
                                if(gp.C==NULL) {
                                        printf("something wrong here\n");
                                        return;
                                }
                                else {
                                        if(N.closest==1) {
                                                cp=gp.C->node;
                                                if((check_intersect(cp,ap)==0)&&(check_intersect(cp,ap)==0))
                                                        found=1;
                                        }
                                        else {
//printf("10.5\n");
                                                cnp=gp.C;
                                                nfound=N.closest-nmax*(ntime-1);
                                                cnp=cnp->previouscnp;
//printf("10.6\n");
                                                for(i=0;i<nfound;i++) {
//printf("10.65\n");
                                                        scnp=cnp->previouscnp;
//printf("10.66\n");
                                                        for(j=0;j<N.closest;j++) {
//printf("N.closest=%d\n",N.closest);
                                                                if(scnp!=cnp) {
//printf("10.75\n");
                                                                        if(get_lamda(fsp,cnp,scnp)<=0.0) {
//printf("10.85\n");
                                                                                cnp->node->reject=1;
                                                                                break;
                                                                        }
//printf("10.95\n");
                                                                }
                                                                scnp=scnp->previouscnp;
                                                        }
                                                        cnp=cnp->previouscnp;
                                                }
//printf("11\n");
                                                cnp=gp.C;
                                                cnp=cnp->previouscnp;
                                                for(i=0;i<nfound;i++) {
                                                        if(cnp->node->reject==0) 
                                                                cnp->alpha=alpha(fsp,cnp->node);
                                                        cnp=cnp->previouscnp;
                                                }
//printf("12\n");
                                                
                                                
                                                for(k=0;k<nfound;k++) {
                                                        fnp=NULL;
                                                        cp=NULL;
                                                        cnp=gp.C;
                                                        cnp=cnp->previouscnp;
                                                        almax=0.0;

                                                        for(i=0;i<nfound;i++) {
                                                                if(cnp->node->reject!=1) {
//printf("cnosidered node =%d\n",cnp->node->n);
                                                                        lm=cnp->firstl;
                                                                        for(j=0;j<cnp->nlamda;j++) {
                                                                                al=lm->l*cnp->alpha;
                                                                                if(al>almax) {
                                                                                        cp=cnp->node;
                                                                                        fnp=cnp;
                                                                                        almax=al;
                                                                                }
                                                                                lm=lm->nextl;
                                                                        }
                                                                }
                                                                cnp=cnp->previouscnp;
                                                        }
                                                        if(cp!=NULL) {
                                                                if((check_intersect(cp,ap)==0)&&(check_intersect(cp,bp)==0)) {
                                                                        found=1;
                                                                        break;
                                                                }
                                                                else
                                                                        fnp->node->reject=1;
                                                        }
                                                        if(cp==NULL)
                                                                break;
                                                }
//printf("13\n");

                                        }
                                }
                        }
                }               
                                        
                                        
//so far we should have found cp


//printf("a=%d\tb=%d\n",ap->n,bp->n);
//printf("point c = %d\n",cp->n);
//if(ap->n==109&&bp->n==81)
//break;

//printf("a=%d\tb=%d\tc=%d\n",ap,bp,cp);
                if(zerop == NULL) {
                        zerop = alloc_tempelement(ap,bp,cp,num);
                        lasttelp = zerop;
                }
                else {
                        elmntp = lasttelp;
                        lasttelp = alloc_tempelement(ap,bp,cp,num);
                        elmntp->nexttelp = lasttelp;
                }
                if(num==0)
                        gp.IE=zerop;
                if(num==1)
                        gp.TE=zerop;
                update_front(ap,bp,cp); 
//printf("a=%d\tb=%d\tc=%d\n",ap->n,bp->n,cp->n);
//printf("element pointer =%d\n",lasttelp);
//temp = gp.FS;
//for(i=0;i<N.frontsegs;i++){
//printf("a=%d\tb=%d\t",temp->nps[0]->n,temp->nps[1]->n);
//temp=temp->nextfsp;
//}
//printf("\n");         
//PlotNodes(0);
//DrawUMesh(2);
                
        }

}

double  get_lamda(fsp,cnp,scnp)
struct  frontseg        *fsp;
struct  closenode       *cnp,*scnp;
{       
        struct  lamda   *lm,*nlm,*alloc_lamda();
        double  lamdav();
        int     i;

//printf("10.96\n");    
        if(cnp->nlamda==0) {
                nlm=alloc_lamda(cnp);
                cnp->firstl=nlm;
        }

        else {
                lm=cnp->firstl;
                for(i=0;i<(cnp->nlamda-1);i++)
                        lm=lm->nextl;
                nlm=alloc_lamda(cnp);
                lm->nextl=nlm;
        }
//printf("10.99\n");    

        nlm->l=lamdav(fsp,cnp->node,scnp->node);
//printf("10.97\n");    

        return(nlm->l);
}

struct  lamda   *alloc_lamda(cnp)       
struct  closenode       *cnp;
{
        struct lamda    *lm;
        int     size;
        
        size = sizeof(struct lamda);
        lm = (struct lamda *) malloc(size);
        cnp->nlamda++;
        return(lm);
}


double  alpha(fsp,np)
struct  frontseg        *fsp;
struct  unode   *np;

{       
        struct  unode   *ap,*bp;
        double  a,area(),lengthsq();
        
        ap=fsp->nps[0];
        bp=fsp->nps[1];
        a = area(ap,bp,np)/(lengthsq(ap,bp)+ lengthsq(bp,np)
                                        +lengthsq(np,ap));
        return(a);
}


double  lamdav(fsp,c1p,c2p)
struct  frontseg        *fsp;
struct  unode           *c1p,*c2p;
{
        double  beta,delta,area(),lengthsq();
        struct  unode   *ap,*bp;
        
        ap=fsp->nps[0];
        bp=fsp->nps[1];
                
        beta=area(c1p,bp,c2p)/(lengthsq(c1p,bp)+ lengthsq(bp,c2p)
                                        +lengthsq(c1p,c2p));
                        
        delta=area(ap,c1p,c2p)/(lengthsq(ap,c1p)+ lengthsq(c1p,c2p)
                                        +lengthsq(c2p,ap));
                        
        if(beta>delta)
                return(beta);
        else
                return(delta);
                
}


void    free_closenodes()
{
        struct  closenode       *cnp,*temp;
        struct  lamda                   *lm;
        int             i,j;
        
        if(gp.C==NULL)
                return;
                
                
        cnp=gp.C;
        for(j=0;j<N.closest;j++) {
                lm=cnp->firstl;
                cnp->node->incnp=0;
                for(i=0;i<cnp->nlamda;i++) {
                        if(lm!=NULL)
                                free(lm);
                        if(lm->nextl!=NULL)
                                lm=lm->nextl;
                }
                cnp=cnp->nextcnp;
        }
                
                
        cnp=gp.C;
        if(cnp->nextcnp==NULL) {
                free(cnp);
                return;
        }
        for(i=0;i<N.closest;i++) {
                temp=cnp->nextcnp;
                if(cnp!=NULL)
                        free(cnp);
                cnp=temp;
        }
        gp.C=NULL;      
        N.closest=0;
        return;
}


int     findclosestneighbours(fsp,np,ntime,nmax)
struct  frontseg        *fsp;
struct  unode           *np;
int             ntime,nmax;

{
        double  lt,length();
        struct  closenode       *cnp,*fnp,*alloc_cnp();
        int     found;
        void    freecnp();
int     i;      
        lt=length(fsp->nps[0],np)+length(fsp->nps[1],np);
//printf("5\n");        
        if(gp.C==NULL) {
                cnp=alloc_cnp(np,lt);
                gp.C = cnp;
        }
        else {
//printf("Nclosest=%d\n",N.closest);
//printf("6\n");
                cnp=gp.C;
                if(N.closest>1)
                        cnp=cnp->previouscnp;
                found=0;
//printf("6.25\n");
                while(cnp->dist>lt) {
                        found=1;
//printf("6.5\n");
                        if(cnp==gp.C) {
//printf("6.75\n");

                                found=2;
                                fnp=alloc_cnp(np,lt);

                                gp.C=fnp;
                                fnp->nextcnp=cnp;
                                if(N.closest==2)
                                        fnp->previouscnp=cnp;
                                else
                                        fnp->previouscnp=cnp->previouscnp;
                                cnp->previouscnp=fnp;
//printf("6.8\n");
                                if(N.closest>(nmax*ntime)) 
                                        freecnp();
                                break;
                        }
                        cnp=cnp->previouscnp;
                }
//printf("7\n");                
                if(found==1) {
// ie found a node to be inserted in the list
//printf("7.5\n");

                        fnp=alloc_cnp(np,lt);
                        fnp->previouscnp=cnp;
                        fnp->nextcnp=cnp->nextcnp;
                        cnp->nextcnp->previouscnp=fnp;
                        cnp->nextcnp=fnp;
                        if(N.closest>(nmax*ntime)) 
                                freecnp();
                }
//printf("8\n");                
                if(found==0) {
// ie the node np is not closer than any node in the list

                        if(N.closest<(nmax*ntime)) {
// append the node to the list
//printf("8.5\n");

                                fnp=alloc_cnp(np,lt);
                                fnp->previouscnp=cnp;
                                cnp->nextcnp=fnp;
                                fnp->nextcnp=gp.C;
                                gp.C->previouscnp=fnp;
                                
                        }
//  might indicate here a logic reject variable for np
                }
//printf("9\n");
        }
//printf("node=%d\n",np->n);
//printf("closest nodes in ascending order are\n");
//cnp=gp.C;
//for(i=0;i<N.closest;i++) {
//printf("i=%d\tcnp=%d\tlt=%g\n",i,cnp->node->n,cnp->dist);
//cnp=cnp->nextcnp;
//}
        return(0);
}       
int     findclosestedgenodes(fsp,np,ntime,nmax)
struct  frontseg        *fsp;
struct  unode           *np;
int             ntime,nmax;

{
        double  la,lb,length();
        static  struct  closenode       *cnpa,*cnpb,*cnp;
        struct  closenode *alloc_cnp();
        int     i;


        la=length(fsp->nps[0],np);
        lb=length(fsp->nps[1],np);
//printf("5\n");        
        if(gp.C==NULL) {
                cnpa=alloc_cnp(np,la);
                cnpb=alloc_cnp(np,lb);
                gp.C=cnpa;
                cnpa->nextcnp=cnpb;
                cnpa->previouscnp=cnpb;
                cnpb->previouscnp=cnpa;
                cnpb->nextcnp=cnpa;
        }
        else {
//printf("i am in else\n");
                if(cnpa->dist>la) {
                        cnpa->node=np;
                        cnpa->dist=la;
                }
                if(cnpb->dist>lb) {
                        cnpb->node=np;
                        cnpb->dist=lb;
                }
        }
                
        if(cnpa->node==cnpb->node) {
                N.closest=1;
                gp.C=cnpa;
        }
        else {
                N.closest=2;
                if(cnpa->dist<cnpb->dist){
                        gp.C=cnpa;
                        cnpa->nextcnp=cnpb;
                        cnpa->previouscnp=cnpb;
                        cnpb->previouscnp=cnpa;
                        cnpb->nextcnp=cnpa;
                }
                else {
                        gp.C=cnpb;
                        cnpb->nextcnp=cnpa;
                        cnpb->previouscnp=cnpa;
                        cnpa->previouscnp=cnpb;
                        cnpa->nextcnp=cnpb;
                }
        }
        
                
                
                
//printf("node=%d\n",np->n);
//printf("closest nodes in ascending order are\n");
//cnp=gp.C;
//for(i=0;i<N.closest;i++) {
//printf("i=%d\tcnp=%d\tlt=%g\n",i,cnp->node->n,cnp->dist);
//cnp=cnp->nextcnp;
//}
        return(0);
}       

void    freecnp()

{
        struct  closenode       *cnp,*tcnp;
        struct  lamda           *lm;
        int     i;
        
        cnp=gp.C;
        cnp=cnp->previouscnp;
        lm=cnp->firstl;
        for(i=0;i<cnp->nlamda;i++) {
                if(lm!=NULL)
                        free(lm);
                if(lm->nextl!=NULL)
                        lm=lm->nextl;
        }
        
        cnp->node->incnp=0;
        tcnp=cnp->previouscnp;
        if(cnp!=NULL)
                free(cnp);
        N.closest--;
        tcnp->nextcnp=gp.C;
        gp.C->previouscnp=tcnp;
}
        





int     check_nodeloc(ap,bp,np)
struct  unode   *ap,*bp,*np;
{

// This function is used to check wether node c1 is inside the triangle 
//  formed by a b c2 , in which case c1 should be chosen
        int     err;
        double  xa,ya,xb,yb,xn,yn,dx,dy,m;
        
        err=0;
        xa=ap->x[0];
        ya=ap->x[1];
        xb=bp->x[0];
        yb=bp->x[1];
        xn=np->x[0];
        yn=np->x[1];
        dx=xb-xa;
        dy=yb-ya;
        
        if(dx==0) {
                if(dy>0) {
                        if(xn<xa) 
                                err=1;
                }
                else {
                        if(xn>xa)
                                err=1;
                }
        }
        else {
                m = dy/dx;
                if(dx>0) {
                        if((yn-m*xn)>(ya-m*xa)) 
                                err=1;
                }
                else {
                        if((yn-m*xn)<(ya-m*xa)) 
                                err=1;
                }
        }
//*****************
err=0;
//*****************
        return(err);
}


int     check_alignment(np,c1p,c2p)
struct  unode   *np,*c1p,*c2p;
{
        int     err;
        double  dx1,dx2,dy1,dy2,m1,m2;
        
        err=0;
        
        dx1=c1p->x[0]-np->x[0];
        dx2=c2p->x[0]-np->x[0];
        
        if((fabs(dx1)<10E-8) && (fabs(dx2)<10E-8)) {
                printf("case 1\n");
                err=1;
        }
        else
                if((fabs(dx1)<10E-8) || (fabs(dx2)<10E-8))
                        err=0;
        
        if((fabs(dx1)>10E-8) && (fabs(dx2)>10E-8)) {
                dy1=c1p->x[1]-np->x[1];
                dy2=c2p->x[1]-np->x[1];
                m1=dy1/dx1;
                m2=dy2/dx2;
                if(fabs(m1-m2)>10E-8)
                        err=0;
                else {
                        printf("case 2\n");             
                        err=1;
                }
        }
//*******************
err=0;  
//******************
        return(err);
}
        



int     check_intersect(ap,bp)
struct  unode   *ap,*bp;

{
        int             i,vert1,vert2;
        double  xa,xb,ya,yb,dx,dy,m,A11,A21,A12,A22,C1,C2,x1max,x1min,x2max,x2min;
        double  y1min,y1max,y2min,y2max,xi,yi,DETA;     
        struct  frontseg        *fsp;
        
        vert1=0;
        vert2=0;
        xa = ap->x[0];
        ya = ap->x[1];
        xb = bp->x[0];
        yb = bp->x[1];
        dx = xb - xa;
        dy = yb - ya;
//printf("xa=%g\tya=%g\txb=%g\tyb=%g\n",xa,ya,xb,yb);   
        if(fabs(dx)>10E-7) {
                m = dy/dx;
                A11=1;
                A12=-m;
                C1=ya-m*xa;
                if(xa<xb) {
                        x1max=xb;
                        x1min=xa;
                }
                else {
                        x1max=xa;
                        x1min=xb;
                }
        }
        else {
                A11=0;
                A12=1;
                C1=xa;
                if(ya<yb) {
                        y1max=yb;
                        y1min=ya;
                }
                else {
                        y1max=ya;
                        y1min=yb;
                }
                vert1=1;
        }
        
        fsp = gp.FS;
        for(i=0;i<N.frontsegs;i++) {
                if(((ap==fsp->nps[0])&&(bp==fsp->nps[1]))||((ap==fsp->nps[1])&&(bp==fsp->nps[0])))
                        return(0);
                vert2=0;
                xa=fsp->nps[0]->x[0];
                ya=fsp->nps[0]->x[1];
                xb=fsp->nps[1]->x[0];
                yb=fsp->nps[1]->x[1];
                
                dx=xb-xa;
                dy=yb-ya;
                if(fabs(dx)>10E-7) {
                        if(xa<xb) {
                                x2max=xb;
                                x2min=xa;
                        }
                        else {
                                x2max=xa;
                                x2min=xb;
                        }
                        m=dy/dx;
                        A21=1;
                        A22=-m;
                        C2=ya-m*xa;
                }
                else {
                        A21 = 0;
                        A22 = 1;
                        C2 = xa;
                        if(ya<yb) {
                                y2max=yb;
                                y2min=ya;
                        }
                        else {
                                y2max=ya;
                                y2min=yb;
                        }
                        vert2=1;
                }
                DETA = A11*A22-A21*A12;

                if(fabs(DETA)>10E-8) {
//printf("detA = %g\n",DETA);
                        x1min+=10E-8;
                        x2min+=10E-8;
                        x1max-=10E-8;
                        x2max-=10E-8;
                        y1min+=10E-8;
                        y2min+=10E-8;
                        y1max-=10E-8;
                        y2max-=10E-8;

                        yi=(A22*C1-A12*C2)/DETA;
                        xi=(-A21*C1+A11*C2)/DETA;
//printf("nodes=%d\t%d%tyi=%g\txi=%g\n",fsp->nps[0]->n,fsp->nps[1]->n,yi,xi);                   

                        if(vert1==1) {
                                if((yi>y1min)&&(yi<y1max)&&(xi>x2min)&&(xi<x2max)) {
printf("case 1\n");
                                        ap->reject=1;
                                        return(1);
                                }
                        }
                        if(vert2==1) {
                                if((yi>y2min)&&(yi<y2max)&&(xi>x1min)&&(xi<x1max)) {
                                        printf("case 2\n");
                                        printf("yi=%g\ny2min=%g\ny2max=%g\n",yi,y2min,y2max);
                                        printf("fsp=%d\t%d\n",fsp->nps[0]->n,fsp->nps[1]->n);
                                        ap->reject=1;
                                        return(1);
                                }
                        }
                        if((vert1==0)&&(vert2==0)) {
                                if((xi>x1min)&&(xi>x2min)&&(xi<x1max)&&(xi<x2max)) {
printf("case 3 for node%d \n",ap->n);
printf("x1min=%g\tx1max=%g\tx2min=%g\tx2max=%g\txi=%g\n",x1min,x1max,x2min,x2max,xi);
printf("front segment =%d\t%d\n",fsp->nps[0]->n,fsp->nps[1]->n);
                                        ap->reject=1;
                                        return(1);
                                }
                        }
                }
                fsp=fsp->nextfsp;
        }
        return(0);
}       
        
                

void    update_front(ap,bp,cp)
struct  unode   *ap,*bp,*cp;

{

        struct frontseg *sfsp,*fspbefore,*fspafter,*newfsp,*fsp,*alloc_frontseg();
        int     i,j,ka,kb;
        
//remove the base ab first
        sfsp = gp.FS;
        sfsp = sfsp->previousfsp;
        fspbefore = sfsp->previousfsp;
        fspafter = sfsp->nextfsp;
        fspbefore->nextfsp = fspafter;
        fspafter->previousfsp = fspbefore;
        if(sfsp != NULL) {
                free(sfsp);
                N.frontsegs--;
        }
//remove ac or bc if they lie on the front allready

        ka=kb=0;
        sfsp=gp.FS;
        for(i=0;i<(N.frontsegs+2);i++) {
                j=0;
                if((sfsp->nps[0]==cp)&&(sfsp->nps[1]==ap)) {
                        fspbefore = sfsp->previousfsp;
                        fspafter = sfsp->nextfsp;
                        fspbefore->nextfsp = fspafter;
                        fspafter->previousfsp = fspbefore;
                        if(sfsp!=NULL) {
                                if(sfsp==gp.FS)
                                        gp.FS=fspafter;
                                free(sfsp);
                                N.frontsegs--;
                        }                       
                        sfsp = fspafter;
                        ka=1;
                        j=1;
                }
                if((sfsp->nps[0]==bp)&&(sfsp->nps[1]==cp)) {
                        fspbefore = sfsp->previousfsp;
                        fspafter = sfsp->nextfsp;
                        fspbefore->nextfsp = fspafter;
                        fspafter->previousfsp = fspbefore;
                        if(sfsp!=NULL) {
                                if(sfsp==gp.FS)
                                        gp.FS=fspafter;
                                free(sfsp);
                                N.frontsegs--;
                        }                       
                        sfsp = fspafter;
                        kb=1;
                        j=1;
                }
                if(j==0)
                        sfsp=sfsp->nextfsp;
        }
        if((ka==1)&&(kb==1))
                return;
        if((kb==1)||(ka==0)) {
                fsp=gp.FS;
                fsp=fsp->previousfsp;
                newfsp=alloc_frontseg();
                fsp->nextfsp=newfsp;
                newfsp->nextfsp=gp.FS;
                newfsp->previousfsp=fsp;
                fsp=gp.FS;
                fsp->previousfsp=newfsp;
                newfsp->nps[0]=ap;
                newfsp->nps[1]=cp;
        }

        if((ka==1)||(kb==0)) {
                fsp=gp.FS;
                fsp=fsp->previousfsp;
                newfsp=alloc_frontseg();
                fsp->nextfsp=newfsp;
                newfsp->nextfsp=gp.FS;
                newfsp->previousfsp=fsp;
                fsp=gp.FS;
                fsp->previousfsp=newfsp;
                newfsp->nps[0]=cp;
                newfsp->nps[1]=bp;
        }

}               

/*
void    findtwomin(fsp,np)
struct  unode   *np;
struct  frontseg        *fsp;
{
        double  lt,length();
        
        lt = length(fsp->nps[0],np) + length(fsp->nps[1],np);

        if(lt<min) { 
                c2p = c1p;
                c1p = np;
                secondmin=min;
                min = lt;
        }
        else    {
                if (lt<secondmin) {
                        c2p = np;
                        secondmin = lt;
                }
        }
        
}
*/
double  length(np1,np2)
struct  unode   *np1,*np2;
{
        double  l;
        l = sqrt((np2->x[0]-np1->x[0])*(np2->x[0]-np1->x[0])+
                                        (np2->x[1]-np1->x[1])*(np2->x[1]-np1->x[1]));
        return(l);
}

double  lengthsq(np1,np2)
struct  unode   *np1,*np2;
{
        double  l;
        l = (np2->x[0]-np1->x[0])*(np2->x[0]-np1->x[0])+
                                        (np2->x[1]-np1->x[1])*(np2->x[1]-np1->x[1]);
        return(l);
}

struct  frontseg        *alloc_frontseg()

{
        struct frontseg *fsp;
        int     size;
        
        size = sizeof(struct frontseg);
        fsp = (struct frontseg *) malloc(size);
        N.frontsegs++;
        fsp->n=N.frontsegs;
        return(fsp);
}
void    place_intnodes(num)
int     num;
{
        double  h,x,y,c;
        void            sortdouble();
        int             i,j,k,ncuts,err,nsegs,get_cuts(),nnodes,inter,nbnodes,ieven;
        struct  unode   *np,*nnp,*alloc_tempnode(),*temp,*zeronp;
        
        h=theMesh->ymin;
        if(num==-1) {
                zeronp=gp.UN;
                inter=1;
                nnodes=N.unodes;
                nbnodes=N.bnodes;
        }
        else {
                zeronp=gp.UIN;
                inter=0;
                nnodes=N.uinodes;
                nbnodes=N.ibnodes;
        }
//printf("4\n");
        ieven=0;
        while((h+=theMesh->dy)<theMesh->ymax) {
                for(i=0;i<20;i++) 
                        xc[i] = 0.0;
                while(((ncuts=get_cuts(h))/2)==((ncuts-1)/2))
                        h+=theMesh->dy/10;
//printf("ncuts=%d\n",ncuts);
                sortdouble(xc,0,(ncuts-1));
//for(i=0;i<10;i++)
//printf("%g\t",xc[i]);
//printf("\n");
//printf("nnodes=%d\n",nnodes);
                nnp=zeronp;
                if(num==-1)
                        nnodes=N.unodes;
                else
                        nnodes=N.uinodes;
                for(i=0;i<(nnodes-1);i++){
                        nnp=nnp->nextunp;
                }
//printf("5\n");
                for(k=0;k<(ncuts/2);k++) {
//printf("k=%d\n",k);
                        y=h;
                        c=0.7*theMesh->dy;
                        i=k*2;
                        x=xc[i];
                        nsegs=(xc[i+1]-xc[i])/(theMesh->dy/0.866);
                        if(((xc[i+1]-xc[i])/(theMesh->dy/0.866)-nsegs)>0.5)
                                nsegs++;
                        theMesh->dx=(xc[i+1]-xc[i])/nsegs;
                        if(ieven%2!=0.0) 
                                x+=theMesh->dx/2.0;
                        while((x+=theMesh->dx)<(xc[i+1]-10E-7)) {
                                np=zeronp;
                                err=0;
                                for(j=0;j<nbnodes;j++) {
                                        if(((x-np->x[0])*(x-np->x[0])+(y-np->x[1])*(y-np->x[1]))<c*c) 
                                                err=1;
                                        np=np->nextunp;
                                }
                                if(err==0) {
                                        lastnodep = alloc_tempnode(0,x,y,inter);
                                        nnp->nextunp=lastnodep;
                                        nnp=lastnodep;
                                }
                        }
                }
                ieven++;
                
        }
        nnp->nextunp = zeronp;
        if(num==-1)
                printf("number of nodes=%d\n",N.unodes);
        else
                printf("number of nodes=%d\n",N.uinodes);
}



void    set_previousunp(num)
int     num;

{       
        struct  unode   *np,*temp,*zeronp;
        int     i,nnodes;

        if(num==-1) {
                zeronp=gp.UN;
                nnodes=N.unodes;
        }
        else {
                zeronp=gp.UIN;
                nnodes=N.uinodes;
        }
        
        np=zeronp;
        temp=np;
        for(i=0;i<nnodes;i++) {
                np=np->nextunp;
                np->previousunp = temp;
                temp = np;
        }
}
void    make_frontsegs()

{
        struct  frontseg        *fsp,*alloc_frontseg(),*lastfsp,*temp;
        struct  tempbelement    *belp;
        struct  unode   *np;
        int     i;
        
        belp = gp.TB;
        gp.FS = alloc_frontseg();
        fsp = gp.FS;
        
        
        for(i=0;i<(N.belms-1);i++) {
                fsp->nps[0] = belp->nps[0];
                fsp->nps[1] = belp->nps[1];
                lastfsp = alloc_frontseg();
                fsp->nextfsp = lastfsp;
                belp = belp->nexttbp;
                fsp = lastfsp;
        }
        fsp->nps[0] = belp->nps[0];
        fsp->nps[1] = belp->nps[1];
        fsp->nextfsp = gp.FS;
        
        fsp = gp.FS;
        for(i=0;i<N.frontsegs;i++) {
                temp = fsp;
                fsp = fsp->nextfsp;
                fsp->previousfsp = temp;
        }       

}

void    make_ubelements()

{
        struct  belement        *belp,*alloc_ubelement(),*lastbelp;
        struct  tempbelement    *tbelp;
        void            put_bcvals();
        int     i;
        
        tbelp = gp.TB;
        gp.B = alloc_ubelement();
        belp = gp.B;
        
        for(i=0;i<N.belms;i++) {
                belp->nps[0] = tbelp->nps[0]->rnp;
                belp->nps[1] = tbelp->nps[1]->rnp;
                put_bcvals(tbelp,belp);
                belp->n=tbelp->n;
                if(theMesh->nnodes==3) {
                        belp->nnds=2;
                        belp->vtype = 111;
                        belp->gtype = 111;
                }
                else    {
                        belp->nnds=3;
                        belp->vtype=121;
                        belp->gtype=121;
                }
                lastbelp = alloc_ubelement();
                belp->nextbelp = lastbelp;
                tbelp = tbelp->nexttbp;
                belp = lastbelp;
        }
//      belp->nps[0] = tbelp->nps[0]->rnp;
//      belp->nps[1] = tbelp->nps[1]->rnp;

/*
belp=gp.B;
for(i=0;i<N.belms;i++) {
printf("belp=%d\tn1 =%d\tn2 =%d\n",belp->n,belp->nps[0]->n,belp->nps[1]->n);
belp=belp->nextbelp;
}
*/
}

void    alloc_nodalbelp(belp)
struct  tempbelement    *belp;
{
        struct  unode   *np;
        
        np=belp->nps[0];
        np->belp[np->nbelms]=belp;
        np->nbelms++;
        
        np=belp->nps[1];
        np->belp[np->nbelms]=belp;
        np->nbelms++;
}


double  area(np1,np2,np3)
struct  unode   *np1,*np2,*np3;
{       
        double  a;
        a = 0.5*((np2->x[0]-np1->x[0])*(np3->x[1]-np1->x[1])
                - (np2->x[1]-np1->x[1])*(np3->x[0]-np1->x[0]));
//printf("nodes=%d\t%d\t%d\tarea = %g\n",np1->n,np2->n,np3->n,a);
        return(a);
}

struct tempelement              *alloc_tempelement(np1,np2,np3,inter)
struct  unode   *np1,*np2,*np3;
int             inter;
{
        int   k,size;
        struct tempelement  *elmntp ;
        struct  unode   *np;
        
        size=sizeof(struct tempelement);

        elmntp = (struct tempelement *) malloc(size);
        if(inter==0){
                N.ielms++;
                elmntp->n = N.ielms ;
        }
        else {
                N.telms++;
                elmntp->n=N.telms;
        }
        elmntp->loc=0;
        
        elmntp->nps[0]=np1;
        elmntp->nps[1]=np2;
        elmntp->nps[2]=np3;
//printf("element =%d\tnode1=%d\tnode2=%d\tnode3=%d\n",elmntp->n,np1->n,np2->n,np3->n); 
        
        np=np1;
        np->neighbour[np->nneighbours]=elmntp;
        np->nneighbours++;
        
        np=np2;
        np->neighbour[np->nneighbours]=elmntp;
        np->nneighbours++;
        
        np=np3;
        np->neighbour[np->nneighbours]=elmntp;
        np->nneighbours++;
        
        
        return(elmntp) ;
}

void    smooth(nsmooth)
int             nsmooth;

{
        struct  tempelement     *elp;
        struct  unode   *np;
        int             i,j,k,l,ninterior;
        double  xtotal,ytotal;
        
        for(l=0;l<nsmooth;l++) {
                np=gp.UN;
                for(i=0;i<N.bnodes;i++)
                        np=np->nextunp;
                ninterior = N.unodes-N.bnodes;
                for(i=0;i<ninterior;i++) {
                        xtotal=0.0;
                        ytotal=0.0;
                        for(j=0;j<np->nneighbours;j++) {
                                elp=np->neighbour[j];
                                for(k=0;k<3;k++) {

                                        xtotal += elp->nps[k]->x[0];
                                        ytotal += elp->nps[k]->x[1];
                                }
                        }
                        xtotal -= (np->nneighbours*np->x[0]);
                        ytotal -= (np->nneighbours*np->x[1]);
                
                        np->x[0] = xtotal/(2*np->nneighbours);
                        np->x[1] = ytotal/(2*np->nneighbours);
                        np = np->nextunp;
                }
        }
}

void    ordermesh()
{
        struct  node    *nnp;
        struct  element *nelp;
        void            make_ubelements(),numberfrontnodes();   
        void            make_ordersegs();
        struct  orderseg        *reorder(),*osp;
        
        make_ordersegs();
        numberfrontnodes();
        while((osp=reorder())!=NULL);
//      osp=reorder();
        make_ubelements();
}


void    make_ordersegs()

{
        int             nsegment;
        struct  orderseg        *osp,*alloc_orderseg(),*newosp;
        struct  tempbelement    *belp;

        belp=gp.TB;
        for(nsegment=0;nsegment<theMesh->nordersegs;) {
                while(theMesh->ordernodes[(2*nsegment)]!=belp->nps[0]->n)
                        belp = belp->nexttbp;
                        
                nsegment++;
                while(theMesh->ordernodes[(2*nsegment-1)]!=belp->nps[1]->n) {
                        if(gp.OR == NULL) {
                                gp.OR = alloc_orderseg();
                                osp = gp.OR;
                                osp->nps[0] = belp->nps[0];
                                osp->nps[1] = belp->nps[1];
                        }
                        else {
                                newosp = alloc_orderseg();
                                newosp->nps[0] = belp->nps[0];
                                newosp->nps[1] = belp->nps[1];
                                osp->nextosp = newosp;
                                osp = newosp;
                        }
                        belp = belp->nexttbp;
                }
//printf("hello -3\n");

                if(gp.OR == NULL) {
                        gp.OR = alloc_orderseg();
                        osp = gp.OR;
                        osp->nps[0] = belp->nps[0];
                        osp->nps[1] = belp->nps[1];
                }
                else {
                        newosp = alloc_orderseg();
                        newosp->nps[0] = belp->nps[0];
                        newosp->nps[1] = belp->nps[1];
                        osp->nextosp = newosp;
                        osp = newosp;
                }
        }

        return;
}

struct  orderseg        *alloc_orderseg()

{
        struct orderseg *osp;
        int     size;
        
        size = sizeof(struct orderseg);
        osp = (struct orderseg *) malloc(size);
        osp->nextosp=NULL;
        N.ordersegs++;
        return(osp);
}

void    numberfrontnodes()

{
        struct  node    *alloc_unode(),*np,*newnp;
        struct  unode   *lastnode;
        struct  orderseg        *osp;
        int             i,j;
        
        osp=gp.OR;
        
        if(gp.N!=NULL)
                free(gp.N);
//note: all other nodes should also be freed and also gp.iptrs, etc

        gp.N = alloc_unode(osp->nps[0]->fxc,osp->nps[0]->x[0],osp->nps[0]->x[1]);
        np = gp.N;

        osp->nps[0]->rnp=np;
                
        newnp = alloc_unode(osp->nps[1]->fxc,osp->nps[1]->x[0],osp->nps[1]->x[1]);      
        osp->nps[1]->rnp = newnp;
        np->nextnp=newnp;
        np=newnp;
        
        lastnode=osp->nps[1];
        if(N.ordersegs==1)
                return;
        else 
                osp=osp->nextosp;
        
        for(i=0;i<(N.ordersegs-1);i++) {
                if(osp->nps[0] == lastnode) {
                        newnp=alloc_unode(osp->nps[1]->fxc,osp->nps[1]->x[0],osp->nps[1]->x[1]);
                        osp->nps[1]->rnp = newnp;
                        np->nextnp=newnp;
                        np=newnp;
                }
                else {
                        for(j=0;j<2;j++) {
                                newnp=alloc_unode(osp->nps[j]->fxc,osp->nps[j]->x[0],osp->nps[j]->x[1]);
                                osp->nps[j]->rnp = newnp;
                                np->nextnp=newnp;
                                np=newnp;
                        }
                }
                lastnode=osp->nps[1];
                osp=osp->nextosp;
        }
        
}       
                

struct element          *alloc_uelm(eltype)

{
        struct element  *elp ;
        long    int     size,k;
        
        size = sizeof(struct element);
        elp = (struct element *) malloc(size) ;

        elp->vtype = eltype;
        elp->gtype = eltype;
        elp->nnds =  nsf(elp->vtype) ;
        
        N.elms++;
        elp->n = N.elms ;
                
        for(k=0;k<N.params;k++)
                elp->p[k] = defelp ;
                
        elp->matrices = NULL ;

//check this one        
        if(eltype == 139) {
                elp->gtype = 111;
                elp->nnds = 6 ;
        }

//check this one
//      elvalues(elp) ;

        return(elp) ;
        
}
        
struct  orderseg        *findnextosp(osp,num)
struct  orderseg        *osp;
int             *num;
{
        struct  orderseg        *newosp;
        struct  unode   *tnp;
        int k,i,j;
        

//      this first check is to account for the case when the last
//      osp is sent
        if(osp==NULL)
                return(osp);


//      check also for the case of islands and left boundaries
//      if nps[1] should be checked
        
        k=0;
        
        while(osp->nextosp!=NULL) {
                for(j=0;j<2;j++) {
                        tnp=osp->nps[j];
                        k=0;
                        for(i=0;i<tnp->nneighbours;i++) {
                                if(tnp->neighbour[i]->loc==0) {
                                        k=1;
                                        break;
                                }
                        }
                        if(k==1)
                                break;
                        else if(j==1) {
                                newosp = osp->nextosp;
                                free(osp);
                                osp=newosp;
                        }
                }
                if(k==1)
                        break;
        }
        if(k==0) {
                if(osp->nextosp==NULL) {
                        for(j=0;j<2;j++) {
                                tnp=osp->nps[j];
                                k=0;
                                for(i=0;i<tnp->nneighbours;i++) {
                                        if(tnp->neighbour[i]->loc==0) {
                                                k=1;
                                                break;
                                        }
                                }
                                if(k==1)
                                        break;
                        }
                        *num=j;
                        if(k==0) {
                                free(osp);
                                return(NULL);
                        }
                        else
                                return(osp);
        
                }
        }
        *num=j;
        if(k==0)
                free(osp);
        return(osp);
}

struct  tempelement     *findfirstnelement(tnp,secondnode,osp)
struct  unode   *tnp,**secondnode;
struct  orderseg        *osp;
{
        struct  unode   *tnp1;
        struct  tempelement     *telp1,*telp2;
        struct  tempbelement    *tbelp;
        int     i,j,k;

//printf("30\n");
        tnp1=NULL;
        for(i=0;i<2;i++) {
                tbelp=tnp->belp[i];
//printf("tbelp=%d\n",tbelp);
                if(tbelp->nps[0]==tnp)
                        k=1;
                else
                        k=0;
                if(tbelp->nps[k]->rnp==NULL) {

                        tnp1=tbelp->nps[k];
                        break;
                }

        }
//printf("31\n");
        
        if(tnp1==NULL) {
//this part is for the cases on page 4/6/93/1
                for(i=0;i<2;i++) {
                        tbelp=tnp->belp[i];
                        if(tbelp->nps[0]==tnp)
                                k=1;
                        else
                                k=0;
//find which belp has a node on the ordering front which lies on an
// osp subsequent to the osp including tnp
                        while(osp!=NULL) {
                                if((tbelp->nps[k]==osp->nps[0])||(tbelp->nps[k]==osp->nps[1])) {
                                        tnp1=tbelp->nps[k];
                                        break;
                                }
                                osp=osp->nextosp;
                        }
                        if(tnp1!=NULL)
                                break;
                }
        }
        
        if(tnp1==NULL) {
                printf("error in finding tnp1 # 1\n");
                return(NULL);
        }
        
//seems to work so far
//now find the element that contains tnp and tnp1
        for(i=0;i<tnp->nneighbours;i++) {
                k=0;
                telp1=tnp->neighbour[i];
                for(j=0;j<tnp1->nneighbours;j++) {
                        telp2=tnp1->neighbour[j];
                        if(telp1==telp2) {
                                k=1;
                                break;
                        }
                }
                if(k==1)
                        break;
        }
//printf("32\n");

        if(k==0) {
                printf("stoooooooooooooooooooooooooooop error 2\n");
                return(NULL);
        }
        else
                telp1->loc=1;

        *secondnode=tnp1;
        return(telp1);
}


struct  orderseg        *reorder()

{       
        int     i,j,k,npsnum,k1,k2,closed;
        struct  unode   *tnp,*tnp1,*tnp2,*tnpn,*findnextelement();
        struct  unode   *closingnode,*savedtnpn,*latestnp;
        struct  orderseg        *osp,*findnextosp(),*newosp,*newosp1;
        struct  tempbelement    *tbelp;
        struct  tempelement     *telp,*findfirstnelement();
        struct  element *elp,*newelp,*alloc_uelm(),*lastelp;
        struct  node    *np,*lastnp;
        
        osp=gp.OR;
        for(i=0;i<N.ordersegs;i++) {
                osp->nps[0]->loc=2;
                osp->nps[1]->loc=2;
                osp=osp->nextosp;
        }
        closed=0;
        closingnode=NULL;
        savedtnpn=NULL;
        latestnp=NULL;
        
        N.ordersegs=0;
        osp = gp.OR;
        newosp=NULL;
        lastelp=NULL;
        lastnp=NULL;
        
        k1=0;

        if(gp.N!=NULL) {
                lastnp=gp.N;
                for(i=0;i<(N.nodes-1);i++)
                        lastnp=lastnp->nextnp;
        }
        if(gp.El!=NULL) {
                lastelp=gp.El;
                for(i=0;i<(N.elms-1);i++)
                        lastelp=lastelp->nextelp;
        }
//printf("before\n");   
        while((osp=findnextosp(osp,&npsnum))!=NULL) {
                k2=0;
//printf("after\n");
                tnp=osp->nps[npsnum];
//printf("tnp=%d\n",tnp->n);
//printf("k1=%d\n",k1);
//if(tnp->n==33)
//return(NULL);
                if(k1==0) {
                        telp=findfirstnelement(tnp,&tnp1,osp);
                        k1=1;
                }
                else {
                        if(closed==1) {
                                if(tnp2->loc==2) {
                                        if(npsnum==0) {
                                                if(tnp2!=osp->nps[1])
                                                        tnpn=tnp2;
                                                else
                                                        tnpn=tnp1;
                                        }
                                        else {
                                                if(tnp2!=osp->nextosp->nps[1])
                                                        tnpn=tnp2;
                                                else
                                                        tnpn=tnp1;
                                        }
                                }
                                else
                                        tnpn=tnp2;
                                tnp1=findnextelement(tnp,tnpn,&telp,closed,osp);
//printf("closed=%d\n",closed);                         
//printf("tnp1=%d\n",tnp1);
                                if(tnp1==NULL) {
                                        if(tnp==closingnode) {
                                                tnpn=savedtnpn;
                                                tnp1=findnextelement(tnp,tnpn,&telp,closed,osp);
                                                if(tnp1==NULL) {
                                                        printf("error #2 in finding tnp1 \n");
                                                        return(NULL);
                                                }
                                                else
                                                        closed=0;
                                        }
                                }
                        }
                        else {
                
                                if(newosp==NULL)
                                        tnpn=tnp1;
                                else
                                        tnpn=latestnp;
//                                      tnpn=newosp->nps[1];
//printf("tnpn=%d\n",tnpn->n);
                                tnp1=findnextelement(tnp,tnpn,&telp,closed,osp);
//printf("chosen tnp1 = %d\n",tnp1->n);
//printf("element=%d\t%d\t%d\n",telp->nps[0]->n,telp->nps[1]->n,telp->nps[2]->n);
                        }
                }
        
//printf("one\n");      

//now find tnp2
        for(i=0;i<3;i++) {
                if((telp->nps[i]!=tnp)&&(telp->nps[i]!=tnp1)) {
                        tnp2=telp->nps[i];
                        break;
                }
        }
//printf("tnp2 =%d\n",tnp2->n);                 
//seems to work up to here
//now check wether tnp2 lies on the old front segment 

        if(tnp1->rnp==NULL) {   
                np=alloc_unode(tnp1->fxc,tnp1->x[0],tnp1->x[1]);
                tnp1->rnp=np;
                lastnp->nextnp=np;
                lastnp=np;
                latestnp=tnp1;
        }

        if(tnp2->rnp==NULL) {
                np=alloc_unode(tnp2->fxc,tnp2->x[0],tnp2->x[1]);
                lastnp->nextnp=np;
                lastnp=np;
                tnp2->rnp=np;
                latestnp=tnp2;
                k2=1;
        }
        elp=alloc_uelm(theMesh->eltype);
//note : the lines below should be modified later to
//       allow for 6-node elements
        if(lastelp==NULL)
                gp.El=elp;
        else
                lastelp->nextelp=elp;
        elp->nps[0]=tnp->rnp;
        elp->nps[1]=tnp2->rnp;
        elp->nps[2]=tnp1->rnp;
        lastelp=elp;
//for(i=0;i<3;i++)
//printf("elment nodes are%d\n",elp->nps[i]->n); 
np=gp.N;
//for(i=0;i<N.nodes;i++){
//printf("i=%d\tnode=%d\tx=%g\n",i,np->n,np->x[0]);
//np=np->nextnp;
//}

                if(k2==1) {
                        if(newosp == NULL) {
                                gp.OR = alloc_orderseg();
                                newosp = gp.OR;
                                newosp->nps[0] = tnp1;
                                newosp->nps[1] = tnp2;
                        }
                        else {
                                newosp1 = alloc_orderseg();
                                newosp1->nps[0] = tnp1;
                                newosp1->nps[1] = tnp2;
                                newosp->nextosp = newosp1;
                                newosp = newosp1;
                        }
                }
                else {
                        if(tnp2->loc==2) {
                                if(npsnum==0) {
                                        if(tnp2!=osp->nps[1]) {
                                                closed=1;
                                                closingnode=tnp2;
                                                savedtnpn=tnp1;
                                        }
                                }
                                else{
                                        if(tnp2!=osp->nextosp->nps[1]) {
                                                closed=1;
                                                closingnode=tnp2;
                                                savedtnpn=tnp1;
                                        }
                                }
                        }
                }
        }
        return(newosp);
        
}       

struct  unode   *findnextelement(tnp,tnpn,telp,closed,osp)
struct  unode   *tnp,*tnpn;
struct  tempelement     **telp;
int             closed;
struct  orderseg        *osp;
{
        int     i,j,k;
        struct  unode   *tnp1;
        struct  tempelement     *telp1,*telp2;
        
        k=0;
        for(i=0;i<tnp->nneighbours;i++) {
                if((tnp->neighbour[i]->loc)!=1) {
                        telp1=tnp->neighbour[i];
                        for(j=0;j<tnpn->nneighbours;j++) {
                                if(tnpn->neighbour[j]->loc!=1) {
                                        telp2=tnpn->neighbour[j];
                                        if(telp1==telp2) {
                                                k=1;
                                                break;
                                        }
                                }
                        }
                }
                if(k==1)
                        break;
        }
        if(k==0) {
                if(closed==0)
                        telp1=findfirstnelement(tnp,&tnp1,osp);
                else
                        return(NULL);
        }
        else {
                tnp1=tnpn;
                telp1->loc=1;           
        }               
        *telp=telp1;
        return(tnp1);
}                                               

/*              
void    free_telements()

{
        struct  tempelement     *elp,*temp;
        int             i;
        
        if(gp.TE==NULL)
                return;
        elp=gp.TE;
//      while(elp!=NULL) {
//              temp=elp->nexttelp;
//              if(elp!=NULL)
//                      free(elp);
//              elp=temp;
//      }
        for(i=0;i<N.telms;i++) {
                temp=elp->nexttelp;
                if(elp!=NULL)
                        free(elp);
                elp=temp;
        }
        gp.TE=NULL;     
        N.telms=0;
        return;
}
*/

void    free_telements(zerop,num)
struct  tempelement     *zerop;
int     num;
{
        struct  tempelement     *elp,*temp;
        int             i;
        
        if(zerop==NULL)
                return;
        elp=zerop;

        for(i=0;i<N.telms;i++) {
                temp=elp->nexttelp;
                if(elp!=NULL)
                        free(elp);
                elp=temp;
        }
        zerop=NULL;
        if(num==0)      
                N.ielms=0;
        if(num==1)
                N.telms=0;
        return;
}


void    free_tbelements()
{
        struct  tempbelement    *elp,*temp;
        int             i;
        
        if(gp.TB==NULL)
                return;
        elp=gp.TB;

        for(i=0;i<N.belms;i++) {
                temp=elp->nexttbp;
                if(elp!=NULL)
                        free(elp);
                elp=temp;
        }
        gp.TB=NULL;     
        N.belms=0;
        return;
}

void    free_bsegs()
{
        struct  bsegment        *bsp,*temp;
        int             i;
        
        if(gp.BS==NULL)
                return;
        bsp=gp.BS;

        for(i=0;i<N.bsegs;i++) {
                temp=bsp->nextbsp;
                if(bsp!=NULL)
                        free(bsp);
                bsp=temp;
        }
        gp.BS=NULL;     
        N.bsegs=0;
        return;
}
        
void    redo_frontsegs()

{
        struct  frontseg        *fsp,*alloc_frontseg(),*lastfsp;
        struct  tempbelement    *belp;
        int     i;
        
        belp = gp.TB;
        gp.FS=alloc_frontseg();
        fsp=gp.FS;
        
        for(i=0;i<(N.belms-1);i++) {
                fsp->nps[0] = belp->nps[0];
                fsp->nps[1] = belp->nps[1];
                lastfsp = alloc_frontseg();
                fsp->nextfsp = lastfsp;
                belp = belp->nexttbp;
                fsp = lastfsp;
        }
        fsp->nps[0] = belp->nps[0];
        fsp->nps[1] = belp->nps[1];
        fsp->nextfsp = gp.FS;
        
        fsp=gp.FS;
        for(i=0;i<N.frontsegs;i++) {
                lastfsp = fsp;
                fsp = fsp->nextfsp;
                fsp->previousfsp = lastfsp;
        }

}

void    get_datapoints(f,num)
FILE    *f;
int     num;

{

// this function reads and allocates memory for the
// interior data points in a mesh
        int     i,fxc,nnodes,inter;
        double  get_dbl(),x,y;
        int     get_int();
        void    get_nodevals();
        struct  unode   *np,*alloc_tempnode(),*nnp,*zeronp;
        
//printf("23\n");
//printf("num=%d\n",num);
        fxc=2;
        if(num==-1){
                zeronp=gp.UN;
                nnodes=N.unodes;
                inter=1;
        }
        else {
                zeronp=gp.UIN;
                nnodes=N.uinodes;
                inter=0;
        }
//printf("23.5\n");

        np=zeronp;
        for(i=0;i<(nnodes-1);i++){
//printf("i=%d\tnp=%d\n",i,np);
                np=np->nextunp; 
        }
//printf("24\n");

        while((i=get_int(f)) != -1) {
                x=get_dbl(f);
                y=get_dbl(f);
                
                nnp=alloc_tempnode(fxc,x,y,inter);
                get_nodevals(f,nnp);
                np->nextunp=nnp;
                np=nnp;
        }
//printf("25\n");

        np->nextunp=zeronp;
}

        
void    get_nodevals(f,np)
FILE    *f;
struct  unode   *np;

{       
        double  get_dbl();
        
        np->p[0]=get_dbl(f);
        np->p[1]=get_dbl(f);
        np->u[0]=get_dbl(f);
        np->u[1]=get_dbl(f);
        np->u[2]=get_dbl(f);
}

void    get_bvals(f,bsp)
FILE    *f;
struct  bsegment        *bsp;
{
        double  get_dbl();
        int     get_int();      
        
        bsp->p[0]=get_dbl(f);
        bsp->p[1]=get_dbl(f);
        bsp->p[2]=get_dbl(f);
        bsp->bcs=get_int(f);
}

void    alloc_bcvals(bsp,belp)
struct  bsegment        *bsp;
struct  tempbelement    *belp;

{
        belp->p[0]=bsp->p[0];
        belp->p[1]=bsp->p[1];
        belp->p[2]=bsp->p[2];
        belp->bcs=bsp->bcs;

}

void    put_bcvals(tbelp,belp)
struct  tempbelement    *tbelp;
struct  belement                        *belp;

{
        belp->p[0]=tbelp->p[0];
        belp->p[1]=tbelp->p[1];
        belp->p[2]=tbelp->p[2];
        belp->bcs[0]=tbelp->bcs;

}

void    find_nodalvalues()

{
        int             i,j,k;
        int             check_unodeloc(),interpolate();
        struct  tempelement     *telp;
        struct  node                    *np;
        struct  unode                   *unp;
        void            set_inoutbcs();
        
        np=gp.N;

//telp=gp.IE;
//for(i=0;i<N.ielms;i++) {
//printf("i=%d\telemn=%d\tnodes=%d\t%d\t%d\n",i,telp->n,telp->nps[0]->n,telp->np//s[1]->n,telp->nps[2]->n);
//telp=telp->nexttelp;
//}
        
        for(i=0;i<N.nodes;i++) {
                telp=gp.IE;
                k=0;
                for(j=0;j<N.ielms;j++) {
                        if((check_unodeloc(telp->nps[0],telp->nps[1],np)==1)&&
                                (check_unodeloc(telp->nps[1],telp->nps[2],np)==1)&&
                                (check_unodeloc(telp->nps[2],telp->nps[0],np)==1)) {
                                k=1;
                                break;
                        }
                telp=telp->nexttelp;
                }
                if(k==0) {
                        printf("error 4 at node number %d\n",np->n);
                        break;
                }
                else
                        interpolate(np,telp);

                np=np->nextnp;
        }
/*
        unp=gp.UN;
        for(i=0;i<N.bnodes;i++) {
                for(j=0;j<2;j++) 
                        unp->rnp->p[j]=unp->p[j];
                for(j=0;j<3;j++)
                        unp->rnp->u[j]=unp->u[j];
                unp=unp->nextunp;
        }
*/      
        set_inoutbcs();         
        
}

void    set_inoutbcs()

{
        struct  belement        *belp;
        int             i,j;
        
        belp=gp.B;
        
        for(i=0;i<N.belms;i++) {
                if(belp->bcs[0] != 0) {
                        for(j=0;j<3;j++) 
                                belp->p[j]=(belp->nps[0]->u[j]+belp->nps[1]->u[j])/2.0;
                }
                belp=belp->nextbelp;
        }
}
                                
        
        
        

int     check_unodeloc(ap,bp,np)
struct  unode   *ap,*bp;
struct  node    *np;
{
        int     err;
        double  xa,ya,xb,yb,xn,yn,dx,dy,m;
        
        err=0;
        xa=ap->x[0];
        ya=ap->x[1];
        xb=bp->x[0];
        yb=bp->x[1];
        xn=np->x[0];
        yn=np->x[1];
        dx=xb-xa;
        dy=yb-ya;
        
        if(fabs(dx)<10E-7) {
                if(dy>0) {
                        if(xn<=xa) 
                                err=1;
                }
                else {
                        if(xn>=xa)
                                err=1;
                }
        }
        else {
                m = dy/dx;
                if(dx>0) {
                        if((yn-m*xn)>=(ya-m*xa)) 
                                err=1;
                }
                else {
                        if((yn-m*xn)<=(ya-m*xa)) 
                                err=1;
                }
        }
        return(err);
}


int             interpolate(np,telp)
struct  node                    *np;
struct  tempelement     *telp;

{       
        int             i;
        double  aa,ab,ac,atotal,area1();
        struct  unode   *ap,*bp,*cp;
        
        ap=telp->nps[0];
        bp=telp->nps[1];
        cp=telp->nps[2];
//printf("a=%d\tb=%d\tc=%d\n",ap->n,bp->n,cp->n);
//for(i=0;i<3;i++)
//printf("i=%d\tn=%d\tu[0]=%g\n",i,telp->nps[i]->n,telp->nps[i]->u[0]);
        aa=area1(bp,cp,np);
        ab=area1(cp,ap,np);
        ac=area1(ap,bp,np);
        
        atotal=aa+ab+ac;
        
        for(i=0;i<2;i++)
                np->p[i]=(ap->p[i]*aa+bp->p[i]*ab+cp->p[i]*ac)/atotal;
        
        for(i=0;i<3;i++)
                np->u[i]=(ap->u[i]*aa+bp->u[i]*ab+cp->u[i]*ac)/atotal;
        return(0);
}

        
double  area1(np1,np2,np3)
struct  unode   *np1,*np2;
struct  node    *np3;
{       
        double  a;
        a = 0.5*((np2->x[0]-np1->x[0])*(np3->x[1]-np1->x[1])
                - (np2->x[1]-np1->x[1])*(np3->x[0]-np1->x[0]));
//printf("nodes=%d\t%d\t%d\tarea = %g\n",np1->n,np2->n,np3->n,a);
        return(a);
}


void    find_tempnodevals(r,bsp,np)
double  r;
struct  bsegment        *bsp;
struct  unode           *np;

{
        int     i;
        
        for(i=0;i<2;i++)
                np->p[i]=bsp->nps[0]->p[i]*(1-r)+bsp->nps[1]->p[i]*r;
                
        for(i=0;i<3;i++) 
                np->u[i]=bsp->nps[0]->u[i]*(1-r)+bsp->nps[1]->u[i]*r;
        
}

