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



int cgns_input(fid)
int fid;
{
	int fid2;
	int ier, ii;
	int nBases, index_base, index_zone, nuser_data;
	int user_iter;
	int celldim, phydim;
	int stat;
    char zonename[33],sectionname[33],username[33];
	char solname[33], cgnsHSFile[250];
	GridLocation_t loc;
	char arrayname1[100], arrayname2[100];
	int iarray,datatype, nndim;
	int dimvals[3];
	int isize[3][1];
	int irmin,irmax;
	int n_pts, n_elems;
	int pflag;
	double *xcord, *ycord;
	double *bndnodex, *bndnodey;
	double *elevation, *roughness;
	double *depth, *vx, *vy;
	int *node, *element;
	int *edgenum, *edgenode1, *edgenode2;
	int *edgenum_s, *edgenode1_s, *edgenode2_s;

    long            bigint ;
	FILE* BndFile;//name of the CDG file
	int bnvert, bndim, bnatt, bnbndmrk;
	int bnseg, segnum;
	int enedges, enbndmrk, enumber,bndedgecount;
	int tint0, tint1, tint2, tint3;
	int orientpos, orientfound;
	int selem, scount;
	char curpath[250], temp[250];
	//for node init
	int i,j;
	struct node *nodep;
	double upstreamInitWSE;
	//for eleminit
	struct element *elmntp;
	struct node *np;
	int jj, nnds, n, rc;
	//for boundary elements init
    struct belement  *belmntp ;
	//for boundary conditions
	struct boundaryseg  *bseg ;
	char c;
	int inflow_count, outflow_count;
	int inflowType, outflowType;
	double fixedStage, fixedQ;
	int bcseg_count;
	int *bndsegs;
	int indexSize;
	int startIndex, endpt1found, endpt2found;
	int checkelem, endseg1, endseg2;
	double dd_k, dd_m;
	// for hotstart
	int hotstart;
	int nsols;

	hotstart = 0;

	orientpos = -1;
	orientfound = 0;

	getcwd(curpath, sizeof(curpath));




	ier = cg_iRIC_Init(fid);
	//ier = cg_iRIC_Read_Integer("m_CalcType", &N.trans);
	ier = cg_iRIC_Read_Integer("initType", &hotstart);
	ier = cg_iRIC_Read_String("hotStartFile", &cgnsHSFile[0]);
	N.trans = 1;
	N.meshtype = 1;
	if(N.trans == 1) {
		N.dims = 2; 
		tvals.nsteps = 1;
		ier = cg_iRIC_Read_Real("m_StartTime", &tvals.t);
		ier = cg_iRIC_Read_Real("m_InitDeltaT", &tvals.dt);
		tvals.theta = 1;
		ier = cg_iRIC_Read_Real("UW", &UW);
		DIF = 0.5;
		sets.latitude = 0;
		ier = cg_iRIC_Read_Real("tvals_S", &tvals.S);
		sets.diffusivewave = 0;
		sets.uwJ = 0;
		sets.plotcode = 2;
		sets.transbcs = 0;
		ier = cg_iRIC_Read_Integer("m_MaxNumIter", &sets.maxitnum);
		sets.smallH = 1;
		sets.JE = 1;
		ier = cg_iRIC_Read_Real("epsilon", &epsilon1);
		ier = cg_iRIC_Read_Real("tvals_gwH", &tvals.gwH);
		ier = cg_iRIC_Read_Real("epsilon3", &epsilon3);
		ier = cg_iRIC_Read_Real("tvals_T", &tvals.T);
		tvals.sg_ice = 0.92;
		tvals.itnum = 0;
		ier = cg_iRIC_Read_Real("m_UpstreamWSE", &upstreamInitWSE);
	}
	N.dims = 2;
	N.vars = 3;
	for(ii = 0; ii < ((N.vars+1)*N.vars); ii++) 
		N.Keqns[ii] = 1;
	N.params = 3;
	N.bparams = 7;

	if(gp.N != NULL) {
		clear_data();
		free(gp.N);
		free(gp.iptrs);
	}
	
    index_base=1;
    index_zone=1;
    stat = cg_zone_read(fid,index_base,index_zone,zonename,isize[0]);
    irmin=1;
    n_pts=isize[0][0];
	n_elems = isize[1][0];
	//xcord = new double[irmax];
	xcord = (double *) malloc(n_pts*sizeof(double));
	ycord = (double *) malloc(n_pts*sizeof(double));
	bndnodex = (double *) malloc(n_pts*sizeof(double));
	bndnodey = (double *) malloc(n_pts*sizeof(double));
	element = (int *) malloc((n_elems*3)*sizeof(int));
	elevation = (double *) malloc(n_pts*sizeof(double));
	roughness = (double *) malloc(n_pts*sizeof(double));

	if(hotstart) {
		depth = (double *) malloc(n_pts*sizeof(double));
		vx = (double *) malloc(n_pts*sizeof(double));
		vy = (double *) malloc(n_pts*sizeof(double));
		stat = cg_open(cgnsHSFile, CG_MODE_READ, &fid2);
		stat = cg_nsols(fid2, index_base, index_zone, &nsols);
		stat = cg_sol_info(fid2,index_base,index_zone,nsols,solname,&loc);
		stat = cg_field_read(fid2,index_base,index_zone,nsols-1,"Depth",
                  RealDouble,&irmin,&n_pts,&depth[0]);
		stat = cg_field_read(fid2,index_base,index_zone,nsols-1,"VelocityX",
                  RealDouble,&irmin,&n_pts,&vx[0]);
		stat = cg_field_read(fid2,index_base,index_zone,nsols-1,"VelocityY",
                  RealDouble,&irmin,&n_pts,&vy[0]);

	}	
/* read grid coordinates */
    stat = cg_coord_read(fid,index_base,index_zone,"CoordinateX",
                  RealDouble,&irmin,&n_pts,&xcord[0]);
    stat = cg_coord_read(fid,index_base,index_zone,"CoordinateY",
                  RealDouble,&irmin,&n_pts,&ycord[0]);

	stat = cg_elements_read(fid, index_base, index_zone, 1, &element[0], &pflag);

	// find grid conditions
	//stat = cg_goto(fid, index_base, "Zone_t", index_zone, "end"); 
	//stat = cg_nuser_data(&nuser_data);
	//stat = cg_user_data_read(nuser_data, username);
	//stat = cg_goto(fid, index_base, "Zone_t", index_zone, "UserDefinedData_t", 1, "end"); 
	//stat = cg_goto(fid, index_base, "Zone_t", index_zone, "UserDefinedData_t", 1, "UserDefinedData_t", 1, "end");
	//stat = cg_narrays(&nuser_data);
	//stat = cg_array_read(1, &elevation[0]);
	//stat = cg_goto(fid, index_base, "Zone_t", index_zone, "UserDefinedData_t", 1, "UserDefinedData_t", 2, "end");
	//stat = cg_array_read(1, &roughness[0]);
	stat = cg_iRIC_Read_Grid_Real_Node("elevation", &elevation[0]);
	stat = cg_iRIC_Read_Grid_Real_Node("roughness", &roughness[0]);

	N.nodes = n_pts;
	Nndso = N.nodes;

    bigint = (5 * N.nodes)/4 * sizeof(struct node *) ;
	if(N.nodes < 2 || (gp.iptrs = (struct node **) malloc(bigint)) == NULL) {
		printf("\n");
		printf("Error in allocating memory for node pointers \n");
		printf("or numnber of nodes is less than 2 \n");
		printf("\n");
		exit(1);
	}
    bigint = (long int) N.nodes * (long int) sizeof(struct node) ;
	if(N.nodes < 2 || (gp.N = (struct node *) malloc(bigint)) == NULL) {
		printf("\n");
		printf("Error in allocating memory for elements \n");
		printf("or numnber of nodes is less than 2 \n");
		printf("\n");
		exit(1);
	}

    N.elms = n_elems ;
    Nelmso = N.elms ;
    bigint = (long int) N.elms * (long int) sizeof(struct element) ;
    if(gp.El != NULL)
            free(gp.El) ;
	if(N.elms < 1 || (gp.El = (struct element *) malloc(bigint)) == NULL) {
		printf("\n");
		printf("Error in allocating memory for elements \n");
		printf("or number of elements is < 1 \n");
		printf("\n");
		exit(1);
	}

	BndFile=fopen("iRICZone.1.edge", "r");
	if (BndFile == NULL) {
	  fprintf(stderr, "Can't open input file iRICZone.1.edge!\n");
	  exit(1);
	}
		fscanf(BndFile, "%d %d", &enedges, &enbndmrk);
		edgenum = (int *) malloc((enedges)*sizeof(int));
		edgenode1 = (int *) malloc((enedges)*sizeof(int));
		edgenode2 = (int *) malloc((enedges)*sizeof(int));
		edgenum_s = (int *) malloc((enedges)*sizeof(int));
		edgenode1_s = (int *) malloc((enedges)*sizeof(int));
		edgenode2_s = (int *) malloc((enedges)*sizeof(int));
		bndedgecount = 0;

		for(enumber = 0; enumber < enedges; enumber++) {
			fscanf(BndFile, "%d %d %d %d", &tint0, &tint1, &tint2, &tint3);
			if(tint3 == 1) {
				++bndedgecount;
				edgenum[bndedgecount-1] = tint0;
				edgenode1[bndedgecount-1] = tint1;
				edgenode2[bndedgecount-1] = tint2;
					if(bndedgecount == 1) {
						edgenum_s[bndedgecount-1] = tint0;
						edgenode1_s[bndedgecount-1] = tint1;
						edgenode2_s[bndedgecount-1] = tint2;
					} 
			}
		}
		N.belms = bndedgecount;
        Nbelmso = N.belms ;
        bigint = (long int) N.belms * (long int) sizeof(struct belement) ;
        if(gp.B != NULL)
                free(gp.B) ;
		if(N.belms < 1 || (gp.B = (struct belement *) malloc( bigint)) == NULL) {
				printf("\n");
				printf("Error in boundary element allocation \n");
				printf("\n");
				exit(1);
		}

		scount = 0;
		do  
		{
			selem = linearSearch(edgenode1, edgenode2_s[scount], bndedgecount);
			if(selem != -1)
			{
				++scount;
			
				edgenum_s[scount] = edgenum[selem];
				edgenode1_s[scount] = edgenode1[selem];
				edgenode2_s[scount] = edgenode2[selem];
				bndnodex[bndedgecount-1] = xcord[tint1];
				bndnodey[bndedgecount-1] = ycord[tint2];

			} else {
				// to do error
				printf("\n");
				printf("Error in ordering boundary elements \n");
				printf("\n");
				exit(1);

			}
		} while( scount <= bndedgecount);
			

	fclose(BndFile);
	
	for( i = 0; i < bndedgecount; i++) {
		bndnodex[i] = xcord[edgenode1_s[i]];
		bndnodey[i] = ycord[edgenode1_s[i]];
	}

	orientpos = ClockWise(bndnodex, bndnodey, bndedgecount);
	//if(orientpos == 0 || orientpos == 1) {
	//	// error - boundary line is oriented 
	//	printf("\n");
	//	printf("Bounding polygon is drawn clockwise - must be anti-clockwise \n");
	//	printf("\n");
	//	exit(1);
	//}

	//Read Nodes 
	ier = cg_iRIC_Read_Real("m_UpstreamWSE", &upstreamInitWSE);

	for( i = 0; i < N.nodes; i++) {
		nodep = gp.N + i;
		nodep->n = (i+1);
		nodep->i = i;
		nodep->fxc = 0;
		gp.iptrs[i] = nodep;
        nodep->x[0] = xcord[i] ;
        nodep->x[1] = ycord[i] ;
        for(j=N.dims;j<NDIMS;j++) 
			nodep->x[j] = 0.0 ;
        nodep->p[0] = elevation[i] ;
        nodep->p[1] = roughness[i] ;
        nodep->p[2] = 10 ;
		if(!hotstart) {
			nodep->u[0] = upstreamInitWSE - elevation[i] ;
			nodep->u[1] = 0.0 ;
			nodep->u[2] = 0.0 ;
		} else {
			nodep->u[0] = depth[i];
			nodep->u[1] = vx[i] ;
			nodep->u[2] = vy[i] ;
		}
		for(j=0;j<N.params;j++)
				nodep->ice[j] = 0.0;
        nodep->nextnp = nodep + 1 ;
	}

	// read elements
	for( i = 0; i < N.elms; i++) {
        elmntp = gp.El + i ;
        elmntp->n = i ;
        elmntp->vtype = 210;
        elmntp->gtype = 210;
        elmntp->nnds =  nsf(elmntp->vtype) ;
        nnds = nsf(elmntp->gtype);
        if(nnds > elmntp->nnds) elmntp->nnds = nnds ;
        if(elmntp->vtype == 229)
                elmntp->nnds = 16 ;
        if(elmntp->vtype == 129)
                elmntp->nnds = 4 ;
        if(elmntp->vtype == 139)
                elmntp->nnds = 6 ;
//        for(j=0;j<elmntp->nnds;j++) {
//                if((n = get_int(f)) >= 0) { 
// /*                      np = gp.N ;
//                        for(jj=0;jj<N.nodes;jj++) {
//                                if(np->n == n)
//                                        break ;
//                                np = np->nextnp ;
//                        } 
//                        if(np == NULL) {
//                                printf(" Non-existent node %d referred to in element %d\n",n,elmntp->n) ;
//                                rc = -1 ;
//                        }
//                        elmntp->nps[j] = np ;
//*/
//						elmntp->nps[j] = gp.N + (n-1);
//                }
//                else
//                        elmntp->nps[j] = NULL ;
//        }
		elmntp->nps[0] = gp.N + (element[(i+1)*3-3]-1);
		elmntp->nps[1] = gp.N + (element[(i+1)*3-2]-1);
		elmntp->nps[2] = gp.N + (element[(i+1)*3-1]-1);

        for(j=0;j<N.params;j++)
                elmntp->p[j] = 0.0 ;
        elmntp->matrices = NULL ;
        if(i < N.elms - 1)
                elmntp->nextelp = elmntp + 1 ;
        else
                elmntp->nextelp = NULL ;
	}

	// read boundary elements
	for(i = 0; i < N.belms; i++) {
		belmntp = gp.B + i ;
        belmntp->n = edgenum_s[i] ;
        belmntp->vtype = 111;
        belmntp->gtype = 111;
        belmntp->nnds =  nsf(belmntp->vtype) ;
        nnds = nsf(belmntp->gtype);
        if(nnds > belmntp->nnds) belmntp->nnds = nnds ;
//        for(j=0;j<belmntp->nnds;j++) {
//                n = get_int(f) ;
// /*               np = gp.N ;
//                for(jj=0;jj<N.nodes;jj++) {
//                        if(np->n == n)
//                                break ;
//                        np = np->nextnp ;
//                } 
//                if(np == NULL) {
//                        printf(" Non-existent node %d referred to in element %d\n",n,belmntp->n) ;
//                        rc = -1 ;
//                }
//                belmntp->nps[j] = np ;
//*/
// 				belmntp->nps[j] = gp.N + (n-1);
//       }
	   belmntp->nps[0] = gp.N + (edgenode1_s[i]-1);
	   belmntp->nps[1] = gp.N + (edgenode2_s[i]-1);
        /*
        np->fxc = 1;
         aaa */
		// for iric case we set the following two loops to zero and let
	   // River2D paramaterize these values correctly according to the 
	   // Boundary Segments
        for(j=0;j<N.bparams;j++)
                belmntp->p[j] = 0.0 ;
        for(j=0;j<N.vars;j++) {
                belmntp->bcs[j] = 0.0 ;
        }
        belmntp->matrices = NULL ;
        belmntp->elp=NULL;
        belmntp->bseg=NULL;
        if(i < N.belms - 1)
                belmntp->nextbelp = belmntp + 1;
        else
                belmntp->nextbelp = gp.B ;
	}

	find_neighbour();
	
	// Get boundary Conditions
    cg_iRIC_Read_BC_Count("inflowCondition", &inflow_count);
    cg_iRIC_Read_BC_Count("OutflowCondition", &outflow_count);

	N.boundarysegs = outflow_count+inflow_count;
	Nbsego = N.boundarysegs;
    if(N.boundarysegs>0) {
            bigint = (long int) N.boundarysegs * (long int) sizeof(struct boundaryseg) ;
            if(gp.BSEG != NULL)
                    free(gp.BSEG) ;
			if(N.boundarysegs < 1 || (gp.BSEG = (struct boundaryseg *) malloc( bigint)) == NULL) {
				printf("\n");
				printf("Error in allocating memory for boundary segments \n");
				printf("or numnber of boundary segments is less than 1 \n");
				printf("\n");
				exit(1);
			}
    }
    //From Peter email on 5/13/2012
	//The eighth column is the type (
	//	1=subcritical inflow(q specified), 
	//	3=subcritical outflow (fixed elevation), 
	//	5=subcritical outflow ( q = f(depth)). 
	//	codes 2 (q and elev fixed) and 4(nothing fixed) are supercritical inflow and outflow. 
	//	These last 2 are not available in the river2D interface. 
	//	The first column is the fixed elevation and the second is the fixed q, as required by the BC type.

	// for now we will assume there are only one of each;
	bcseg_count = 0;
	for(i = 1; i <= outflow_count; i++) {
		++bcseg_count;
		cg_iRIC_Read_BC_Integer("OutflowCondition", i, "Type", &outflowType);
		//<Enumeration value="0" caption="Fixed Elevation" />
        //<Enumeration value="1" caption="Time Varying Elevation" />
        //<Enumeration value="2" caption="Stage-Total Discharge Relationship (Rating Curve)" />
        //<Enumeration value="3" caption="Depth-Unit Discharge Relationship q = Kh^m" />

        cg_iRIC_Read_BC_IndicesSize("OutflowCondition", i, &indexSize);
			
		bndsegs = (int *) malloc((indexSize)*sizeof(int));

		cg_iRIC_Read_BC_Indices("OutflowCondition", i, &bndsegs[0]);

		//find boundary segment endpoints

		startIndex = bndsegs[0];
		endpt1found = 0;
		endpt2found = 0;
		//forward search to find endpoint1
		selem = linearSearch(edgenode1_s, startIndex, N.belms);
		do {
		
		//check if edgenode2_s is in bndsegs
		checkelem = linearSearch(bndsegs, edgenode2_s[selem], indexSize);
		if(checkelem != -1) {
			// not at endpoint
			++selem;
			if(selem >= N.belms) {
				selem = 1;
			}

		} else {
			endpt1found = 1;
			endseg1 = edgenode1_s[selem];
		}
		} while (endpt1found == 0);

		//backward search to find endpoint1
		selem = linearSearch(edgenode1_s, startIndex, N.belms);
		--selem;
		if(selem <= 0) {
			selem = N.belms;
		}
		do {
		
		//check if edgenode1_s is in bndsegs
		checkelem = linearSearch(bndsegs, edgenode1_s[selem], indexSize);
		if(checkelem != -1) {
			// not at endpoint
			--selem;
			if(selem <= 0) {
				selem = N.belms;
			}

		} else {
			endpt2found = 1;
			endseg2 = edgenode2_s[selem];
		}
		} while (endpt2found == 0);

		bseg = gp.BSEG + bcseg_count - 1 ;
        bseg->n = bcseg_count ;
		if(outflowType == 0) {
			bseg->bcs[0] = 3;
			cg_iRIC_Read_BC_Real("OutflowCondition", i, "fixedElevation", &fixedStage);
			bseg->p[0] = fixedStage;
			bseg->p[1] = 0;
		} else if(outflowType == 1 || outflowType == 2) {
			printf("\n");
			printf("Transient flow solver has not been implemented yet \n");
			printf("\n");
			exit(1);
		} else {
			bseg->bcs[0] = 5;
			cg_iRIC_Read_BC_Real("OutflowCondition", i, "depthdischarge_k", &dd_k);
			cg_iRIC_Read_BC_Real("OutflowCondition", i, "depthdischarge_m", &dd_m);
			bseg->p[0] = dd_k;
			bseg->p[1] = dd_m;
		}
        
        bseg->nstart = endseg2;
        bseg->nend = endseg1;

		// no variable bc for now so
		c = "";
		//c = getc(f);
		//if(isdigit(c))
		//	ungetc(c,f);
		//else
		//{
		//while (c != '\n')
		//	c = getc(f);
		//}
		    
        //if(i < N.boundarysegs - 1)
                bseg->nextbseg = bseg + 1;
        //else
        //        bseg->nextbseg = NULL ;	
	free(bndsegs);
	}
	

	for(i = 1; i <= inflow_count; i++) {
		++bcseg_count;
		cg_iRIC_Read_BC_Integer("inflowCondition", i, "Type", &inflowType);
		cg_iRIC_Read_BC_Real("inflowCondition", i, "fixedDischarge", &fixedQ);
        cg_iRIC_Read_BC_IndicesSize("inflowCondition", i, &indexSize);
			
		bndsegs = (int *) malloc((indexSize)*sizeof(int));

		cg_iRIC_Read_BC_Indices("inflowCondition", i, &bndsegs[0]);

		//find boundary segment endpoints

		startIndex = bndsegs[0];
		endpt1found = 0;
		endpt2found = 0;
		//forward search to find endpoint1
		selem = linearSearch(edgenode1_s, startIndex, N.belms);
		do {
		
		//check if edgenode2_s is in bndsegs
		checkelem = linearSearch(bndsegs, edgenode2_s[selem], indexSize);
		if(checkelem != -1) {
			// not at endpoint
			++selem;
			if(selem >= N.belms) {
				selem = 1;
			}
		} else {
			endpt1found = 1;
			endseg1 = edgenode1_s[selem];
		}
		} while (endpt1found == 0);

		//backward search to find endpoint1
		selem = linearSearch(edgenode1_s, startIndex, N.belms);
		--selem;
		if(selem <= 0) {
			selem = N.belms;
		}

		do {
		
		//check if edgenode1_s is in bndsegs
		checkelem = linearSearch(bndsegs, edgenode1_s[selem], indexSize);
		if(checkelem != -1) {
			// not at endpoint
			--selem;
			if(selem <= 0) {
				selem = N.belms;
			}

		} else {
			endpt2found = 1;
			endseg2 = edgenode2_s[selem];
		}
		} while (endpt2found == 0);

		bseg = gp.BSEG + bcseg_count - 1 ;
        bseg->n = bcseg_count ;
        bseg->bcs[0] = 1;
        bseg->p[0] = 0.0;
        bseg->p[1] = fixedQ;
        bseg->nstart = endseg2;
        bseg->nend = endseg1;

		// no variable bc for now so
		c = "";
		//c = getc(f);
		//if(isdigit(c))
		//	ungetc(c,f);
		//else
		//{
		//while (c != '\n')
		//	c = getc(f);
		//}
		    
        //if(i < N.boundarysegs - 1)
        //        bseg->nextbseg = bseg + 1;
        //else
                bseg->nextbseg = bseg + 1 ;	
	free(bndsegs);
	}
	bseg->nextbseg = NULL;
	set_bvalues();

	free(xcord);
	free(ycord);
	free(bndnodex);
	free(bndnodey);
	free(element);
	free(elevation);
	free(roughness);
	free(edgenum);
	free(edgenode1);
	free(edgenode2);
	free(edgenum_s);
	free(edgenode1_s);
	free(edgenode2_s);	

	if(hotstart) {
		free(depth);
		free(vx);
		free(vy);
	}

		return(0);

}

int linearSearch( const int array[], int key, int size )
{
   int n;
   for ( n = 0; n < size; ++n ) {
      if ( array[ n ] == key ) { 
         return n;
      } 
   } 
   return -1;

}
/*
Return the clockwise status of a curve, clockwise or counterclockwise
n vertices making up curve p
return 0 for incomputables eg: colinear points
CLOCKWISE == 1
COUNTERCLOCKWISE == -1
It is assumed that
- the polygon is closed
- the last point is not repeated.
- the polygon is simple (does not intersect itself or have holes)
*/
int ClockWise(double *x, double *y, int n)
{
int i,j,k;
int count = 0;
double z;
if (n < 3)
return(0);
for (i=0;i<n;i++) {
j = (i + 1) % n;
k = (i + 2) % n;
z = (x[j] - x[i]) * (y[k] - y[j]);
z -= (y[j] - y[i]) * (x[k] - x[j]);
//z = (p[j].x - p[i].x) * (p[k].y - p[j].y);
//z -= (p[j].y - p[i].y) * (p[k].x - p[j].x);
if (z < 0)
count--;
else if (z > 0)
count++;
}
if (count > 0)
return(1);
else if (count < 0)
return(-1);
else
return(0);
}

int				cgns_output(FID, SolID, nin)
int				FID, SolID, nin;
{
         int             err,i ;
        struct node     *np ;
        struct element  *elp ;
        struct belement *belp ;
        struct  boundaryseg     *bseg;
		double *elevation, *wse, *depth, *vx, *vy, *tbx, *tby;
		double *chezy;
		int *ibc;
	    double Cs;
	    double ks;
	    double tau_bx;
	    double U,V,H;
	    double water_elevation,bed_elevation;
	    double rho;
	    double g;
	    // for cgns output
	    char buffer[16];
		int BID, ZID, SolIndex;
		BID = ZID = 1;


		rho = 1000.;
		g = 9.81;
  


		elevation = (double *) malloc(N.nodes*sizeof(double));
		wse = (double *) malloc(N.nodes*sizeof(double));
		depth = (double *) malloc(N.nodes*sizeof(double));
		vx = (double *) malloc(N.nodes*sizeof(double));
		vy = (double *) malloc(N.nodes*sizeof(double));
		tbx = (double *) malloc(N.nodes*sizeof(double));
		tby = (double *) malloc(N.nodes*sizeof(double));
		ibc = (int *) malloc(N.nodes*sizeof(int));
		chezy = (double *) malloc(N.nodes*sizeof(double));

		np=gp.N;
		for(i = 0; i < N.nodes; i++) {
			ks = np->p[1];
			chezy[i] = ks;

			U=np->ud[1];
			vx[i] = U;

			V=np->ud[2];
			vy[i] = V;

			H=np->ud[0];
			depth[i] = H;

			bed_elevation=np->p[0];
			elevation[i] = bed_elevation;

			wse[i] = bed_elevation+H;

			Cs=5.75*log (12* H / ks);
			if(H>0)
			{
				ibc[i] = 1;
				chezy[i] = Cs;
				tbx[i]=sqrt(U*U +V*V) * U / (Cs *Cs) * rho ; //calculates the bed shear stress;
				tby[i]=sqrt(U*U +V*V) * V / (Cs *Cs) * rho ; //calculates the bed shear stress;
			} else {
				ibc[i] = 0;
				chezy[i] = 0.0;
				tbx[i]=0.0; 
				tby[i]=0.0 ; 
			}
			np = np->nextnp;
		}
    
		err = cg_field_write(FID,BID,ZID,SolID,
		 RealDouble,"Elevation",elevation,&SolIndex);
		err = cg_field_write(FID,BID,ZID,SolID,
		 RealDouble,"WaterSurfaceElevation",wse,&SolIndex);
		err = cg_field_write(FID,BID,ZID,SolID,
		 Integer,"IBC",ibc,&SolIndex);
		err = cg_field_write(FID,BID,ZID,SolID,
		 RealDouble,"Depth",depth,&SolIndex);
		err = cg_field_write(FID,BID,ZID,SolID,
		 RealDouble,"VelocityX",vx,&SolIndex);
		err = cg_field_write(FID,BID,ZID,SolID,
		 RealDouble,"VelocityY",vy,&SolIndex);
		err = cg_field_write(FID,BID,ZID,SolID,
		 RealDouble,"Chezy",chezy,&SolIndex);

		free(elevation);
		free(wse);
		free(depth);
		free(vx);
		free(vy);
		free(tbx);
		free(tby);
		free(ibc);
		free(chezy);

}