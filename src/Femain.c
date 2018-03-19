//#pragma options cpluscmt	
#include	"Fe_PS.h"

#if defined(__MWERKS__)
#include <Events.h>
#include <OSEvents.h>
#include <OSUtils.h>
#define SLICE 1
EventRecord gmyEvent;
long lastTick;
extern Boolean SIOUXQuitting; 

#endif

/*
main(argc,argv)
int	argc;
char	**argv;
{
	char	ain ;
	int	nin, errcode, resfilenum, i ;
	double x, y;

	initXwindow(argc,argv);	
	init_gp() ;
	printf(" CDG2D Depth Averaged Hydrodynamic Simulation\n\n");
	printf(" Development Version, Intended for research and evaluation purposes\n");
	printf(" Copyright: A. Ghanem and P. Steffler\n\n");
	printf(" Input command letter and integer code, for example:\n\n");
	printf(" i 1 = input a mesh (.cdg file),\n");
	printf(" o 1 = output a mesh (.cdg file),\n");
	printf(" s 0 = run to steady state (fixed simulation time actually),\n");
	printf(" d -1 = run transient simulation,\n");
	printf(" q 0 = quit application. \n\n") ;
	ain = 't' ;
	
	while(ain != 'q'){
		printf("> ");
		scanf(" %c %d",&ain,&nin) ;
		switch (ain) {
			
			case 't' :
				errcode = test_2Delm(nin) ;
		
				break ;
			case 'i' :
				errcode = input(nin) ;
				ScaleXwindow(0); 
				break ;
			case 'R' :
				ScaleXwindow(nin);
				break;
			case 'Q' :
				SetAuxPlot(nin);
				ScaleXwindow(2);
				break;
			case 'n' :
				errcode = test_mesh(nin) ;
				break ;
			case 'm' :
				errcode = test_map(nin) ;
				break ;
			case 'l' :
				errcode = list_vars(nin) ;
				break ;
			case 'o' :
				errcode = output(nin) ;
				break ;
			case 'e' :
				errcode = test_elK(nin,-1) ;// -1 to use 2 vars at once 
				break ;
			case 'g' :
				errcode = test_belK(nin,-1) ;
				break ;
			case 'G' :
				errcode = test_globalK(nin,-1) ;
				break ;
			case 'A' :
				errcode = test_globalK1(nin,-1) ;
				break ;
			case 'B' :
//				errcode = switch_2D(nin,-1) ;
				errcode = matrixform(nin);
				break ;

			case 'j' :
//				errcode = find_neighbour() ;
				break ;
			case 'f' :
				errcode = assemble(nin,0,0,2,0) ;
				//errcode = solve(nin,1,0) ;
				
				//errcode = get_L2norm(nin) ;
				break ;
			case 'x' :
				for(i=0;i<10;i++) {
					printf(" iteration # %d\n",i+1) ;
					update_eta() ;
					assemble(0,0,0,2) ; solve(0,1,0) ;
					set_etabc() ;
					assemble(1,0,0,2) ;
					solve(1,1,0) ;
				}
				break ;
			case 'a' :
				//errcode = assemble(nin,0,0,1) ;
				errcode = assemble(nin,0,0,2,0) ;// FAYE 88-08-03 
				printf(" assemble code = %d\n",errcode);
				break ;
			case 's' :
				errcode = test_steady(nin) ;
				break ;
			case 'p' :
				errcode = perimeter(nin) ;
				break ;
			case 'd' :
				errcode = test_trans(nin) ;
				break ;
			case 'b' :
				errcode = test_bvalues(nin) ;
				break ;
         case 'S' :
            errcode = smooth(nin) ;
            break ;
         case 'T' :
            errcode = test_time(nin) ;
           	break ; 
         case 'c' :
				errcode = test_cons(nin);
            break ;
			case 'M' : 
				errcode = DrawMesh(nin);
				break;
			case 'E' : 
				errcode = DrawUMesh(nin);
				break;
			case 'W' : 
				errcode = DrawUIMesh(nin);
				break;

			case 'N' :
				errcode = PlotNodes(nin);
				break;
			case 'U' :
				errcode = PlotUNodes(nin);
				break;

			case 'X' :
				errcode = DrawXsec(nin);
				break;	

			case 'D' :
				errcode = Draw3D(nin) ;
				break ;
			case 'C' :
				errcode = DrawContour(nin);
				break ;
			case 'V' :
				errcode = PlotNodeVels(nin);
				break;
			case 'z' :
				errcode = set_param(nin) ;
				break ;
			case 'u' :
				errcode = set_ptou(nin) ;
				break ;
			case 'y' :
				errcode = assemble(1,0,0) ;
				errcode = solve(1,1,0) ;
				errcode = fix_psi(1) ;
				errcode = assemble(0,0,0) ;
				errcode = solve(0,1,0) ;
				break ;
			case 'v' :
				errcode = Velocities(nin) ;
				break ;
			case 'P'	:
				errcode = test_UMesh(nin) ;
				break ;
			case 'H'	:
				errcode = make_elements(nin) ;
				break ;
			case 'F'	:
				errcode = free_telements() ;
				break ;
			case 'I'	:
//				errcode = redo_frontsegs() ;
//				errcode = input_itriangles();
				errcode = find_nodalvalues();
				break ;
			case 'O'  :
				errcode = ordermesh(nin);
				break;
			case 'J'  :
				errcode = get_IFG4xsecs(nin);
				break;
			case 'K'  :
				errcode = get_xsecvel(nin);
				break;
			case 'Y'  :
				errcode = get_IFG4dat(nin);
				break;

			case 'Z'  :
				errcode = DrawBounds(nin);
				break;
			case 'L'  :
//				errcode = test_mixedgps(nin);
//				errcode = test_boundgps(nin);
				errcode = test_contour(nin);
				break;
			case 'h'  :
				errcode = Nodal_Values(nin);
				break;
			case 'k'  :
				errcode = outxsec(nin);
				break;

			
		}
	}
	killXwindow();
	exit(0);
}
*/

FILE	*get_fptr(code)
int		code ;

{
	FILE	*fp, *fopen() ;
	char	s[63] ;
	
	if(code == 1) {
		printf("\n Input Data File Name  ") ;
		scanf(" %s",s) ;
		fp = fopen(s,"r") ;
	}
	else
		fp = stdin ;
	return( fp ) ;
}

FILE	*put_fptr(code)
int		code ;

{
	FILE	*fp, *fopen() ;
	char	s[63] ;
	
	if(code == 1) {
	  printf("\n Input Output File Name  ") ;
		scanf(" %s",s) ;
		fp = fopen(s,"a") ;//The "w" attribute is changed to "a" by C.Frias June 16 2010
	}
	else {
		if(code == 2)
			fp = fopen("fe.out","a") ;//The "w" attribute is changed to "a" by C.Frias June 16 2010
		else
			fp = stdin ;
	}
	return( fp ) ;
}

FILE	*log_fptr()
{
	FILE	*fp, *fopen() ;
	char	s[63] ;
	
	printf("\n Input Log File Name  ") ;
	scanf(" %s",s) ;
	fp = fopen(s,"a") ;
	return( fp ) ;
}

int		StartMes(astr,bstr,cstr,dstr)
char		*astr, *bstr, *cstr, *dstr ;

{
 
	printf("%s\n%s\n%s\n%s\n",astr,bstr,cstr,dstr) ;
 
	return(0);
}

int		UpDateMes(astr,bstr, cstr, dstr)
char		*astr, *bstr, *cstr, *dstr ;

{
	printf("%s\n",dstr) ;
	return(0);
}

int		EndMes()

{
	return(0);
}

void DoMultiTask(long sTime)
{
#if defined(__MWERKS__)
	long currentTick;
	
	currentTick = LMGetTicks();
	if((currentTick - lastTick) > SLICE){
		if(WaitNextEvent(everyEvent,&gmyEvent, sTime, NULL)) {
			/* put any event handling here */
		}
	lastTick = LMGetTicks();
	}
#endif
}


