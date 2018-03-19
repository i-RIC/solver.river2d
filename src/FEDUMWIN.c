//#pragma options cpluscmt	
#include "Fe_PS.h"
   
int	scr;
unsigned long	fg, bg;

char name[] = {"Two Dimensional Free Surface Flow"};

extern struct control N;
extern struct pointers gp;
extern struct transient tvals;
extern struct settings	sets;
struct element		*theElement;
struct node		*theNode;
struct belement		*theBoundEl;

double	winfactor, horigin, hscale, vorigin, vscale, v1origin, v1scale;
double	ymin, ymax, xmin, xmax ;
int		wtp, ICOUNT;
double	y2min, y2max, v2origin, v2scale,scale,h2origin,h2scale;
double	h3origin, v3origin, h3scale, v3scale, gama, beta, varscale, xo, yo;

initXwindow(argc,argv)
int	argc;
char	**argv;
{
	return(0);
}

killXwindow()

{
	return(0);
}


ScaleUXwindow(scalecode)
int	scalecode;

{
	return(0);
}

ScaleUIXwindow(scalecode)
int	scalecode;

{
	return(0);
}

ScaleXwindow(scalecode,crdnum)
int	scalecode,crdnum;

{
	return(0);
}

PlotNodeVels(code)
int	code;

{
}

PlotNodeVels_trans(code)
int	code;

{
}

DrawBElement(belp)

struct belement *belp;

{
	return(0);
}


DrawElement(elmntp)
struct element	*elmntp;

{
	return(0);
}

DrawUBElement(belp)

struct tempbelement *belp;

{
	return(0);
}


DrawUElement(elmntp)
struct tempelement	*elmntp;

{
	return(0);
}


DrawMesh(code)
int code;
{
}

DrawUMesh(code)
int code;
{
}

DrawUIMesh(code)
int code;
{
}

PlotNodes(code)
int	code;

{
}

PlotUNodes(code)
int	code;

{
}

DrawUNodeNum(code)
int 	code;

{
}

DrawNodeNum(code)
int 	code;

{
}


DrawNode(np)
struct node *np;

{
}

DrawUNode(np)
struct unode *np;

{
	return(0);
}

DrawDoubleNum(d,x,y)
double	d;
short	x, y;

{
	return(0);
}

DrawIntNum(d,x,y)

int d;
short x,y;

{
	return(0);
}


DrawContour(code)
int	code;
{
}

DrawContour_trans(code)
int	code;
{
}

	
DrawElementContour(elp, nvar, niv, cval)

struct	element *elp;
int	nvar, niv;
double	cval;
{
	return(0);
}

DrawElementContour1(elp, nvar, niv, cval, ncoord, ndraw)
struct	element *elp;
int	nvar, niv, ncoord, ndraw;
double	cval;
{
	return(0);
}


DrawXsec(code)
int	code;
{
}

DrawXsec_trans(code)
int	code;
{
}	
			

Draw3D(scalecode)

{
	return(0);
}		 

Draw3D_trans(scalecode)

{
	return(0);
}

	
DrawElementContour3(elp, nvar, niv, cval, ncoord, ndraw)
struct	element *elp;
int	nvar, niv, ncoord, ndraw;
double	cval;
{
	return(0);
}

DrawBounds(code)
int code;
{
}

DrawBSegment(bsp)

struct bsegment *bsp;

{
	return(0);
}


SetAuxPlot(n)
int		n;

{ 
	return(0);
}
