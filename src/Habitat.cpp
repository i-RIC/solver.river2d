
//		Habitat.cpp

#include "Habitat.h"
#include "math.h"
#include "string.h"



HabitatNode::HabitatNode(Habitat* hab,int nm, double xc, double yc, double zc, double rf,
						double chi,	double depth, double xDis, double yDis)
				:ShallowNode(nm,xc,yc,zc)
{
//	setn(nm);
//	assign(xc, yc, zc);
	theHabitat = hab;
	ks = rf;
	d = depth;
	qx = xDis;
	qy = yDis;
	ci = chi;
	dsi = 0.0;
	vsi = 0.0;
	ssi = 0.0;
	csi = 0.0;
	WUA = 0.0;
	psi = 0.0;
	tice = 0.0;
	kice = 0.0;
	vorticity = 0.0;
	setVandF();
}

void HabitatNode::setVandF()
{
	if((d-theHabitat->getSg_ice()*tice) >= theHabitat->getMinDepth()){
		V = sqrt(qx*qx + qy*qy)/(d-theHabitat->getSg_ice()*tice);
		Fr = V / sqrt(G*(d-theHabitat->getSg_ice()*tice));
	}
	else {
		V = 0.0;
		Fr = 0.0;
	}
}

HabitatNode::HabitatNode(Habitat* hab,istream& is)
				:ShallowNode(is)
{
	theHabitat = hab;
/*	int num = -990;
	double xcoord=0.0, ycoord=0.0, zcoord=0.0;
	double rough = 0.05;
	while(is) {
		if(!(is >> num)){
			num = -991;
			break;
		}
		if(!(is >> xcoord)){
			num = -992;
			break;
		}
		if(!(is >> ycoord)){
			num = -993;
			break;
		}
		if(!(is >> zcoord)){
			num = -994;
			break;
		}
		if(!(is >> rough)){
			num = -995;
			break;
		}
		break;
	}

	is >> num;
	is >> xcoord;
	is >> ycoord;
	is >> zcoord;
	is >> rough;
*/
//	setn(num);
//	assign(xcoord, ycoord, zcoord);
//	ks = rough;
//	d = 1.0;
//	qx = 0.0;
//	qy = 0.0;
	V = 0.0;
	Fr = 0.0;
	ci = ks;
	dsi = 0.0;
	vsi = 0.0;
	ssi = 0.0;
	csi = 0.0;
	WUA = 0.0;
	psi = 0.0;
	tice = 0.0;
	kice = 0.0;
	vorticity = 0.0;
}

void HabitatNode::setSI(double d, double v, double s, double c)
{
	dsi = d;
	vsi = v;
	ssi = s;
	csi = c;
}

void HabitatNode::setWUA(double wua)
{
	WUA = wua;
}

void HabitatNode::incWUA(double wua)
{
	WUA += wua;
}

double HabitatNode::getPar(int i)
{
	switch(i){
		case 0 :
			return z;
		case 1 :
			return z;
		case 2 :
			return ks;
		case 3 :
			return (d-(theHabitat->getSg_ice()*tice));
		case 4 :
			return qx;
		case 5 :
			return qy;
		case 6 :
			return z+d;     // water surface elevation
		case 7 :
			if((d-theHabitat->getSg_ice()*tice) > theHabitat->getMinDepth())
				return V;  // Velocity
			else
				return 0.0;
		case 8 :			
			if((d-theHabitat->getSg_ice()*tice) > theHabitat->getMinDepth())
				return Fr;  //  Froude Number
			else
				return 0.0;
		case 9 :
			return ci;
		case 10 :
			return dsi;
		case 11 :
			return vsi;
		case 12 :
			return ssi;
		case 13 :
			return csi;
		case 14 :
			return WUA;
		case 15:
			return psi;
		case 16:
			if((d-theHabitat->getSg_ice()*tice) > theHabitat->getMinDepth())
			{
				if(((d-theHabitat->getSg_ice()*tice)/ks) > (2.718/12))
					return V/(2.5*log(12.0*(d-theHabitat->getSg_ice()*tice)/ks));
				else
					return V/((30.0*(d-theHabitat->getSg_ice()*tice))/(2.718*ks));
			}
			else
				return 0.0;
		case 17:
			return tice;
		case 18:
			return kice;
		case 19:
			return vorticity;
	}
	return 0.0;
}

double HabitatNode::getVar(int i)
{
	switch(i){
		case 1 :
			return d;
		case 2 :
			return qx;
		case 3 :
			return qy;
	}
	return 0.0;
}

void HabitatNode::getParName(int i, char* parName)
{
	switch(i){
		case 0 :
			strcpy(parName,"Bed Elevation");
			break;
		case 1 :
			strcpy(parName,"Bed Elevation");
			break;
		case 2 :
			strcpy(parName,"Bed Roughness");
			break;
		case 3 :
			strcpy(parName,"Depth");
			break;
		case 4 :
			strcpy(parName,"qx");
			break;
		case 5 :
			strcpy(parName,"qy");
			break;
		case 6 :
			strcpy(parName,"Water Surface Elev");
			break;
		case 7 :
			strcpy(parName,"Velocity");
			break;
		case 8 :
			strcpy(parName,"Froude #");
			break;
		case 9 :
			strcpy(parName,"Channel Index");
			break;
		case 10 :
			strcpy(parName,"Depth Suitability");
			break;
		case 11 :
			strcpy(parName,"Velocity Suitability");
			break;
		case 12 :
			strcpy(parName,"Substrate Suitability");
			break;
		case 13 :
			strcpy(parName,"Combined Suitability");
			break;
		case 14 :
			strcpy(parName,"Weighted Useable Area");
			break;
		case 15:
			strcpy(parName,"Cumulative Discharge");
			break;
		case 16:
			strcpy(parName,"Shear Velocity Magnitude");
			break;
		case 17:
			strcpy(parName,"Ice Thickness");
			break;
		case 18:
			strcpy(parName,"Ice Roughness");
			break;
		case 19:
			strcpy(parName,"Vorticity");

	}
}

void HabitatNode::interp(int n, Node** nPtrs, double* wts)
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
	ks = 0.0;
	ci = 0.0;
	d = 0.0;
	qx = 0.0;
	qy = 0.0;
	V = 0.0;
	Fr = 0.0;
	ci = 0.0;
	dsi = 0.0;
	vsi = 0.0;
	ssi = 0.0;
	csi = 0.0;
	psi = 0.0;
	tice = 0.0;
	kice = 0.0;
	vorticity = 0.0;
	for(int i=0;i<n;i++) {
		if(nPtrs[i] != NULL) {
			x += wts[i] * nPtrs[i]->getXo();
			y += wts[i] * nPtrs[i]->getYo();
			z += wts[i] * nPtrs[i]->getZc();
			ks += wts[i] * nPtrs[i]->getPar(2);
			ci += wts[i] * nPtrs[i]->getPar(9);
			d += wts[i] * nPtrs[i]->getPar(3);
			qx += wts[i] * nPtrs[i]->getPar(4);
			qy += wts[i] * nPtrs[i]->getPar(5);
			V += wts[i] * nPtrs[i]->getPar(7);
			Fr += wts[i] * nPtrs[i]->getPar(8);
			dsi += wts[i] * nPtrs[i]->getPar(10);
			vsi += wts[i] * nPtrs[i]->getPar(11);
			ssi += wts[i] * nPtrs[i]->getPar(12);
			csi += wts[i] * nPtrs[i]->getPar(13);
			psi += wts[i] * nPtrs[i]->getPar(15);
			tice +=wts[i] * nPtrs[i]->getPar(17);
			kice +=wts[i] * nPtrs[i]->getPar(18);
			vorticity +=wts[i] * nPtrs[i]->getPar(19);
		}
	}
	d = d + theHabitat->getSg_ice()*tice;
	if(theHabitat->getCImode() == discrete){
		ssi = nPtrs[0]->getPar(12);
		ci = nPtrs[0]->getPar(9);
		double maxwt = wts[0];
		for(int ii=1;ii<n;ii++){
			if(wts[ii] > maxwt){
				ssi = nPtrs[ii]->getPar(12);
				ci = nPtrs[ii]->getPar(9);
				maxwt = wts[ii];
			}
		}
	}
	saveLocation();

}

Habitat::Habitat() : Shallow()
{
	dPref = NULL;
	vPref = NULL;
	sPref = NULL;
	calcWUA = product;
	interpCI = continuous;
	sg_ice = 0.92;
}

Node* Habitat::makeNewNode(int n, double x, double y, double z)
{
	return new HabitatNode(this,n,x,y,z,0.1);
}

Node* Habitat::readInNode(istream& is)
{
	return new HabitatNode(this,is);
}

void Habitat::loadPrefs(ifstream &f)
{
	delete vPref;
	delete dPref;
	delete sPref;

	vPref = new InterpLst(f);
	dPref = new InterpLst(f);
	sPref = new InterpLst(f);
}

void Habitat::calcSI(Node *nP)
{
	HabitatNode *hNP = (HabitatNode *) nP;
	double csi;
	double dsi = dPref->value(hNP->getPar(3));
	double vsi = vPref->value(hNP->getPar(7));
	double ssi = sPref->value(hNP->getPar(9));

	csi = dsi * vsi * ssi;
	if(calcWUA == geomean)
		csi = pow(csi,0.33333333);
	if(calcWUA == minimum){
		csi = dsi;
		csi = (vsi < csi)? vsi : csi;
		csi = (ssi < csi)? ssi : csi;
	}
	hNP->setSI(dsi,vsi,ssi,csi);
	hNP->setWUA(0.0);
}

