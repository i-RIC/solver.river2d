
//			Habitat.h


//		Habitat is a subclaas of Shallow is a subclass of Physics which describes
//		river bed topography description.
//		The nodes produced are of the Habitat
//		Node class with 3 coordinates plus a bed roughness height.

//		Initiated:	Jan 11, 1997
//
//		Last Revised: Jan 11, 1997

#ifndef HABITAT_H
#define HABITAT_H

#include "Shallow.h"
#include "Intrplst.h"
#include <fstream>
using namespace std;

enum WUAmode {
	product,
	geomean,
	minimum
};

enum interpCImode {
	continuous,
	discrete
};


class Habitat : public Shallow
{
	protected:
		InterpLst *dPref;
		InterpLst *vPref;
		InterpLst *sPref;
		WUAmode calcWUA;	//	how to calc csi
		interpCImode interpCI;		//	how to interpolate channel index
		double sg_ice;

	public:
		Habitat();
		void setCImode(interpCImode newMode) {interpCI = newMode;}
		interpCImode getCImode() {return interpCI;}
		void setWUAmode(WUAmode newMode) {calcWUA = newMode;}
		WUAmode getWUAcalc() {return calcWUA;}
		virtual Node* makeNewNode(int n=1,double x= 100.0,
				double y= 100.0,double z= 100.0);
		virtual Node* readInNode(istream& is);
		virtual int nVars() {return 3;}
		virtual int nParams() {return 7;}
		virtual char* typeName() {return "Habitat";}
		void loadPrefs(ifstream &f);
		void calcSI(Node * hNP);
		void setSg_ice(double specific_grav_ice){sg_ice = specific_grav_ice;}
		double getSg_ice(){return sg_ice;}

};


class HabitatNode : public ShallowNode
{
	protected:
		double V;			//	Velocity magnitude
		double Fr;			//	Froude Number
		double ci;			//	channel index
		double dsi;			//	depth suitability index
		double vsi;			//	velocity suitability index
		double ssi;			//	substrate suitability index
		double csi;			//	combined suitability index
		double WUA;			//  weighted useable area in vicinity
		double psi;			//  stream function
		double tice;		//  ice thickness
		double kice;		//  ice roughness
		double vorticity;	//	vorticity
		Habitat* theHabitat; //	physics object which governs behaviour

	public:
		void getParName(int i, char* parName);
		HabitatNode(	
					Habitat* hab,
					int nm=1,
					double xc=1000.0,
					double yc=1000.0,
					double zc=100.0,
					double rf = 0.05,
					double chi = 10.0,
					double depth = 1.0,
					double xDis = 0.0,
					double yDis = 0.0
					);
		HabitatNode(Habitat* hab,istream& is);
		void setVandF();
		void setCi(double chIndex) {ci = chIndex;}
		double getCi() const {return ci;}
		void setVorticity(double vVal) {vorticity = vVal;}
		double getVorticity() const {return vorticity;}
		void setSI(double d, double v, double s, double c);
		void setWUA(double wua);
		void incWUA(double wua);
		virtual int getNParams() {return 3;}
		virtual int getNVars() {return 3;}
		virtual double getPar(int i);
		virtual double getVar(int i);
		virtual void interp(int n, Node** nPtrs, double* wts);
		void setPsi(double strmFunc) {psi = strmFunc;}
		double getPsi() const {return psi;}
		void setTice(double thickness) {tice = thickness;}
		double getTice() const {return tice;}
		void setKice(double rough) {kice = rough;}
		double getKice() const {return kice;}
};

#endif

