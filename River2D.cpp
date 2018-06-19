// River2DWin.cpp : Defines the entry point for the console application.
//

#include <tchar.h>
#include <string>
#include "src\Tin.h"
#include "src\HabitatTIN.h"
#include "src\Shallow.h"
#include "src\Outcdg2d.h"
#include "src\Fe_PS.h"
#include "iriclib.h"
#include "cgnslib.h"



extern "C" struct transient tvals; //call the tvals struct in the Fe files
extern "C" int steady(int a, int b, double tol, double tMax, FILE* lfp, int i);///call the steady function in the Fe files
extern "C" int steadynew(int a, int b, double tol, double tMax, FILE* lfp, int i, int m, int k, double gmrestol);///call the steady function in the Fe files
extern "C" int set_ptou(int num);//calls the set_ptou in the Fe files (initializes the velocities?
extern "C" void updateVelocities();
extern "C" int assembleFlag;
extern "C" int input(FILE* CDGFile);
extern "C" int cgns_input(int fid);
extern "C" int init_gp(); //the initialization of all the variables of the Finite element method
extern "C" double test_outflow(int code);//declares the inflow and outlow values
extern "C" double uchange;//declares the change in velocity or Solution change
extern "C" int tecplot_output(int nin, char s[63]); //what this does?
extern "C" int cgns_output(int FID, int SolID, int nin); //what this does?
extern "C" int output(int nin); //what this does?

int main(int argc, char* argv[])
{
	HabitatTIN *meshP;
	int FID, error;
	int solvType;
	int m, k;
	double gmrestol;
	//if (argc != 2) {
	//    fprintf (stderr, "open_cgns CGNSfile\n");
	//    exit (1);
	//}

	/************************************************************************
	RUN THE STEADY SIMULATION WITH A DIRECT SOLVER
	River2D hast 2 solver options: Direct Solver and Iterative Solver
	************************************************************************/
	//For this time I'm using the file from the tutorial
	init_gp();
	int tsnum; //number of step
			   // int i=0;
	int err; // residual, solution change or error?
	double tolerance = 0.0001; // tolerance 
	double dtMaximum = 100; //Maximum time increment
	FILE* logFile; //name of the log file
	FILE* CDGFile;//name of the CDG file
	string CDGFileName;
	char OutputFileName[63];
	int WriteTimeStep;
	int itcount;
	static char** namearray = 0;
	static char* solpointers;
	char solname[32];
	char sn[10000][32];
	int namelen;
	double *tstep;
	char   *sol_names = NULL;

	// for cgns output
	char buffer[32];
	int BID, ZID, SolID, SolIndex;
	BID = ZID = 1;
	itcount = 0;



	logFile = fopen("logfile.txt", "a");//open the log file for writing

	if (argc > 1) {
		char *pFilename = argv[1];
		error = cg_open(pFilename, CG_MODE_MODIFY, &FID);
		cg_error_print();
		cgns_input(FID);
	}
	else {
		cout << "Enter the name of the CDG file: ";
		cin >> CDGFileName;
		CDGFile = fopen(CDGFileName.c_str(), "r");//open the CDG file for reading
		input(CDGFile);//inputs the data in the CDGFile in the structures of the Feio.c source file
	}
	//  fclose(CDGFile); this causes an error of memory allocation

	//cout<<"Enter the name of the output file:";
	//scanf("%s",OutputFileName);

	//cout<<"Every how many iterations do you want to output to the file: ";
	//cin>>WriteTimeStep;

	set_ptou(0);/* function to set old values (uo, uoo) */ //copies the values at the nodes of the current time step to the space allocated for the next time space
	updateVelocities();//	April, 2003  -  P. Steffler
					   //	Calculates nodal velocities on the basis of a projection of the 
					   //	discharge intensity and depth distribution as an alternative
					   //	to velocities based simply on nodal qx, qy, and H.
					   //	Also performs a projection on depth which has a smoothing effect.

	err = cg_iRIC_Read_Real("m_MaxTimedt", &dtMaximum);
	err = cg_iRIC_Read_Real("m_GoalSolChange", &tolerance);
	err = cg_iRIC_Read_Integer("m_PlotInc", &WriteTimeStep);
	err = cg_iRIC_Read_Real("m_FinalTime", &tvals.tfinal);
	err = cg_iRIC_Read_Integer("m_SolverType", &solvType);
	err = cg_iRIC_Read_Integer("m_NumStepsBfrRestart", &m);
	err = cg_iRIC_Read_Integer("m_MaxNumIter", &k);
	err = cg_iRIC_Read_Real("m_ConvergenceTol", &gmrestol);

	int numIter = (tvals.tfinal - tvals.t / WriteTimeStep);
	tstep = new double[numIter];

	solpointers = (char*)malloc(sizeof(char) * (numIter+1) * 32);
	for (int j = 0; j <32 * numIter; j++)
	{
		solpointers[j] = ' ';
	}
	assembleFlag = 0;

	cout << endl;
	cout << "***********STARTING THE SIMULATION**************" << endl;
	//cout<<"The initial time is: "<<tvals.t<<endl;
	//cout<<"Enter the final time for the simulation: ";
	//cin>>tvals.tfinal;

	tsnum = 0;
	while (tvals.t <tvals.tfinal)
	{
		++itcount;
		for (int i = 0; i<WriteTimeStep; i++)
			//i++;
		{
			tsnum = tsnum++;
			if (solvType == 1) {
				err = steady(-1, 2, tolerance, dtMaximum, logFile, tsnum);
			}
			else {
				err = steadynew(-1, 2, tolerance, dtMaximum, logFile, tsnum, m, k, gmrestol);
			}
			cout << endl << "Number of Iteration: " << tsnum << endl <<
				"Current Time: " << tvals.t << endl <<
				"Time Difference: " << tvals.dt << endl;
			cout << " Solution Change: " << uchange << endl;
			cout << " Total Inflow: " << test_outflow(1) << endl;
			cout << " Total Outflow: " << test_outflow(2) << endl;
			cout << "------------------------------------------------------------------------" << endl << endl;
			if (err<0) cout << "Memory Allocation Failure" << endl;
		}

		sprintf(buffer, "FlowSolution%d", itcount);
		err = cg_sol_write(FID, BID, ZID, buffer, Vertex, &SolID);

		tstep[SolID - 1] = tvals.t;

		err = cg_biter_write(FID, BID, "BaseIterativeData", SolID);
		err = cg_goto(FID, BID, "BaseIterativeData_t", 1, "end");
		//note made change below for iric v3.  cgsize_t was int
		cgsize_t nuse[1];
		//int nuse[1];
		nuse[0] = SolID;
		err = cg_array_write("TimeValues", RealDouble, 1, nuse, tstep);

		cgsize_t idata[2];
		idata[0] = 32;
		idata[1] = SolID;

		sol_names = (char *)malloc(idata[0] * idata[1] * sizeof(char));

		for (int j = 0; j < idata[0] * idata[1]; j++)
		{
			sol_names[j] = ' ';
		}
		for (int si = 0; si < SolID; si++) {
			sprintf(solname, "FlowSolution%d", si + 1);
			namelen = strlen(solname);
			strncpy(sol_names + si * 32, solname, namelen);
		}

		cgns_output(FID, SolID, itcount);
		err = cg_goto(FID, BID, "Zone_t", ZID, "ZoneIterativeData_t", 1, "end");
		
		err = cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, idata, sol_names);
		
		//cg_error_print();
		free(sol_names);
	}
	cout << "************END OF SIMULATION*****************" << endl;
	err = cg_close(FID);
	//cg_error_print();
	return 0;
}

