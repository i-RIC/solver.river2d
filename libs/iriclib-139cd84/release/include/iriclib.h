#include "iriclib_global.h"
#include <cgnslib.h>

#define IRIC_GEO_UNKNOWN        0
#define IRIC_GEO_POLYGON        1
#define IRIC_GEO_RIVERSURVEY    2
#define IRIC_GEO_POINTMAP       3

#define IRIC_OPTION_UNKNOWN         0
#define IRIC_OPTION_CANCEL          1
#define IRIC_OPTION_DIVIDESOLUTIONS 2
#define IRIC_OPTION_STDSOLUTION     3

#define IRIC_LOCKED 1

#define IRIC_CANCELED 1

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************/
/* Initialization                             */
/**********************************************/

int IRICLIBDLL cg_iRIC_Init(int fid);

int IRICLIBDLL cg_iRIC_InitRead_Base(int fid, const char* basename);

int IRICLIBDLL cg_iRIC_InitRead(int fid);

void IRICLIBDLL cg_iRIC_SetFilename(int fid, const char* fname);

int IRICLIBDLL iRIC_InitOption(int option);

int IRICLIBDLL cg_iRIC_Flush(const char* name, int* fid);

int IRICLIBDLL cg_iRIC_GotoBase(int fid, int* B);

int IRICLIBDLL cg_iRIC_GotoCC(int fid);

int IRICLIBDLL cg_iRIC_GotoRawDataTop(int fid);

int IRICLIBDLL cg_iRIC_Set_ZoneId_Mul(int fid, int zoneid);

/**********************************************/
/* iRIC GUI coorporation                      */
/**********************************************/

int IRICLIBDLL iRIC_Write_Sol_Start(const char* fname);

int IRICLIBDLL iRIC_Write_Sol_End(const char* fname);

int IRICLIBDLL iRIC_Check_Lock(const char* fname);

int IRICLIBDLL iRIC_Check_Cancel();

/**********************************************/
/* Reading Calculation Conditions             */
/**********************************************/

int IRICLIBDLL cg_iRIC_Read_Integer_Mul(int fid, const char* name, int* intvalue);

int IRICLIBDLL cg_iRIC_Read_Real_Mul(int fid, const char* name, double* realvalue);

int IRICLIBDLL cg_iRIC_Read_RealSingle_Mul(int fid, const char* name, float* realvalue);

int IRICLIBDLL cg_iRIC_Read_StringLen_Mul(int fid, const char* name, int* length);

int IRICLIBDLL cg_iRIC_Read_String_Mul(int fid, const char* name, char* strvalue);

int IRICLIBDLL cg_iRIC_Read_FunctionalSize_Mul(int fid, const char* name, cgsize_t* size);

int IRICLIBDLL cg_iRIC_Read_Functional_Mul(int fid, const char* name, double* x, double* y);

int IRICLIBDLL cg_iRIC_Read_FunctionalWithName_Mul(int fid, const char* name, const char* paramname, double* data);

int IRICLIBDLL cg_iRIC_Read_Functional_RealSingle_Mul(int fid, const char* name, float* x, float* y);

int IRICLIBDLL cg_iRIC_Read_FunctionalWithName_RealSingle_Mul(int fid, const char* name, const char* paramname, float* data);

/**********************************************/
/* Writing Calculation Conditions             */
/**********************************************/

int IRICLIBDLL cg_iRIC_Write_Integer_Mul(int fid, const char* name, int intvalue);

int IRICLIBDLL cg_iRIC_Write_Real_Mul(int fid, const char* name, double realvalue);

int IRICLIBDLL cg_iRIC_Write_String_Mul(int fid, const char* name, const char* strvalue);

int IRICLIBDLL cg_iRIC_Write_Functional_Mul(int fid, const char* name, int length, double* realarray_x, double* realarray_y);

int IRICLIBDLL cg_iRIC_Write_FunctionalWithName_Mul(int fid, const char* name, const char* paramname, int length, double* data);

/**********************************************/
/* Reading Complex Calculation Conditions     */
/**********************************************/

int IRICLIBDLL cg_iRIC_Read_Complex_Count_Mul(int fid, const char* groupname, int* num);

int IRICLIBDLL cg_iRIC_Read_Complex_Integer_Mul(int fid, const char* groupname, int num, const char* name, int* intvalue);

int IRICLIBDLL cg_iRIC_Read_Complex_Real_Mul(int fid, const char* groupname, int num, const char* name, double* realvalue);

int IRICLIBDLL cg_iRIC_Read_Complex_RealSingle_Mul(int fid, const char* groupname, int num, const char* name, float* realvalue);

int IRICLIBDLL cg_iRIC_Read_Complex_StringLen_Mul(int fid, const char* groupname, int num, const char* name, int* length);

int IRICLIBDLL cg_iRIC_Read_Complex_String_Mul(int fid, const char* groupname, int num, const char* name, char* strvalue);

int IRICLIBDLL cg_iRIC_Read_Complex_FunctionalSize_Mul(int fid, const char* groupname, int num, const char* name, cgsize_t* size);

int IRICLIBDLL cg_iRIC_Read_Complex_Functional_Mul(int fid, const char* groupname, int num, const char* name, double* x, double* y);

int IRICLIBDLL cg_iRIC_Read_Complex_FunctionalWithName_Mul(int fid, const char* groupname, int num, const char* name, const char* paramname, double* data);

int IRICLIBDLL cg_iRIC_Read_Complex_Functional_RealSingle_Mul(int fid, const char* groupname, int num, const char* name, float* x, float* y);

int IRICLIBDLL cg_iRIC_Read_Complex_FunctionalWithName_RealSingle_Mul(int fid, const char* groupname, int num, const char* name, const char* paramname, float* data);

int IRICLIBDLL cg_iRIC_Read_Grid_Complex_Node_Mul(int fid, const char* groupname, int* values);

int IRICLIBDLL cg_iRIC_Read_Grid_Complex_Cell_Mul(int fid, const char *groupname, int* values);

/**********************************************/
/* Writing Complex Calculation Conditions     */
/**********************************************/

int IRICLIBDLL cg_iRIC_Clear_Complex_Mul(int fid);

int IRICLIBDLL cg_iRIC_Write_Complex_Integer_Mul(int fid, const char* groupname, int num, const char* name, int intvalue);

int IRICLIBDLL cg_iRIC_Write_Complex_Real_Mul(int fid, const char* groupname, int num, const char* name, double realvalue);

int IRICLIBDLL cg_iRIC_Write_Complex_String_Mul(int fid, const char* groupname, int num, const char* name, const char* strvalue);

int IRICLIBDLL cg_iRIC_Write_Complex_Functional_Mul(int fid, const char* groupname, int num, const char* name, int length, double* realarray_x, double* realarray_y);

int IRICLIBDLL cg_iRIC_Write_Complex_FunctionalWithName_Mul(int fid, const char* groupname, int num, const char* name, const char* paramname, int length, double* data);

int IRICLIBDLL cg_iRIC_Write_Grid_Complex_Node_Mul(int fid, const char* groupname, int* values);

int IRICLIBDLL cg_iRIC_Write_Grid_Complex_Cell_Mul(int fid, const char *groupname, int* values);

/**********************************************/
/* Reading Grids                              */
/**********************************************/

int IRICLIBDLL cg_iRIC_GotoGridCoord2d_Mul(int fid, cgsize_t* isize, cgsize_t* jsize);

int IRICLIBDLL cg_iRIC_GetGridCoord2d_Mul(int fid, double *x, double *y);

int IRICLIBDLL cg_iRIC_GotoGridCoord3d_Mul(int fid, cgsize_t* isize, cgsize_t* jsize, cgsize_t* ksize);

int IRICLIBDLL cg_iRIC_GetGridCoord3d_Mul(int fid, double *x, double *y, double *z);

int IRICLIBDLL cg_iRIC_Read_Grid_Real_Node_Mul(int fid, const char* name, double* values);

int IRICLIBDLL cg_iRIC_Read_Grid_Integer_Node_Mul(int fid, const char* name, int* values);

int IRICLIBDLL cg_iRIC_Read_Grid_Real_Cell_Mul(int fid, const char *name, double* values);

int IRICLIBDLL cg_iRIC_Read_Grid_Integer_Cell_Mul(int fid, const char *name, int* values);

int IRICLIBDLL cg_iRIC_Read_Grid_FunctionalDimensionSize_Mul(int fid, const char* name, const char* dimname, cgsize_t* count);

int IRICLIBDLL cg_iRIC_Read_Grid_FunctionalDimension_Integer_Mul(int fid, const char* name, const char* dimname, int* value);

int IRICLIBDLL cg_iRIC_Read_Grid_FunctionalDimension_Real_Mul(int fid, const char* name, const char* dimname, double* value);

int IRICLIBDLL cg_iRIC_Read_Grid_FunctionalTimeSize_Mul(int fid, const char* name, cgsize_t* count);

int IRICLIBDLL cg_iRIC_Read_Grid_FunctionalTime_Mul(int fid, const char* name, double* time);

int IRICLIBDLL cg_iRIC_Read_Grid_Functional_Integer_Node_Mul(int fid, const char* name, int dimid, int* value);

int IRICLIBDLL cg_iRIC_Read_Grid_Functional_Real_Node_Mul(int fid, const char* name, int dimid, double* value);

int IRICLIBDLL cg_iRIC_Read_Grid_Functional_Integer_Cell_Mul(int fid, const char* name, int dimid, int* value);

int IRICLIBDLL cg_iRIC_Read_Grid_Functional_Real_Cell_Mul(int fid, const char* name, int dimid, double* value);

/**********************************************/
/* Writing Grids                              */
/**********************************************/

int IRICLIBDLL cg_iRIC_WriteGridCoord1d_Mul(int fid, cgsize_t isize, double* x);

int IRICLIBDLL cg_iRIC_WriteGridCoord2d_Mul(int fid, cgsize_t isize, cgsize_t jsize, double* x, double* y);

int IRICLIBDLL cg_iRIC_WriteGridCoord3d_Mul(int fid, cgsize_t isize, cgsize_t jsize, cgsize_t ksize, double* x, double* y, double* z);

int IRICLIBDLL cg_iRIC_InitGrid_Mul(int fid, int zoneId);

int IRICLIBDLL cg_iRIC_Write_Grid_Real_Node_Mul(int fid, const char* name, double* values);

int IRICLIBDLL cg_iRIC_Write_Grid_Integer_Node_Mul(int fid, const char* name, int* values);

int IRICLIBDLL cg_iRIC_Write_Grid_Real_Cell_Mul(int fid, const char *name, double* values);

int IRICLIBDLL cg_iRIC_Write_Grid_Integer_Cell_Mul(int fid, const char *name, int* values);

/**********************************************/
/* Reading Solutions                          */
/**********************************************/

int IRICLIBDLL cg_iRIC_Read_Sol_Count_Mul(int fid, int* count);

int IRICLIBDLL cg_iRIC_Read_Sol_Time_Mul(int fid, int step, double* time);

int IRICLIBDLL cg_iRIC_Read_Sol_Iteration_Mul(int fid, int step, int* index);

int IRICLIBDLL cg_iRIC_Read_Sol_BaseIterative_Integer_Mul(int fid, int step, const char* name, int* value);

int IRICLIBDLL cg_iRIC_Read_Sol_BaseIterative_Real_Mul(int fid, int step, const char* name, double* value);

int IRICLIBDLL cg_iRIC_Read_Sol_GridCoord2d_Mul(int fid, int step, double* x, double* y);

int IRICLIBDLL cg_iRIC_Read_Sol_GridCoord3d_Mul(int fid, int step, double* x, double* y, double* z);

int IRICLIBDLL cg_iRIC_Read_Sol_Integer_Mul(int fid, int step, const char *name, int* data);

int IRICLIBDLL cg_iRIC_Read_Sol_Real_Mul(int fid, int step, const char *name, double* data);

int IRICLIBDLL cg_iRIC_Read_Sol_Cell_Integer_Mul(int fid, int step, const char *name, int* data);

int IRICLIBDLL cg_iRIC_Read_Sol_Cell_Real_Mul(int fid, int step, const char *name, double* data);

/**********************************************/
/* Writing Solutions                          */
/**********************************************/

int IRICLIBDLL cg_iRIC_Write_Sol_Time_Mul(int fid, double time);

int IRICLIBDLL cg_iRIC_Write_Sol_Iteration_Mul(int fid, int index);

int IRICLIBDLL cg_iRIC_Write_Sol_BaseIterative_Integer_Mul(int fid, const char *name, int value);

int IRICLIBDLL cg_iRIC_Write_Sol_BaseIterative_Real_Mul(int fid, const char *name, double value);

int IRICLIBDLL cg_iRIC_Write_Sol_GridCoord2d_Mul(int fid, double *x, double *y);

int IRICLIBDLL cg_iRIC_Write_Sol_GridCoord3d_Mul(int fid, double *x, double *y, double *z);

int IRICLIBDLL cg_iRIC_Write_Sol_Integer_Mul(int fid, const char *name, int* data);

int IRICLIBDLL cg_iRIC_Write_Sol_Cell_Integer_Mul(int fid, const char *name, int* data);

int IRICLIBDLL cg_iRIC_Write_Sol_Real_Mul(int fid, const char *name, double* data);

int IRICLIBDLL cg_iRIC_Write_Sol_Cell_Real_Mul(int fid, const char *name, double* data);

int IRICLIBDLL cg_iRIC_Write_ErrorCode_Mul(int fid, int errorcode);

/**********************************************/
/* Reading Boundary Conditions                */
/**********************************************/

void IRICLIBDLL cg_iRIC_Init_BC_Names_Mul(int fid);

void IRICLIBDLL cg_iRIC_Read_BC_Count_Mul(int fid, const char* type, int* num);

int IRICLIBDLL cg_iRIC_Read_BC_IndicesSize_Mul(int fid, const char* type, int num, cgsize_t* size);

int IRICLIBDLL cg_iRIC_Read_BC_Indices_Mul(int fid, const char* type, int num, cgsize_t* indices);

int IRICLIBDLL cg_iRIC_Read_BC_Integer_Mul(int fid, const char* type, int num, const char* name, int* intvalue);

int IRICLIBDLL cg_iRIC_Read_BC_Real_Mul(int fid, const char* type, int num, const char* name, double* realvalue);

int IRICLIBDLL cg_iRIC_Read_BC_RealSingle_Mul(int fid, const char* type, int num, const char* name, float* realvalue);

int IRICLIBDLL cg_iRIC_Read_BC_StringLen_Mul(int fid, const char* type, int num, const char* name, int* length);

int IRICLIBDLL cg_iRIC_Read_BC_String_Mul(int fid, const char* type, int num, const char* name, char* strvalue);

int IRICLIBDLL cg_iRIC_Read_BC_FunctionalSize_Mul(int fid, const char* type, int num, const char* name, cgsize_t* size);

int IRICLIBDLL cg_iRIC_Read_BC_Functional_Mul(int fid, const char* type, int num, const char* name, double* x, double* y);

int IRICLIBDLL cg_iRIC_Read_BC_FunctionalWithName_Mul(int fid, const char* type, int num, const char* name, const char* paramname, double* data);

int IRICLIBDLL cg_iRIC_Read_BC_Functional_RealSingle_Mul(int fid, const char* type, int num, const char* name, float* x, float* y);

int IRICLIBDLL cg_iRIC_Read_BC_FunctionalWithName_RealSingle_Mul(int fid, const char* type, int num, const char* name, const char* paramname, float* data);

/**********************************************/
/* Writing Boundary Conditions                */
/**********************************************/

int IRICLIBDLL cg_iRIC_Clear_BC_Mul(int fid);

int IRICLIBDLL cg_iRIC_Write_BC_Indices_Mul(int fid, const char* type, int num, cgsize_t size, cgsize_t* indices);

int IRICLIBDLL cg_iRIC_Write_BC_Integer_Mul(int fid, const char* type, int num, const char* name, int intvalue);

int IRICLIBDLL cg_iRIC_Write_BC_Real_Mul(int fid, const char* type, int num, const char* name, double realvalue);

int IRICLIBDLL cg_iRIC_Write_BC_String_Mul(int fid, const char* type, int num, const char* name, const char* strvalue);

int IRICLIBDLL cg_iRIC_Write_BC_Functional_Mul(int fid, const char* type, int num, const char* name, int length, double* realarray_x, double* realarray_y);

int IRICLIBDLL cg_iRIC_Write_BC_FunctionalWithName_Mul(int fid, const char* type, int num, const char* name, const char* paramname, int length, double* data);

/**********************************************/
/* Writing Particles                          */
/**********************************************/

int IRICLIBDLL cg_iRIC_Write_Sol_Particle_Pos2d_Mul(int fid, cgsize_t count, double* x, double* y);
int IRICLIBDLL cg_iRIC_Write_Sol_Particle_Pos3d_Mul(int fid, cgsize_t count, double* x, double* y, double* z);
int IRICLIBDLL cg_iRIC_Write_Sol_Particle_Real_Mul(int fid, const char* name, double* value);
int IRICLIBDLL cg_iRIC_Write_Sol_Particle_Integer_Mul(int fid, const char* name, int* value);

/**********************************************/
/* Reading Geographic Data                    */
/**********************************************/

int IRICLIBDLL cg_iRIC_Read_Geo_Count_Mul(int fid, const char* name, int* count);

int IRICLIBDLL cg_iRIC_Read_Geo_Filename_Mul(int fid, const char* name, int geoid, char* filename, int* type);

/***********************************************/
/* Reading Geographic Data (Polygon)           */
/***********************************************/

int IRICLIBDLL iRIC_Geo_Polygon_Open(const char* filename, int *id);

int IRICLIBDLL iRIC_Geo_Polygon_Read_IntegerValue(int id, int* value);

int IRICLIBDLL iRIC_Geo_Polygon_Read_RealValue(int id, double* value);

int IRICLIBDLL iRIC_Geo_Polygon_Read_PointCount(int id, int* size);

int IRICLIBDLL iRIC_Geo_Polygon_Read_Points(int id, double* x, double* y);

int IRICLIBDLL iRIC_Geo_Polygon_Read_HoleCount(int id, int* count);

int IRICLIBDLL iRIC_Geo_Polygon_Read_HolePointCount(int id, int holeid, int* count);

int IRICLIBDLL iRIC_Geo_Polygon_Read_HolePoints(int id, int holeid, double* x, double* y);

int IRICLIBDLL iRIC_Geo_Polygon_Close(int id);

/***********************************************/
/* Reading Geographic Data (RiverSurvey)       */
/***********************************************/

int IRICLIBDLL iRIC_Geo_RiverSurvey_Open(const char* filename, int* id);

int IRICLIBDLL iRIC_Geo_RiverSurvey_Read_Count(int id, int* count);

int IRICLIBDLL iRIC_Geo_RiverSurvey_Read_Position(int id, int pointid, double* x, double* y);

int IRICLIBDLL iRIC_Geo_RiverSurvey_Read_Direction(int id, int pointid, double* vx, double* vy);

int IRICLIBDLL iRIC_Geo_RiverSurvey_Read_Name(int id, int pointid, const char** name);

int IRICLIBDLL iRIC_Geo_RiverSurvey_Read_RealName(int id, int pointid, double* name);

int IRICLIBDLL iRIC_Geo_RiverSurvey_Read_LeftShift(int id, int pointid, double* shift);

int IRICLIBDLL iRIC_Geo_RiverSurvey_Read_AltitudeCount(int id, int pointid, int* count);

int IRICLIBDLL iRIC_Geo_RiverSurvey_Read_Altitudes(int id, int pointid, double* position, double* height, int* active);

int IRICLIBDLL iRIC_Geo_RiverSurvey_Read_FixedPointL(int id, int pointid, int* set, double* directionX, double* directionY, int* index);

int IRICLIBDLL iRIC_Geo_RiverSurvey_Read_FixedPointR(int id, int pointid, int* set, double* directionX, double* directionY, int* index);

int IRICLIBDLL iRIC_Geo_RiverSurvey_Read_WaterSurfaceElevation(int id, int pointid, int *set, double* value);

int IRICLIBDLL iRIC_Geo_RiverSurvey_Close(int id);

/************************************************/
/* Functions for handling single CGNS files     */
/*                                              */
/* These functions are created using            */
/* make_iriclib_singlefuncs_header.php.         */
/************************************************/

int IRICLIBDLL cg_iRIC_Set_ZoneId(int zoneid);

int IRICLIBDLL cg_iRIC_Read_Integer(const char* name, int* intvalue);

int IRICLIBDLL cg_iRIC_Read_Real(const char* name, double* realvalue);

int IRICLIBDLL cg_iRIC_Read_RealSingle(const char* name, float* realvalue);

int IRICLIBDLL cg_iRIC_Read_StringLen(const char* name, int* length);

int IRICLIBDLL cg_iRIC_Read_String(const char* name, char* strvalue);

int IRICLIBDLL cg_iRIC_Read_FunctionalSize(const char* name, cgsize_t* size);

int IRICLIBDLL cg_iRIC_Read_Functional(const char* name, double* x, double* y);

int IRICLIBDLL cg_iRIC_Read_FunctionalWithName(const char* name, const char* paramname, double* data);

int IRICLIBDLL cg_iRIC_Read_Functional_RealSingle(const char* name, float* x, float* y);

int IRICLIBDLL cg_iRIC_Read_FunctionalWithName_RealSingle(const char* name, const char* paramname, float* data);

int IRICLIBDLL cg_iRIC_Write_Integer(const char* name, int intvalue);

int IRICLIBDLL cg_iRIC_Write_Real(const char* name, double realvalue);

int IRICLIBDLL cg_iRIC_Write_String(const char* name, char* strvalue);

int IRICLIBDLL cg_iRIC_Write_Functional(const char* name, int length, double* realarray_x, double* realarray_y);

int IRICLIBDLL cg_iRIC_Write_FunctionalWithName(const char* name, const char* paramname, int length, double* data);

int IRICLIBDLL cg_iRIC_Read_Complex_Count(const char* groupname, int* num);

int IRICLIBDLL cg_iRIC_Read_Complex_Integer(const char* groupname, int num, const char* name, int* intvalue);

int IRICLIBDLL cg_iRIC_Read_Complex_Real(const char* groupname, int num, const char* name, double* realvalue);

int IRICLIBDLL cg_iRIC_Read_Complex_RealSingle(const char* groupname, int num, const char* name, float* realvalue);

int IRICLIBDLL cg_iRIC_Read_Complex_StringLen(const char* groupname, int num, const char* name, int* length);

int IRICLIBDLL cg_iRIC_Read_Complex_String(const char* groupname, int num, const char* name, char* strvalue);

int IRICLIBDLL cg_iRIC_Read_Complex_FunctionalSize(const char* groupname, int num, const char* name, cgsize_t* size);

int IRICLIBDLL cg_iRIC_Read_Complex_Functional(const char* groupname, int num, const char* name, double* x, double* y);

int IRICLIBDLL cg_iRIC_Read_Complex_FunctionalWithName(const char* groupname, int num, const char* name, const char* paramname, double* data);

int IRICLIBDLL cg_iRIC_Read_Complex_Functional_RealSingle(const char* groupname, int num, const char* name, float* x, float* y);

int IRICLIBDLL cg_iRIC_Read_Complex_FunctionalWithName_RealSingle(const char* groupname, int num, const char* name, const char* paramname, float* data);

int IRICLIBDLL cg_iRIC_Read_Grid_Complex_Node(const char* groupname, int* values);

int IRICLIBDLL cg_iRIC_Read_Grid_Complex_Cell(const char *groupname, int* values);

int IRICLIBDLL cg_iRIC_Clear_Complex();

int IRICLIBDLL cg_iRIC_Write_Complex_Integer(const char* groupname, int num, const char* name, int intvalue);

int IRICLIBDLL cg_iRIC_Write_Complex_Real(const char* groupname, int num, const char* name, double realvalue);

int IRICLIBDLL cg_iRIC_Write_Complex_String(const char* groupname, int num, const char* name, char* strvalue);

int IRICLIBDLL cg_iRIC_Write_Complex_Functional(const char* groupname, int num, const char* name, int length, double* realarray_x, double* realarray_y);

int IRICLIBDLL cg_iRIC_Write_Complex_FunctionalWithName(const char* groupname, int num, const char* name, const char* paramname, int length, double* data);

int IRICLIBDLL cg_iRIC_Write_Grid_Complex_Node(const char* groupname, int* values);

int IRICLIBDLL cg_iRIC_Write_Grid_Complex_Cell(const char *groupname, int* values);

int IRICLIBDLL cg_iRIC_GotoGridCoord2d(cgsize_t* isize, cgsize_t* jsize);

int IRICLIBDLL cg_iRIC_GetGridCoord2d(double *x, double *y);

int IRICLIBDLL cg_iRIC_GotoGridCoord3d(cgsize_t* isize, cgsize_t* jsize, cgsize_t* ksize);

int IRICLIBDLL cg_iRIC_GetGridCoord3d(double *x, double *y, double *z);

int IRICLIBDLL cg_iRIC_Read_Grid_Real_Node(const char* name, double* values);

int IRICLIBDLL cg_iRIC_Read_Grid_Integer_Node(const char* name, int* values);

int IRICLIBDLL cg_iRIC_Read_Grid_Real_Cell(const char *name, double* values);

int IRICLIBDLL cg_iRIC_Read_Grid_Integer_Cell(const char *name, int* values);

int IRICLIBDLL cg_iRIC_Read_Grid_FunctionalDimensionSize(const char* name, const char* dimname, cgsize_t* count);

int IRICLIBDLL cg_iRIC_Read_Grid_FunctionalDimension_Integer(const char* name, const char* dimname, int* value);

int IRICLIBDLL cg_iRIC_Read_Grid_FunctionalDimension_Real(const char* name, const char* dimname, double* value);

int IRICLIBDLL cg_iRIC_Read_Grid_FunctionalTimeSize(const char* name, cgsize_t* count);

int IRICLIBDLL cg_iRIC_Read_Grid_FunctionalTime(const char* name, double* time);

int IRICLIBDLL cg_iRIC_Read_Grid_Functional_Integer_Node(const char* name, int dimid, int* value);

int IRICLIBDLL cg_iRIC_Read_Grid_Functional_Real_Node(const char* name, int dimid, double* value);

int IRICLIBDLL cg_iRIC_Read_Grid_Functional_Integer_Cell(const char* name, int dimid, int* value);

int IRICLIBDLL cg_iRIC_Read_Grid_Functional_Real_Cell(const char* name, int dimid, double* value);

int IRICLIBDLL cg_iRIC_WriteGridCoord1d(cgsize_t isize, double* x);

int IRICLIBDLL cg_iRIC_WriteGridCoord2d(cgsize_t isize, cgsize_t jsize, double* x, double* y);

int IRICLIBDLL cg_iRIC_WriteGridCoord3d(cgsize_t isize, cgsize_t jsize, cgsize_t ksize, double* x, double* y, double* z);

int IRICLIBDLL cg_iRIC_InitGrid(int zoneId);

int IRICLIBDLL cg_iRIC_Write_Grid_Real_Node(const char* name, double* values);

int IRICLIBDLL cg_iRIC_Write_Grid_Integer_Node(const char* name, int* values);

int IRICLIBDLL cg_iRIC_Write_Grid_Real_Cell(const char *name, double* values);

int IRICLIBDLL cg_iRIC_Write_Grid_Integer_Cell(const char *name, int* values);

int IRICLIBDLL cg_iRIC_Read_Sol_Count(int* count);

int IRICLIBDLL cg_iRIC_Read_Sol_Time(int step, double* time);

int IRICLIBDLL cg_iRIC_Read_Sol_Iteration(int step, int* index);

int IRICLIBDLL cg_iRIC_Read_Sol_BaseIterative_Integer(int step, const char* name, int* value);

int IRICLIBDLL cg_iRIC_Read_Sol_BaseIterative_Real(int step, const char* name, double* value);

int IRICLIBDLL cg_iRIC_Read_Sol_GridCoord2d(int step, double* x, double* y);

int IRICLIBDLL cg_iRIC_Read_Sol_GridCoord3d(int step, double* x, double* y, double* z);

int IRICLIBDLL cg_iRIC_Read_Sol_Integer(int step, const char *name, int* data);

int IRICLIBDLL cg_iRIC_Read_Sol_Real(int step, const char *name, double* data);

int IRICLIBDLL cg_iRIC_Read_Sol_Cell_Integer(int step, const char *name, int* data);

int IRICLIBDLL cg_iRIC_Read_Sol_Cell_Real(int step, const char *name, double* data);

int IRICLIBDLL cg_iRIC_Write_Sol_Time(double time);

int IRICLIBDLL cg_iRIC_Write_Sol_Iteration(int index);

int IRICLIBDLL cg_iRIC_Write_Sol_BaseIterative_Integer(const char *name, int value);

int IRICLIBDLL cg_iRIC_Write_Sol_BaseIterative_Real(const char *name, double value);

int IRICLIBDLL cg_iRIC_Write_Sol_GridCoord2d(double *x, double *y);

int IRICLIBDLL cg_iRIC_Write_Sol_GridCoord3d(double *x, double *y, double *z);

int IRICLIBDLL cg_iRIC_Write_Sol_Integer(const char *name, int* data);

int IRICLIBDLL cg_iRIC_Write_Sol_Cell_Integer(const char *name, int* data);

int IRICLIBDLL cg_iRIC_Write_Sol_Real(const char *name, double* data);

int IRICLIBDLL cg_iRIC_Write_Sol_Cell_Real(const char *name, double* data);

int IRICLIBDLL cg_iRIC_Write_ErrorCode(int errorcode);

void IRICLIBDLL cg_iRIC_Init_BC_Names();

void IRICLIBDLL cg_iRIC_Read_BC_Count(const char* type, int* num);

int IRICLIBDLL cg_iRIC_Read_BC_IndicesSize(const char* type, int num, cgsize_t* size);

int IRICLIBDLL cg_iRIC_Read_BC_Indices(const char* type, int num, cgsize_t* indices);

int IRICLIBDLL cg_iRIC_Read_BC_Integer(const char* type, int num, const char* name, int* intvalue);

int IRICLIBDLL cg_iRIC_Read_BC_Real(const char* type, int num, const char* name, double* realvalue);

int IRICLIBDLL cg_iRIC_Read_BC_RealSingle(const char* type, int num, const char* name, float* realvalue);

int IRICLIBDLL cg_iRIC_Read_BC_StringLen(const char* type, int num, const char* name, int* length);

int IRICLIBDLL cg_iRIC_Read_BC_String(const char* type, int num, const char* name, char* strvalue);

int IRICLIBDLL cg_iRIC_Read_BC_FunctionalSize(const char* type, int num, const char* name, cgsize_t* size);

int IRICLIBDLL cg_iRIC_Read_BC_Functional(const char* type, int num, const char* name, double* x, double* y);

int IRICLIBDLL cg_iRIC_Read_BC_FunctionalWithName(const char* type, int num, const char* name, const char* paramname, double* data);

int IRICLIBDLL cg_iRIC_Read_BC_Functional_RealSingle(const char* type, int num, const char* name, float* x, float* y);

int IRICLIBDLL cg_iRIC_Read_BC_FunctionalWithName_RealSingle(const char* type, int num, const char* name, const char* paramname, float* data);

int IRICLIBDLL cg_iRIC_Clear_BC();

int IRICLIBDLL cg_iRIC_Write_BC_Indices(const char* type, int num, cgsize_t size, cgsize_t* indices);

int IRICLIBDLL cg_iRIC_Write_BC_Integer(const char* type, int num, const char* name, int intvalue);

int IRICLIBDLL cg_iRIC_Write_BC_Real(const char* type, int num, const char* name, double realvalue);

int IRICLIBDLL cg_iRIC_Write_BC_String(const char* type, int num, const char* name, const char* strvalue);

int IRICLIBDLL cg_iRIC_Write_BC_Functional(const char* type, int num, const char* name, int length, double* realarray_x, double* realarray_y);

int IRICLIBDLL cg_iRIC_Write_BC_FunctionalWithName(const char* type, int num, const char* name, const char* paramname, int length, double* data);

int IRICLIBDLL cg_iRIC_Write_Sol_Particle_Pos2d(cgsize_t count, double* x, double* y);

int IRICLIBDLL cg_iRIC_Write_Sol_Particle_Pos3d(cgsize_t count, double* x, double* y, double* z);

int IRICLIBDLL cg_iRIC_Write_Sol_Particle_Real(const char* name, double* value);

int IRICLIBDLL cg_iRIC_Write_Sol_Particle_Integer(const char* name, int* value);

int IRICLIBDLL cg_iRIC_Read_Geo_Count(const char* name, int* count);

int IRICLIBDLL cg_iRIC_Read_Geo_Filename(const char* name, int geoid, char* filename, int* type);

#ifdef __cplusplus
}
#endif
