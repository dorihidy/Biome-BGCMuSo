/* 
groundwater_calculations.c
Summarize groundwater routines

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2025, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ini.h"
#include "bgc_struct.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"
#include "bgc_func.h"


int groundwater_calculations(control_struct* ctrl, const siteconst_struct* sitec, const groundwaterINIT_struct* GWS, soilprop_struct* sprop, soilInfo_struct* soilInfo, epvar_struct* epv,
	                     wstate_struct* ws, wflux_struct* wf, cstate_struct* cs, nstate_struct* ns)
{
	int errorCode = 0;



	/* ***************************************************************************************************** */
	/* Groundwater preproc: if  */

	if (ctrl->spinup != 1)
	{
		sprop->GWD_pre = sprop->GWD;
		if (!errorCode && groundwater_preproc(ctrl, GWS, sitec, sprop, soilInfo, ws, wf, cs, ns))
		{
			printf("\n");
			printf("ERROR in groundwater_preproc.c for groundwater_calculation.c\n");
			errorCode = 1;
		}
	}
	else
		sprop->GWD = DATA_GAP;

	/* ***************************************************************************************************** */
    /* calculation of GWD for VWCeq calculation (in case of GWD-simulation: input data, else: calculated data */

	if (!errorCode && calc_VWCeq(ctrl, sitec, sprop))
	{
		printf("\n");
		printf("ERROR in calc_VWCeq.c for groundwater_calculation.c\n");
		errorCode = 1;
	}

	/* ***************************************************************************************************** */
	/* calculation of effect of GW-moving: water content and material pools */
	
	if (sprop->GWlayer != DATA_GAP && sprop->GWD != sprop->GWD_pre)
	{

		if (!errorCode && groundwater_movement(sitec, sprop, epv, soilInfo, ws, wf, cs, ns))
		{
			printf("\n");
			printf("ERROR in groundwater_movement.c for groundwater_calculation.c\n");
			errorCode = 1;
		}
	}
	

	return (errorCode);
}
