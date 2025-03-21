/* 
conductLimit_calculations.c
calculate the limitation factors of conductance
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2022, D. Hidy [dori.hidy@gmail.com]
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


int conductLimit_information(file logfile, const control_struct* ctrl, const epconst_struct* epc,  epvar_struct* epv)
{
	int errorCode=0;

	if (ctrl->spinup < 2) fprintf(logfile.ptr, "CRITICAL VALUES OF VWC and PSI FOR LIMITATION OF WATER IN 10 SOIL LAYERS (-9999: no limitation defined by the User)\n");  
	if (epc->VWCratio_WScrit1 != DATA_GAP)
	{
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "VWC [m3/m3] at start of drought limitation:      %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n", epv->VWC_WScrit1[0], epv->VWC_WScrit1[1], epv->VWC_WScrit1[2], epv->VWC_WScrit1[3], epv->VWC_WScrit1[4], epv->VWC_WScrit1[5], epv->VWC_WScrit1[6], epv->VWC_WScrit1[7], epv->VWC_WScrit1[8], epv->VWC_WScrit1[9]);
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "PSI [MPa]   at start of drought limitation:      %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n", epv->PSI_WScrit1[0], epv->PSI_WScrit1[1], epv->PSI_WScrit1[2], epv->PSI_WScrit1[3], epv->PSI_WScrit1[4], epv->PSI_WScrit1[5], epv->PSI_WScrit1[6], epv->PSI_WScrit1[7], epv->PSI_WScrit1[8], epv->PSI_WScrit1[9]);
	}
	else 
	{
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "VWC [m3/m3] at start of drought limitation:      %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f\n", epv->VWC_WScrit1[0], epv->VWC_WScrit1[1], epv->VWC_WScrit1[2], epv->VWC_WScrit1[3], epv->VWC_WScrit1[4], epv->VWC_WScrit1[5], epv->VWC_WScrit1[6], epv->VWC_WScrit1[7], epv->VWC_WScrit1[8], epv->VWC_WScrit1[9]);
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "PSI [MPa]   at start of drought limitation:      %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f\n", epv->PSI_WScrit1[0], epv->PSI_WScrit1[1], epv->PSI_WScrit1[2], epv->PSI_WScrit1[3], epv->PSI_WScrit1[4], epv->PSI_WScrit1[5], epv->PSI_WScrit1[6], epv->PSI_WScrit1[7], epv->PSI_WScrit1[8], epv->PSI_WScrit1[9]);

	}

	if (epc->VWCratio_WScrit2 != DATA_GAP)
	{
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "VWC [m3/m3] at start of anoxic  limitation:      %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n", epv->VWC_WScrit2[0], epv->VWC_WScrit2[1], epv->VWC_WScrit2[2], epv->VWC_WScrit2[3], epv->VWC_WScrit2[4], epv->VWC_WScrit2[5], epv->VWC_WScrit2[6], epv->VWC_WScrit2[7], epv->VWC_WScrit2[8], epv->VWC_WScrit2[9]);
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "PSI [MPa]   at start of anoxic  limitation:      %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n", epv->PSI_WScrit2[0], epv->PSI_WScrit2[1], epv->PSI_WScrit2[2], epv->PSI_WScrit2[3], epv->PSI_WScrit2[4], epv->PSI_WScrit2[5], epv->PSI_WScrit2[6], epv->PSI_WScrit2[7], epv->PSI_WScrit2[8], epv->PSI_WScrit2[9]);
	}
	else
	{
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "VWC [m3/m3] at start of anoxic  limitation:      %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f\n", epv->VWC_WScrit2[0], epv->VWC_WScrit2[1], epv->VWC_WScrit2[2], epv->VWC_WScrit2[3], epv->VWC_WScrit2[4], epv->VWC_WScrit2[5], epv->VWC_WScrit2[6], epv->VWC_WScrit2[7], epv->VWC_WScrit2[8], epv->VWC_WScrit2[9]);
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "PSI [MPa]   at start of anoxic  limitation:      %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f\n", epv->PSI_WScrit2[0], epv->PSI_WScrit2[1], epv->PSI_WScrit2[2], epv->PSI_WScrit2[3], epv->PSI_WScrit2[4], epv->PSI_WScrit2[5], epv->PSI_WScrit2[6], epv->PSI_WScrit2[7], epv->PSI_WScrit2[8], epv->PSI_WScrit2[9]);
	}
	if (ctrl->spinup < 2) fprintf(logfile.ptr, " \n");
		
	return (errorCode);
}
