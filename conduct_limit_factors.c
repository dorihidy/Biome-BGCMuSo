/* 
conduct_limit_factors.c
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


int conduct_limit_factors(file logfile, const control_struct* ctrl, const soilprop_struct* sprop, const epconst_struct* epc, 
						  epvar_struct* epv)
{
	int errorCode=0;
	int layer;
	double VWCsat,VWCfc,VWCwp, VWC_WScrit1, VWC_WScrit2;
    double PSI_WScrit1[N_SOILLAYERS];
	double PSI_WScrit2[N_SOILLAYERS];

	VWCsat=VWCfc=VWCwp=VWC_WScrit1=VWC_WScrit2=0;

	

	/* calculations layer by layer (due to different soil properties) */
	for (layer=0; layer < N_SOILLAYERS; layer++)
	{
		VWCsat = sprop->VWCsat[layer]; 
		VWCfc  = sprop->VWCfc[layer]; 
		VWCwp  = sprop->VWCwp[layer];

		/* VWCratio_WScrit1 = -9999 - no drought limitation defined by the User */
		if (epc->VWCratio_WScrit1 != DATA_GAP)
		{

			VWC_WScrit1 = VWCwp + epc->VWCratio_WScrit1 * (VWCfc - VWCwp);
			PSI_WScrit1[layer] = sprop->PSIsat[layer] * pow((VWC_WScrit1 / VWCsat), -1 * sprop->soilB[layer]);
	
		}
		else
		{ 
			VWC_WScrit1 = DATA_GAP;
			PSI_WScrit1[layer] = DATA_GAP;
		}
			
		/* VWCratio_WScrit2 = -9999 - no anoxic limitation defined by the User */
		if (epc->VWCratio_WScrit2 != DATA_GAP)
		{
			VWC_WScrit2 = VWCfc + epc->VWCratio_WScrit2 * (VWCsat - VWCfc);
			PSI_WScrit2[layer] = sprop->PSIsat[layer] * pow((VWC_WScrit2 / VWCsat), -1 * sprop->soilB[layer]);
		}		
		else
		{
			VWC_WScrit2 = DATA_GAP;
			PSI_WScrit2[layer] = DATA_GAP;
		}
			


		/* CONTROL */
		if (!errorCode && VWC_WScrit2 > VWCsat)
		{
			if (VWC_WScrit2 - VWCsat > CRIT_PRECwater)
			{
				printf("\n");
				printf("ERROR in conduct_limit_factors.c: VWC_WScrit2 data is greater than saturation value in layer:%i\n", layer);
				errorCode=1; 
			}
			else
				VWC_WScrit2 = VWCsat;
		}

		if (!errorCode && VWC_WScrit2 < VWC_WScrit1)
		{
			if (VWC_WScrit1 - VWC_WScrit2 > CRIT_PRECwater)
			{
				printf("\n");
				printf("ERROR in conduct_limit_factors.c: VWC_WScrit1 data is greater then VWC_WScrit2 data in layer:%i\n", layer);
				errorCode=1; 
			}
			else
				VWC_WScrit1 = VWC_WScrit2;
		}
	
		epv->VWC_WScrit1[layer]		= VWC_WScrit1;
		epv->VWC_WScrit2[layer]		= VWC_WScrit2;


	}
	if (ctrl->spinup < 2) fprintf(logfile.ptr, "CRITICAL VALUES OF VWC and PSI FOR LIMITATION OF WATER IN 10 SOIL LAYERS (-9999: no limitation defined by the User)\n");  
	if (epc->VWCratio_WScrit1 != DATA_GAP)
	{
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "VWC [m3/m3] at start of drought limitation:      %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n", epv->VWC_WScrit1[0], epv->VWC_WScrit1[1], epv->VWC_WScrit1[2], epv->VWC_WScrit1[3], epv->VWC_WScrit1[4], epv->VWC_WScrit1[5], epv->VWC_WScrit1[6], epv->VWC_WScrit1[7], epv->VWC_WScrit1[8], epv->VWC_WScrit1[9]);
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "PSI [MPa]   at start of drought limitation:      %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n", PSI_WScrit1[0], PSI_WScrit1[1], PSI_WScrit1[2], PSI_WScrit1[3], PSI_WScrit1[4], PSI_WScrit1[5], PSI_WScrit1[6], PSI_WScrit1[7], PSI_WScrit1[8], PSI_WScrit1[9]);
	}
	else 
	{
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "VWC [m3/m3] at start of drought limitation:      %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f\n", epv->VWC_WScrit1[0], epv->VWC_WScrit1[1], epv->VWC_WScrit1[2], epv->VWC_WScrit1[3], epv->VWC_WScrit1[4], epv->VWC_WScrit1[5], epv->VWC_WScrit1[6], epv->VWC_WScrit1[7], epv->VWC_WScrit1[8], epv->VWC_WScrit1[9]);
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "PSI [MPa]   at start of drought limitation:      %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f\n", PSI_WScrit1[0], PSI_WScrit1[1], PSI_WScrit1[2], PSI_WScrit1[3], PSI_WScrit1[4], PSI_WScrit1[5], PSI_WScrit1[6], PSI_WScrit1[7], PSI_WScrit1[8], PSI_WScrit1[9]);

	}

	if (epc->VWCratio_WScrit2 != DATA_GAP)
	{
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "VWC [m3/m3] at start of anoxic  limitation:      %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n", epv->VWC_WScrit2[0], epv->VWC_WScrit2[1], epv->VWC_WScrit2[2], epv->VWC_WScrit2[3], epv->VWC_WScrit2[4], epv->VWC_WScrit2[5], epv->VWC_WScrit2[6], epv->VWC_WScrit2[7], epv->VWC_WScrit2[8], epv->VWC_WScrit2[9]);
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "PSI [MPa]   at start of anoxic  limitation:      %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n", PSI_WScrit2[0], PSI_WScrit2[1], PSI_WScrit2[2], PSI_WScrit2[3], PSI_WScrit2[4], PSI_WScrit2[5], PSI_WScrit2[6], PSI_WScrit2[7], PSI_WScrit2[8], PSI_WScrit2[9]);
	}
	else
	{
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "VWC [m3/m3] at start of anoxic  limitation:      %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f\n", epv->VWC_WScrit2[0], epv->VWC_WScrit2[1], epv->VWC_WScrit2[2], epv->VWC_WScrit2[3], epv->VWC_WScrit2[4], epv->VWC_WScrit2[5], epv->VWC_WScrit2[6], epv->VWC_WScrit2[7], epv->VWC_WScrit2[8], epv->VWC_WScrit2[9]);
		if (ctrl->spinup < 2) fprintf(logfile.ptr, "PSI [MPa]   at start of anoxic  limitation:      %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f\n", PSI_WScrit2[0], PSI_WScrit2[1], PSI_WScrit2[2], PSI_WScrit2[3], PSI_WScrit2[4], PSI_WScrit2[5], PSI_WScrit2[6], PSI_WScrit2[7], PSI_WScrit2[8], PSI_WScrit2[9]);
	}
	if (ctrl->spinup < 2) fprintf(logfile.ptr, " \n");
		
	return (errorCode);
}
