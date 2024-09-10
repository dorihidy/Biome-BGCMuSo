/* 
multilayer_hydrolprocess.c
calculation of soil water content layer by layer taking into account soil hydrological processes 
(precipitation, evaporation, runoff, percolation, diffusion)

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2022, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_constants.h"
#include "bgc_func.h"    

int multilayer_hydrolprocess(control_struct* ctrl, siteconst_struct* sitec, soilprop_struct* sprop, soilInfo_struct* soilInfo, const epconst_struct* epc, epvar_struct* epv,
	                         wstate_struct* ws, wflux_struct* wf, nstate_struct* ns, nflux_struct* nf, cstate_struct* cs, cflux_struct* cf, 
	                         groundwaterINIT_struct* GWS, flooding_struct* FLS, int* mondays)
{
	/* given a list of site constants and the soil water mass (kg/m2),
	this function returns the soil water potential (MPa)
	inputs:


	For further discussion see:
	Cosby, B.J., G.M. Hornberger, R.B. Clapp, and T.R. Ginn, 1984.     

	Balsamo et al 2009 - 
	A Revised Hydrology for the ECMWF Model - Verification from Field Site to Water Storage IFS - JHydromet.pdf

	Chen and Dudhia 2001 - 
	Coupling an Advanced Land Surface-Hydrology Model with the PMM5 Modeling System Part I - MonWRev.pdf*/
	

	/* internal variables */



	int layer;
	int errorCode=0;
	double soilw_before = 0;


	/* update of hydrolparams */
	if (!errorCode && multilayer_hydrolparams(sitec, sprop, ws, epv))
	{
		printf("\n");
		printf("ERROR in multilayer_hydrolparams.c from multilayer_hydrolprocess.c\n");
		errorCode = 52610;
	}

	/* GROUNDWATER preprocess: calculation depth of GW and CF, GW-movchange */	
	if (!errorCode && groundwater_preproc(ctrl, sitec, GWS, sprop, epv, soilInfo, ws, wf, cs, cf, ns, nf))
	{
		printf("ERROR in groundwater_preproc.c from bgc.c\n");
		errorCode=52601;
	}

	/* INFILTRATION AND PONDW FORMATION */
	if (!errorCode && infiltANDpond(sitec,sprop, epv,ws, wf))
	{
		printf("\n");
		printf("ERROR in infiltANDpond.c from multilayer_hydrolprocess.c\n");
		errorCode=52602; 
	} 

	/* PERCOLATION  in layers without groundwater */	
	if (!errorCode && tipping(sitec, sprop, epv, ws, wf))
	{
		printf("\n");
		printf("ERROR in tipping.c from multilayer_hydrolprocess.c\n");
		errorCode=52603;
	} 

	/* PERCOLATION  in layers with groundwater */
	if (sprop->GWlayer != DATA_GAP)
	{
		if (!errorCode && groundwater_tipping(sitec, sprop, epv, ws, wf))
		{
			printf("\n");
			printf("ERROR in tipping.c from multilayer_hydrolprocess.c\n");
			errorCode = 52604;
		}
	}

	/* SOIL EVAPORATION */
	if (!errorCode && soilEVP_calc(ctrl,sitec,sprop, epv, ws,wf))
	{
		printf("ERROR in soilEVP_calc.c from multilayer_hydrolprocess.c\n");
		errorCode=52605;
	}

	/* TRANSPIRATION */
	if (!errorCode && multilayer_transpiration(ctrl, sitec, sprop, epv, ws, wf))
	{
		printf("ERROR in multilayer_transpiration.c from multilayer_hydrolprocess.c\n");
		errorCode=52606;
	}

	/* POND AND RUNOFF */	
	if (!errorCode && pondANDrunoffD(ctrl,sitec,sprop, epv,ws, wf))
	{
		printf("\n");
		printf("ERROR in pondANDrunoffD.c from multilayer_hydrolprocess.c\n");
		errorCode=52607; 
	} 

	/* FLOODING : fills the soil layers and pondw */
	if (!errorCode && flooding(ctrl, sitec, FLS, sprop, epv, ws, wf, cs, cf, ns, nf, soilInfo, mondays))
	{
		printf("ERROR in flooding.c from multilayer_hydrolprocess.c\n");
		errorCode = 52608;
	}

	/* diffusion in GWlayers */
	if (sprop->GWlayer != DATA_GAP)
	{
		if (!errorCode && groundwater_diffusion(sitec, sprop, epv, ws, wf))
		{
			printf("\n");
			printf("ERROR in groundwater_diffusion.c from multilayer_hydrolprocess.c\n");
			errorCode = 52609;
		}
	}

	/* diffusion in non-GWlayers */
	if (!errorCode && diffusion(sitec, sprop, epv, ws, wf))
	{
		printf("\n");
		printf("ERROR in diffusion.c from multilayer_hydrolprocess.c\n");
		errorCode = 52610;
	}

	

	/* ET-calculation */	
	wf->ET         = wf->EVPcanopyw + wf->EVPsoilw + wf->TRPsoilw_SUM + wf->EVPpondw + wf->SUBLsnoww;
	wf->EVPsurface = wf->EVPsoilw + wf->EVPpondw;

	
	/* BOTTOM LAYER IS SPECIAL: percolated water is net loss for the system, water content does not change 	*/
	
	if (sprop->GWD == DATA_GAP || (sprop->GWlayer == DATA_GAP))
	{
		soilw_before              = ws->soilw[N_SOILLAYERS-1];
		epv->VWC[N_SOILLAYERS-1]  = sprop->VWCfc[N_SOILLAYERS-1];
		ws->soilw[N_SOILLAYERS-1] = sprop->VWCfc[N_SOILLAYERS-1] * (sitec->soillayer_thickness[N_SOILLAYERS-1]) * water_density;
	
		wf->soilwFlux[N_SOILLAYERS-1] += soilw_before - ws->soilw[N_SOILLAYERS-1];

	}
	
	/* Soilstress calculation based on VWC or transpiration demand-possibitiy */
	if (!errorCode && soilstress_calculation(ctrl, epc, sprop, epv, ws, wf))
	{
		printf("\n");
		printf("ERROR in soilstress_calculation.c from multilayer_hydrolprocess.c\n");
		errorCode=52611; 
	} 

	/* State update of hydroparams */	
	if (!errorCode && multilayer_hydrolparams(sitec, sprop, ws, epv))
	{
		printf("\n");
		printf("ERROR in multilayer_hydrolparams.c from multilayer_hydrolprocess.c\n");
		errorCode = 52612;
	}
	
	/* control to avoid irrealistic PET */
	if (wf->ET - wf->PET > CRIT_PRECwater)
	{
		printf("ERROR in potential evaporation calculation (ET > PET) in multilayer_hydrolprocess.c\n");
		errorCode = 52613;
	}


	/* CONTROL and calculating averages - unrealistic VWC content (higher than saturation value or less then hygroscopic) - GWtest */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		/* control */
		if (wf->GWrecharge && wf->GWdischarge && !errorCode)
		{
			printf("ERROR in multilayer_hydrolprocess.c: GWrecharge AND GWdischarge in the same time and space\n");
			errorCode = 52613;
		}

		if (epv->VWC[layer] < sprop->VWChw[layer])
		{
			if (sprop->VWChw[layer] - epv->VWC[layer] < 1e-3)
			{
				wf->soilwFlux[N_SOILLAYERS - 1] -= (sprop->VWChw[layer] - epv->VWC[layer]) * water_density * sitec->soillayer_thickness[layer];
				epv->VWC[layer] = sprop->VWChw[layer];
				ws->soilw[layer] = epv->VWC[layer] * water_density * sitec->soillayer_thickness[layer];
			}
			else
			{
				if (!errorCode)
				{ 
					printf("\n");
					printf("ERROR in soil water content calculation (multilayer_hydrolprocess.c) - actual is less than hygroscopic value\n");
					errorCode = 52613;
				}
			}

		}

		if (epv->VWC[layer] - sprop->VWCsat[layer] > CRIT_PREC)
		{
			if (!errorCode)
			{
				printf("\n");
				printf("ERROR in soil water content calculation - actual is higher than saturation value (multilayer_hydrolprocess)\n");
				errorCode = 52613;
			}
		}

	}

	return (errorCode);
}