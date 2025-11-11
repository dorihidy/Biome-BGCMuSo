/* 
multilayer_hydrolprocess.c
calculation of hydrological processes routine by routne

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2025, D. Hidy [dori.hidy@gmail.com]
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
	wstate_struct* ws, wflux_struct* wf, nstate_struct* ns, nflux_struct* nf, cstate_struct* cs, cflux_struct* cf, flooding_struct* FLS, int* mondays)
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
	int errorCode = 0;
	double soilw_before = 0;

	/* to check balance in check_virtualLayer_balance */
	for (layer = 0; layer < N_SOILLAYERS; layer++) ws->soilw_pre[layer] = ws->soilw[layer];
	sprop->soilw_NORMcf_pre = sprop->soilw_NORMcf;
	sprop->soilw_CAPILcf_pre = sprop->soilw_CAPILcf;
	sprop->soilw_NORMgw_pre = sprop->soilw_NORMgw;
	sprop->soilw_CAPILgw_pre = sprop->soilw_CAPILgw;


	/* update of hydrolparams */
	if (!errorCode && multilayer_hydrolparams(sitec, sprop, ws, epv))
	{
		printf("\n");
		printf("ERROR in multilayer_hydrolparams.c from multilayer_hydrolprocess.c\n");
		errorCode = 52601;
	}


	/* INFILTRATION AND PONDW FORMATION */
	if (!errorCode && infiltANDpond(ctrl, ws, wf))
	{
		printf("\n");
		printf("ERROR in infiltANDpond.c from multilayer_hydrolprocess.c\n");
		errorCode = 52603;
	}


	/* PERCOLATION  in layers without groundwater */	
	if (!errorCode && tipping(ctrl, sitec, sprop, epv, ws, wf))
	{
		printf("\n");
		printf("ERROR in tipping.c from multilayer_hydrolprocess.c\n");
		errorCode=52604;
	} 

	/* PERCOLATION  in presence of GW */
	if (sprop->GWlayer != DATA_GAP)
	{
		/* in layers in capillary zone(without GW) */
		if (sprop->dz_NORMcf + sprop->dz_CAPILcf)
		{ 
			if (!errorCode && capillary_tipping(ctrl, sitec, sprop, epv, ws, wf))
			{
				printf("\n");
				printf("ERROR in capillary_tipping.c from multilayer_hydrolprocess.c\n");
				errorCode = 52605;
			}
		}
		/* in layers in capillary zone (with GW) */		
		if (!errorCode && groundwater_tipping(ctrl, sitec, sprop, epv, ws, wf))
		{
			printf("\n");
			printf("ERROR in groundwater_tipping.c from multilayer_hydrolprocess.c\n");
			errorCode = 52606;
		}

	
	}

	
	/* SOIL EVAPORATION */
	if (!errorCode && soilEVP_calc(ctrl, sitec, sprop, epv, ws, wf))
	{
		printf("ERROR in soilEVP_calc.c from multilayer_hydrolprocess.c\n");
		errorCode = 52607;
	}


	/* TRANSPIRATION */
	if (!errorCode && multilayer_transpiration(ctrl, sitec, sprop, epv, ws, wf))
	{
		printf("ERROR in multilayer_transpiration.c from multilayer_hydrolprocess.c\n");
		errorCode = 52608;
	}

	/* POND AND RUNOFF */	
	if (!errorCode && pondANDrunoffD(ctrl,sitec,sprop, epv,ws, wf))
	{
		printf("\n");
		printf("ERROR in pondANDrunoffD.c from multilayer_hydrolprocess.c\n");
		errorCode=52609; 
	} 

	/* FLOODING : fills the soil layers and pondw */
	if (!errorCode && ctrl->spinup != 1 && flooding(ctrl, sitec, FLS, sprop, epv, ws, wf, cs, cf, ns, nf, soilInfo, mondays))
	{
		printf("ERROR in flooding.c from multilayer_hydrolprocess.c\n");
		errorCode = 52610;
	}

	/* diffusion in GWlayers */
	if (sprop->GWlayer != DATA_GAP)
	{
		if (!errorCode && groundwater_diffusion(sitec, sprop, epv, ws, wf))
		{
			printf("\n");
			printf("ERROR in groundwater_diffusion.c from multilayer_hydrolprocess.c\n");
			errorCode = 52611;
		}
		if (!errorCode && capillary_diffusion(sitec, sprop, epv, ws, wf))
		{
			printf("\n");
 			printf("ERROR in capillary_diffusion.c from multilayer_hydrolprocess.c\n");
			errorCode = 52612;
		}
	}

	/* diffusion in non-GWlayers */
	if (!errorCode && diffusion(sitec, sprop, epv, ws, wf))
	{
		printf("\n");
		printf("ERROR in diffusion.c from multilayer_hydrolprocess.c\n");
		errorCode = 52613;
	}

	/* calculation of net water transport */
	for (layer = 0; layer < N_SOILLAYERS; layer++) 	wf->soilwFlux[layer] = wf->soilwPercol[layer] + wf->soilwDiffus[layer];
	wf->soilwFlux_NORMvsCAPILgw       = wf->soilwPercol_NORMvsCAPILgw       + wf->soilwDiffus_NORMvsCAPILgw;
	wf->soilwFlux_NORMvsCAPILcf = wf->soilwPercol_NORMvsCAPILcf + wf->soilwDiffus_NORMvsCAPILcf;

	/* top soil flux from downward fo ET-limitation calculation on the next day */
	if (wf->soilwDiffus[0] < 0)
		wf->soilwFLuxFromBelow = -1 * wf->soilwDiffus[0] + wf->GWdischarge[0];
	else
		wf->soilwFLuxFromBelow = wf->GWdischarge[0];

	/* ET-calculation */	
	wf->ET         = wf->EVPcanopyw + wf->EVPsoilw + wf->TRPsoilw_SUM + wf->EVPpondw + wf->SUBLsnoww;
	wf->EVPsurface = wf->EVPsoilw + wf->EVPpondw;

	
	/* BOTTOM LAYER IS SPECIAL: percolated water is net loss for the system, water content does not change 	*/	
	if (sprop->GWD == DATA_GAP || (sprop->GWlayer == DATA_GAP))
	{
		soilw_before              = ws->soilw[N_SOILLAYERS-1];
		if (epv->VWC[N_SOILLAYERS - 1] < sprop->VWCeq[N_SOILLAYERS - 1])
		{ 
			epv->VWC[N_SOILLAYERS - 1] = sprop->VWCeq[N_SOILLAYERS - 1];
			ws->soilw[N_SOILLAYERS-1] = epv->VWC[N_SOILLAYERS - 1] * (sitec->soillayer_thickness[N_SOILLAYERS-1]) * water_density;


			wf->soilwFlux[N_SOILLAYERS-1] += soilw_before - ws->soilw[N_SOILLAYERS-1];
		}

	}

	/* Soilstress calculation based on VWC or transpiration demand-possibitiy */
	if (!errorCode && soilstress_calculation(ctrl, epc, sprop, epv, ws, wf))
	{
		printf("\n");
		printf("ERROR in soilstress_calculation.c from multilayer_hydrolprocess.c\n");
		errorCode=52614; 
	} 

	/* State update of hydroparams */	
	if (!errorCode && multilayer_hydrolparams(sitec, sprop, ws, epv))
	{
		printf("\n");
		printf("ERROR in multilayer_hydrolparams.c from multilayer_hydrolprocess.c\n");
		errorCode = 52601;
	}
	
	
	if (!errorCode && hydrol_control(sitec, sprop, ws, wf, epv))
	{
		printf("\n");
		printf("ERROR in hydrol_control.c from multilayer_hydrolprocess.c\n");
		errorCode = 52615;
	}




	return (errorCode);
}