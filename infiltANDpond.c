/* 
infiltANDpond.c
calculation of waterFromAbove, pond water accumulation and potential infiltration

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


int infiltANDpond(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{

	
	/* internal variables */
	int errorCode, layer, flagEXTRA;

	double soilwEXTRA, ratio;

	 errorCode=layer=flagEXTRA=0;
     
	/*------------------------------------------*/
	/* 1. calculation of water from above */ 
	
	 wf->waterFromAbove = (wf->prcp_to_soilSurface + wf->snoww_to_soilw + wf->canopyw_to_soilw + wf->IRG_to_soilSurface);

	
	 /*------------------------------------------*/
	/* 2.calculation of the amount of water which can still fits into the soil */ 
	
	

	soilwEXTRA = (sprop->VWCsat[0] - epv->VWC[0]) * sitec->soillayer_thickness[0] * water_density;


	
	layer = 1;
	while (flagEXTRA == 0 && layer < N_SOILLAYERS-1)
	{
		if (sprop->infiltDepth_max > sitec->soillayer_depth[layer - 1])
		{
			if (sprop->infiltDepth_max > sitec->soillayer_depth[layer])
				soilwEXTRA += (sprop->VWCsat[layer] - epv->VWC[layer]) * sitec->soillayer_thickness[layer] * water_density;
			else
			{
				ratio = ((sprop->infiltDepth_max - sitec->soillayer_depth[layer - 1]) / sitec->soillayer_thickness[layer]);
				soilwEXTRA += ((sprop->VWCsat[layer] - epv->VWC[layer]) * sitec->soillayer_thickness[layer] * water_density) * ratio;
			}
		}
		else
			flagEXTRA = 1;

		layer += 1;
	}
	if (soilwEXTRA < CRIT_PRECwater) soilwEXTRA = 0;
	
	/* ---------------------------------------*/
	/* 3. if there is pond water in the area: pondw_to_soilw, infiltPOT */ 

	ws->pondw += wf->GW_to_pondw;

	if (ws->pondw)
	{
		ws->pondw += wf->waterFromAbove;
	

		/* PRCP to pond water (limitaion: pondmax)  */
		if (ws->pondw > soilwEXTRA)
			wf->pondw_to_soilw = soilwEXTRA;
		else
			wf->pondw_to_soilw = ws->pondw;

		ws->pondw         -= wf->pondw_to_soilw;
		
		wf->infiltPOT      = wf->pondw_to_soilw;



	/* ---------------------------------------*/
	/* 4. if there is pond water in the area: soilw_to_pondw, infiltPOT */ 
   	}
	else
	{
		
		/* if empty space in soil is greater than infiltration -> pond water formation (and runoff - if pond water is small)*/
		if (wf->waterFromAbove > soilwEXTRA)
		{
			wf->infiltPOT       = soilwEXTRA;				
			wf->prcp_to_pondw = wf->waterFromAbove - soilwEXTRA;

		}
		else
			wf->infiltPOT = wf->waterFromAbove; 


	} 


	return (errorCode);
}
