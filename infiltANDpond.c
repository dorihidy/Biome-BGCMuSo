/* 
infiltANDpond.c
calculation of waterFromAbove, pond water accumulation and potential infiltration

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


int infiltANDpond(wstate_struct* ws, wflux_struct* wf)
{

	
	/* internal variables */
	int errorCode, layer;

	 errorCode=layer=0;
     
	/*------------------------------------------*/
	/* 1. calculation of water from above */ 
	
	 wf->waterFromAbove = (wf->prcp_to_soilSurface + wf->snoww_to_soilw + wf->canopyw_to_soilw + wf->IRG_to_soilSurface);

	
	/* ---------------------------------------*/
	/* 2. if there is pond water in the area: pondw_to_soilw, infiltPOT */ 

	ws->pondw += wf->GW_to_pondw;

	if (ws->pondw)
	{
		ws->pondw         += wf->waterFromAbove;
		wf->pondw_to_soilw = ws->pondw;
		ws->pondw         -= wf->pondw_to_soilw;	
		wf->infiltPOT      = wf->pondw_to_soilw;

   	}
	else
		wf->infiltPOT = wf->waterFromAbove;

	/* rain flag for tipping calculation */
	if (wf->infiltPOT)
		wf->flagRAIN = 1;
	else
		wf->flagRAIN = 0;


	return (errorCode);
}
