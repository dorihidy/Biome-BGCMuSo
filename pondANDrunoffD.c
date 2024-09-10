/* 
pondANDrunoffD.c
state update of pond water, calculation of Dunnian runoff, pond water limitation

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


int pondANDrunoffD(control_struct* ctrl, siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{

	/* internal variables */
	int errorCode, layer, flagEXTRA;
	double soilw_diff, soilw_sat, pondmax, soilwEXTRA, ratio;

	errorCode=layer=flagEXTRA=0;
	soilw_diff=soilw_sat=0;


	/*--------------------------------------*/
	/* Water flux from soil to pond */
	

	ws->pondw += wf->soilw_to_pondw + wf->prcp_to_pondw;

	
	/*--------------------------------------*/
	/* Pond water evaporation: water stored on surface  */

	if (ws->pondw)
	{
		if (wf->potEVPsurface < ws->pondw)
			wf->EVPpondw = wf->potEVPsurface;
		else
			wf->EVPpondw = ws->pondw;
	}


	if (ws->GW_waterlogging > 0)
		ws->GWsrc_W += wf->EVPpondw;
	else
		ws->pondw    -= wf->EVPpondw;

	/* control to avoid exceed the potential evaporation level 	*/
	soilw_diff = wf->EVPsoilw + wf->EVPpondw - wf->potEVPsurface;

	if (soilw_diff > 0)
	{
		wf->EVPsoilw -= soilw_diff;

		/* in case of GW in layer 0, the source of soil evaporatiaton is groundwater */
		if ((int)sprop->GWlayer == 0)
			wf->GWevap -= soilw_diff;
		else
		{
			ws->soilw[0] += soilw_diff;
			epv->VWC[0] = ws->soilw[0] / (sitec->soillayer_thickness[0] * water_density);
		}
	}
	
	/*--------------------------------------*/
	/* saturation of top soil layers: correction of water flux from pond to soil */


	layer = 0;
	flagEXTRA = 0;
	while (flagEXTRA == 0 && layer < N_SOILLAYERS - 1)
	{
		if ((layer == 0 || sprop->infiltDepth_max > sitec->soillayer_depth[layer-1]) && ws->pondw > 0)
		{
			if (sprop->infiltDepth_max > sitec->soillayer_depth[layer] || layer == 0)
				soilwEXTRA = (sprop->VWCsat[layer] - epv->VWC[layer]) * sitec->soillayer_thickness[layer] * water_density;
			else
			{
				ratio = ((sprop->infiltDepth_max - sitec->soillayer_depth[layer - 1]) / sitec->soillayer_thickness[layer]);
				soilwEXTRA = ((sprop->VWCsat[layer] - epv->VWC[layer]) * sitec->soillayer_thickness[layer] * water_density) * ratio;
			}
		
			if (soilwEXTRA)
			{
				if (soilwEXTRA > ws->pondw)
					soilwEXTRA = ws->pondw;

				ws->pondw          -= soilwEXTRA;
				wf->pondw_to_soilw += soilwEXTRA;

				ws->soilw[layer] += soilwEXTRA;
				epv->VWC[layer] = ws->soilw[layer] / (sitec->soillayer_thickness[layer] * water_density);
				if (epv->VWC[layer] - sprop->VWCsat[layer] > CRIT_PRECwater && ws->pondw)
				{
					printf("\n");
					printf("ERROR in pondwANDrunoffD.c for tipping.c\n");
					errorCode = 1;
				}
			}
		}
			
		else
			flagEXTRA = 1;

		layer += 1;
	}

	/* Dunnian runoff */
	if (ws->GW_waterlogging > sprop->pondmax)
		pondmax = ws->GW_waterlogging;
	else
		pondmax = sprop->pondmax;

	if (ws->pondw > pondmax)
	{
		wf->pondw_to_runoff  = ws->pondw - pondmax;
		ws->pondw            = pondmax;	
	}


	/* pond_flag: flag of WARNING writing (only at first time) */
	if (!ctrl->pond_flag && ws->pondw) ctrl->pond_flag = 1;


	return (errorCode);
}
