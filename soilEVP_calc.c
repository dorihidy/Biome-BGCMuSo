/* 
Elimit_and_PET.c
Calculate the limitation of soil evaporation (actual EVP), GWevap and update top soil water content

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
#include "bgc_func.h"
#include "bgc_constants.h"

int soilEVP_calc(control_struct* ctrl, const siteconst_struct* sitec,const soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{
	int errorCode=0;
	double EVP_lack, soilw_hw0, soilw_diff;
	
	/* SOILEVAP UPDATE: if GW-table is in the top soil layer, the source of evaporation is the GW-table */
	

	if (ws->pondw == 0)
	{
		if ((int)sprop->GWlayer == 0)
		{
			soilw_diff = wf->EVPsoilw - wf->potEVPsurface;
			if (soilw_diff) wf->EVPsoilw -= soilw_diff;
			wf->GWevap = wf->EVPsoilw;
		}
		else
		{
			/* soilw evaporation - limitation: hygroscopic water */
			soilw_hw0 = sprop->VWChw[0] * sitec->soillayer_thickness[0] * water_density;
			EVP_lack = wf->EVPsoilw - (ws->soilw[0] - soilw_hw0);

			/* theoretical lower limit of water content: hygroscopic water content. */
			if (EVP_lack > 0)
			{
				if (ws->EVPsurface2cum >= sprop->soilEVPlim)
				{
					if (ws->EVPsurface2cum > wf->EVPsoilw)
					{
						ws->EVPsurface2cum = ws->EVPsurface2cum - EVP_lack;
						epv->DSR = pow((ws->EVPsurface2cum / sprop->coeff_EVPcum), 2);
					}
					else
					{
						ws->EVPsurface1cum = ws->EVPsurface1cum - (wf->EVPsoilw - ws->EVPsurface2cum);

						ws->EVPsurface2cum = (ws->soilw[0] - soilw_hw0) - sprop->soilEVPlim;
						if (ws->EVPsurface2cum < 0) ws->EVPsurface2cum = 0;

						ws->EVPsurface1cum = ws->EVPsurface1cum + (ws->soilw[0] - soilw_hw0);
						if (ws->EVPsurface1cum < sprop->soilEVPlim) ws->EVPsurface2cum = sprop->soilEVPlim;

						epv->DSR = pow((ws->EVPsurface2cum / sprop->coeff_EVPcum), 2);
					}
				}
				else
				{
					ws->EVPsurface1cum = ws->EVPsurface1cum - wf->EVPsoilw + (ws->soilw[0] - soilw_hw0);
				}

				wf->EVPsoilw = (ws->soilw[0] - soilw_hw0);
				/* limitEVP_flag: flag of WARNING writing in log file (only at first time) */
				if (fabs(EVP_lack) > CRIT_PREC && !ctrl->limitEVP_flag) ctrl->limitEVP_flag = 1;
			}

			soilw_diff = wf->EVPsoilw - wf->potEVPsurface;

			if (soilw_diff > 0) wf->EVPsoilw -= soilw_diff;

			ws->soilw[0] -= (wf->EVPsoilw);	
			epv->VWC[0] = ws->soilw[0] / water_density / sitec->soillayer_thickness[0];


		}
	}
	else
		wf->EVPsoilw = 0;



		
	return(errorCode);
}

