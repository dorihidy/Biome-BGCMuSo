/* 
CWDextract.c
CWD-extract management type - remove CWD from forest

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
#include "pointbgc_struct.h"
#include "bgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"

int CWDextract(control_struct* ctrl, const CWDextract_struct* CWE, cstate_struct* cs, nstate_struct*ns, cflux_struct* cf, nflux_struct* nf)
{

	int errorCode=0;
	int md, year, layer;
	double removePROP, cwdc_to_CWE;			            

	removePROP = 0;
	year = ctrl->simstartyear + ctrl->simyr;
	md = CWE->mgmdCWE-1;

	cf->cwdc_to_CWE = 0;
	nf->cwdn_to_CWE = 0;
	layer = 0;

	/* 2. we assume that the transfer pools contain the palnt material of seeds. Therefore planting increase the transfer pools */ 
	if (CWE->CWE_num && md >= 0)
	{
		if (year == CWE->CWEyear_array[md] && ctrl->month == CWE->CWEmonth_array[md] && ctrl->day == CWE->CWEday_array[md])
		{

			removePROP = CWE->removePROP_CWE[md] / 100.;
			cwdc_to_CWE = cs->cwdCabove_total * removePROP;

			/* if removed material is greater than 0 */
			if (cwdc_to_CWE)
			{
				/* if the first layer covers CWE */
				while (cwdc_to_CWE > CRIT_PREC && layer < N_SOILLAYERS)
				{
					if (cs->cwdCabove[layer] >= cwdc_to_CWE)
					{
						cf->cwdc_to_CWElayer[layer] = cwdc_to_CWE;
						cwdc_to_CWE = 0;
					}
					else
					{
						cf->cwdc_to_CWElayer[layer] = cs->cwdCabove[layer];
						cwdc_to_CWE -= cf->cwdc_to_CWElayer[layer];

					}
					if (cs->cwdc[layer]) nf->cwdn_to_CWElayer[layer] = cf->cwdc_to_CWElayer[layer] * ns->cwdn[layer] / cs->cwdc[layer];

					/* state update */
					cs->cwdc[layer] -= cf->cwdc_to_CWElayer[layer];
					ns->cwdn[layer] -= nf->cwdn_to_CWElayer[layer];

					if (cs->cwdc[layer] < CRIT_PREC) cs->cwdc[layer] = 0;
					if (ns->cwdn[layer] < CRIT_PREC) ns->cwdn[layer] = 0;

					cf->cwdc_to_CWE += cf->cwdc_to_CWElayer[layer];
					nf->cwdn_to_CWE += nf->cwdn_to_CWElayer[layer];


					cs->cwdCabove[layer] -= cf->cwdc_to_CWElayer[layer];
					cs->cwdCabove_total -= cf->cwdc_to_CWElayer[layer];


					layer += 1;
				}
			}


			cs->CWEsnk_C += cf->cwdc_to_CWE;
			ns->CWEsnk_N += nf->cwdn_to_CWE;

		}
		else
			for (layer = 0; layer < N_SOILLAYERS; layer++) cf->cwdc_to_CWElayer[layer] = 0;

	}

	/* nitrogen fluxes from the ratio of CWD in given layer */

	

	
	return (errorCode);
}