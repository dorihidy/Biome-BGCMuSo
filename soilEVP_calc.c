/* 
soilEVP.c
Calculate bare soil evaporation

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
#include "bgc_func.h"
#include "bgc_constants.h"
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

int soilEVP_calc(control_struct* ctrl, const siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{
	int errorCode=0;
	int CFlayer, GWlayer;
	double EVP_lack, soilw_hw0, soilw_diff;
	double ratioNORM, ratioCAPIL, diffNORM, diffCAPIL, soilwAVAIL_NORMcf, soilwAVAIL_CAPILcf;
	CFlayer = (int)sprop->CFlayer;
	GWlayer = (int)sprop->GWlayer;
	soilwAVAIL_NORMcf = soilwAVAIL_CAPILcf = 0;
	/* SOILEVAP UPDATE: if GW-table is in the top soil layer, the source of evaporation is the GW-table */
	

	if (ws->pondw == 0 && wf->EVPsoilw)
	{
		if (GWlayer == 0)
		{
			soilw_diff = wf->EVPsoilw - wf->potEVPsurface;
			if (soilw_diff) wf->EVPsoilw -= soilw_diff;
			wf->EVPfromGW = wf->EVPsoilw;
		}
		else
		{
			/* soilw evaporation - limitation: hygroscopic water */
			soilw_hw0 = sprop->VWChw[0] * sitec->soillayer_thickness[0] * water_density;
			EVP_lack = wf->EVPsoilw - (ws->soilw[0] - soilw_hw0);

			/* theoretical lower limit of water content: hygroscopic water content. */
			if (EVP_lack > 0)
			{
				if (ws->EVPsurface2cum >= sprop->soilEVPcrit)
				{
					if (ws->EVPsurface2cum > wf->EVPsoilw)
					{
						ws->EVPsurface2cum = ws->EVPsurface2cum - EVP_lack;
						epv->DSR = pow((ws->EVPsurface2cum / sprop->coeff_EVPcum), 2);
					}
					else
					{
						ws->EVPsurface1cum = ws->EVPsurface1cum - (wf->EVPsoilw - ws->EVPsurface2cum);

						ws->EVPsurface2cum = (ws->soilw[0] - soilw_hw0) - sprop->soilEVPcrit;
						if (ws->EVPsurface2cum < 0) ws->EVPsurface2cum = 0;

						ws->EVPsurface1cum = ws->EVPsurface1cum + (ws->soilw[0] - soilw_hw0);
						if (ws->EVPsurface1cum < sprop->soilEVPcrit) ws->EVPsurface2cum = sprop->soilEVPcrit;

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

			/* if capillary zone exists in unsaturated zone (not in GWlayer) and capillary zone is in the top soil layer */
			if ((sprop->dz_CAPILcf + sprop->dz_NORMcf) && CFlayer == 0)
			{
				if (sprop->soilw_NORMcf) soilwAVAIL_NORMcf = MAX(0, sprop->soilw_NORMcf - sprop->VWChw[CFlayer] * sprop->dz_NORMcf * water_density);
				if (sprop->dz_CAPILcf) soilwAVAIL_CAPILcf = MAX(0, sprop->soilw_CAPILcf - sprop->VWChw[CFlayer] * sprop->dz_CAPILcf * water_density);
				if (soilwAVAIL_CAPILcf + soilwAVAIL_NORMcf)
				{
					ratioNORM = soilwAVAIL_NORMcf / (soilwAVAIL_NORMcf + soilwAVAIL_CAPILcf);
					ratioCAPIL = soilwAVAIL_CAPILcf / (soilwAVAIL_NORMcf + soilwAVAIL_CAPILcf);
				}
				else
				{
					printf("\n");
					printf("ERROR in ratio calculation in multilayer_transpiration.c\n");
					errorCode = 1;
				}
				if (fabs(1 - ratioNORM - ratioCAPIL) > CRIT_PREC)
				{
					printf("\n");
					printf("ERROR in ratio calculation in multilayer_transpiration.c\n");
					errorCode = 1;
				}

				diffNORM = wf->EVPsoilw * ratioNORM;
				diffCAPIL = wf->EVPsoilw * ratioCAPIL;

				/* NORM zone: minimum value: hygroscopic water */
				soilw_hw0 = sprop->VWChw[0] * sprop->dz_NORMcf * water_density;
				if (sprop->soilw_NORMcf - diffNORM - soilw_hw0 < 0)
				{
					if (fabs(sprop->soilw_NORMcf - diffNORM - soilw_hw0) > CRIT_PREC)
					{
						printf("ERROR in content_NORMgw calculation in soilEVP.c\n");
						errorCode = 1;
					}
					else
						diffNORM = sprop->soilw_NORMcf - soilw_hw0;
				}

				sprop->soilw_NORMcf -= diffNORM;
				if (sprop->dz_NORMcf) sprop->VWC_NORMcf = sprop->soilw_NORMcf / (sprop->dz_NORMcf * water_density);


				/* CAPIL zone: minimum value: hygroscopic water */
				soilw_hw0 = sprop->VWChw[0] * sprop->dz_CAPILcf * water_density;
				if (sprop->soilw_CAPILcf - diffCAPIL - soilw_hw0 < 0)
				{
					if (fabs(sprop->soilw_CAPILcf - diffCAPIL - soilw_hw0) > CRIT_PREC)
					{
						printf("ERROR in content_CAPILgw calculation in soilEVP.c\n");
						errorCode = 1;
					}
					else
						diffCAPIL = sprop->soilw_CAPILcf - soilw_hw0;
				}
				sprop->soilw_CAPILcf -= diffCAPIL;
				if (sprop->dz_CAPILcf) sprop->VWC_CAPILcf = sprop->soilw_CAPILcf / (sprop->dz_CAPILcf * water_density);

				wf->EVPsoilwNORMcf = diffNORM;
				wf->EVPsoilwCAPILcf = diffCAPIL;
			}


		}
	}
	else
		wf->EVPsoilw = 0;



		
	return(errorCode);
}

