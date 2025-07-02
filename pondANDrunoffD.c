/* 
pondANDrunoffD.c
state update of pond water, calculation of Dunnian runoff, pond water limitation

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


int pondANDrunoffD(control_struct* ctrl, siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{

	/* internal variables */
	int errorCode, layer, flagEXTRA;
	double pondmax, soilwEXTRA;
	double soilw_NORMgw_full, soilw_CAPILgw_full, ratioNORM, ratioCAPIL;

	errorCode=layer=flagEXTRA=0;


	/*--------------------------------------*/
	/* Water flux from soil to pond */
	

	ws->pondw += wf->soilw_to_pondw + wf->prcp_to_pondw;

	/*--------------------------------------*/
	/* Dunnian runoff */
	if (ws->GW_waterlogging > sprop->pondmax)
		pondmax = ws->GW_waterlogging;
	else
		pondmax = sprop->pondmax;

	if (sprop->FLD > sprop->pondmax)
		pondmax = ws->GW_waterlogging;
	else
		pondmax = sprop->pondmax;

	

	if (ws->pondw > pondmax)
	{
		wf->pondw_to_runoff = ws->pondw - pondmax;
		ws->pondw = pondmax;
	}



	
	/*--------------------------------------*/
	/* saturation of top soil layers: correction of water flux from pond to soil */

	if (ws->pondw > 0)
	{
		flagEXTRA = 0;
		soilwEXTRA = (sprop->VWCsat[0] - epv->VWC[0]) * sitec->soillayer_thickness[0] * water_density;
		if (soilwEXTRA > sprop->hydrCONDUCTsat[0] * nSEC_IN_DAY * m_to_mm) soilwEXTRA = sprop->hydrCONDUCTsat[0] * nSEC_IN_DAY* m_to_mm;
		if (fabs(soilwEXTRA) < CRIT_PREC_lenient) soilwEXTRA = 0;
		if (soilwEXTRA)
		{
			if (soilwEXTRA > ws->pondw) soilwEXTRA = ws->pondw;

			ws->pondw -= soilwEXTRA;
			wf->pondw_to_soilw += soilwEXTRA;

			ws->soilw[0] += soilwEXTRA;
			epv->VWC[0] = ws->soilw[0] / (sitec->soillayer_thickness[0] * water_density);

			if ((int)sprop->CFlayer == 0)
			{
				soilw_NORMgw_full = sprop->VWCsat[0] * sprop->dz_NORMcf * water_density;
				soilw_CAPILgw_full = sprop->VWCsat[0] * sprop->dz_CAPILcf * water_density;
				ratioNORM = soilw_NORMgw_full/ (soilw_NORMgw_full+ soilw_CAPILgw_full);
				ratioCAPIL = soilw_CAPILgw_full / (soilw_NORMgw_full + soilw_CAPILgw_full);
				sprop->soilw_NORMcf += soilwEXTRA * ratioNORM;
				sprop->soilw_CAPILcf += soilwEXTRA * ratioCAPIL;
				if (sprop->dz_NORMcf) sprop->VWC_NORMcf = sprop->soilw_NORMcf / sprop->dz_NORMcf / water_density;
				if (sprop->dz_CAPILcf) sprop->VWC_CAPILcf = sprop->soilw_CAPILcf / sprop->dz_CAPILcf / water_density;
			}
			if (epv->VWC[0] - sprop->VWCsat[0] > CRIT_PREC_lenient || sprop->VWC_CAPILcf - sprop->VWCsat[0] > CRIT_PREC_lenient || sprop->VWC_NORMcf - sprop->VWCsat[0] > CRIT_PREC_lenient)
			{
				printf("\n");
				printf("ERROR in pondwANDrunoffD.c for tipping.c\n");
				errorCode = 1;
			}

			if ((int)sprop->GWlayer==0)
			{
				soilw_NORMgw_full = sprop->VWCsat[0] * sprop->dz_NORMgw * water_density;
				soilw_CAPILgw_full = sprop->VWCsat[0] * sprop->dz_CAPILgw * water_density;
				ratioNORM = (soilw_NORMgw_full - sprop->soilw_NORMgw) / soilwEXTRA;
				ratioCAPIL = (soilw_CAPILgw_full - sprop->soilw_CAPILgw) / soilwEXTRA;
				sprop->soilw_NORMgw += soilwEXTRA * ratioNORM;
				sprop->soilw_CAPILgw += soilwEXTRA * ratioCAPIL;
				if (sprop->dz_NORMgw) sprop->VWC_NORMgw = sprop->soilw_NORMgw / sprop->dz_NORMgw / water_density;
				if (sprop->dz_CAPILgw) sprop->VWC_CAPILgw = sprop->soilw_CAPILgw / sprop->dz_CAPILgw / water_density;
			}
			if (epv->VWC[0] - sprop->VWCsat[0] > CRIT_PREC_lenient || sprop->VWC_CAPILgw - sprop->VWCsat[0] > CRIT_PREC_lenient || sprop->VWC_NORMgw - sprop->VWCsat[0] > CRIT_PREC_lenient)
			{
				printf("\n");
				printf("ERROR in pondwANDrunoffD.c for tipping.c\n");
				errorCode = 1;
			}
		}
	}

	/*--------------------------------------*/
	/* Pond water evaporation: water stored on surface  */

	if (ws->pondw)
	{
		if (wf->potEVPsurface - wf->EVPsoilw < ws->pondw)
			wf->EVPpondw = wf->potEVPsurface - wf->EVPsoilw;
		else
			wf->EVPpondw = ws->pondw;
	}




	/* Water sink: pond water evaporation */
	if (ws->GW_waterlogging > 0)
		ws->GWsrc_W += wf->EVPpondw;
	else
		ws->pondw -= wf->EVPpondw;



	/* pond_flag: flag of WARNING writing (only at first time) */
	if (!ctrl->pond_flag && ws->pondw) ctrl->pond_flag = 1;


	return (errorCode);
}
