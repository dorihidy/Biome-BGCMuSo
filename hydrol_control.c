/* 
hydrol_control.c
control of VWC values

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

int hydrol_control(siteconst_struct* sitec, soilprop_struct* sprop, wstate_struct* ws, wflux_struct* wf, epvar_struct* epv)
{


	int errorCode=0;
	int layer, GWlayer, CFlayer;


	GWlayer = (int)sprop->GWlayer;
	CFlayer = (int)sprop->CFlayer;


	/*--------------------------------------------------------------------------*/
	/* 1. control to avoid irrealistic PET */
	if (wf->ET - wf->PET > CRIT_PREC_lenient)
	{
		printf("ERROR in potential evaporation calculation (ET > PET) in hydrol_control.c of multilayer_hydrolprocess.c\n");
		errorCode = 1;
	}

	/*--------------------------------------------------------------------------*/
	/* 2. CONTROL and calculating averages - unrealistic VWC content (higher than saturation value or less then hygroscopic) - GWtest */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
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
					printf("ERROR in soil water content calculation (hydrol_control.c of multilayer_hydrolprocess.c) - actual is less than hygroscopic value\n");
					errorCode = 1;
				}
			}

		}

		if (epv->VWC[layer] - sprop->VWCsat[layer] > CRIT_PREC)
		{
			if (!errorCode)
			{
				printf("\n");
				printf("ERROR in soil water content calculation - actual is higher than saturation value (hydrol_control.c of multilayer_hydrolprocess)\n");
				errorCode = 1;
			}
		}

		if (epv->VWC[layer] && epv->VWC[layer] / epv->VWC[layer] != 1)
		{
			if (!errorCode)
			{
				printf("\n");
				printf("ERROR in soil water content calculation - invalid VWC value (hydrol_control.c of multilayer_hydrolprocess)\n");
				errorCode = 1;
			}
		}


	}

	/*--------------------------------------------------------------------------*/
	/* 3. in case of groundwater - control of virtaul layers */
	if (sprop->GWlayer != DATA_GAP)
	{

		if ((sprop->dz_NORMcf && sprop->dz_NORMgw))
		{
			printf("\n");
			printf("ERROR in tipping calculation in case of GW in hydrol_control.c of multilayer_hydrolprocess.c\n");
			errorCode = 1;
		}

		if (sprop->VWC_NORMgw && (sprop->VWC_NORMgw - sprop->VWCsat[GWlayer] > CRIT_PREC_lenient || sprop->VWChw[GWlayer] - sprop->VWC_NORMgw > CRIT_PREC_lenient || sprop->VWC_NORMgw / sprop->VWC_NORMgw != 1))
		{
			printf("\n");
			printf("ERROR in soil water content calculation - invalid VWC value (hydrol_control.c of multilayer_hydrolprocess)\n");
			errorCode = 1;
		}

		if (sprop->VWC_CAPILgw && (sprop->VWC_CAPILgw - sprop->VWCsat[GWlayer] > CRIT_PREC_lenient || sprop->VWChw[GWlayer] - sprop->VWC_CAPILgw > CRIT_PREC_lenient || sprop->VWC_CAPILgw / sprop->VWC_CAPILgw != 1))
		{
			printf("\n");
			printf("ERROR in soil water content calculation - invalid VWC value (hydrol_control.c of multilayer_hydrolprocess)\n");
			errorCode = 1;
		}

		if (sprop->VWC_NORMcf && (sprop->VWC_NORMcf - sprop->VWCsat[CFlayer] > CRIT_PREC_lenient || sprop->VWChw[CFlayer] - sprop->VWC_NORMcf > CRIT_PREC_lenient || sprop->VWC_NORMcf / sprop->VWC_NORMcf != 1))
		{
			printf("\n");
			printf("ERROR in soil water content calculation - invalid VWC value (hydrol_control.c of multilayer_hydrolprocess)\n");
			errorCode = 1;
		}

		if (sprop->VWC_CAPILcf && (sprop->VWC_CAPILcf - sprop->VWCsat[CFlayer] > CRIT_PREC_lenient || sprop->VWC_CAPILcf < CRIT_PREC_lenient || sprop->VWC_CAPILcf / sprop->VWC_CAPILcf != 1))
		{
			printf("\n");
			printf("ERROR in soil water content calculation - invalid VWC value (hydrol_control.c of multilayer_hydrolprocess)\n");
			errorCode = 1;
		}

		if (CFlayer != GWlayer && (wf->GWdischarge[CFlayer] - wf->GWdischargeCAPILcf - wf->GWdischargeNORMcf > CRIT_PREC_lenient))
		{
			printf("\n");
			printf("ERROR in soil GWdischarge partitioning (hydrol_control.c of multilayer_hydrolprocess.c)\n");
			errorCode = 1;
		}

		if (wf->GWdischarge[GWlayer] - wf->GWdischargeCAPILgw - wf->GWdischargeNORMgw > CRIT_PREC_lenient)
		{
			printf("\n");
			printf("ERROR in soil GWdischarge partitioning (hydrol_control.c of multilayer_hydrolprocess.c)\n");
			errorCode = 1;
		}

		if (CFlayer != GWlayer && (wf->soilwDiffus[CFlayer - 1] - wf->soilwDiffus_aboveCFlayer_vs_CAPILcf - wf->soilwDiffus_aboveCFlayer_vs_NORMcf > CRIT_PREC_lenient))
		{
			printf("\n");
			printf("ERROR in  soilwDiffus partitioning in soilwDiffus_aboveCFlayer_vs_CAPILcf and wf->soilwDiffus_aboveCFlayer_vs_NORMcf (hydrol_control.c of multilayer_hydrolprocess.c)\n");
			errorCode = 1;
		}

		if (CFlayer != GWlayer && (wf->soilwDiffus[CFlayer] - wf->soilwDiffusCAPILcf - wf->soilwDiffusNORMcf > CRIT_PREC_lenient))
		{
			printf("\n");
			printf("ERROR in  soilwDiffus partitioning  in soilwDiffusCAPILcf and wf->soilwDiffusNORMcf (hydrol_control.c of multilayer_hydrolprocess.c)\n");
			errorCode = 1;
		}

		if (CFlayer == GWlayer && wf->soilwDiffus[GWlayer - 1] - wf->soilwDiffus_aboveGWlayer_vs_CAPILgw - wf->soilwDiffus_aboveGWlayer_vs_NORMgw > CRIT_PREC_lenient)
		{
			printf("\n");
			printf("ERROR in soilwDiffus partitioning in wf->soilwDiffus_aboveGWlayer_vs_CAPILgw ansd wf->soilwDiffus_aboveGWlayer_vs_NORMgw (hydrol_control.c of multilayer_hydrolprocess.c)\n");
			errorCode = 1;
		}


	}


	return(errorCode);
}

