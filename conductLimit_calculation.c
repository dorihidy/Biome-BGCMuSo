/* 
conductLimit_calculations.c
calculate the limitation factors of conductance
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2022, D. Hidy [dori.hidy@gmail.com]
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
#include "bgc_struct.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"


int conductLimit_calculations(const soilprop_struct* sprop, const epconst_struct* epc, epvar_struct* epv)
{
	int errorCode = 0;
	int layer;
	double VWCsat, VWCfc, VWCwp, VWC_WScrit1, VWC_WScrit2, PSI_WScrit1, PSI_WScrit2;


	VWCsat = VWCfc = VWCwp = VWC_WScrit1 = VWC_WScrit2 = PSI_WScrit1 = PSI_WScrit2 = 0;



	/* calculations layer by layer (due to different soil properties) */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		VWCsat = sprop->VWCsat[layer];
		VWCfc = sprop->VWCfc[layer];
		VWCwp = sprop->VWCwp[layer];

		/* VWCratio_WScrit1 = -9999 - no drought limitation defined by the User */
		if (epc->VWCratio_WScrit1 != DATA_GAP)
		{

			VWC_WScrit1 = VWCwp + epc->VWCratio_WScrit1 * (VWCfc - VWCwp);
			PSI_WScrit1 = sprop->PSIsat[layer] * pow((VWC_WScrit1 / VWCsat), -1 * sprop->soilB[layer]);

		}
		else
		{
			VWC_WScrit1 = DATA_GAP;
			PSI_WScrit1 = DATA_GAP;
		}

		/* VWCratio_WScrit2 = -9999 - no anoxic limitation defined by the User */
		if (epc->VWCratio_WScrit2 != DATA_GAP)
		{
			VWC_WScrit2 = VWCfc + epc->VWCratio_WScrit2 * (VWCsat - VWCfc);
			PSI_WScrit2 = sprop->PSIsat[layer] * pow((VWC_WScrit2 / VWCsat), -1 * sprop->soilB[layer]);
		}
		else
		{
			VWC_WScrit2 = DATA_GAP;
			PSI_WScrit2 = DATA_GAP;
		}



		/* CONTROL */
		if (!errorCode && VWC_WScrit2 > VWCsat)
		{
			if (VWC_WScrit2 - VWCsat > CRIT_PREC_lenient)
			{
				printf("\n");
				printf("ERROR in conductLimit_calculations.c: VWC_WScrit2 data is greater than saturation value in layer:%i\n", layer);
			}
			else
				VWC_WScrit2 = VWCsat;
		}

		if (!errorCode && VWC_WScrit2 < VWC_WScrit1)
		{
			if (VWC_WScrit1 - VWC_WScrit2 > CRIT_PREC_lenient)
			{
				printf("\n");
				printf("ERROR in conductLimit_calculations.c: VWC_WScrit1 data is greater then VWC_WScrit2 data in layer:%i\n", layer);
				errorCode = 1;
			}
			else
				VWC_WScrit1 = VWC_WScrit2;
		}

		epv->VWC_WScrit1[layer] = VWC_WScrit1;
		epv->VWC_WScrit2[layer] = VWC_WScrit2;
		epv->PSI_WScrit1[layer] = PSI_WScrit1;
		epv->PSI_WScrit2[layer] = PSI_WScrit2;


	}
		
	return (errorCode);
}
