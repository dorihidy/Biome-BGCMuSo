 /*
tipping.c
Calculation of percolation and diffusion fluxes (optionally GWdischarge)

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
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_func.h"
#include "bgc_constants.h"
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

int tipping(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{

	int errorCode = 0;
	int layer, ll, N_NAGlayers;

	double VWC, soilw_sat1, soilw1;
	double INFILT, conduct_cmday, conductSAT_cmday;
	double VWCsat, VWCfc, dz0, dz1, dz0_cm, HOLD;

	double DC, DRN, EXCESS, VWCnew, soilw0;

	/* tipping is used only for layers without GW */
	if (sprop->GWD == DATA_GAP)
		N_NAGlayers = N_SOILLAYERS;
	else
		N_NAGlayers = (int) sprop->CFlayer;

	
	/* --------------------------------------------------------------------------------------------------------------------*/
	/* 1.PERCOLATION */

	INFILT = wf->infiltPOT * mm_to_cm;

	
	/* BEGIN LOOP: layer from top to CFlayer (is GW is exists) */
	for (layer = 0; layer < N_NAGlayers; layer++)
	{

		VWC    = epv->VWC[layer];
		VWCsat = sprop->VWCsat[layer];
		VWCfc  = MAX(sprop->VWCfc[layer], sprop->VWCeq[layer]);
		soilw0 = ws->soilw[layer];

		/* unit change from m to cm */
		dz0 = sitec->soillayer_thickness[layer];
		dz0_cm = dz0 * m_to_cm;

		DC = sprop->drainCoeff[layer];

		/* hydraulic conductivity in actual layer (cm/day = m/s * 100 * sec/day) */
		conductSAT_cmday = sprop->hydrCONDUCTsat[layer] * m_to_cm * nSEC_IN_DAY;
		conduct_cmday = conductSAT_cmday;


		if (!errorCode && calc_drainage(wf->flagRAIN, INFILT, VWC, VWCsat, VWCfc, dz0_cm, DC, conduct_cmday, &DRN, &EXCESS, &VWCnew))
		{
			printf("\n");
			printf("ERROR calc_drainage.c for tipping.c\n");
			errorCode = 1;
		}

		/* water flux: cm/day to kg/(m2*day) */
		wf->soilwPercol[layer] = DRN / mm_to_cm;


		/* state update: with new VWC calcualte soilw */
		epv->VWC[layer] = VWCnew;
		ws->soilw[layer] = epv->VWC[layer] * sitec->soillayer_thickness[layer] * water_density;
	
		/* control */
		if (fabs(VWCnew * dz0 * water_density - soilw0 - ((INFILT - DRN - EXCESS) / mm_to_cm)) > CRIT_PREC_lenient)
		{
			printf("\n");
			printf("ERROR in tipping.c: soilw update\n");
			errorCode = 1;
		}

		INFILT = DRN;

		/* if there is excess water, redistribute it in layers above */
		if (EXCESS > CRIT_PREC_lenient)
		{
			for (ll = layer - 1; ll >= 0; ll--)
			{
				dz1 = sitec->soillayer_thickness[ll];
				soilw_sat1 = sprop->VWCsat[ll] * dz1 * water_density;
				soilw1 = epv->VWC[ll] * dz1 * water_density;
				HOLD = MIN((soilw_sat1 - soilw1) * mm_to_cm, EXCESS);
				soilw1 += HOLD / mm_to_cm;
				if (soilw1 > soilw_sat1)
				{
					printf("\n");
					printf("ERROR in tipping.c: EXCESS calculation\n");
					errorCode = 1;
				}
				ws->soilw[ll] = soilw1;
				epv->VWC[ll] = ws->soilw[ll] / dz1 / water_density;

				wf->soilwPercol[ll] -= EXCESS / mm_to_cm;

				EXCESS = EXCESS - HOLD;

			}
			/* if too much pondwater -> runoff */
			wf->soilw_to_pondw += EXCESS / mm_to_cm;
		}



	} /* END FOR (layer) */




	return (errorCode);

}
