/*
capillary_tipping.c
Calculation of percolation fluxes in capillary layers

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
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_func.h"
#include "bgc_constants.h"
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

int capillary_tipping(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{

	int errorCode = 0;
	int  ll, layer, flagRAIN;

	double VWC, soilw_sat1, soilw1, soilB;
	double INFILT, conductSAT_cmday;
	double VWCsat, VWCcrit, dz0, dz0_cm, dz1, HOLD, soilw0;

	double DC, DRN, EXCESS, VWCnew;

	int CFlayer = (int)sprop->CFlayer;
	int GWlayer = (int)sprop->GWlayer;


	/* infiltration: percolation from CFlayer above */
	if (CFlayer == 0)
		INFILT = wf->infiltPOT * mm_to_cm;
	else
		INFILT = wf->soilwPercol[CFlayer - 1] * mm_to_cm;


	/* -----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* I. normZone in CFlayer without GW */

	if (sprop->dz_NORMcf)
	{

		VWC = sprop->VWC_NORMcf;
		VWCsat = sprop->VWCsat[CFlayer];
		VWCcrit = MAX(sprop->VWCfc[CFlayer], sprop->VWCeq[CFlayer]);
		soilw0 = sprop->soilw_NORMcf;


		/* unit change from m to cm */
		dz0 = sprop->dz_NORMcf;
		dz0_cm = dz0 * m_to_cm;

		DC = sprop->drainCoeff[CFlayer];
		soilB = sprop->soilB[CFlayer];
		flagRAIN = wf->flagRAIN[CFlayer];

		/* saturated hydraulic conductivity in actual CFlayer (cm/day = m/s * 100 * sec/day) */
		conductSAT_cmday = sprop->hydrCONDUCTsat[CFlayer] * m_to_cm * nSEC_IN_DAY;

		if (!errorCode && calc_drainage(flagRAIN, soilB, INFILT, VWC, VWCsat, VWCcrit, dz0_cm, DC, conductSAT_cmday, &DRN, &EXCESS, &VWCnew))
		{
			printf("\n");
			printf("ERROR calc_drainage.c for capillary_tipping.c\n");
			errorCode = 1;
		}

		/* water flux: cm/day to kg/(m2*day) */
		wf->soilwPercol_NORMvsCAPILcf = DRN / mm_to_cm;

		/* state update: with new VWC calcualte soilw */
		sprop->VWC_NORMcf = VWCnew;
		sprop->soilw_NORMcf = sprop->VWC_NORMcf * dz0 * water_density;

		/* control */
		if (fabs(VWCnew * dz0 * water_density - soilw0 - ((INFILT - DRN - EXCESS) / mm_to_cm)) > CRIT_PREC_lenient)
		{
			printf("\n");
			printf("ERROR in capillary_tipping.c: soilw udpate\n");
			errorCode = 1;
		}

		INFILT = DRN;




		/* if there is excess water, redistribute it in layers above */
		if (EXCESS > 0)
		{
			for (ll = CFlayer - 1; ll >= 0; ll--)
			{
				dz1 = sitec->soillayer_thickness[ll];
				soilw_sat1 = sprop->VWCsat[ll] * dz1 * water_density;
				soilw1 = epv->VWC[ll] * dz1 * water_density;
				HOLD = MIN((soilw_sat1 - soilw1) * mm_to_cm, EXCESS);
				soilw1 += HOLD / mm_to_cm;
				if (soilw1 > soilw_sat1)
				{
					printf("\n");
					printf("ERROR in capillary_tipping.c: EXCESS calculation\n");
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
	}

	/* -----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* 2. capillZone in CFlayer without GW  */

	if (sprop->dz_CAPILcf)
	{
		VWC = sprop->VWC_CAPILcf;
		VWCsat = sprop->VWCsat[CFlayer];
		VWCcrit = sprop->VWCsat[CFlayer];
		soilw0 = sprop->soilw_CAPILcf;

		/* unit change from m to cm */
		dz0 = sprop->dz_CAPILcf;
		dz0_cm = dz0 * m_to_cm;

		DC = sprop->drainCoeff[CFlayer];
		soilB = sprop->soilB[CFlayer];
		flagRAIN = wf->flagRAIN[CFlayer];

		/* hydraulic conductivity in actual CFlayer (cm/day = m/s * 100 * sec/day) */
		conductSAT_cmday = sprop->hydrCONDUCTsat[CFlayer] * m_to_cm * nSEC_IN_DAY;

		if (!errorCode && calc_drainage(flagRAIN, soilB, INFILT, VWC, VWCsat, VWCcrit, dz0_cm, DC, conductSAT_cmday, &DRN, &EXCESS, &VWCnew))
		{
			printf("\n");
			printf("ERROR calc_drainage.c for capillary_tipping.c\n");
			errorCode = 1;
		}

		/* water flux: cm/day to kg/(m2*day) */
		wf->soilwPercol[CFlayer] = DRN / mm_to_cm;
		

		/* state update: with new VWC calcualte soilw */
		sprop->VWC_CAPILcf = VWCnew;
		sprop->soilw_CAPILcf = sprop->VWC_CAPILcf * dz0 * water_density;

		/* control */
		if (fabs(VWCnew * dz0 * water_density - soilw0 - ((INFILT - DRN - EXCESS) / mm_to_cm)) > CRIT_PREC_lenient)
		{
			printf("\n");
			printf("ERROR in capillary_tipping.c: soilw update\n");
			errorCode = 1;
		}


		INFILT = DRN;

		/* if there is excess water, redistribute it in layers above */
		if (EXCESS > CRIT_PREC_lenient)
		{
			/* first of all: normal layer is saturated */
			if (sprop->dz_NORMcf)
			{
				dz1 = sprop->dz_NORMcf;
				soilw_sat1 = sprop->VWCsat[CFlayer] * dz1 * water_density;
				soilw1 = sprop->VWC_NORMcf * dz1 * water_density;
				HOLD = MIN((soilw_sat1 - soilw1) * mm_to_cm, EXCESS);
				soilw1 += HOLD / mm_to_cm;
				if (soilw1 > soilw_sat1)
				{
					printf("\n");
					printf("ERROR in tipping.c: EXCESS calculation\n");
					errorCode = 1;
				}
				sprop->soilw_NORMcf = soilw1;
				sprop->VWC_NORMcf = sprop->soilw_NORMcf / dz1 / water_density;

				wf->soilwPercol_NORMvsCAPILcf -= EXCESS / mm_to_cm;

				EXCESS = EXCESS - HOLD;
			}
			if (EXCESS > CRIT_PREC_lenient)
			{
				for (ll = CFlayer - 1; ll >= 0; ll--)
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
					epv->VWC[ll] = ws->soilw[ll] / sitec->soillayer_thickness[ll] / water_density;

					wf->soilwPercol[ll] -= EXCESS / mm_to_cm;

					EXCESS = EXCESS - HOLD;

				}
			}
			/* if too much pondwater -> runoff */
			wf->soilw_to_pondw += EXCESS / mm_to_cm;
		}

	}
	else
	{ 
		wf->soilwPercol[CFlayer] = wf->soilwPercol_NORMvsCAPILcf;
		wf->soilwPercol_NORMvsCAPILcf = 0;

	}

	/* state update: with new VWC calcualte soilw */
	ws->soilw[CFlayer] = sprop->soilw_NORMcf + sprop->soilw_CAPILcf;
	epv->VWC[CFlayer] = ws->soilw[CFlayer] / (sitec->soillayer_thickness[CFlayer] * water_density);

	/* -----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* 3. Capillary layers betwen CFlayer and GWlayer */

	INFILT = wf->soilwPercol[CFlayer] * mm_to_cm;

	for (layer = CFlayer+1; layer < GWlayer; layer++)
	{

		VWC = epv->VWC[layer];
		VWCsat = sprop->VWCsat[layer];
		VWCcrit = sprop->VWCsat[layer];
		soilw0 = ws->soilw[layer];

		/* unit change from m to cm */
		dz0 = sitec->soillayer_thickness[layer];
		dz0_cm = dz0 * m_to_cm;

		DC = sprop->drainCoeff[layer];
		soilB = sprop->soilB[layer];
		flagRAIN = wf->flagRAIN[layer];

		/* hydraulic conductivity in actual layer (cm/day = m/s * 100 * sec/day) */
		conductSAT_cmday = sprop->hydrCONDUCTsat[layer] * m_to_cm * nSEC_IN_DAY;

		if (!errorCode && calc_drainage(flagRAIN, soilB, INFILT, VWC, VWCsat, VWCcrit, dz0_cm, DC, conductSAT_cmday, &DRN, &EXCESS, &VWCnew))
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
			printf("ERROR in capillary_tipping.c: soilw update\n");
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
					printf("ERROR in capillary_tipping.c: EXCESS calculation\n");
					errorCode = 1;
				}
				
				ws->soilw[ll] = soilw1;
				epv->VWC[ll] = ws->soilw[ll] / dz1 / water_density;
				if (ll == CFlayer)
				{
					sprop->soilw_NORMcf = ws->soilw[ll] * sprop->ratioNORMcf;
					sprop->soilw_CAPILcf = ws->soilw[ll] * sprop->ratioCAPILcf;
					if (sprop->dz_NORMcf) sprop->VWC_NORMcf = sprop->soilw_NORMcf / sprop->dz_NORMcf / water_density;
					if (sprop->dz_CAPILcf) sprop->VWC_CAPILcf = sprop->soilw_CAPILcf / sprop->dz_CAPILcf / water_density;
				}

				wf->soilwPercol[ll] -= EXCESS / mm_to_cm;

				EXCESS = EXCESS - HOLD;

			}
			/* if too much pondwater -> runoff */
			wf->soilw_to_pondw += EXCESS / mm_to_cm;
		}



	} /* END FOR (layer) */



	return (errorCode);

}
