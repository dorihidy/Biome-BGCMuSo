/*
groundwater_tipping.c
Calculation of percolation fluxes in groundwater layers

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

int groundwater_tipping(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{

	int errorCode = 0;
	int  ll;

	double VWC, soilw_sat1, soilw1;
	double INFILT, conduct_cmday, conductSAT_cmday;
	double VWCsat, VWCcrit, dz0, dz1, dz0_cm, HOLD, soilw0;

	double DC, DRN, EXCESS, VWCnew;

	int GWlayer = (int)sprop->GWlayer;
	int CFlayer = (int)sprop->CFlayer;


	/* infiltration: percolation from layer above */
	if (GWlayer == 0)
		INFILT = wf->infiltPOT * mm_to_cm;
	else
		INFILT = wf->soilwPercol[GWlayer - 1] * mm_to_cm;

	/* hydraulic conductivity in actual GWlayer (cm/day = m/s * 100 * sec/day) */
	conductSAT_cmday = sprop->hydrCONDUCTsat[GWlayer] * m_to_cm * nSEC_IN_DAY;
	conduct_cmday = conductSAT_cmday;

	/* -----------------------------*/
	/* 1. normZone */

	if (sprop->dz_NORMgw)
	{

		VWC = sprop->VWC_NORMgw;
		VWCsat = sprop->VWCsat[GWlayer];
		VWCcrit = MAX(sprop->VWCfc[GWlayer], sprop->VWCeq[GWlayer]);
		soilw0 = sprop->soilw_NORMgw;

		/* unit change from m to cm */
		dz0 = sprop->dz_NORMgw;
		dz0_cm = dz0 * m_to_cm;

		DC = sprop->drainCoeff[GWlayer];


		if (!errorCode && calc_drainage(wf->flagRAIN, INFILT, VWC, VWCsat, VWCcrit, dz0_cm, DC, conduct_cmday, &DRN, &EXCESS, &VWCnew))
		{
			printf("\n");
			printf("ERROR calc_drainage.c for tipping.c\n");
			errorCode = 1;
		}

		/* water flux: cm/day to kg/(m2*day) */
		wf->soilwPercol_NORMvsCAPILgw = DRN / mm_to_cm;

		/* state update: with new VWC calcualte soilw */
		sprop->VWC_NORMgw = VWCnew;
		sprop->soilw_NORMgw = sprop->VWC_NORMgw * dz0 * water_density;

		/* control */
		if (fabs(VWCnew * dz0 * water_density - soilw0 - ((INFILT - DRN - EXCESS) / mm_to_cm)) > CRIT_PREC_lenient)
		{
			printf("\n");
			printf("ERROR in tipping.c: soilw udpate\n");
			errorCode = 1;
		}

		INFILT = DRN;




		/* if there is excess water, redistribute it in layers above */
		if (EXCESS > CRIT_PREC_lenient)
		{
			for (ll = GWlayer - 1; ll >= 0; ll--)
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
				/* udapte of NORMcf and CAPILcf*/
				if (ll == CFlayer)
				{
					sprop->soilw_CAPILcf += soilw1 - ws->soilw[ll];
					soilw_sat1 = sprop->VWCsat[ll] * sprop->dz_CAPILcf * water_density;
					if (sprop->soilw_CAPILcf > soilw_sat1)
					{
						sprop->soilw_NORMcf += sprop->soilw_CAPILcf - soilw_sat1;
						sprop->soilw_CAPILcf = soilw_sat1;
					}
					if (sprop->dz_CAPILcf) sprop->VWC_CAPILcf = sprop->soilw_CAPILcf / sprop->dz_CAPILcf / water_density;
					if (sprop->dz_NORMcf) sprop->VWC_NORMcf = sprop->soilw_NORMcf / sprop->dz_NORMcf / water_density;
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

	/* -----------------------------*/
	/* 2. capillZone */

	if (sprop->dz_CAPILgw)
	{
		VWC = sprop->VWC_CAPILgw;
		VWCsat = sprop->VWCsat[GWlayer];
		VWCcrit = sprop->VWCsat[GWlayer];
		soilw0 = sprop->soilw_CAPILgw;

		/* unit change from m to cm */
		dz0 = sprop->dz_CAPILgw;
		dz0_cm = dz0 * m_to_cm;

		DC = sprop->drainCoeff[GWlayer];

		/* hydraulic conductivity in actual GWlayer (cm/day = m/s * 100 * sec/day) */
		conductSAT_cmday = sprop->hydrCONDUCTsat[GWlayer] * m_to_cm * nSEC_IN_DAY;
		conduct_cmday = conductSAT_cmday;


		if (!errorCode && calc_drainage(wf->flagRAIN, INFILT, VWC, VWCsat, VWCcrit, dz0_cm, DC, conduct_cmday, &DRN, &EXCESS, &VWCnew))
		{
			printf("\n");
			printf("ERROR calc_drainage.c for tipping.c\n");
			errorCode = 1;
		}

		/* water flux: cm/day to kg/(m2*day) */
		wf->GWrecharge_CAPILgw = DRN / mm_to_cm;


		/* state update: with new VWC calcualte soilw */
		sprop->VWC_CAPILgw = VWCnew;
		sprop->soilw_CAPILgw = sprop->VWC_CAPILgw * dz0 * water_density;
		/* control */
		if (fabs(VWCnew * dz0 * water_density - soilw0 - ((INFILT - DRN - EXCESS) / mm_to_cm)) > CRIT_PREC_lenient)
		{
			printf("\n");
			printf("ERROR in groundwater_tipping.c: soilw update\n");
			errorCode = 1;
		}


		INFILT = DRN;

		/* if there is excess water, redistribute it in layers above */
		if (EXCESS > 0)
		{
			/* first of all: normal layer is saturated */
			if (sprop->dz_NORMgw)
			{
				dz1 = sprop->dz_NORMgw;
				soilw_sat1 = sprop->VWCsat[GWlayer] * dz1 * water_density;
				soilw1 = sprop->VWC_NORMgw * dz1 * water_density;
				HOLD = MIN((soilw_sat1 - soilw1) * mm_to_cm, EXCESS);
				soilw1 += HOLD / mm_to_cm;
				if (soilw1 > soilw_sat1)
				{
					printf("\n");
					printf("ERROR in tipping.c: EXCESS calculation\n");
					errorCode = 1;
				}
	
				sprop->soilw_NORMgw = soilw1;
				sprop->VWC_NORMgw = sprop->soilw_NORMgw / dz1 / water_density;

				wf->soilwPercol_NORMvsCAPILgw -= EXCESS / mm_to_cm;
				EXCESS = EXCESS - HOLD;
			}
			if (EXCESS > CRIT_PREC_lenient)
			{ 
				for (ll = GWlayer - 1; ll >= 0; ll--)
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

					/* udapte of NORMcf and CAPILcf*/
					if (ll == CFlayer)
					{
						sprop->soilw_CAPILcf += soilw1 - ws->soilw[ll];
						soilw_sat1 = sprop->VWCsat[ll] * sprop->dz_CAPILcf * water_density;
						if (sprop->soilw_CAPILcf > soilw_sat1)
						{
							sprop->soilw_NORMcf += sprop->soilw_CAPILcf - soilw_sat1;
							sprop->soilw_CAPILcf = soilw_sat1;
						}
						if (sprop->dz_CAPILcf) sprop->VWC_CAPILcf = sprop->soilw_CAPILcf / sprop->dz_CAPILcf / water_density;
						if (sprop->dz_NORMcf) sprop->VWC_NORMcf = sprop->soilw_NORMcf / sprop->dz_NORMcf / water_density;
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


		/* state update: with new VWC calcualte soilw */
		ws->soilw[GWlayer] = sprop->soilw_NORMgw + sprop->soilw_CAPILgw + sprop->soilw_SATgw;
		epv->VWC[GWlayer] = ws->soilw[GWlayer] / (sitec->soillayer_thickness[GWlayer] * water_density);

	}
	else
	{
		/* no capillary zone in GW-layer */
		if (GWlayer - 1 == CFlayer)
		{
			wf->GWrecharge_CAPILcf = INFILT / mm_to_cm;
			wf->soilwPercol[CFlayer] = 0;
		}
		else
			wf->GWrecharge_lastCAPIL = INFILT / mm_to_cm;
	}


	return (errorCode);

}
