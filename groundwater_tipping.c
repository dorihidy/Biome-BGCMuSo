/*
groundwater_tipping.c
Calculation of percolation and diffusion fluxes 

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

int groundwater_tipping(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{

	int errorCode = 0;
	int ll;
	int GWlayer;
	double INFILT, DRN, GWR;

	double VWC, soilw_sat1, soilw1;
	double EXCESS, HOLD, DRAIN, DRMX, conduct_sat;

	double VWC_zoneNORM, VWC_zoneCAPIL, soilw_zoneNORM, soilw_zoneCAPIL, changeNORM1, changeCAPIL1, changeNORM2, changeCAPIL2;

	double dz0, dz1, VWCfc;



	GWlayer = (int)sprop->GWlayer;
	DRN = GWR = DRMX = 0;

	/* saturated hydraulic conductivity in actual layer (cm/day = m/s * 100 * sec/day) */
	conduct_sat = sprop->hydrCONDUCTsat[GWlayer] * m_to_cm * nSEC_IN_DAY;

	/* infiltration: percolation from layer above */
	if (GWlayer == 0)
		INFILT = wf->infiltPOT * mm_to_cm;
	else
		INFILT = wf->soilwPercol[GWlayer - 1] * mm_to_cm;

	/* -----------------------------*/
	/* 1. normZone */

	if (sprop->dz_zoneNORM > 0)
	{

		VWC = sprop->VWC_zoneNORM;
		dz0 = sprop->dz_zoneNORM * m_to_cm;
		VWCfc = sprop->VWCfc_base[GWlayer];


		/* [cm = m3/m3 * cm */
		HOLD = (sprop->VWCsat[GWlayer] - VWC) * dz0;


		/* 1.1.  IF: INFILT > HOLD */
		if (INFILT > 0.0 && INFILT > HOLD)
		{
			/* drainage from soil profile [cm = m3/m3 * cm ] */
			DRAIN = sprop->drainCoeff[GWlayer] * (sprop->VWCsat[GWlayer] - VWCfc) * dz0;

			/* drainage rate throug soil layer (cm/day) */
			DRN = INFILT - HOLD + DRAIN;


			/* drainage is limited: cm/h * h/day */
			if ((DRN - conduct_sat) > 0.0)
			{
				DRN = conduct_sat;
				DRAIN = DRN + HOLD - INFILT;
			}


			/* state update temporal varialbe */
			VWC = VWC + (INFILT - DRN) / dz0;

			/* above saturation - */
			if (VWC >= sprop->VWCsat[GWlayer])
			{
				EXCESS = (VWC - sprop->VWCsat[GWlayer]) * dz0;
				VWC = sprop->VWCsat[GWlayer];

	
				for (ll = GWlayer - 1; ll >= 0; ll--)
				{
					dz1 = sitec->soillayer_thickness[ll] * mm_to_cm;
					soilw_sat1 = sprop->VWCsat[ll] * sitec->soillayer_thickness[ll] * water_density * mm_to_cm;
					soilw1 = epv->VWC[ll] * sitec->soillayer_thickness[ll] * water_density * mm_to_cm;
					HOLD = MIN(soilw_sat1 - soilw1, EXCESS);
					ws->soilw[ll] += HOLD / mm_to_cm;
					epv->VWC[ll] = ws->soilw[ll] / sitec->soillayer_thickness[ll] / water_density;

					wf->soilwPercol[ll] = MAX(((wf->soilwPercol[ll] - EXCESS) / m_to_cm) * water_density, 0);

					EXCESS = EXCESS - HOLD;

				}
				/* if too much pondwater -> runoff */
				wf->soilw_to_pondw += EXCESS / mm_to_cm;
			}

			INFILT = DRN;
			

		} /* END IF: INFILT > HOLD */
		else
		/* 1.2 BEGIN ELSE: INFILT < HOLD */
		{ 

			VWC = VWC + INFILT / dz0;


			/* BEGIN IF-ELSE: VWC > FC */
			if (VWC >= VWCfc)
			{

				DRAIN = (VWC - VWCfc) * sprop->drainCoeff[GWlayer] * dz0;

				/* drainage rate throug soil layer (cm/day) */
				if (DRAIN > CRIT_PRECwater)
					DRN = DRAIN;
				else
					DRN = 0;


				/* drainage is limited */
				if ((DRN - conduct_sat) > 0.0)
				{
					DRN = conduct_sat;
					DRAIN = DRN;
				}


				VWC = VWC - DRN / dz0;
				INFILT = DRN;
			}
			else
			{
				INFILT = 0.0;
				DRN = 0.0;

			} /* END IF-ELSE: VWC > FC */

			

		} /* END ELSE: INFILT < HOLD */


		VWC_zoneNORM = VWC;


	}
	else
		VWC_zoneNORM = 0;

	/* -----------------------------*/
	/* 2. capillZone */

	VWC = sprop->VWC_zoneCAPIL;
	dz0 = sprop->dz_zoneCAPIL * m_to_cm;
	VWCfc = sprop->VWCsat[GWlayer];


	/* [cm = m3/m3 * cm */
	HOLD = (sprop->VWCsat[GWlayer] - VWC) * dz0;


	/* 2.1.  IF: INFILT > HOLD */
	if (INFILT > 0.0 && INFILT > HOLD)
	{
		/* drainage from soil profile [cm = m3/m3 * cm ] */
		DRAIN = sprop->drainCoeff[GWlayer] * (sprop->VWCsat[GWlayer] - VWCfc) * dz0;

		/* drainage rate throug soil layer (cm/day) */
		GWR = INFILT - HOLD + DRAIN;


		/* drainage is limited: cm/h * h/day */
		if ((GWR - conduct_sat) > 0.0)
		{
			GWR = conduct_sat;
			DRAIN = GWR + HOLD - INFILT;
		}


		/* state update temporal varialbe */
		VWC = VWC + (INFILT - GWR) / dz0;

		/* above saturation - */
		if (VWC >= sprop->VWCsat[GWlayer])
		{
			EXCESS = (VWC - sprop->VWCsat[GWlayer]) * dz0;
			VWC = sprop->VWCsat[GWlayer];
			GWR += EXCESS;
		}


	} /* END IF: INFILT > HOLD */
	else
	/* 2.2 BEGIN ELSE: INFILT < HOLD */
	{

		VWC = VWC + INFILT / dz0;


		/* BEGIN IF-ELSE: VWC > FC */
		if (VWC >= VWCfc)
		{

			DRAIN = (VWC - VWCfc) * sprop->drainCoeff[GWlayer] * dz0;

			/* drainage rate throug soil layer (cm/day) */
			if (DRAIN > CRIT_PRECwater)
				GWR = DRAIN;
			else
				GWR = 0;



			VWC = VWC - GWR / dz0;
		}
		else
			GWR = 0.0;

	} /* END ELSE: INFILT < HOLD */

	VWC_zoneCAPIL = VWC;


	wf->soilwPercol_NORMvsCAPIL = DRN / mm_to_cm;
	wf->GWrecharge              = GWR / mm_to_cm;

	sprop->VWC_zoneNORM = VWC_zoneNORM;
	sprop->VWC_zoneCAPIL = VWC_zoneCAPIL;


	/* control*/
	soilw_zoneNORM = VWC_zoneNORM * (sprop->dz_zoneNORM * water_density);
	soilw_zoneCAPIL = VWC_zoneCAPIL * (sprop->dz_zoneCAPIL * water_density);


	if (GWlayer == 0)
		INFILT = wf->infiltPOT;
	else
		INFILT = wf->soilwPercol[GWlayer - 1];

	if (sprop->dz_zoneNORM)
	{
		changeNORM1 = INFILT - DRN / mm_to_cm;
		changeCAPIL1 = DRN / mm_to_cm  - GWR / mm_to_cm;
	}
	else
	{
		changeNORM1 = 0;
		changeCAPIL1 = INFILT - GWR / mm_to_cm;
	}


	changeNORM2 = soilw_zoneNORM - sprop->soilw_zoneNORM;
	changeCAPIL2 = soilw_zoneCAPIL - sprop->soilw_zoneCAPIL;


	if (fabs(changeNORM1 - changeNORM2) > CRIT_PRECwater || fabs(changeCAPIL1 - changeCAPIL2) > CRIT_PRECwater)
	{
		printf("\n");
		printf("ERROR in groundwater_tipping in multilayer_hydrolprocess.c\n");
		errorCode = 1;
	}

	sprop->soilw_zoneNORM = soilw_zoneNORM;
	sprop->soilw_zoneCAPIL = soilw_zoneCAPIL;


	/* state update: with new VWC calcualte soilw */
	ws->soilw[GWlayer] = sprop->soilw_zoneNORM + sprop->soilw_zoneCAPIL + sprop->soilw_zoneSAT;
	epv->VWC[GWlayer] = ws->soilw[GWlayer] / (sitec->soillayer_thickness[GWlayer] * water_density);






	return (errorCode);

}
