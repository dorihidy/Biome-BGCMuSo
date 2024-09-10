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
	int layer, ll, N_noGWlayers;

	double VWC, soilw_sat1, soilw1;
	double INFILT, EXCESS, HOLD, DRAIN, DRMX, conduct_sat;


	double dz0, dz1;

	double DRN[N_SOILLAYERS]; /* drainage rate throug soil layer (cm/day) */

	for (layer = 0; layer < N_SOILLAYERS; layer++) DRN[layer] = 0;


	/* tipping is used only for layers without GW */
	if (sprop->GWD == DATA_GAP)
		N_noGWlayers = N_SOILLAYERS;
	else
		N_noGWlayers = (int) sprop->GWlayer;

	/* --------------------------------------------------------------------------------------------------------------------*/
	/* 1.PERCOLATION */

	INFILT = wf->infiltPOT * mm_to_cm;

	/* -----------------------------*/
	/* 1.1. rainy days */
	if (INFILT > 0)
	{

		/* 1.1.1. BEGIN LOOP: layer */
		for (layer = 0; layer < N_noGWlayers; layer++)
		{

			VWC = epv->VWC[layer];
			dz0 = sitec->soillayer_thickness[layer] * m_to_cm;

			/* saturated hydraulic conductivity in actual layer (cm/day = m/s * 100 * sec/day) */
			conduct_sat = sprop->hydrCONDUCTsat[layer] * m_to_cm * nSEC_IN_DAY;


			/* [cm = m3/m3 * cm */
			HOLD = (sprop->VWCsat[layer] - VWC) * dz0;


			/* 1.1.2.  IF: INFILT > HOLD */
			if (INFILT > 0.0 && INFILT > HOLD)
			{
				/* drainage from soil profile [cm = m3/m3 * cm ] */
				DRAIN = sprop->drainCoeff[layer] * (sprop->VWCsat[layer] - sprop->VWCfc[layer]) * dz0;

				/* drainage rate throug soil layer (cm/day) */
				DRN[layer] = INFILT - HOLD + DRAIN;


				/* drainage is limited: cm/h * h/day */
				if ((DRN[layer] - conduct_sat) > 0.0)
				{
					DRN[layer] = conduct_sat;
					DRAIN = DRN[layer] + HOLD - INFILT;
				}


				/* state update temporal varialbe */
				VWC = VWC + (INFILT - DRN[layer]) / dz0;

				/* above saturation - */
				if (VWC >= sprop->VWCsat[layer])
				{

					EXCESS = (VWC - sprop->VWCsat[layer]) * dz0;
					VWC = sprop->VWCsat[layer];

					/* if there is excess water, redistribute it in layers above */
					if (EXCESS > 0)
					{
						for (ll = layer - 1; ll >= 0; ll--)
						{
							dz1 = sitec->soillayer_thickness[ll] * mm_to_cm;
							soilw_sat1 = sprop->VWCsat[ll] * dz1 * water_density;
							soilw1 = epv->VWC[ll] * dz1 * water_density;
							HOLD = MIN(soilw_sat1 - soilw1, EXCESS);
							ws->soilw[ll] += HOLD / mm_to_cm;
							epv->VWC[ll] = ws->soilw[ll] / sitec->soillayer_thickness[ll] / water_density;

							DRN[ll] = MAX(DRN[ll] - EXCESS, 0.0);
							wf->soilwPercol[ll] = (DRN[ll] / m_to_cm) * water_density;
							
							EXCESS = EXCESS - HOLD;

						}
						/* if too much pondwater -> runoff */
						wf->soilw_to_pondw += EXCESS / mm_to_cm;;
					}
				}

				INFILT = DRN[layer];

			} /* END IF: INFILT > HOLD */
			else
			{ /* 1.1.3. BEGIN ELSE: INFILT < HOLD */

				VWC = VWC + INFILT / dz0;


				/* BEGIN IF-ELSE: VWC > FC */
				if (VWC >= sprop->VWCfc[layer])
				{

					DRAIN = (VWC - sprop->VWCfc[layer]) * sprop->drainCoeff[layer] * dz0;

					/* drainage rate throug soil layer (cm/day) */
					if (DRAIN > CRIT_PRECwater)
						DRN[layer] = DRAIN;
					else
						DRN[layer] = 0;


					/* drainage is limited */
					if ((DRN[layer] - conduct_sat) > 0.0)
					{
						DRN[layer] = conduct_sat;
						DRAIN = DRN[layer];
					}

					VWC = VWC - DRN[layer] / dz0;
					INFILT = DRN[layer];
				}
				else
				{
					INFILT = 0.0;
					DRN[layer] = 0.0;

				} /* END IF-ELSE: VWC > FC */

			} /* END ELSE: INFILT < HOLD */

			/* water flux: cm/day to kg/(m2*day) */
			wf->soilwPercol[layer] = (DRN[layer] / m_to_cm) * water_density;
	

			/* state update: with new VWC calcualte soilw */
			epv->VWC[layer] = VWC;
			ws->soilw[layer] = epv->VWC[layer] * sitec->soillayer_thickness[layer] * water_density;


		} /* END FOR (layer) */


	}
	/* -----------------------------*/
	else /* 1.2. rainless days */
	{

		/* BEGIN LOOP: VWCsat flow */
		for (layer = 0; layer < N_noGWlayers; layer++)
		{

			VWC = epv->VWC[layer];
			dz0 = sitec->soillayer_thickness[layer] * m_to_cm;

			/* saturated hydraulic conductivity in actual layer (cm/day = m/s * 100 * sec/day) */
			conduct_sat = sprop->hydrCONDUCTsat[layer] * m_to_cm * nSEC_IN_DAY;


			if (VWC > sprop->VWCfc[layer])
			{
				DRMX = (VWC - sprop->VWCfc[layer]) * sprop->drainCoeff[layer] * dz0;
				DRMX = MAX(0.0, DRMX);

				DRMX = MAX((VWC - sprop->VWCfc[layer]) * sprop->drainCoeff[layer] * dz0, 0);

			}
			else
				DRMX = 0;

			/* BEGIN IF-ELSE: layer == 0 */
			if (layer == 0)
			{
				DRN[layer] = DRMX;
			}
			else
			{

				if (epv->VWC[layer] < sprop->VWCfc[layer])
					HOLD = (sprop->VWCfc[layer] - epv->VWC[layer]) * dz0;
				else
					HOLD = 0.0;

				DRN[layer] = MAX(DRN[layer - 1] + DRMX - HOLD, 0.0);


			} 	/* BEGIN IF-ELSE: layer == 0 */


			/* limitation of drainage: saturation conductivity */
			if ((DRN[layer] - conduct_sat) > 0.0) DRN[layer] = conduct_sat;
			if (DRN[layer] < CRIT_PRECwater) DRN[layer] = 0;


		} /* END LOOP: VWCsat flow */

		for (layer = N_noGWlayers - 1; layer >= 0; layer--)
		{

			VWC = epv->VWC[layer];
			dz0 = sitec->soillayer_thickness[layer] * m_to_cm;


			if (layer > 0)
			{
				VWC = epv->VWC[layer] + DRN[layer - 1]/dz0 - DRN[layer]/dz0;
				if (VWC - sprop->VWCsat[layer] > CRIT_PREC)
				{
					DRN[layer - 1] = (sprop->VWCsat[layer] - epv->VWC[layer]) * dz0 + DRN[layer];
					VWC = sprop->VWCsat[layer];
				}
			}
			else
			{
				VWC = epv->VWC[layer] - (DRN[layer] / dz0);
			}

			if (DRN[layer] < CRIT_PRECwater) DRN[layer] = 0;

			/* water flux: cm/day to kg/(m2*day) */
			wf->soilwPercol[layer] = (DRN[layer] / m_to_cm) * water_density;

	
			/* state update: with new VWC calcualte soilw */
			epv->VWC[layer] = VWC;
			ws->soilw[layer] = epv->VWC[layer] * sitec->soillayer_thickness[layer] * water_density;


		}

	}



	return (errorCode);

}
