/*
capillary_diffusion.c
Calculate diffusion fluxes in  capillary layer

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2025, D. Hidy [dori.hidy@gmail.com]GWdischar
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

int capillary_diffusion(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{

	int errorCode = 0;
	int GWlayer, CFlayer, layer, realCAPILlayer, ll;

	double soilwDiffus_act, DBAR;
	double VWC0, VWC0_sat, VWC0_eq, VWC0_wp, VWC0_hw, VWC0_EqFC, VWC1, VWC1_sat, VWC1_eq, VWC1_wp, VWC1_hw, VWC1_EqFC;
	double rVWC0, rVWC1, EXCESS, soilwPercolDiffus_fromNORM_act, fluxLimit, fl0, fl1;

	double dz0, dz1, dLk;

	double diffusGW_sum, soilw_satNORM, soilw_satCAPIL, diffNORM, soilw_hwNORM, diffCAPIL;
	double soilwSAT, diff, extra, soilwNORM_pre;
	soilwDiffus_act = diffusGW_sum = 0;
	GWlayer = (int)sprop->GWlayer;
	CFlayer = (int)sprop->CFlayer;
	realCAPILlayer = (int)sprop->realCAPILlayer;


	/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
    /* 1. diffusion fluxes between unsaturated capillay layers  */

	for (layer = realCAPILlayer - 1; layer >= CFlayer+1; layer--)
	{
		dz0 = sitec->soillayer_thickness[layer];
		VWC0 = epv->VWC[layer];
		VWC0_sat = sprop->VWCsat[layer];
		VWC0_eq = sprop->VWCsat[layer];
		VWC0_wp = sprop->VWCwp[layer];
		VWC0_hw = sprop->VWChw[layer];
		VWC0_EqFC = MAX(VWC0_eq, sprop->VWCfc[layer]);

		/* if bottom neigbourg is GWlayer - diffusion between CAPIL zone of GWlayer */
		if (layer + 1 == GWlayer)
		{	
			dz1 = sprop->dz_CAPILgw;
			VWC1 = sprop->VWC_CAPILgw;
			if (!sprop->dz_CAPILgw) errorCode = 1;
		}
		else
		{ 
			dz1 = sitec->soillayer_thickness[layer + 1];
			VWC1 = epv->VWC[layer + 1];
		}
		VWC1_sat = sprop->VWCsat[layer + 1];
		VWC1_eq = sprop->VWCsat[layer + 1];
		VWC1_wp = sprop->VWCwp[layer + 1];
		VWC1_hw = sprop->VWChw[layer + 1];
		VWC1_EqFC = MAX(VWC1_eq, sprop->VWCfc[layer + 1]);

		rVWC0 = (VWC0 - VWC0_wp) / (VWC0_EqFC - VWC0_wp);
		rVWC1 = (VWC1 - VWC1_wp) / (VWC1_EqFC - VWC1_wp);

		fl0 = 1.0 / (dz0 * (VWC0_EqFC - VWC0_wp));
		fl1 = 1.0 / (dz1 * (VWC1_EqFC - VWC1_wp));
		fluxLimit = (rVWC1 - rVWC0) / (fl0 + fl1) * 1000;


		dLk = DATA_GAP; // only interpreted in GW and CF layers



		if (!errorCode && calc_diffus(layer, sprop, dz0, VWC0, VWC0_sat, VWC0_EqFC, VWC0_wp, VWC0_hw, dz1, VWC1, VWC1_sat, VWC1_EqFC, VWC1_wp, VWC1_hw, fluxLimit, dLk, &DBAR, &soilwDiffus_act))
		{
			printf("\n");
			printf("ERROR in calc_diffus.c for capillary_diffusion.c\n");
			errorCode = 1;
		}

		sprop->DBARarray[layer] = DBAR;

		/* udpate of pools and fluxes */
		if (fabs(soilwDiffus_act) > CRIT_PREC_lenient)
		{
			wf->soilwDiffus[layer] = soilwDiffus_act;

			ws->soilw[layer] -= wf->soilwDiffus[layer];
			epv->VWC[layer] = ws->soilw[layer] / dz0 / water_density;

			/* if bottom neigbourg is GWlayer - diffusion between CAPIL zone of GWlayer */
			if (layer + 1 == GWlayer)
			{
				sprop->soilw_CAPILgw += wf->soilwDiffus[layer];
				sprop->VWC_CAPILgw = sprop->soilw_CAPILgw / dz1 / water_density;
				wf->soilwDiffus_aboveGWlayer_vs_CAPILgw = wf->soilwDiffus[layer];
			}
			ws->soilw[layer + 1] += wf->soilwDiffus[layer];
			epv->VWC[layer + 1] = ws->soilw[layer + 1] / sitec->soillayer_thickness[layer + 1] / water_density;
		}

	}



	/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* 2. between last fully capillary layer and CAPILzone of CF layer  */

	if (sprop->dz_CAPILcf && sprop->realCAPILlayer > CFlayer)
	{
		dz0 = sitec->soillayer_thickness[CFlayer];
		VWC0 = epv->VWC[CFlayer];
		VWC0_sat = sprop->VWCsat[CFlayer];
		VWC0_eq = sprop->VWCeq[CFlayer] * sprop->ratioNORMcf + sprop->VWCsat[CFlayer] * sprop->ratioCAPILcf;
		VWC0_wp = sprop->VWCwp[CFlayer];
		VWC0_hw = sprop->VWChw[CFlayer];
		VWC0_EqFC = MAX(VWC0_eq, sprop->VWCfc[CFlayer]);

		if (CFlayer + 1 == GWlayer)
		{
			dz1 = sprop->dz_CAPILgw + sprop->dz_NORMgw;
			VWC1 = sprop->VWC_CAPILgw * sprop->ratioCAPILgw + sprop->VWC_NORMgw * sprop->ratioNORMgw;		
			VWC1_eq = sprop->VWCsat[CFlayer + 1] * sprop->ratioCAPILgw + sprop->VWCeq[CFlayer+1] * sprop->ratioNORMgw;
		}
		else
		{ 
			dz1 = sitec->soillayer_thickness[CFlayer + 1];
			VWC1 = epv->VWC[CFlayer+1];
			VWC1_eq = sprop->VWCsat[CFlayer+1];
		}

		VWC1_sat = sprop->VWCsat[CFlayer + 1];
		VWC1_wp = sprop->VWCwp[CFlayer + 1];
		VWC1_hw = sprop->VWChw[CFlayer + 1];
		VWC1_EqFC = MAX(VWC1_eq, sprop->VWCfc[CFlayer + 1]);

		rVWC0 = (VWC0 - VWC0_wp) / (VWC0_EqFC - VWC0_wp);
		rVWC1 = (VWC1 - VWC1_wp) / (VWC1_EqFC - VWC1_wp);

		fl0 = 1.0 / (dz0 * (VWC0_EqFC - VWC0_wp));
		fl1 = 1.0 / (dz1 * (VWC1_EqFC - VWC1_wp));
		fluxLimit = (rVWC1 - rVWC0) / (fl0 + fl1) * 1000;

		dLk = DATA_GAP;


		if (!errorCode && calc_diffus(CFlayer, sprop, dz0, VWC0, VWC0_sat, VWC0_EqFC, VWC0_wp, VWC0_hw, dz1, VWC1, VWC1_sat, VWC1_EqFC, VWC1_wp, VWC1_hw, fluxLimit, dLk, &DBAR, &soilwDiffus_act))
		{
			printf("\n");
			printf("ERROR in calc_diffus.c for capillary_diffusion.c\n");
			errorCode = 1;
		}

		sprop->DBARarray[CFlayer] = DBAR;

		/* udpate of pools and fluxes */
		if (fabs(soilwDiffus_act) > CRIT_PREC_lenient)
		{
			/* if there is space for diffusion flux in capillary layer (else: NORMcf is filling */
			if (soilwDiffus_act < 0)
			{
				soilwNORM_pre = sprop->soilw_NORMcf;
				soilwSAT = sprop->VWCsat[CFlayer] * sprop->dz_CAPILcf * water_density;
				extra = -1 * soilwDiffus_act;
				diff = soilwSAT - sprop->soilw_CAPILcf;

				if (diff > 0)
				{
					if (diff < extra)
					{
						extra -= diff;
						sprop->soilw_CAPILcf = soilwSAT;
						sprop->soilw_NORMcf += extra;
					}
					else
						sprop->soilw_CAPILcf -= soilwDiffus_act;


				}
				else
					sprop->soilw_NORMcf -= soilwDiffus_act;

 				wf->soilwDiffusNORMcf = soilwNORM_pre - sprop->soilw_NORMcf;
				wf->soilwDiffusCAPILcf = soilwDiffus_act - wf->soilwDiffusNORMcf;
			}
			else
			{
				wf->soilwDiffusNORMcf = sprop->ratioNORMcf * soilwDiffus_act;
				wf->soilwDiffusCAPILcf = sprop->ratioCAPILcf * soilwDiffus_act;
				sprop->soilw_NORMcf -= wf->soilwDiffusNORMcf;
				sprop->soilw_CAPILcf -= wf->soilwDiffusCAPILcf;
			}

			if (sprop->dz_CAPILcf) sprop->VWC_CAPILcf = sprop->soilw_CAPILcf / sprop->dz_CAPILcf / water_density;
			if (sprop->dz_NORMcf) sprop->VWC_NORMcf = sprop->soilw_NORMcf / sprop->dz_NORMcf / water_density;

		
			ws->soilw[CFlayer] -= soilwDiffus_act;
			epv->VWC[CFlayer] = ws->soilw[CFlayer] / sitec->soillayer_thickness[CFlayer] / water_density;

			wf->soilwDiffus[CFlayer] = soilwDiffus_act;
			ws->soilw[CFlayer + 1] += soilwDiffus_act;
			epv->VWC[CFlayer + 1] = ws->soilw[CFlayer + 1] / sitec->soillayer_thickness[CFlayer + 1] / water_density;


			if (fabs(ws->soilw[CFlayer] - sprop->soilw_CAPILcf - sprop->soilw_NORMcf) > CRIT_PREC_lenient)
			{
				printf("\n");
				printf("ERROR in capillary_diffusion: CAPILcf soilw calculation\n");
				errorCode = 1;
			}
			/* if  bottom neighbour is GWlayer */
			if (CFlayer + 1 == GWlayer)
			{
				wf->soilwDiffus_aboveGWlayer_vs_CAPILgw = soilwDiffus_act * sprop->ratioCAPILgw;
				wf->soilwDiffus_aboveGWlayer_vs_NORMgw = soilwDiffus_act * sprop->ratioNORMgw;
				sprop->soilw_CAPILgw += wf->soilwDiffus_aboveGWlayer_vs_CAPILgw;
				sprop->soilw_NORMgw += wf->soilwDiffus_aboveGWlayer_vs_NORMgw;
				if (sprop->dz_CAPILgw) sprop->VWC_CAPILgw = sprop->soilw_CAPILgw / sprop->dz_CAPILgw / water_density;
				if (sprop->dz_NORMgw) sprop->VWC_NORMgw = sprop->soilw_NORMgw / sprop->dz_NORMgw / water_density;
			}
		}

	}

	/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* 3. NORMcf and CAPILcf  */

	if (sprop->dz_NORMcf && sprop->dz_CAPILcf)
	{
		/* control */
		if (GWlayer == CFlayer)
		{
			printf("\n");
			printf("ERROR in sublayer differenciation in groundwater_diffusion.c\n");
			errorCode = 1;
		}

		/* --------------------------------------------------------------------------------------------------------*/
		/* II./2. between CAPILcf vs. NORMcf  */

		dz0 = sprop->dz_NORMcf;
		VWC0 = sprop->VWC_NORMcf;
		VWC0_sat = sprop->VWCsat[CFlayer];
		VWC0_eq = sprop->VWCeq[CFlayer]; 
		VWC0_wp = sprop->VWCwp[CFlayer];
		VWC0_hw = sprop->VWChw[CFlayer];
		VWC0_EqFC = MAX(VWC0_eq, sprop->VWCfc[CFlayer]);

		dz1 = sprop->dz_CAPILcf;
		VWC1 = sprop->VWC_CAPILcf;
		VWC1_sat = sprop->VWCsat[CFlayer];
		VWC1_eq = sprop->VWCsat[CFlayer]; 
		VWC1_wp = sprop->VWCwp[CFlayer];
		VWC1_hw = sprop->VWChw[CFlayer];
		VWC1_EqFC = MAX(VWC1_eq, sprop->VWCfc[CFlayer]);

		rVWC0 = (VWC0 - VWC0_wp) / (VWC0_EqFC - VWC0_wp);
		rVWC1 = (VWC1 - VWC1_wp) / (VWC1_EqFC - VWC1_wp);
	
		fl0 = 1.0 / (dz0 * (VWC0_EqFC - VWC0_wp));
		fl1 = 1.0 / (dz1 * (VWC1_EqFC - VWC1_wp));
		fluxLimit = (rVWC1 - rVWC0) / (fl0 + fl1) * 1000;

		dLk = DATA_GAP;
		

		if (!errorCode && calc_diffus(CFlayer, sprop, dz0, VWC0, VWC0_sat, VWC0_EqFC, VWC0_wp, VWC0_hw, dz1, VWC1, VWC1_sat, VWC1_EqFC, VWC1_wp, VWC1_hw, fluxLimit, dLk, &DBAR, &soilwDiffus_act))
		{
			printf("\n");
			printf("ERROR in calc_diffus.c for capillary_diffusion.c\n");
			errorCode = 1;
		}
		sprop->DBARarray[CFlayer] = DBAR;

		/* update of pools and fluxes */
		if (fabs(soilwDiffus_act) > CRIT_PREC_lenient)
		{
			/* soilwPercolDiffus_fromNORM: if rVWC_NORM > rVWC_CAPIL, but CAPIL is saturated: filling of capillary layers below and GWrecharge_NORMcf*/
			EXCESS = (sprop->soilw_CAPILcf + soilwDiffus_act) - sprop->VWCsat[CFlayer] * sprop->dz_CAPILcf * water_density;
			if (soilwDiffus_act > 0 && EXCESS > 0)
			{
				sprop->soilw_NORMcf -= soilwDiffus_act;
				wf->soilwDiffus_NORMvsCAPILcf += soilwDiffus_act - EXCESS;
				sprop->soilw_CAPILcf += soilwDiffus_act - EXCESS;

				for (ll = CFlayer + 1; ll <= GWlayer; ll++)
				{
					/* if CAPIL is in GWlayer - only capillary zone */
					if (ll == GWlayer && sprop->dz_CAPILgw)
					{
						dz1 = sprop->dz_CAPILgw;
						diff = (sprop->VWCsat[ll] - sprop->VWC_CAPILgw) * dz1 * water_density;
					}
					else
					{
						dz1 = sitec->soillayer_thickness[ll];
						diff = (sprop->VWCsat[ll] - epv->VWC[ll]) * dz1 * water_density;
					}
					if (diff > 0)
					{
						if (diff > EXCESS)
						{
							soilwPercolDiffus_fromNORM_act = EXCESS;		
							EXCESS = 0;
						}
						else
						{
							soilwPercolDiffus_fromNORM_act = diff;
							EXCESS -= diff;

						}
						wf->soilwPercolDiffus_fromNORM[ll] += soilwPercolDiffus_fromNORM_act;
						/* if CAPIL is in GWlayer - update of capillary zone */
						if (ll == GWlayer && sprop->dz_CAPILgw)
						{
							sprop->soilw_CAPILgw += soilwPercolDiffus_fromNORM_act;
							sprop->VWC_CAPILgw = sprop->soilw_CAPILgw / dz1 / water_density;
						}
						ws->soilw[ll] += soilwPercolDiffus_fromNORM_act;
						epv->VWC[ll] = ws->soilw[ll] / sitec->soillayer_thickness[ll] / water_density;
					}
			

				}
				/* Excess water goes into groundwater (GWrechrage) */
				if (EXCESS) wf->GWrecharge_NORMcf += EXCESS;

			}
			else
			{ 
				wf->soilwDiffus_NORMvsCAPILcf = soilwDiffus_act;
				sprop->soilw_NORMcf -= wf->soilwDiffus_NORMvsCAPILcf;
				sprop->soilw_CAPILcf += wf->soilwDiffus_NORMvsCAPILcf;
			}



			sprop->VWC_NORMcf = sprop->soilw_NORMcf / sprop->dz_NORMcf / water_density;
			sprop->VWC_CAPILcf = sprop->soilw_CAPILcf / sprop->dz_CAPILcf / water_density;
		}
	}

	/* --------------------------------------------------------------------------------------------------------*/
	/* --------------------------------------------------------------------------------------------------------*/
	/* IV. between first normal layer and CFlayer */

	if (CFlayer != 0 && sprop->dz_CAPILcf)
	{
		dz0 = sitec->soillayer_thickness[CFlayer - 1];
		VWC0 = epv->VWC[CFlayer - 1];
		VWC0_sat = sprop->VWCsat[CFlayer - 1];
		VWC0_eq = sprop->VWCeq[CFlayer - 1];
		VWC0_wp = sprop->VWCwp[CFlayer - 1];
		VWC0_hw = sprop->VWChw[CFlayer - 1];
		VWC0_EqFC =MAX(VWC0_eq, sprop->VWCfc[CFlayer-1]);
	
		dz1 = sitec->soillayer_thickness[CFlayer];
		VWC1 = epv->VWC[CFlayer];
		VWC1_sat = sprop->VWCsat[CFlayer];
		VWC1_eq = sprop->VWCeq[CFlayer] * sprop->ratioNORMcf + sprop->VWCsat[CFlayer] * sprop->ratioCAPILcf;
		VWC1_wp = sprop->VWCwp[CFlayer];
		VWC1_hw = sprop->VWChw[CFlayer];
		VWC1_EqFC = MAX(VWC1_eq, sprop->VWCfc[CFlayer]);

		soilwNORM_pre = sprop->soilw_NORMcf;

		rVWC0 = (VWC0 - VWC0_wp) / (VWC0_EqFC - VWC0_wp);
		rVWC1 = (VWC1 - VWC1_wp) / (VWC1_EqFC - VWC1_wp);
	
		fl0 = 1.0 / (dz0 * (VWC0_EqFC - VWC0_wp));
		fl1 = 1.0 / (dz1 * (VWC1_EqFC - VWC1_wp));
		fluxLimit = (rVWC1 - rVWC0) / (fl0 + fl1) * 1000;

		dLk = DATA_GAP;


		if (!errorCode && calc_diffus(CFlayer, sprop, dz0, VWC0, VWC0_sat, VWC0_EqFC, VWC0_wp, VWC0_hw, dz1, VWC1, VWC1_sat, VWC1_EqFC, VWC1_wp, VWC1_hw, fluxLimit, dLk, &DBAR, &soilwDiffus_act))
		{
			printf("\n");
			printf("ERROR in calc_diffus.c for capillary_diffusion.c\n");
			errorCode = 1;
		}

		/* element 11: top = noGW */
		sprop->DBARarray[CFlayer - 1] = DBAR;

		if (fabs(soilwDiffus_act) > CRIT_PREC_lenient)
		{
			/*  avoiding downward percolation into NORMgw from above layer and upward diffusion from CAPILgw into above layer in the same time - source of diffusion is NORMgw primaraly  */
			if (soilwDiffus_act < 0 && wf->soilwPercol[CFlayer - 1] > fabs(soilwDiffus_act))
			{
				soilw_hwNORM = sprop->VWChw[CFlayer] * sprop->dz_NORMcf * water_density;
				if (sprop->soilw_NORMcf - soilw_hwNORM + soilwDiffus_act > 0)
					wf->soilwDiffus_aboveCFlayer_vs_NORMcf = soilwDiffus_act;
				else
					wf->soilwDiffus_aboveCFlayer_vs_NORMcf = -1 * (sprop->soilw_NORMcf - soilw_hwNORM);

				wf->soilwDiffus_aboveCFlayer_vs_CAPILcf = soilwDiffus_act - wf->soilwDiffus_aboveCFlayer_vs_NORMcf;

			}
			else
			{
				wf->soilwDiffus_aboveCFlayer_vs_NORMcf = soilwDiffus_act * sprop->ratioNORMcf;
				wf->soilwDiffus_aboveCFlayer_vs_CAPILcf = soilwDiffus_act * sprop->ratioCAPILcf;

				/* avoid overdrainage - limitation: VWC_hw */
				if (wf->soilwDiffus_aboveCFlayer_vs_NORMcf && (sprop->soilw_NORMcf + wf->soilwDiffus_aboveCFlayer_vs_NORMcf - VWC1_hw * sprop->dz_NORMcf * water_density < CRIT_PREC))
				{
					wf->soilwDiffus_aboveCFlayer_vs_NORMcf = sprop->soilw_NORMcf - VWC1_hw * sprop->dz_NORMcf * water_density;
					wf->soilwDiffus_aboveCFlayer_vs_CAPILcf = soilwDiffus_act - wf->soilwDiffus_aboveCFlayer_vs_NORMcf;

				}
				else
				{ 
					if (wf->soilwDiffus_aboveCFlayer_vs_CAPILcf && (sprop->soilw_CAPILcf + wf->soilwDiffus_aboveCFlayer_vs_CAPILcf - VWC1_hw * sprop->dz_CAPILcf * water_density < CRIT_PREC))
					{
						wf->soilwDiffus_aboveCFlayer_vs_CAPILcf = sprop->soilw_CAPILcf - VWC1_hw * sprop->dz_CAPILcf * water_density;
						wf->soilwDiffus_aboveCFlayer_vs_NORMcf = soilwDiffus_act - wf->soilwDiffus_aboveCFlayer_vs_CAPILcf;

					}
				}
			}

			/* avoiding oversaturation of NORM and CAPIL */
			if (soilwDiffus_act > 0)
			{
				soilw_satNORM = sprop->VWCsat[CFlayer] * sprop->dz_NORMcf * water_density;
				soilw_satCAPIL = sprop->VWCsat[CFlayer] * sprop->dz_CAPILcf * water_density;
				diffNORM = sprop->soilw_NORMcf + wf->soilwDiffus_aboveCFlayer_vs_NORMcf - soilw_satNORM;
				diffCAPIL = sprop->soilw_CAPILcf+ wf->soilwDiffus_aboveCFlayer_vs_CAPILcf - soilw_satCAPIL;

				if (diffNORM > CRIT_PREC_lenient)
				{
					if (diffCAPIL + diffNORM < CRIT_PREC_lenient)
					{
						wf->soilwDiffus_aboveCFlayer_vs_NORMcf = soilw_satNORM - sprop->soilw_NORMcf;
						wf->soilwDiffus_aboveCFlayer_vs_CAPILcf = soilwDiffus_act - wf->soilwDiffus_aboveCFlayer_vs_NORMcf;
					}
					else
					{
						printf("\n");
						printf("ERROR in diffusion calulation between first normal layer and CFlayer in capillary_diffusion.c\n");
						errorCode = 1;
					}
				}

				if (diffCAPIL > 0)
				{
					if (diffCAPIL + diffNORM < CRIT_PREC_lenient)
					{
						wf->soilwDiffus_aboveCFlayer_vs_CAPILcf = soilw_satCAPIL - sprop->soilw_CAPILcf;
						wf->soilwDiffus_aboveCFlayer_vs_NORMcf = soilwDiffus_act - wf->soilwDiffus_aboveCFlayer_vs_CAPILcf;
					}
					else
					{
						printf("\n");
						printf("ERROR in diffusion calulation between first normal layer and CFlayer in capillary_diffusion.c\n");
						errorCode = 1;
					}
				}
			}

			sprop->soilw_NORMcf += wf->soilwDiffus_aboveCFlayer_vs_NORMcf;
			sprop->soilw_CAPILcf += wf->soilwDiffus_aboveCFlayer_vs_CAPILcf;

			if (sprop->dz_NORMcf)  sprop->VWC_NORMcf = sprop->soilw_NORMcf / sprop->dz_NORMcf / water_density;
			if (sprop->dz_CAPILcf) sprop->VWC_CAPILcf = sprop->soilw_CAPILcf / sprop->dz_CAPILcf / water_density;

			ws->soilw[CFlayer - 1] -= soilwDiffus_act;
			epv->VWC[CFlayer - 1] = ws->soilw[CFlayer - 1] / dz0 / water_density;


			wf->soilwDiffus[CFlayer - 1] = soilwDiffus_act;

			ws->soilw[CFlayer] = sprop->soilw_NORMcf + sprop->soilw_CAPILcf;
			epv->VWC[CFlayer] = ws->soilw[CFlayer] / sitec->soillayer_thickness[CFlayer] / water_density;
		}

	
	}
	




	return (errorCode);

}
