/*
groundwater_diffusion.c
Calculate diffusion in groundwater layers

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

int groundwater_diffusion(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{

	int errorCode = 0;
	int GWlayer, CFlayer, layer, flagSAT, realCAPILlayer, ll;

	double soilwDiffus_act, DBAR, EXCESS;
	double VWC0, VWC0_sat, VWC0_eq, VWC0_wp, VWC0_hw, VWC0_EqFC, VWC1, VWC1_sat, VWC1_eq, VWC1_wp, VWC1_hw, VWC1_EqFC;
	double rVWC0, rVWC1, fluxLimit, fl0, fl1;

	double dz0, dz1, dLk;

	double K_sum, K_mean, K_eff,  diffusGW_sum, soilw_satNORM, soilw_satCAPIL, diffNORM, diffCAPIL;
	double diffusGW[N_SOILLAYERS];

	double soilwSAT, diff, extra, soilwNORM_pre, soilw_NORMgw, soilw_CAPILgw, soilw_hwNORM;

	for (layer = 0; layer < N_SOILLAYERS; layer++) diffusGW[layer] = 0;

	soilwDiffus_act = K_sum = K_mean = diffusGW_sum = soilw_satNORM = soilw_satCAPIL = 0.0;
	
	realCAPILlayer = 0;

	GWlayer = (int)sprop->GWlayer;
	CFlayer = (int)sprop->CFlayer;
 


	/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* 1. between SATgw vs. zoneCAP -  */

	if (sprop->dz_CAPILgw)
	{
		dz0 = sprop->dz_CAPILgw + sprop->dz_NORMgw;
		VWC0 = sprop->VWC_NORMgw * sprop->ratioNORMgw + sprop->VWC_CAPILgw * sprop->ratioCAPILgw;
		VWC0_sat = sprop->VWCsat[GWlayer];
		VWC0_eq = sprop->VWCeq[GWlayer] * sprop->ratioNORMgw + sprop->VWCsat[GWlayer] * sprop->ratioCAPILgw;
		VWC0_wp = sprop->VWCwp[GWlayer];
		VWC0_hw = sprop->VWChw[GWlayer];
		VWC0_EqFC = MAX(VWC0_eq, sprop->VWCfc[GWlayer]);

		soilwNORM_pre = sprop->soilw_NORMgw;

		dz1 = sprop->dz_SATgw;
		VWC1 = sprop->VWCsat[GWlayer];
		VWC1_sat = sprop->VWCsat[GWlayer];
		VWC1_eq = sprop->VWCsat[GWlayer];
		VWC1_wp = sprop->VWCwp[GWlayer];
		VWC1_hw = sprop->VWCsat[GWlayer]; // special for GWlayer
		VWC1_EqFC = MAX(VWC1_eq, sprop->VWCfc[GWlayer]);

		/* no limitation of relative VWC because of groundwater source, only VWC0_eq */
		rVWC0 = (VWC0 - VWC0_wp) / (VWC0_EqFC - VWC0_wp);
		rVWC1 = (VWC1 - VWC1_wp) / (VWC1_EqFC - VWC1_wp);

		fluxLimit = (VWC0_eq - VWC0) * dz0 * water_density;


		dLk = DATA_GAP;

		if (!errorCode && calc_diffus(GWlayer, sprop, dz0, VWC0, VWC0_sat, VWC0_EqFC, VWC0_wp, VWC0_hw, dz1, VWC1, VWC1_sat, VWC1_EqFC, VWC1_wp, VWC1_hw, fluxLimit, dLk, &DBAR, &soilwDiffus_act))
		{
			printf("\n");
			printf("ERROR in calc_diffus.c for groundwater_diffusion.c\n");
			errorCode = 1;
		}

		sprop->DBARarray[GWlayer] = DBAR;

		/* udpate of pools and fluxes */
		if (fabs(soilwDiffus_act) > CRIT_PREC)
		{
			/* upward flux - discharge from groudnwater */
			if (soilwDiffus_act < 0)
			{
				/* calculation of limitation from Ksat to direct filling */
				K_eff = MAX(0, sprop->hydrCONDUCTsat[GWlayer] * nSEC_IN_DAY * m_to_mm);

				/* limitation  */
				if (fabs(soilwDiffus_act) > K_eff)
					wf->GWdischarge[GWlayer] = K_eff;
				else
					wf->GWdischarge[GWlayer] = -1 * soilwDiffus_act;

				ws->soilw[GWlayer] -= soilwDiffus_act;

				/* if NORMgw */
				if (sprop->dz_NORMgw)
				{
					soilwSAT = sprop->VWCsat[GWlayer] * sprop->dz_CAPILgw * water_density;
					extra = -1 * soilwDiffus_act;
					diff = soilwSAT - sprop->soilw_CAPILgw;

					/* if there is space for diffusion flux in capillary layer (else: NORMcf is filling */
					if (diff > CRIT_PREC_lenient)
					{
						/* if not enough space in capillary zone: extra goes into NORMgw */
						if (extra - diff > CRIT_PREC_lenient)
						{
							extra -= diff;
							sprop->soilw_CAPILgw = soilwSAT;
							sprop->soilw_NORMgw += extra;

						}
						else
							sprop->soilw_CAPILgw -= soilwDiffus_act;


					}
					else
						sprop->soilw_NORMgw -= soilwDiffus_act;
				}
				else
					sprop->soilw_CAPILgw -= soilwDiffus_act;

				wf->GWdischargeNORMgw = sprop->soilw_NORMgw - soilwNORM_pre;
				wf->GWdischargeCAPILgw = wf->GWdischarge[GWlayer] - wf->GWdischargeNORMgw;

				if (fabs(wf->GWdischargeCAPILgw) < CRIT_PREC) wf->GWdischargeCAPILgw = 0;

			}
			/* downward flux: filling of capillary zone and recharge */
			else
			{
				soilwSAT = sprop->VWCsat[GWlayer] * sprop->dz_CAPILgw * water_density;
				extra = soilwDiffus_act;
				diff = soilwSAT - sprop->soilw_CAPILgw;

				/* if there is space for diffusion flux in capillary layer */
				if (diff > CRIT_PREC_lenient)
				{
					/* if not enough space in capillary zone: extra goes into NORMgw */
					if (extra - diff > CRIT_PREC_lenient)
					{
						extra -= diff;
						wf->soilwDiffus_NORMvsCAPILgw = soilwSAT - sprop->soilw_CAPILgw;
						wf->GWrecharge_NORMgw = extra;
						sprop->soilw_CAPILgw = soilwSAT;
						sprop->soilw_NORMgw -= extra + wf->soilwDiffus_NORMvsCAPILgw;

					}
					else
					{
						wf->soilwDiffus_NORMvsCAPILgw = soilwDiffus_act;
						sprop->soilw_CAPILgw += soilwDiffus_act;
						sprop->soilw_NORMgw -= soilwDiffus_act;
					}
				}
				else
				{ 
					wf->GWrecharge_NORMgw = soilwDiffus_act;
					sprop->soilw_NORMgw -= soilwDiffus_act;
				}

				ws->soilw[GWlayer] -= wf->GWrecharge_NORMgw;

			}


			epv->VWC[GWlayer] = ws->soilw[GWlayer] / sitec->soillayer_thickness[GWlayer] / water_density;

			if (sprop->dz_CAPILgw) sprop->VWC_CAPILgw = sprop->soilw_CAPILgw / sprop->dz_CAPILgw / water_density;
			if (sprop->dz_NORMgw) sprop->VWC_NORMgw = sprop->soilw_NORMgw / sprop->dz_NORMgw / water_density;
		}





		/* flagSAT: saturated or not?: if saturated, the layer above is realCAPIL layer */
		if (sprop->VWCsat[GWlayer] - sprop->VWC_CAPILgw > CRIT_PREC_lenient)
		{
			flagSAT = 0;
			sprop->realCAPILlayer = GWlayer;
		}
		else
		{
			flagSAT = 1;
			sprop->realCAPILlayer = GWlayer - 1;
		}
	}
	/* if there  is no CAPIL zone in GWlayer (GWD is on soillayer boundary) - the layer above is realCAPILlayer is */
	else
	{ 
		flagSAT = 1;
		sprop->realCAPILlayer = GWlayer - 1;
	}




	/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* 2. If GWlayer has three sub-layers: SAT, CAPIL and NORM */

	if (sprop->dz_NORMgw)
	{
		/* control */
		if (GWlayer != CFlayer || !sprop->dz_CAPILgw)
		{
			printf("\n");
			printf("ERROR in sublayer differenciation in groundwater_diffusion.c\n");
			errorCode = 1;
		}

		/* --------------------------------------------------------------------------------------------------------*/
		/* 2. between CAPILgw vs. NORMgw  */
	
		dz0 = sprop->dz_NORMgw;
		VWC0 = sprop->VWC_NORMgw;
		VWC0_sat = sprop->VWCsat[GWlayer];
		VWC0_eq = sprop->VWCeq[GWlayer];
		VWC0_wp = sprop->VWCwp[GWlayer];
		VWC0_hw = sprop->VWChw[GWlayer];
		VWC0_EqFC = MAX(VWC0_eq, sprop->VWCfc[GWlayer]);

		dz1 = sprop->dz_CAPILgw;
		VWC1 = sprop->VWC_CAPILgw;
		VWC1_sat = sprop->VWCsat[GWlayer];
		VWC1_eq = sprop->VWCsat[GWlayer]; 
		VWC1_wp = sprop->VWCwp[GWlayer];
		VWC1_hw = sprop->VWChw[GWlayer];
		VWC1_EqFC = MAX(VWC1_eq, sprop->VWCfc[GWlayer]);

		rVWC0 = (VWC0 - VWC0_wp) / (VWC0_EqFC - VWC0_wp);
		rVWC1 = (VWC1 - VWC1_wp) / (VWC1_EqFC - VWC1_wp);

		fl0 = 1.0 / (dz0 * (VWC0_EqFC - VWC0_wp));
		fl1 = 1.0 / (dz1 * (VWC1_EqFC - VWC1_wp));
		fluxLimit = (rVWC1 - rVWC0) / (fl0 + fl1) * 1000;
		
		dLk = DATA_GAP; // only interpreted in GW and CF layers

		if (!errorCode && calc_diffus(GWlayer, sprop, dz0, VWC0, VWC0_sat, VWC0_EqFC, VWC0_wp, VWC0_hw, dz1, VWC1, VWC1_sat, VWC1_EqFC, VWC1_wp, VWC1_hw, fluxLimit, dLk, &DBAR, &soilwDiffus_act))
		{
			printf("\n");
			printf("ERROR in calc_diffus.c for groundwater_diffusion.c\n");
			errorCode = 1;
		}

		sprop->DBARarray[GWlayer] = DBAR;

		/* udpate of pools and fluxes */
		if (fabs(soilwDiffus_act) > CRIT_PREC)
		{
			/* GWrecharge_NORMgw: if rVWC_NORMgw > rVWC_CAPILgw, but CAPILgw is saturated: GWrecharge_NORMgw*/
			EXCESS = (sprop->soilw_CAPILgw + soilwDiffus_act) - sprop->VWCsat[GWlayer] * sprop->dz_CAPILgw * water_density;
			if (soilwDiffus_act > 0 && EXCESS > 0)
			{
				sprop->soilw_NORMgw -= soilwDiffus_act;
				wf->soilwDiffus_NORMvsCAPILgw += soilwDiffus_act - EXCESS;
				sprop->soilw_CAPILgw += soilwDiffus_act - EXCESS;
				wf->GWrecharge_NORMgw += EXCESS;
			}
			else
			{
				wf->soilwDiffus_NORMvsCAPILgw += soilwDiffus_act;
				sprop->soilw_NORMgw -= soilwDiffus_act;
				sprop->soilw_CAPILgw += soilwDiffus_act;
			}



			sprop->VWC_NORMgw = sprop->soilw_NORMgw / dz0 / water_density;
			sprop->VWC_CAPILgw = sprop->soilw_CAPILgw / dz1 / water_density;
		}

		/* --------------------------------------------------------------------------------------------------------*/
		/* 3. between CAPILgw+NORMgw  vs noGWlayer*/

		if (GWlayer != 0)
		{
			dz0 = sitec->soillayer_thickness[GWlayer - 1];
			VWC0 = epv->VWC[GWlayer - 1];
			VWC0_sat = sprop->VWCsat[GWlayer - 1];
			VWC0_eq = sprop->VWCeq[GWlayer - 1];
			VWC0_wp = sprop->VWCwp[GWlayer - 1];
			VWC0_hw = sprop->VWChw[GWlayer - 1];
			VWC0_EqFC = MAX(VWC0_eq, sprop->VWCfc[GWlayer-1]);

			dz1 = sprop->dz_NORMgw + sprop->dz_CAPILgw;	
			VWC1 = sprop->VWC_NORMgw * sprop->ratioNORMgw + sprop->VWC_CAPILgw * sprop->ratioCAPILgw;
			VWC1_sat = sprop->VWCsat[GWlayer];
			VWC1_eq = sprop->VWCeq[GWlayer] * sprop->ratioNORMgw + sprop->VWCsat[GWlayer] * sprop->ratioCAPILgw;
			VWC1_wp = sprop->VWCwp[GWlayer];
			VWC1_hw = sprop->VWChw[GWlayer];
			VWC1_EqFC = MAX(VWC1_eq, sprop->VWCfc[GWlayer]);

			rVWC0 = (VWC0 - VWC0_wp) / (VWC0_EqFC - VWC0_wp);
			rVWC1 = (VWC1 - VWC1_wp) / (VWC1_EqFC - VWC1_wp);
	
			fl0 = 1.0 / (dz0 * (VWC0_EqFC - VWC0_wp));
			fl1 = 1.0 / (dz1 * (VWC1_EqFC - VWC1_wp));
			fluxLimit = (rVWC1 - rVWC0) / (fl0 + fl1) * 1000;

			soilwNORM_pre = sprop->soilw_NORMgw;

			if (!errorCode && calc_diffus(GWlayer, sprop, dz0, VWC0, VWC0_sat, VWC0_EqFC, VWC0_wp, VWC0_hw, dz1, VWC1, VWC1_sat, VWC1_EqFC, VWC1_wp, VWC1_hw, fluxLimit, dLk, &DBAR, &soilwDiffus_act))
			{
				printf("\n");
				printf("ERROR in calc_diffus.c for groundwater_diffusion.c\n");
				errorCode = 1;
			}


			sprop->DBARarray[GWlayer - 1] = DBAR;

			/* udpate of pools and fluxes */
			if (fabs(soilwDiffus_act) > CRIT_PREC_lenient)
			{

				/* avoiding oversaturation of NORM and CAPIL */
				soilw_NORMgw   = sprop->soilw_NORMgw;
				soilw_CAPILgw  = sprop->soilw_CAPILgw;
	
				soilw_satNORM = sprop->VWCsat[GWlayer] * sprop->dz_NORMgw * water_density;
				soilw_satCAPIL = sprop->VWCsat[GWlayer] * sprop->dz_CAPILgw * water_density;

				/*  avoiding downward percolation into NORMgw from above layer and upward diffusion from CAPILgw into above layer in the same time - source of diffusion is NORMgw primaraly */
				if (soilwDiffus_act < 0 && wf->soilwPercol[GWlayer - 1] > fabs(soilwDiffus_act))
				{
					soilw_hwNORM = sprop->VWChw[GWlayer - 1] * sitec->soillayer_thickness[GWlayer - 1]* water_density;
					if (soilw_NORMgw - soilw_hwNORM + soilwDiffus_act > 0)
						soilw_NORMgw += soilwDiffus_act;
					else
					{	
						soilw_CAPILgw += (soilwDiffus_act + soilw_NORMgw - soilw_hwNORM);
						soilw_NORMgw = soilw_hwNORM;		
					}
				}
				else
				{ 
					soilw_NORMgw += soilwDiffus_act * sprop->ratioNORMgw;
					soilw_CAPILgw += soilwDiffus_act * sprop->ratioCAPILgw;
				}
				diffNORM = soilw_NORMgw - soilw_satNORM;
				diffCAPIL = soilw_CAPILgw - soilw_satCAPIL;
		
				if (diffNORM > 0)
				{
					if (fabs(diffCAPIL - diffNORM) > CRIT_PREC_lenient)
					{							
						soilw_NORMgw = soilw_satNORM;
						soilw_CAPILgw += diffNORM;
					}
					else
						errorCode = 1;
				}

				if (diffCAPIL > 0)
				{
					if (fabs(diffCAPIL - diffNORM) > CRIT_PREC_lenient)
					{
						soilw_NORMgw += diffCAPIL;
						soilw_CAPILgw = soilw_satCAPIL;
					}
					else
						errorCode = 1;
				}


				ws->soilw[GWlayer - 1] -= soilwDiffus_act;
				epv->VWC[GWlayer - 1] = ws->soilw[GWlayer - 1] / dz0 / water_density;
				wf->soilwDiffus[GWlayer - 1] = soilwDiffus_act;

				wf->soilwDiffus_aboveGWlayer_vs_NORMgw = soilw_NORMgw - sprop->soilw_NORMgw;
				wf->soilwDiffus_aboveGWlayer_vs_CAPILgw = soilw_CAPILgw - sprop->soilw_CAPILgw;

				sprop->soilw_NORMgw = soilw_NORMgw;
				sprop->soilw_CAPILgw = soilw_CAPILgw;

				if (sprop->dz_NORMgw)  sprop->VWC_NORMgw = sprop->soilw_NORMgw / sprop->dz_NORMgw / water_density;
				if (sprop->dz_CAPILgw) sprop->VWC_CAPILgw = sprop->soilw_CAPILgw / sprop->dz_CAPILgw / water_density;


				ws->soilw[GWlayer] += soilwDiffus_act;
				epv->VWC[GWlayer] = ws->soilw[GWlayer] / sitec->soillayer_thickness[GWlayer] / water_density;
			}
	
		}
	
		sprop->realCAPILlayer = GWlayer - 1;
	}
	else
	{

		/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
		/* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
		/* 4. If the capillary zone spans several layers: direct filling from GW layer only when capillary layers below are satured */

		realCAPILlayer = (int)sprop->realCAPILlayer;

		while (flagSAT == 1 && realCAPILlayer >= CFlayer)
		{
		
			dz0 = sitec->soillayer_thickness[realCAPILlayer];
			VWC0 = epv->VWC[realCAPILlayer];
			VWC0_sat = sprop->VWCsat[realCAPILlayer];
			if (realCAPILlayer == CFlayer)
				VWC0_eq = sprop->VWCeq[realCAPILlayer] * sprop->ratioNORMcf + sprop->VWCsat[realCAPILlayer] * sprop->ratioCAPILcf;
			else
				VWC0_eq = sprop->VWCsat[realCAPILlayer];
			VWC0_wp = sprop->VWCwp[realCAPILlayer];
			VWC0_hw = sprop->VWChw[realCAPILlayer];
			VWC0_EqFC = MAX(VWC0_eq, sprop->VWCfc[realCAPILlayer]);

			dz1 = sprop->dz_SATgw;
			VWC1 = sprop->VWCsat[GWlayer];
			VWC1_sat = sprop->VWCsat[GWlayer];
			VWC1_eq = sprop->VWCsat[GWlayer];
			VWC1_wp = sprop->VWCwp[GWlayer]; 
			VWC1_hw = sprop->VWCsat[GWlayer];  // special for GWlayer
			VWC1_EqFC = MAX(VWC1_eq, sprop->VWCfc[GWlayer]);

			rVWC0 = (VWC0 - VWC0_wp) / (VWC0_EqFC - VWC0_wp);
			rVWC1 = (VWC1 - VWC1_wp) / (VWC1_EqFC - VWC1_wp);

			fluxLimit = (VWC0_eq - VWC0) * dz0 * water_density;

			soilwNORM_pre = sprop->soilw_NORMcf;

			dLk = (sprop->GWD - sitec->soillayer_midpoint[realCAPILlayer]) * m_to_cm;


			if (!errorCode && calc_diffus(realCAPILlayer, sprop, dz0, VWC0, VWC0_sat, VWC0_EqFC, VWC0_wp, VWC0_hw, dz1, VWC1, VWC1_sat, VWC1_EqFC, VWC1_wp, VWC1_hw, fluxLimit, dLk, &DBAR, &soilwDiffus_act))
			{
				printf("\n");
				printf("ERROR in calc_diffus.c for groundwater_diffusion.c\n");
				errorCode = 1;
			}

			sprop->DBARarray[realCAPILlayer] = DBAR;

			/* update of pools and fluxes */
			if (fabs(soilwDiffus_act) > CRIT_PREC_lenient)
			{
				/* upward flux: discharge from groundwater */
				if (soilwDiffus_act < 0)
				{
					/* calculation of limitation from Ksat to direct filling */
					for (layer = GWlayer; layer > realCAPILlayer; layer--)
					{
						K_sum += sprop->hydrCONDUCTsat[layer];
						diffusGW_sum += diffusGW[layer];
					}
					K_mean = K_sum / (double)(GWlayer - realCAPILlayer);

					K_eff = MAX(0, K_mean * nSEC_IN_DAY * m_to_mm - fabs(diffusGW_sum));

					/* limitation only if there is intermediate capillary layer : if (fabs(soilwDiffus_act) > K_eff && sprop->CFeff[realCAPILlayer + 1] != 0) */
					if (fabs(soilwDiffus_act) > K_eff)
						diffusGW[realCAPILlayer] = -1 * K_eff;
					else
						diffusGW[realCAPILlayer] = soilwDiffus_act;


					/* update of pools */
					ws->soilw[realCAPILlayer] -= diffusGW[realCAPILlayer];
					epv->VWC[realCAPILlayer] = ws->soilw[realCAPILlayer] / sitec->soillayer_thickness[realCAPILlayer] / water_density;

					/* if upper layer is CFlayer - mixed layer is fillin in order to avoid underestimation if upward diffusion: first capillary layer is filled */
					if (realCAPILlayer == CFlayer)
					{
						/* if capillary zone is too thin */
						if (sprop->dz_NORMcf != 0)
						{
							soilwSAT = sprop->VWCsat[realCAPILlayer] * sprop->dz_CAPILcf * water_density;
							extra = -1 * diffusGW[realCAPILlayer];
							diff = soilwSAT - sprop->soilw_CAPILcf;

							/* if there is space for diffusion flux in capillary layer (else: NORMcf is filling */
							if (diff > 0)
							{
								if (extra - diff > CRIT_PREC_lenient)
								{
									extra -= diff;
									sprop->soilw_CAPILcf = soilwSAT;
									if (extra > CRIT_PREC) sprop->soilw_NORMcf += extra;
								}
								else
									sprop->soilw_CAPILcf -= diffusGW[realCAPILlayer];

							}
							else
								sprop->soilw_NORMcf -= diffusGW[realCAPILlayer];
						}
						else
							sprop->soilw_CAPILcf -= diffusGW[realCAPILlayer];

						if (sprop->dz_CAPILcf) sprop->VWC_CAPILcf = sprop->soilw_CAPILcf / sprop->dz_CAPILcf / water_density;
						if (sprop->dz_NORMcf) sprop->VWC_NORMcf = sprop->soilw_NORMcf / sprop->dz_NORMcf / water_density;

						if (fabs(ws->soilw[realCAPILlayer] - sprop->soilw_CAPILcf - sprop->soilw_NORMcf) > CRIT_PREC_lenient)
						{
							printf("\n");
							printf("ERROR in groundwater_diffusion: CAPILcf soilw calculation\n");
							errorCode = 1;
						}


						wf->GWdischargeNORMcf = sprop->soilw_NORMcf - soilwNORM_pre;
						wf->GWdischargeCAPILcf = -1 * diffusGW[realCAPILlayer] - wf->GWdischargeNORMcf;

					}

					wf->GWdischarge[realCAPILlayer] = -1 * diffusGW[realCAPILlayer];

				}
				/* downward flux: filling of capillary zone below and recharge */
				else
				{
				
					EXCESS = (sprop->soilw_CAPILcf + diffusGW[realCAPILlayer]) - sprop->VWCsat[realCAPILlayer] * sprop->dz_CAPILcf * water_density;
					if (EXCESS > 0)
					{
						sprop->soilw_NORMcf -= diffusGW[realCAPILlayer];
						wf->soilwDiffus_NORMvsCAPILcf = diffusGW[realCAPILlayer] - EXCESS;
						sprop->soilw_CAPILcf += wf->soilwDiffus_NORMvsCAPILcf;

						for (ll = realCAPILlayer + 1; ll <= GWlayer; ll++)
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
									wf->soilwPercolDiffus_fromNORM[ll] = EXCESS;
									EXCESS = 0;
								}
								else
								{
									wf->soilwPercolDiffus_fromNORM[ll] = diff;
									EXCESS -= diff;

								}

								/* if CAPIL is in GWlayer - update of capillary zone */
								if (ll == GWlayer && sprop->dz_CAPILgw)
								{
									sprop->soilw_CAPILgw += wf->soilwPercolDiffus_fromNORM[ll];
									sprop->VWC_CAPILgw = sprop->soilw_CAPILgw / sprop->dz_CAPILgw / water_density;
								}
								ws->soilw[ll] += wf->soilwPercolDiffus_fromNORM[ll];
								epv->VWC[ll] = ws->soilw[ll] / sitec->soillayer_thickness[ll] / water_density;
							}


						}
						/* Excess water goes into groundwater (GWrechrage) */
						if (EXCESS) wf->GWrecharge_NORMcf += EXCESS;

					}
				}

	

			}
			/* flagSAT: saturated or not? */
			if (sprop->VWCsat[realCAPILlayer] - epv->VWC[realCAPILlayer] < CRIT_PREC_lenient)
			{
				realCAPILlayer -= 1;
				flagSAT = 1;
			}
			else
				flagSAT = 0;



			if (realCAPILlayer >= CFlayer)
				sprop->realCAPILlayer = realCAPILlayer;
			else
				sprop->realCAPILlayer = CFlayer;



		}
	}


	return (errorCode);

}
