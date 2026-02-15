/* 
multilayer_transpiration.c
Calculation of part-transpiration (regarding to the different layers of the soil) calculation based on the layer's soil water content and TRPfromGW

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

int multilayer_transpiration(control_struct* ctrl, const siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv,
	                         wstate_struct* ws, wflux_struct* wf)
{
	/* given a list of site constants and the soil water mass (kg/m2),
	this function returns the soil water potential (MPa)
	inputs:


	For further discussion see:
	Cosby, B.J., G.M. Hornberger, R.B. Clapp, and T.R. Ginn, 1984.     

	Balsamo et al 2009 - 
	A Revised Hydrology for the ECMWF Model - Verification from Field Site to Water Storage IFS - JHydromet.pdf

	Chen and Dudhia 2001 - 
	Coupling an Advanced Land Surface-Hydrology Model with the PMM5 Modeling System Part I - MonWRev.pdf*/
	

	#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
	#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

	/* internal variables */
	int layer, n_ROOTlayers, GWlayer, CFlayer;
	double TRPsoilw_SUM, soilw_wp, ratio, upperBoundary;
	double TRPinGWlayer, depthNORM, depthCF;
	double TRPdemandNORM, TRPdemandCF, TRPdemandSAT, soilwAVAILNORM, soilwAVAILCF;
	double ratioNORM, diffNORM, diffCAPIL, ratioCAPIL, soilwAVAIL_NORMcf, soilwAVAIL_CAPILcf;
	
	int errorCode=0;

	n_ROOTlayers = (int)epv->n_rootlayers;
	GWlayer      = (int)sprop->GWlayer;
	CFlayer      = (int)sprop->CFlayer;

	TRPsoilw_SUM=soilw_wp= TRPinGWlayer= depthNORM= depthCF=0;

	/* determination of GW variables */
	if (GWlayer != DATA_GAP)
	{ 
		if (GWlayer > 0)
			depthNORM = sitec->soillayer_depth[GWlayer - 1] + sprop->dz_NORMgw;
		else
			depthNORM = sprop->dz_NORMgw;

		depthCF = depthNORM + sprop->dz_CAPILgw;
	}

	/* *****************************************************************************************************************/
	/* 0. calculation of actual transpiration based on demand */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		/* actual soil water content at theoretical lower limit of water content: hygroscopic water point */
		soilw_wp = sprop->VWCwp[layer] * sitec->soillayer_thickness[layer] * water_density;

		/* TRP_lack: control parameter to avoid negative soil water content (due to overestimated transpiration + dry soil) */
		ws->soilwAVAIL[layer] = MAX(0, (ws->soilw[layer] - soilw_wp));
	}
	
	
	/* 1. PART-transpiration: first approximation tanspiration from every soil layer equally */

	for (layer = epv->germ_layer; layer < n_ROOTlayers; layer++)
	{		
	
		
		/* transpiration based on rootlenght proportion */
		wf->TRPsoilw_demand[layer] = wf->potTRPsoilw * epv->rootlengthProp[layer]; 

		/* upper boundary of the given layer */
		if (layer == 0)
			upperBoundary = 0;
		else
			upperBoundary = sitec->soillayer_depth[layer - 1];

		/* layers without groundwater */
		if (sprop->GWlayer == DATA_GAP || layer < sprop->GWlayer)
		{
			/* soilwAVAIL in last rooting layer: only proportion of rooting depth*/
			if (layer < epv->n_rootlayers - 1 || epv->n_rootlayers == 1)
				ratio = 1;
			else
				ratio = (epv->rootDepth - sitec->soillayer_depth[n_ROOTlayers - 2]) / sitec->soillayer_thickness[n_ROOTlayers];

			/* control */
			if (ratio < 0)
			{
				if (fabs(ratio) > CRIT_PREC)
				{
					printf("\n");
					printf("ERROR in transpiration calculation in multilayer_transpiration.c:\n");
					errorCode = 1;
				}
				else
					ratio = 0;
			}
			if (ratio > 1)
			{
				if (fabs(1-ratio) > CRIT_PREC)
				{
					printf("\n");
					printf("ERROR in transpiration calculation in multilayer_transpiration.c:\n");
					errorCode = 1;
				}
				else
					ratio = 1;
			}

			
			/* if transpiration demand is greater than theoretical lower limit of water content: wilting point -> limited transpiration flux)  */
			if (wf->TRPsoilw_demand[layer] > ws->soilwAVAIL[layer]*ratio)
			{
				/* theoretical limit */
				if (ws->soilwAVAIL[layer] > CRIT_PREC)
					wf->TRPsoilw[layer] = ws->soilwAVAIL[layer] * ratio;
				else
					wf->TRPsoilw[layer] = 0;


				/* limitTRP_flag: writing in log file (only at first time) */
				if (wf->TRPsoilw_demand[layer] - ws->soilwAVAIL[layer] > CRIT_PREC && !ctrl->limitTRP_flag) ctrl->limitTRP_flag = 1;
			}
			else
				wf->TRPsoilw[layer] = wf->TRPsoilw_demand[layer];

			ws->soilw[layer] -= wf->TRPsoilw[layer];
			epv->VWC[layer] = ws->soilw[layer] / sitec->soillayer_thickness[layer] / water_density;
		}
		/* layers with goundwater */
		else
		{
			/* GWlayer: change of normal zone, capillary zone */
			if (layer == GWlayer)
			{
				TRPinGWlayer = wf->TRPsoilw_demand[layer];

				/* rootDepth in normZone */
				if (epv->rootDepth < depthNORM)
				{
					TRPdemandNORM = wf->TRPsoilw_demand[layer];
					TRPdemandCF   = 0;
					TRPdemandSAT  = 0;
				}
				else
				{
					/* rootDepth in capillZone */
					if (epv->rootDepth < depthCF)
					{ 
						TRPdemandNORM = wf->TRPsoilw_demand[layer] * sprop->dz_NORMgw / (epv->rootDepth - upperBoundary);
						TRPdemandCF   = wf->TRPsoilw_demand[layer] - TRPdemandNORM;
						TRPdemandSAT  = 0;
						/* control */
						if (sprop->dz_NORMgw / (epv->rootDepth - upperBoundary) < 0 || sprop->dz_NORMgw / (epv->rootDepth - upperBoundary) > 1)
						{
							printf("\n");
							printf("ERROR in transpiration calculation in multilayer_transpiration.c:\n");
							errorCode = 1;
						}
					}
					else
					{
						/* rootDepth in satZone */
						if (epv->rootDepth < sitec->soillayer_depth[GWlayer])
						{
							TRPdemandNORM = wf->TRPsoilw_demand[layer] * sprop->dz_NORMgw / (epv->rootDepth - upperBoundary);
							TRPdemandCF   = wf->TRPsoilw_demand[layer] * sprop->dz_CAPILgw / (epv->rootDepth - upperBoundary);
							TRPdemandSAT  = wf->TRPsoilw_demand[layer] - TRPdemandNORM - TRPdemandCF;
						}
						/* rootDepth below GWlayer */
						else
						{
							TRPdemandNORM = wf->TRPsoilw_demand[layer] * sprop->dz_NORMgw / sitec->soillayer_thickness[layer];
							TRPdemandCF   = wf->TRPsoilw_demand[layer] * sprop->dz_CAPILgw / sitec->soillayer_thickness[layer];
							TRPdemandSAT  = wf->TRPsoilw_demand[layer] * sprop->dz_SATgw  / sitec->soillayer_thickness[layer];
							/* control */
							if (fabs(TRPdemandNORM+ TRPdemandCF+ TRPdemandSAT - wf->TRPsoilw_demand[layer]) > CRIT_PREC_lenient)
							{
								printf("\n");
								printf("ERROR in transpiration calculation in multilayer_transpiration.c:\n");
								errorCode = 1;
							}
						}
					}
						TRPdemandNORM = wf->TRPsoilw_demand[layer] * sprop->dz_NORMgw / sitec->soillayer_thickness[GWlayer];
				}

				/* calculation of transpiration fluxes and water content in normZone */
				soilwAVAILNORM = sprop->soilw_NORMgw - sprop->VWCwp[GWlayer] * sprop->dz_NORMgw * water_density;
				if (TRPdemandNORM > soilwAVAILNORM)
					wf->TRPsoilw_NORMgw = soilwAVAILNORM;
				else
					wf->TRPsoilw_NORMgw = TRPdemandNORM;

				sprop->soilw_NORMgw -= wf->TRPsoilw_NORMgw;
				if (sprop->dz_NORMgw) sprop->VWC_NORMgw = sprop->soilw_NORMgw / (sprop->dz_NORMgw * water_density);


				/* calculation of transpiration fluxes and water content in capillone */
				soilwAVAILCF   = sprop->soilw_CAPILgw - sprop->VWCwp[GWlayer] * sprop->dz_CAPILgw * water_density;
				if (TRPdemandCF > soilwAVAILCF)
					wf->TRPsoilw_CAPILgw = soilwAVAILCF;
				else
					wf->TRPsoilw_CAPILgw = TRPdemandCF;

				sprop->soilw_CAPILgw -= wf->TRPsoilw_CAPILgw;
				if (sprop->soilw_CAPILgw) sprop->VWC_CAPILgw = sprop->soilw_CAPILgw / (sprop->dz_CAPILgw * water_density);

				/* calculation of transpiration fluxes and water content in satZone */
				wf->TRPfromGW[layer] = TRPdemandSAT;
			
				/* udpate water content of GWlayer */
				wf->TRPsoilw[layer] = wf->TRPsoilw_CAPILgw + wf->TRPsoilw_NORMgw;
				ws->soilw[layer]    = sprop->soilw_NORMgw + sprop->soilw_CAPILgw + sprop->soilw_SATgw;
				epv->VWC[layer]     = ws->soilw[layer] / sitec->soillayer_thickness[layer] / water_density;
			}
			else
				wf->TRPfromGW[layer] = wf->TRPsoilw_demand[layer];
		}



		TRPsoilw_SUM += wf->TRPsoilw[layer] + wf->TRPfromGW[layer];
	}

	wf->TRPsoilw_SUM = TRPsoilw_SUM;

	/* if capillary zone exists in unsaturated zone (not in GWlayer) and capillary zone is in the top soil layer */
	if ((sprop->dz_CAPILcf+sprop->dz_NORMcf) && wf->TRPsoilw[CFlayer] > 0)
	{
		soilwAVAIL_NORMcf = MAX(0, sprop->soilw_NORMcf - sprop->VWCwp[CFlayer] * sprop->dz_NORMcf * water_density);
		soilwAVAIL_CAPILcf = MAX(0, sprop->soilw_CAPILcf - sprop->VWCwp[CFlayer] * sprop->dz_CAPILcf * water_density);

		if (soilwAVAIL_NORMcf+soilwAVAIL_CAPILcf)
		{
			ratioNORM = soilwAVAIL_NORMcf / (soilwAVAIL_NORMcf + soilwAVAIL_CAPILcf);
			ratioCAPIL = soilwAVAIL_CAPILcf / (soilwAVAIL_NORMcf + soilwAVAIL_CAPILcf);
		}
		else
		{
			ratioNORM = sprop->dz_NORMcf / (sprop->dz_CAPILcf + sprop->dz_NORMcf);
			ratioCAPIL = sprop->dz_CAPILcf / (sprop->dz_CAPILcf + sprop->dz_NORMcf);
		}
		if (fabs(1 - ratioNORM - ratioCAPIL) > CRIT_PREC)
		{
			printf("\n");
			printf("ERROR in ratio calculation in multilayer_transpiration.c\n");
			errorCode = 1;
		}


		diffNORM = wf->TRPsoilw[CFlayer] * ratioNORM;

		/* NORM zone: minimum value: wilting point */
		if (soilwAVAIL_NORMcf)
		{ 
			soilw_wp = sprop->VWCwp[CFlayer] * sprop->dz_NORMcf * water_density;
			if (sprop->soilw_NORMcf - diffNORM - soilw_wp < 0)
			{
				if (fabs(sprop->soilw_NORMcf - diffNORM - soilw_wp) > CRIT_PREC)
				{
					printf("ERROR in content_NORMgw calculation in multilayer_transpiration.c\n");
					errorCode = 1;
				}
				else
					diffNORM = sprop->soilw_NORMcf - soilw_wp;
			}
			sprop->soilw_NORMcf -= diffNORM;
			if (sprop->dz_NORMcf) sprop->VWC_NORMcf = sprop->soilw_NORMcf / (sprop->dz_NORMcf * water_density);
		}



		diffCAPIL = wf->TRPsoilw[CFlayer] - diffNORM;

		/* CAPIL zone: minimum value: wilting point  */
		if (soilwAVAIL_CAPILcf)
		{
			soilw_wp = sprop->VWCwp[CFlayer] * sprop->dz_CAPILcf * water_density;
			if (sprop->soilw_CAPILcf - diffCAPIL - soilw_wp < 0)
			{
				if (fabs(sprop->soilw_CAPILcf - diffCAPIL - soilw_wp) > CRIT_PREC)
				{
					printf("ERROR in content_CAPILgw calculation in multilayer_transpiration.c\n");
					errorCode = 1;
				}
				else
					diffCAPIL = sprop->soilw_CAPILcf - soilw_wp;
			}
			sprop->soilw_CAPILcf -= diffCAPIL;
			if (sprop->dz_CAPILcf) sprop->VWC_CAPILcf = sprop->soilw_CAPILcf / (sprop->dz_CAPILcf * water_density);

		}
		wf->TRPsoilwNORMcf = diffNORM;
		wf->TRPsoilwCAPILcf = diffCAPIL;

	}

	/* control */
	if (wf->TRPsoilw_SUM - wf->potTRPsoilw > CRIT_PREC)
	{
		printf("\n");
		printf("ERROR in transpiration calculation in multilayer_transpiration.c:\n");
		errorCode=1;
	}

	/* extreme dry soil - no transpiration occurs */
	if (TRPsoilw_SUM == 0 && wf->TRPsoilw_SUM != 0)
	{
		wf->TRPsoilw_SUM = 0;
		/* noTRP_flag: flag of WARNING writing in log file (only at first time) */
		if (!ctrl->noTRP_flag) ctrl->noTRP_flag = 1;
	}





	return (errorCode);
}



