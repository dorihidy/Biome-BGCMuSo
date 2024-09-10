/* 
multilayer_transpiration.c
Calculation of part-transpiration (regarding to the different layers of the soil) calculation based on the layer's soil water content and GWtransp

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2022, D. Hidy [dori.hidy@gmail.com]
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
	

	/* internal variables */
	int layer, n_ROOTlayers, GWlayer;
	double TRPsoilw_SUM, soilw_wp, ratio, upperBoundary;
	double TRPinGWlayer, depthNORM, depthCF;
	double TRPdemandNORM, TRPdemandCF, TRPdemandSAT, soilw_availNORM, soilw_availCF;
	
	int errorCode=0;

	n_ROOTlayers = (int)epv->n_rootlayers;
	GWlayer      = (int)sprop->GWlayer;

	TRPsoilw_SUM=soilw_wp= TRPinGWlayer=0;

	/* determination of GW variables */
	if (GWlayer > 0)
		depthNORM = sitec->soillayer_depth[GWlayer - 1] + sprop->dz_zoneNORM;
	else
		depthNORM = sprop->dz_zoneNORM;

	depthCF = depthNORM + sprop->dz_zoneCAPIL;



	/* calculation of actual transpiration based on demand */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		/* actual soil water content at theoretical lower limit of water content: hygroscopic water point */
		soilw_wp = sprop->VWCwp[layer] * sitec->soillayer_thickness[layer] * water_density;

		/* TRP_lack: control parameter to avoid negative soil water content (due to overestimated transpiration + dry soil) */
		ws->soilw_avail[layer] = (ws->soilw[layer] - soilw_wp);
	}
	
		/* *****************************************************************************************************************/
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
			/* soilw_avail in last rooting layer: only proportion of rooting depth*/
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
			if (wf->TRPsoilw_demand[layer] > ws->soilw_avail[layer])
			{
				/* theoretical limit */
				if (ws->soilw_avail[layer] > CRIT_PREC)
					wf->TRPsoilw[layer] = ws->soilw_avail[layer] * ratio;
				else
					wf->TRPsoilw[layer] = 0;


				/* limitTRP_flag: writing in log file (only at first time) */
				if (wf->TRPsoilw_demand[layer] - ws->soilw_avail[layer] > CRIT_PREC && !ctrl->limitTRP_flag) ctrl->limitTRP_flag = 1;
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
						TRPdemandNORM = wf->TRPsoilw_demand[layer] * sprop->dz_zoneNORM / (epv->rootDepth - upperBoundary);
						TRPdemandCF   = wf->TRPsoilw_demand[layer] - TRPdemandNORM;
						TRPdemandSAT  = 0;
						/* control */
						if (sprop->dz_zoneNORM / (epv->rootDepth - upperBoundary) < 0 || sprop->dz_zoneNORM / (epv->rootDepth - upperBoundary) > 1)
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
							TRPdemandNORM = wf->TRPsoilw_demand[layer] * sprop->dz_zoneNORM / (epv->rootDepth - upperBoundary);
							TRPdemandCF   = wf->TRPsoilw_demand[layer] * sprop->dz_zoneCAPIL / (epv->rootDepth - upperBoundary);
							TRPdemandSAT  = wf->TRPsoilw_demand[layer] - TRPdemandNORM - TRPdemandCF;
						}
						/* rootDepth below GWlayer */
						else
						{
							TRPdemandNORM = wf->TRPsoilw_demand[layer] * sprop->dz_zoneNORM / sitec->soillayer_thickness[layer];
							TRPdemandCF   = wf->TRPsoilw_demand[layer] * sprop->dz_zoneCAPIL   / sitec->soillayer_thickness[layer];
							TRPdemandSAT  = wf->TRPsoilw_demand[layer] * sprop->dz_zoneSAT  / sitec->soillayer_thickness[layer];
							/* control */
							if (fabs(TRPdemandNORM+ TRPdemandCF+ TRPdemandSAT - wf->TRPsoilw_demand[layer]) > CRIT_PRECwater)
							{
								printf("\n");
								printf("ERROR in transpiration calculation in multilayer_transpiration.c:\n");
								errorCode = 1;
							}
						}
					}
						TRPdemandNORM = wf->TRPsoilw_demand[layer] * sprop->dz_zoneNORM / sitec->soillayer_thickness[GWlayer];
				}

				/* calculation of transpiration fluxes and water content in normZone */
				soilw_availNORM = sprop->soilw_zoneNORM - sprop->VWCwp[GWlayer] * sprop->dz_zoneNORM * water_density;
				if (TRPdemandNORM > soilw_availNORM)
					wf->GWtransp_zoneNORM = soilw_availNORM;
				else
					wf->GWtransp_zoneNORM = TRPdemandNORM;

				sprop->soilw_zoneNORM -= wf->GWtransp_zoneNORM;
				sprop->VWC_zoneNORM    = sprop->soilw_zoneNORM / (sprop->dz_zoneNORM * water_density);

				/* calculation of transpiration fluxes and water content in capillone */
				soilw_availCF   = sprop->soilw_zoneCAPIL   - sprop->VWCwp[GWlayer] * sprop->dz_zoneCAPIL * water_density;
				if (TRPdemandCF > soilw_availCF)
					wf->GWtransp_zoneCAPIL = soilw_availCF;
				else
					wf->GWtransp_zoneCAPIL = TRPdemandCF;

				sprop->soilw_zoneCAPIL -= wf->GWtransp_zoneCAPIL;
				sprop->VWC_zoneNORM = sprop->soilw_zoneCAPIL / (sprop->dz_zoneCAPIL* water_density);

				/* calculation of transpiration fluxes and water content in satZone */
				wf->GWtransp[layer] = TRPdemandSAT;
			
				/* udpate water content of GWlayer */
				wf->TRPsoilw[layer] = wf->GWtransp_zoneCAPIL + wf->GWtransp_zoneNORM;
				ws->soilw[layer]   -= wf->TRPsoilw[layer];
				epv->VWC[layer]     = ws->soilw[layer] / sitec->soillayer_thickness[layer] / water_density;
			}
			else
				wf->GWtransp[layer] = wf->TRPsoilw_demand[layer];
		}



		TRPsoilw_SUM += wf->TRPsoilw[layer];
	}

	wf->TRPsoilw_SUM = TRPsoilw_SUM;

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



