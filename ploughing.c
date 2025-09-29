 /* 
ploughing.c
do ploughing  - decrease the plant material (leafc, leafn, canopy water)

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
#include <time.h>
#include <math.h>
#include "ini.h"     
#include "pointbgc_struct.h"
#include "bgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_func.h"
#include "bgc_constants.h"

int ploughing(const control_struct* ctrl, const epconst_struct* epc, siteconst_struct* sitec, soilprop_struct* sprop, soilInfo_struct* soilInfo, metvar_struct* metv,  epvar_struct* epv,
			  ploughing_struct* PLG, cstate_struct* cs, nstate_struct* ns, wstate_struct* ws, cflux_struct* cf, nflux_struct* nf, wflux_struct* wf)
{

	/* ploughing parameters */
	int PLGlayer, layer, dm;
	double PLGdepth, PLGcoeff, soilw_SUM, TsoilSUM, sand_SUM, silt_SUM;	 
	double litr1c_SUM, litr2c_SUM, litr3c_SUM, litr4c_SUM, litr1n_SUM, litr2n_SUM, litr3n_SUM, litr4n_SUM;
	double content_soilSUM[N_DISSOLVMATER], content_soilPRE[N_DISSOLVMATER][N_SOILLAYERS];
	double ratioNORM, ratioCAPIL, diffNORM, diffCAPIL;
	int md, year, flag;
	
	int errorCode = 0;
	flag = 0;

	year = ctrl->simstartyear + ctrl->simyr;
	md = PLG->mgmdPLG-1;

	errorCode=0;
	PLGdepth=0;
	PLGcoeff=soilw_SUM=TsoilSUM=sand_SUM=silt_SUM=0;	 
	litr1c_SUM=litr2c_SUM=litr3c_SUM=litr4c_SUM=litr1n_SUM=litr2n_SUM=litr3n_SUM=litr4n_SUM=0;

	PLGcoeff = 0;

	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		content_soilSUM[dm] = 0;
		for (layer = 0; layer < N_SOILLAYERS; layer++) content_soilPRE[dm][layer] = soilInfo->content_soil[dm][layer];
	}

	/**********************************************************************************************/
	/* I. CALCULATING PLGcoeff AND PLGdepth */

	if (PLG->PLG_num && md >= 0)
	{
		if (year == PLG->PLGyear_array[md] && ctrl->month == PLG->PLGmonth_array[md] && ctrl->day == PLG->PLGday_array[md])
		{
			/* decrease of plant material caused by ploughing: difference between plant material before and after harvesting */
			PLGcoeff      = 1.0;  
			PLGdepth = PLG->PLGdepths_array[md];
			
			/* ploughing layer from depth */
			layer = 0;
			while (layer < N_SOILLAYERS && flag == 0)
			{
				if (PLGdepth <= sitec->soillayer_depth[layer])
				{

					flag = 1;
					PLGlayer = layer;
				}
				layer += 1;
				
			}

			/* ploughing is not possible in layers which contain GW */
			if (sprop->GWlayer != DATA_GAP)
			{
				if (PLGlayer >= sprop->GWlayer) PLGlayer = (int) sprop->GWlayer - 1;
			}

		}
	}


	/**********************************************************************************************/
	/* II. UNIFORM DISTRIBUTION OF SMINN, VWC litterC, litterN and soilC and soilN AFTER PLOUGHING */


	if (PLGcoeff > 0)
	{

		/* *****************************************************************************************************/
		/* II./1. firsttime_flag=0 (before calculation note initial values, partlyORtotal_flag=1 (TOTAL (BOUND+DISSOLV) is affected  */

		if (!errorCode && calc_DISSOLVandBOUND(0, 1, sprop, soilInfo))
		{
			printf("ERROR in calc_DISSOLVandBOUND.c for ploughing.c\n");
			errorCode = 1;
		}

		/* *****************************************************************************************************/
		/* II./2. summarizing content data of affected layers */
		for (layer = 0; layer <= PLGlayer; layer++)
		{

			TsoilSUM += metv->Tsoil[layer] * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[PLGlayer];

			soilw_SUM += ws->soilw[layer];
			sand_SUM += sprop->sand[layer];
			silt_SUM += sprop->silt[layer];
			litr1c_SUM += cs->litr1c[layer];
			litr2c_SUM += cs->litr2c[layer];
			litr3c_SUM += cs->litr3c[layer];
			litr4c_SUM += cs->litr4c[layer];
			litr1n_SUM += ns->litr1n[layer];
			litr2n_SUM += ns->litr2n[layer];
			litr3n_SUM += ns->litr3n[layer];
			litr4n_SUM += ns->litr4n[layer];

			for (dm = 0; dm < N_DISSOLVMATER; dm++) content_soilSUM[dm] += soilInfo->content_soil[dm][layer];


		}

		/* *****************************************************************************************************/
		/* II./3. uniforming content data of affected layers */

		for (layer = 0; layer<=PLGlayer; layer++)
		{
			metv->Tsoil[layer] = TsoilSUM;

    		ws->soilw[layer]   = soilw_SUM * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[PLGlayer];
			epv->VWC[layer]    = ws->soilw[layer] / (water_density * sitec->soillayer_thickness[layer]);

			sprop->sand[layer] = sand_SUM/PLGdepth;
			sprop->silt[layer] = silt_SUM/PLGdepth;
			sprop->clay[layer] = 100-sprop->sand[layer]-sprop->silt[layer];

			cs->litr1c[layer] = litr1c_SUM * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[PLGlayer];
			cs->litr2c[layer] = litr2c_SUM * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[PLGlayer];
			cs->litr3c[layer] = litr3c_SUM * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[PLGlayer];
			cs->litr4c[layer] = litr4c_SUM * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[PLGlayer];

			ns->litr1n[layer] = litr1n_SUM * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[PLGlayer];
			ns->litr2n[layer] = litr2n_SUM * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[PLGlayer];
			ns->litr3n[layer] = litr3n_SUM * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[PLGlayer];
			ns->litr4n[layer] = litr4n_SUM * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[PLGlayer];


			for (dm = 0; dm < N_DISSOLVMATER; dm++) soilInfo->content_soil[dm][layer] = content_soilSUM[dm] * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[PLGlayer] ;

	
			/* if capillary zone exists in unsaturated zone (not in GWlayer) */
			if (layer == sprop->CFlayer && sprop->dz_CAPILcf + sprop->dz_NORMcf)
			{

				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					/* calculation of NORM and CAPIL ratio */
					if (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm])
					{
						ratioNORM = soilInfo->content_NORMcf[dm] / (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm]);
						ratioCAPIL = soilInfo->content_CAPILcf[dm] / (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm]);
					}
					else
					{
						ratioNORM = sprop->dz_NORMcf / (sprop->dz_NORMcf + sprop->dz_CAPILcf);
						ratioCAPIL = sprop->dz_CAPILcf / (sprop->dz_NORMcf + sprop->dz_CAPILcf);
					}
					if (fabs(1 - ratioNORM - ratioCAPIL) > CRIT_PREC)
					{
						printf("\n");
						printf("ERROR in ratio calculation in multilayer_sminn.c\n");
						errorCode = 1;
					}

					diffNORM = (soilInfo->content_soil[dm][layer] - content_soilPRE[dm][layer]) * ratioNORM;
					diffCAPIL = (soilInfo->content_soil[dm][layer] - content_soilPRE[dm][layer]) * ratioCAPIL;

					/* NORM zone: negative storage value is not possible - covered by GW */
					if (diffNORM + soilInfo->content_NORMcf[dm] < 0) diffNORM = -1 * soilInfo->content_NORMcf[dm];

					soilInfo->content_NORMcf[dm] += diffNORM;

					/* CAPIL zone: negative storage value is not possible - covered by GW */
					diffCAPIL = (soilInfo->content_soil[dm][layer] - content_soilPRE[dm][layer]) - diffNORM;


					if (soilInfo->content_CAPILcf[dm] + diffCAPIL < 0)
					{
						if (fabs(soilInfo->content_CAPILcf[dm] + diffCAPIL) > CRIT_PREC)
						{
							printf("ERROR in content_CAPILgw calculation in multilayer_sminn.c\n");
							errorCode = 1;
						}
						else
							diffCAPIL = soilInfo->content_CAPILcf[dm];
					}
					soilInfo->content_CAPILcf[dm] += diffCAPIL;

					if (fabs(soilInfo->content_soil[dm][layer] - (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm])) > CRIT_PREC_lenient * 10)
					{
						printf("ERROR in content_CAPILcf calculation in multilayer_sminn.c\n");
						errorCode = 1;
					}

				}
			}


	

		}

		/* update of litrCabove and litrCbelow */
		cs->litrCabove_total = 0;
		cs->litrCbelow_total = 0;
		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{		
			cs->litrCabove[layer] = 0;
			cs->litrCbelow[layer] = cs->litr1c[layer] + cs->litr2c[layer] + cs->litr3c[layer] + cs->litr4c[layer];
			cs->litrCbelow_total += cs->litrCbelow[layer];
		}

		/* update Tsoil values */
		metv->Tsoil_surface     = metv->Tsoil[0];


		/**********************************************************************************************/
		/* III. CALCULATING FLUXES */
	
	
		/* III.1. leaf, froot, yield, sofstem, gresp*/

		if (epc->leaf_cn)
		{
			cf->leafc_to_PLG				= cs->leafc * PLGcoeff;
			cf->leafc_storage_to_PLG		= cs->leafc_storage * PLGcoeff;
			cf->leafc_transfer_to_PLG		= cs->leafc_transfer * PLGcoeff;
		
			nf->leafn_to_PLG				= ns->leafn * PLGcoeff;
			nf->leafn_storage_to_PLG		= ns->leafn_storage * PLGcoeff;
			nf->leafn_transfer_to_PLG		= ns->leafn_transfer * PLGcoeff;
		}

		if (epc->froot_cn)
		{
			cf->frootc_to_PLG				= cs->frootc * PLGcoeff;
			cf->frootc_storage_to_PLG		= cs->frootc_storage * PLGcoeff;
			cf->frootc_transfer_to_PLG		= cs->frootc_transfer * PLGcoeff;

			nf->frootn_to_PLG				= ns->frootn * PLGcoeff;
			nf->frootn_storage_to_PLG		= ns->frootn_storage * PLGcoeff;
			nf->frootn_transfer_to_PLG		= ns->frootn_transfer * PLGcoeff;
		}

		if (epc->yield_cn)
		{
			cf->yieldc_to_PLG				= cs->yieldc * PLGcoeff;
			cf->yieldc_storage_to_PLG		= cs->yieldc_storage * PLGcoeff;
			cf->yieldc_transfer_to_PLG		= cs->yieldc_transfer * PLGcoeff;

			nf->yieldn_to_PLG			    = ns->yieldn * PLGcoeff;
			nf->yieldn_storage_to_PLG		= ns->yieldn_storage * PLGcoeff;
			nf->yieldn_transfer_to_PLG		= ns->yieldn_transfer * PLGcoeff;
		}

		if (epc->softstem_cn)
		{
			cf->softstemc_to_PLG				= cs->softstemc * PLGcoeff;
			cf->softstemc_storage_to_PLG		= cs->softstemc_storage * PLGcoeff;
			cf->softstemc_transfer_to_PLG		= cs->softstemc_transfer * PLGcoeff;

			nf->softstemn_to_PLG				= ns->softstemn * PLGcoeff;
			nf->softstemn_storage_to_PLG		= ns->softstemn_storage * PLGcoeff;
			nf->softstemn_transfer_to_PLG		= ns->softstemn_transfer * PLGcoeff;
		}

		cf->gresp_storage_to_PLG		= cs->gresp_storage  * PLGcoeff;
		cf->gresp_transfer_to_PLG		= cs->gresp_transfer * PLGcoeff;

		nf->retransn_to_PLG              = ns->retransn * PLGcoeff;

		/* III.2. standing dead biome to cut-down belowground materail: PLGcoeff part of aboveground and whole belowground */

		cf->STDBc_leaf_to_PLG	         = cs->STDBc_leaf     * PLGcoeff;
		cf->STDBc_froot_to_PLG	         = cs->STDBc_froot    * PLGcoeff;
		cf->STDBc_yield_to_PLG	         = cs->STDBc_yield    * PLGcoeff;
		cf->STDBc_softstem_to_PLG        = cs->STDBc_softstem * PLGcoeff;

		nf->STDBn_leaf_to_PLG	         = ns->STDBn_leaf     * PLGcoeff;
		nf->STDBn_froot_to_PLG	         = ns->STDBn_froot    * PLGcoeff;
		nf->STDBn_yield_to_PLG	         = ns->STDBn_yield    * PLGcoeff;
		nf->STDBn_softstem_to_PLG        = ns->STDBn_softstem * PLGcoeff;


		 /* 3. cut-down dead biome: aboveground to belowground  */
		cf->CTDBc_leaf_to_PLG	         = cs->CTDBc_leaf     * PLGcoeff;
		cf->CTDBc_yield_to_PLG	         = cs->CTDBc_yield    * PLGcoeff;
		cf->CTDBc_softstem_to_PLG        = cs->CTDBc_softstem * PLGcoeff;

		nf->CTDBn_leaf_to_PLG	         = ns->CTDBn_leaf     * PLGcoeff;
		nf->CTDBn_yield_to_PLG	         = ns->CTDBn_yield    * PLGcoeff;
		nf->CTDBn_softstem_to_PLG        = ns->CTDBn_softstem * PLGcoeff;

		/* 4. WATER */ 
		wf->canopyw_to_PLG              = ws->canopyw * PLGcoeff;


		/**********************************************************************************************/
		/* IV. STATE UPDATE */

		/* IV. 1. OUT */
		/* IV. 1.1. leaf, froot, yield, sofstem, gresp*/
		cs->leafc				-= cf->leafc_to_PLG;
		cs->leafc_transfer		-= cf->leafc_transfer_to_PLG;
		cs->leafc_storage		-= cf->leafc_storage_to_PLG;
		cs->gresp_transfer		-= cf->gresp_transfer_to_PLG;
		cs->gresp_storage		-= cf->gresp_storage_to_PLG;
		cs->frootc				-= cf->frootc_to_PLG;
		cs->frootc_transfer		-= cf->frootc_transfer_to_PLG;
		cs->frootc_storage		-= cf->frootc_storage_to_PLG;
		cs->yieldc				-= cf->yieldc_to_PLG;
		cs->yieldc_transfer		-= cf->yieldc_transfer_to_PLG;
		cs->yieldc_storage		-= cf->yieldc_storage_to_PLG;
		cs->softstemc			-= cf->softstemc_to_PLG;
		cs->softstemc_transfer  -= cf->softstemc_transfer_to_PLG;
		cs->softstemc_storage   -= cf->softstemc_storage_to_PLG;

		ns->leafn				-= nf->leafn_to_PLG;
		ns->leafn_transfer		-= nf->leafn_transfer_to_PLG;
		ns->leafn_storage		-= nf->leafn_storage_to_PLG;
		ns->frootn				-= nf->frootn_to_PLG;
		ns->frootn_transfer		-= nf->frootn_transfer_to_PLG;
		ns->frootn_storage		-= nf->frootn_storage_to_PLG;
		ns->yieldn				-= nf->yieldn_to_PLG;
		ns->yieldn_transfer		-= nf->yieldn_transfer_to_PLG;
		ns->yieldn_storage		-= nf->yieldn_storage_to_PLG;
		ns->softstemn			-= nf->softstemn_to_PLG;
		ns->softstemn_transfer  -= nf->softstemn_transfer_to_PLG;
		ns->softstemn_storage	-= nf->softstemn_storage_to_PLG;
		ns->retransn			-= nf->retransn_to_PLG;
   

		/* IV.1.2. standing dead biome */
		cs->STDBc_leaf     -= cf->STDBc_leaf_to_PLG;
		cs->STDBc_froot    -= cf->STDBc_froot_to_PLG;
		cs->STDBc_yield    -= cf->STDBc_yield_to_PLG;
		cs->STDBc_softstem -= cf->STDBc_softstem_to_PLG;

		ns->STDBn_leaf     -= nf->STDBn_leaf_to_PLG;
		ns->STDBn_froot    -= nf->STDBn_froot_to_PLG;
		ns->STDBn_yield    -= nf->STDBn_yield_to_PLG;
		ns->STDBn_softstem -= nf->STDBn_softstem_to_PLG;

		/* IV.1.3. cut-down dead biome: aboveground to belowground */
		cs->CTDBc_leaf     -= cf->CTDBc_leaf_to_PLG;
		cs->CTDBc_yield    -= cf->CTDBc_yield_to_PLG;
		cs->CTDBc_softstem -= cf->CTDBc_softstem_to_PLG;

		ns->CTDBn_leaf     -= nf->CTDBn_leaf_to_PLG;
		ns->CTDBn_yield    -= nf->CTDBn_yield_to_PLG;
		ns->CTDBn_softstem -= nf->CTDBn_softstem_to_PLG;
	
	
		/* IV.1.4 water*/
		ws->canopyw        -= wf->canopyw_to_PLG;

		/*--------------------------------------------------------------------*/
		/* IV.2. IN: fluxes to belowground cut-down biomass */
	 
		cs->CTDBc_froot += cf->leafc_to_PLG     + cf->leafc_storage_to_PLG     + cf->leafc_transfer_to_PLG  + 
						   cf->yieldc_to_PLG     + cf->yieldc_storage_to_PLG     + cf->yieldc_transfer_to_PLG + 
						   cf->softstemc_to_PLG + cf->softstemc_storage_to_PLG + cf->softstemc_transfer_to_PLG + 
						   cf->frootc_to_PLG    + cf->frootc_storage_to_PLG    + cf->frootc_transfer_to_PLG +
						   cf->gresp_storage_to_PLG + cf->gresp_transfer_to_PLG +
						   cf->STDBc_leaf_to_PLG     + cf->CTDBc_leaf_to_PLG +  
						   cf->STDBc_yield_to_PLG    + cf->CTDBc_yield_to_PLG + 
						   cf->STDBc_softstem_to_PLG + cf->CTDBc_softstem_to_PLG +
						   cf->STDBc_froot_to_PLG;

		ns->CTDBn_froot += nf->leafn_to_PLG     + nf->leafn_storage_to_PLG     + nf->leafn_transfer_to_PLG  + 
							nf->yieldn_to_PLG     + nf->yieldn_storage_to_PLG     + nf->yieldn_transfer_to_PLG + 
							nf->softstemn_to_PLG + nf->softstemn_storage_to_PLG + nf->softstemn_transfer_to_PLG + 
							nf->frootn_to_PLG    + nf->frootn_storage_to_PLG    + nf->frootn_transfer_to_PLG +
							nf->retransn_to_PLG +
							nf->STDBn_leaf_to_PLG     + nf->CTDBn_leaf_to_PLG +  
							nf->STDBn_yield_to_PLG    + nf->CTDBn_yield_to_PLG + 
							nf->STDBn_softstem_to_PLG + nf->CTDBn_softstem_to_PLG +
							nf->STDBn_froot_to_PLG;

		/* *****************************************************************************************************/
		/* V. udpate pools */


		/* transfer value: NH4, NO3, DOC, DON  to content_array etc.*/
		if (!errorCode && check_soilcontent(-1, 1, sprop, cs, ns, soilInfo))
		{
			printf("ERROR in check_soilcontent.c for multilayer_sminn.c\n");
			errorCode = 1;
		}

		/* firsttime_flag=1 (after calculation note initial values, int partlyORtotal_flag=1 (TOTAL (BOUND+DISSOLV) is affected  */
		if (!errorCode && calc_DISSOLVandBOUND(1, 1, sprop, soilInfo))
		{
			printf("ERROR in calc_DISSOLVandBOUND.c for multilayer_sminn.c\n");
			errorCode = 1;
		}


	}

   return (errorCode);
}
	