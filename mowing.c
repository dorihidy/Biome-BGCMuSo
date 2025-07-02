/* 
mowing.c
do mowing  - decrease the plant material (leafc, leafn, canopy water)

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
#include "bgc_constants.h"

int mowing(const control_struct* ctrl, const epconst_struct* epc, const mowing_struct* MOW, epvar_struct* epv, 
		   cstate_struct* cs, nstate_struct* ns, wstate_struct* ws, cflux_struct* cf, nflux_struct* nf, wflux_struct* wf)
{

	/* mowing parameters */
	double LAI_limit, outc, outn, inc, inn, MOW_to_transpC, MOW_to_transpN;
	double remained_prop;						/* remained proportion of plabnt material is calculated from transport coefficient */

	/* local parameters */
	double befgrass_LAI;						/* value of LAI before mowing */
	double MOWcoeff;							/* coefficient determining the decrease of plant material caused by mowing */
	int errorCode=0;
	int md, year,ap, condMOW_startyr;

	year = ctrl->simstartyear + ctrl->simyr;
	md = MOW->mgmdMOW-1;

	/* actual phenological phase */
	ap = (int) epv->n_actphen-1; 

	LAI_limit=MOW_to_transpC=MOW_to_transpN=outc=outn=inc=inn=MOWcoeff=0;
	condMOW_startyr=0;

	/**********************************************************************************************/
	/* I. CALCULATING: MOWcoeff */
	
	/* 1. mowing calculation based on LAI_limit: we assume that due to the mowing LAI (leafc) reduces to the value of LAI_limit  */

	befgrass_LAI = cs->leafc * epc->avg_proj_sla[ap];

	/* 2. mowing type: fixday or condMGM */
	if (MOW->condMOW_flag == 0)
	{	
		if (MOW->MOW_num && md >= 0)
		{
			if (year == MOW->MOWyear_array[md] && ctrl->month == MOW->MOWmonth_array[md] && ctrl->day == MOW->MOWday_array[md]) 
			{
				remained_prop = (100 - MOW->transportMOW_array[md])/100.;
				LAI_limit = MOW->LAI_limit_array[md];
			}
			else 
			{
				remained_prop = 0;
				LAI_limit = befgrass_LAI;
			}
		}
	}
	else
	{
		/* type of conditional mowing: 0 -no, 1: - based on LAIratio, 12002 - based on LAIratio but only after 2002 */
		if (MOW->condMOW_flag > 10000 && MOW->condMOW_flag < 20000) condMOW_startyr = MOW->condMOW_flag - 10000;

		/* mowing if LAI greater than a given value (fixLAIbef_condMOW) */
		if (ctrl->simstartyear + ctrl->simyr >= condMOW_startyr && cs->leafc * epc->avg_proj_sla[ap] > MOW->fixLAIbef_condMOW)
		{
			remained_prop = (100 - MOW->transpCOEFF_condMOW)/100.;
			LAI_limit = MOW->fixLAIaft_condMOW;
			if (ctrl->onscreen && (ctrl->spinup == 0 || ctrl->simyr == 0)) printf("conditional MOWING on %i%s%i\n", ctrl->month, "/", ctrl->day);
		}
		else
		{
			remained_prop = 0;
			LAI_limit = befgrass_LAI;
		}
	}

	/* 3. effect of mowing: MOWcoeff */
	if ((MOW->condMOW_flag || (MOW->MOW_num && md >= 0)) && befgrass_LAI > LAI_limit)
	{	
		MOWcoeff = (1- LAI_limit/befgrass_LAI);
	}


	if (MOWcoeff > 0)
	{
		/**********************************************************************************************/
		/* II. CALCULATING FLUXES */

		/*----------------------------------------------------------*/
		/* 1. daily loss due to MOW */ 

		/* 1.1. actual and transfer plant pools*/
		if (epc->leaf_cn)
		{
			cf->leafc_to_MOW              = cs->leafc * MOWcoeff;
			cf->leafc_transfer_to_MOW     = cs->leafc_transfer* MOWcoeff * 0.1;
			cf->leafc_storage_to_MOW      = cs->leafc_storage * MOWcoeff * 0.1;

			nf->leafn_to_MOW              = cf->leafc_to_MOW          / epc->leaf_cn;
			nf->leafn_transfer_to_MOW     = cf->leafc_transfer_to_MOW / epc->leaf_cn;
			nf->leafn_storage_to_MOW      = cf->leafc_storage_to_MOW  / epc->leaf_cn;
		}
	
		if (epc->yield_cn)
		{
			cf->yieldc_to_MOW              = cs->yieldc * MOWcoeff;
			cf->yieldc_transfer_to_MOW     = cs->yieldc_transfer * MOWcoeff * 0.1;
			cf->yieldc_storage_to_MOW      = cs->yieldc_storage  * MOWcoeff * 0.1;

			nf->yieldn_to_MOW              = cf->yieldc_to_MOW          / epc->yield_cn;
			nf->yieldn_transfer_to_MOW     = cf->yieldc_transfer_to_MOW / epc->yield_cn;
			nf->yieldn_storage_to_MOW      = cf->yieldc_storage_to_MOW  / epc->yield_cn;
		}

		if (epc->softstem_cn)
		{
			cf->softstemc_to_MOW              = cs->softstemc * MOWcoeff;
			cf->softstemc_transfer_to_MOW     = cs->softstemc_transfer * MOWcoeff * 0.1;
			cf->softstemc_storage_to_MOW      = cs->softstemc_storage  * MOWcoeff * 0.1;

			nf->softstemn_to_MOW              = cf->softstemc_to_MOW          / epc->softstem_cn;
			nf->softstemn_transfer_to_MOW     = cf->softstemc_transfer_to_MOW / epc->softstem_cn;
			nf->softstemn_storage_to_MOW      = cf->softstemc_storage_to_MOW  / epc->softstem_cn;
		}
	
		cf->gresp_transfer_to_MOW     = cs->gresp_transfer * MOWcoeff * 0.1;
		cf->gresp_storage_to_MOW      = cs->gresp_storage  * MOWcoeff * 0.1;



		/* 1.2.  standing dead biome */
		cf->STDBc_leaf_to_MOW     = cs->STDBc_leaf     * MOWcoeff;
		cf->STDBc_yield_to_MOW    = cs->STDBc_yield    * MOWcoeff; 
		cf->STDBc_softstem_to_MOW = cs->STDBc_softstem * MOWcoeff;

		nf->STDBn_leaf_to_MOW     = ns->STDBn_leaf     * MOWcoeff; 
		nf->STDBn_yield_to_MOW    = ns->STDBn_yield    * MOWcoeff; 
		nf->STDBn_softstem_to_MOW = ns->STDBn_softstem * MOWcoeff;



		/* 1.3 water */
		wf->canopyw_to_MOW = ws->canopyw * MOWcoeff;


		/*----------------------------------------------------------*/
		/* 2. transport: part of the plant material is transported (MOW_to_transpC and MOW_to_transpN; 1-remained_prop),*/

		MOW_to_transpC = (cf->leafc_to_MOW     + cf->leafc_transfer_to_MOW     + cf->leafc_storage_to_MOW     + cf->STDBc_leaf_to_MOW +
			              cf->gresp_storage_to_MOW      + cf->gresp_transfer_to_MOW +
			              cf->yieldc_to_MOW     + cf->yieldc_transfer_to_MOW     + cf->yieldc_storage_to_MOW     + cf->STDBc_yield_to_MOW + 
						  cf->softstemc_to_MOW + cf->softstemc_transfer_to_MOW + cf->softstemc_storage_to_MOW + cf->STDBc_softstem_to_MOW)  * (1-remained_prop);
		MOW_to_transpN = (nf->leafn_to_MOW     + nf->leafn_transfer_to_MOW     + nf->leafn_storage_to_MOW     + nf->STDBn_leaf_to_MOW +
			              nf->yieldn_to_MOW    + nf->yieldn_transfer_to_MOW   + nf->yieldn_storage_to_MOW    + nf->STDBn_yield_to_MOW + 
						  nf->softstemn_to_MOW + nf->softstemn_transfer_to_MOW + nf->softstemn_storage_to_MOW + nf->STDBn_softstem_to_MOW)  * (1-remained_prop);

		cf->MOW_to_transpC = MOW_to_transpC;
		/*----------------------------------------------------------*/
		/* 3. cut-down biomass: the rest remains at the site (MOW_to_litrc_strg, MOW_to_litrn_strg)*/
	
		cf->MOW_to_CTDBc_leaf     = (cf->leafc_to_MOW     + cf->leafc_transfer_to_MOW    + cf->leafc_storage_to_MOW     + cf->STDBc_leaf_to_MOW +
			                         cf->gresp_storage_to_MOW      + cf->gresp_transfer_to_MOW) * remained_prop;

		cf->MOW_to_CTDBc_yield    = (cf->yieldc_to_MOW     + cf->yieldc_transfer_to_MOW     + cf->yieldc_storage_to_MOW     + cf->STDBc_yield_to_MOW)  * remained_prop;

		cf->MOW_to_CTDBc_softstem = (cf->softstemc_to_MOW + cf->softstemc_transfer_to_MOW + cf->softstemc_storage_to_MOW + cf->STDBc_softstem_to_MOW)  * remained_prop;



		nf->MOW_to_CTDBn_leaf     = (nf->leafn_to_MOW     + nf->leafn_transfer_to_MOW     + nf->leafn_storage_to_MOW     + nf->STDBn_leaf_to_MOW) * remained_prop;

		nf->MOW_to_CTDBn_yield    = (nf->yieldn_to_MOW    + nf->yieldn_transfer_to_MOW    + nf->yieldn_storage_to_MOW    + nf->STDBn_yield_to_MOW)  * remained_prop;

		nf->MOW_to_CTDBn_softstem = (nf->softstemn_to_MOW + nf->softstemn_transfer_to_MOW + nf->softstemn_storage_to_MOW + nf->STDBn_softstem_to_MOW) * remained_prop;

		

		/**********************************************************************************************/
		/* III. STATE UPDATE */

		/* 1. OUT */
		/* 1.1.actual and transfer plant pools */	
		cs->leafc				-= cf->leafc_to_MOW;
		cs->yieldc				-= cf->yieldc_to_MOW;
		cs->softstemc			-= cf->softstemc_to_MOW;

		cs->leafc_storage		-= cf->leafc_storage_to_MOW;
		cs->yieldc_storage		-= cf->yieldc_storage_to_MOW;
		cs->softstemc_storage	-= cf->softstemc_storage_to_MOW;

		cs->leafc_transfer		-= cf->leafc_transfer_to_MOW;
		cs->yieldc_transfer		-= cf->yieldc_transfer_to_MOW;
		cs->softstemc_transfer	-= cf->softstemc_transfer_to_MOW;

		cs->gresp_transfer      -= cf->gresp_transfer_to_MOW;
		cs->gresp_storage       -= cf->gresp_storage_to_MOW;


		ns->leafn				-= nf->leafn_to_MOW;
		ns->yieldn				-= nf->yieldn_to_MOW;
		ns->softstemn			-= nf->softstemn_to_MOW;

		ns->leafn_storage		-= nf->leafn_storage_to_MOW;
		ns->yieldn_storage		-= nf->yieldn_storage_to_MOW;
		ns->softstemn_storage	-= nf->softstemn_storage_to_MOW;

		ns->leafn_transfer		-= nf->leafn_transfer_to_MOW;
		ns->yieldn_transfer		-= nf->yieldn_transfer_to_MOW;
		ns->softstemn_transfer	-= nf->softstemn_transfer_to_MOW;


		/* 1.2. dead standing biomass */
		cs->STDBc_leaf     -= cf->STDBc_leaf_to_MOW;
		cs->STDBc_yield    -= cf->STDBc_yield_to_MOW;
		cs->STDBc_softstem -= cf->STDBc_softstem_to_MOW;

		ns->STDBn_leaf     -= nf->STDBn_leaf_to_MOW;
		ns->STDBn_yield    -= nf->STDBn_yield_to_MOW;
		ns->STDBn_softstem -= nf->STDBn_softstem_to_MOW;


		/* 1.3. water */
		ws->canopyw_MOWsnk += wf->canopyw_to_MOW;
		ws->canopyw        -= wf->canopyw_to_MOW;

	
		/* 2. TRANSPORT*/
		cs->MOW_snkC  += MOW_to_transpC;
		ns->MOWsnk_N  += MOW_to_transpN;
	
	
		/* 3. IN: cut-down biome  */
		cs->CTDBc_leaf     += cf->MOW_to_CTDBc_leaf;
		cs->CTDBc_yield    += cf->MOW_to_CTDBc_yield;
		cs->CTDBc_softstem += cf->MOW_to_CTDBc_softstem;

		ns->CTDBn_leaf     += nf->MOW_to_CTDBn_leaf;
		ns->CTDBn_yield    += nf->MOW_to_CTDBn_yield;
		ns->CTDBn_softstem += nf->MOW_to_CTDBn_softstem;
	

		/**********************************************************************************************/
		/* IV. CONTROL */

		outc = cf->leafc_to_MOW + cf->leafc_transfer_to_MOW + cf->leafc_storage_to_MOW +
			   cf->yieldc_to_MOW + cf->yieldc_transfer_to_MOW + cf->yieldc_storage_to_MOW +
			   cf->softstemc_to_MOW + cf->softstemc_transfer_to_MOW + cf->softstemc_storage_to_MOW +
			   cf->gresp_storage_to_MOW + cf->gresp_transfer_to_MOW + 
			   cf->STDBc_leaf_to_MOW + cf->STDBc_yield_to_MOW + cf->STDBc_softstem_to_MOW;

		outn = nf->leafn_to_MOW + nf->leafn_transfer_to_MOW + nf->leafn_storage_to_MOW +
			   nf->yieldn_to_MOW + nf->yieldn_transfer_to_MOW + nf->yieldn_storage_to_MOW +
			   nf->softstemn_to_MOW + nf->softstemn_transfer_to_MOW + nf->softstemn_storage_to_MOW +
			   nf->STDBn_leaf_to_MOW + nf->STDBn_yield_to_MOW + nf->STDBn_softstem_to_MOW;

		inc = cf->MOW_to_CTDBc_leaf + cf->MOW_to_CTDBc_yield  + cf->MOW_to_CTDBc_softstem;
		inn = nf->MOW_to_CTDBn_leaf + nf->MOW_to_CTDBn_yield  + nf->MOW_to_CTDBn_softstem;

		if (fabs(inc + MOW_to_transpC - outc) > CRIT_PREC || fabs(inn + MOW_to_transpN - outn) > CRIT_PREC )
		{
			printf("\n");
			printf("BALANCE ERROR in mowing calculation in mowing.c\n");
			errorCode=1;
		}

}

  return (errorCode);
}	

