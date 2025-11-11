/* 
precision_control.c
Detects very low values in state variable structures, and forces them to
0.0, in order to avoid rounding and overflow errors.

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2025, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

Updated:
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_func.h"
#include "bgc_constants.h"

int precision_control(control_struct* ctrl, soilprop_struct* sprop, wstate_struct* ws, cstate_struct* cs, nstate_struct* ns,  soilInfo_struct* soilInfo)
{
	int errorCode = 0;
	int layer;
	double cpool, npool;


	
	/***********************************************************************************************************************************/
	/* I. CARBON AND NITROGEN STATE VARIABLES */
	/* force very low softstem C to 0.0, to avoid roundoff error in canopy radiation routines. Send excess to Cprec_SNK and Nprec_SNK
	Fine root C and N follow softstem C and N. This control is triggered at a higher value than the others because leafc is used in exp.c in radtrans, and so can cause rounding error at larger values. 
	From MuSo7: precision_reset function is used to reset the very low values and control C:N ration (first variable: litter flag - modification of litrCabove and litrCbelow) */

	/************************/
	/* leaf */

	cpool = cs->leafc;
	npool = ns->leafn;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->leafc = cpool;
	ns->leafn = npool;

	cpool = cs->leafc_transfer;
	npool = ns->leafn_transfer;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->leafc_transfer = cpool;
	ns->leafn_transfer = npool;

	cpool = cs->leafc_storage;
	npool = ns->leafn_storage;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->leafc_storage = cpool;
	ns->leafn_storage = npool;

	/************************/
	/* froot */

	cpool = cs->frootc;
	npool = ns->frootn;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->frootc = cpool;
	ns->frootn = npool;

	cpool = cs->frootc_transfer;
	npool = ns->frootn_transfer;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->frootc_transfer = cpool;
	ns->frootn_transfer = npool;

	cpool = cs->frootc_storage;
	npool = ns->frootn_storage;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->frootc_storage = cpool;
	ns->frootn_storage = npool;

	/************************/
	/* yield */

	cpool = cs->yieldc;
	npool = ns->yieldn;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->yieldc = cpool;
	ns->yieldn = npool;

	cpool = cs->yieldc_transfer;
	npool = ns->yieldn_transfer;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->yieldc_transfer = cpool;
	ns->yieldn_transfer = npool;

	cpool = cs->yieldc_storage;
	npool = ns->yieldn_storage;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->yieldc_storage = cpool;
	ns->yieldn_storage = npool;

	/************************/
	/* softstem */

	cpool = cs->softstemc;
	npool = ns->softstemn;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->softstemc = cpool;
	ns->softstemn = npool;

	cpool = cs->softstemc_transfer;
	npool = ns->softstemn_transfer;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->softstemc_transfer = cpool;
	ns->softstemn_transfer = npool;

	cpool = cs->softstemc_storage;
	npool = ns->softstemn_storage;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->softstemc_storage = cpool;
	ns->softstemn_storage = npool;


	/************************/
	/* livestem */

	cpool = cs->livestemc;
	npool = ns->livestemn;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->livestemc = cpool;
	ns->livestemn = npool;

	cpool = cs->livestemc_transfer;
	npool = ns->livestemn_transfer;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->livestemc_transfer = cpool;
	ns->livestemn_transfer = npool;

	cpool = cs->livestemc_storage;
	npool = ns->livestemn_storage;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->livestemc_storage = cpool;
	ns->livestemn_storage = npool;

	/************************/
	/* deadstem */

	cpool = cs->deadstemc;
	npool = ns->deadstemn;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->deadstemc = cpool;
	ns->deadstemn = npool;

	cpool = cs->deadstemc_transfer;
	npool = ns->deadstemn_transfer;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->deadstemc_transfer = cpool;
	ns->deadstemn_transfer = npool;

	cpool = cs->deadstemc_storage;
	npool = ns->deadstemn_storage;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->deadstemc_storage = cpool;
	ns->deadstemn_storage = npool;

	/************************/
	/* livecroot */

	cpool = cs->livecrootc;
	npool = ns->livecrootn;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->livecrootc = cpool;
	ns->livecrootn = npool;

	cpool = cs->livecrootc_transfer;
	npool = ns->livecrootn_transfer;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->livecrootc_transfer = cpool;
	ns->livecrootn_transfer = npool;

	cpool = cs->livecrootc_storage;
	npool = ns->livecrootn_storage;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->livecrootc_storage = cpool;
	ns->livecrootn_storage = npool;

	/************************/
	/* deadcroot */

	cpool = cs->deadcrootc;
	npool = ns->deadcrootn;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->deadcrootc = cpool;
	ns->deadcrootn = npool;

	cpool = cs->deadcrootc_transfer;
	npool = ns->deadcrootn_transfer;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->deadcrootc_transfer = cpool;
	ns->deadcrootn_transfer = npool;

	cpool = cs->deadcrootc_storage;
	npool = ns->deadcrootn_storage;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->deadcrootc_storage = cpool;
	ns->deadcrootn_storage = npool;

	
	/************************/
	/* STDB   */

	cpool = cs->STDBc_leaf;
	npool = ns->STDBn_leaf;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->STDBc_leaf = cpool;
	ns->STDBn_leaf = npool;

	cpool = cs->STDBc_froot;
	npool = ns->STDBn_froot;
 	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->STDBc_froot = cpool;
	ns->STDBn_froot = npool;

	cpool = cs->STDBc_yield;
	npool = ns->STDBn_yield;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->STDBc_yield = cpool;
	ns->STDBn_yield = npool;

	cpool = cs->STDBc_softstem;
	npool = ns->STDBn_softstem;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->STDBc_softstem = cpool;
	ns->STDBn_softstem = npool;

	


	/************************/
    /* CTDB */

	cpool = cs->CTDBc_softstem;
	npool = ns->CTDBn_softstem;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->CTDBc_softstem = cpool;
	ns->CTDBn_softstem = npool;

	cpool = cs->CTDBc_froot;
	npool = ns->CTDBn_froot;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->CTDBc_froot = cpool;
	ns->CTDBn_froot = npool;

	cpool = cs->CTDBc_yield;
	npool = ns->CTDBn_yield;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->CTDBc_yield = cpool;
	ns->CTDBn_yield = npool;

	cpool = cs->CTDBc_softstem;
	npool = ns->CTDBn_softstem;
	precision_reset(-1, -1, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
	cs->CTDBc_softstem = cpool;
	ns->CTDBn_softstem = npool;


	/*********************************************************************************************************************/
	/* II. Soil and litter variables  - layer by layer  */

	 
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{   	

		cpool = cs->soil1c[layer];
		npool = ns->soil1n[layer];
		precision_reset(2, layer, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
		cs->soil1c[layer] = cpool;
		ns->soil1n[layer] = npool;

		cpool = cs->soil2c[layer];
		npool = ns->soil2n[layer];
		precision_reset(3, layer, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
		cs->soil2c[layer] = cpool;
		ns->soil2n[layer] = npool;

		cpool = cs->soil3c[layer];
		npool = ns->soil3n[layer];
		precision_reset(4, layer, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
		cs->soil3c[layer] = cpool;
		ns->soil3n[layer] = npool; 

		cpool = cs->soil4c[layer];
		npool = ns->soil4n[layer];
		precision_reset(5, layer, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
		cs->soil4c[layer] = cpool;
		ns->soil4n[layer] = npool; 
 
		/**/
		cpool = cs->litr1c[layer];
		npool = ns->litr1n[layer];
		precision_reset(0, layer, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
		cs->litr1c[layer] = cpool;
		ns->litr1n[layer] = npool;


		cpool = cs->litr2c[layer];
		npool = ns->litr2n[layer];
		precision_reset(0, layer, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
		cs->litr2c[layer] = cpool;
		ns->litr2n[layer] = npool;


		cpool = cs->litr3c[layer];
		npool = ns->litr3n[layer];
		precision_reset(0, layer, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
		cs->litr3c[layer] = cpool;
		ns->litr3n[layer] = npool;


		cpool = cs->litr4c[layer];
		npool = ns->litr4n[layer];
		precision_reset(0, layer, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
		cs->litr4c[layer] = cpool;
		ns->litr4n[layer] = npool;

		cpool = cs->cwdc[layer];
		npool = ns->cwdn[layer];
		precision_reset(1, layer, ctrl, sprop, cs, ns, soilInfo, &cpool, &npool, &errorCode);
		cs->cwdc[layer] = cpool;
		ns->cwdn[layer] = npool;
		

	}
		
	/***********************************************************************************************************************************/
	/* II: CARBON OR NITROGEN VARIABLES - no paired variables */

	if (fabs(cs->gresp_transfer) < CRIT_PREC && cs->gresp_transfer != 0)
	{
		cs->Cprec_snk += cs->gresp_transfer;
		cs->gresp_transfer = 0.0;
	}

	if (fabs(cs->gresp_storage) < CRIT_PREC && cs->gresp_storage != 0)
	{
		cs->Cprec_snk += cs->gresp_storage;
		cs->gresp_storage = 0.0;
	}

	if (fabs(cs->cpool) < CRIT_PREC && cs->cpool != 0)
	{
		cs->Cprec_snk += cs->cpool;
		cs->cpool = 0.0;

	}

	if (fabs(cs->NSCnw) < CRIT_PREC && cs->NSCnw != 0)
	{
		cs->Cprec_snk += cs->NSCnw;
		cs->NSCnw = 0.0;

	}

	if (fabs(cs->NSCw) < CRIT_PREC && cs->NSCw != 0)
	{
		cs->Cprec_snk += cs->NSCw;
		cs->NSCw = 0.0;

	}

	if (fabs(cs->SCnw) < CRIT_PREC && cs->SCnw != 0)
	{
		cs->Cprec_snk += cs->SCnw;
		cs->SCnw = 0.0;

	}

	if (fabs(cs->SCw) < CRIT_PREC && cs->SCw != 0)
	{
		cs->Cprec_snk += cs->SCw;
		cs->SCw = 0.0;

	}

	if (fabs(ns->npool) < CRIT_PREC && ns->npool != 0)
	{
		ns->Nprec_snk += ns->npool;
		ns->npool = 0.0;

	}
	if (fabs(ns->retransn) < CRIT_PREC && ns->retransn != 0)
	{
		ns->Nprec_snk += ns->retransn;
		ns->retransn = 0.0;
	}

	/***********************************************************************************************************************************/
    /* III: SOIL MINERAL VALUES  */

	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (fabs(ns->NH4[layer]) < CRIT_PREC_RIG && ns->NH4[layer] != 0)
		{
			ns->Nprec_snk += ns->NH4[layer];
			ns->NH4[layer] = 0.0;
		}

		if (fabs(ns->NO3[layer]) < CRIT_PREC_RIG && ns->NO3[layer] != 0)
		{
			ns->Nprec_snk += ns->NO3[layer];
			ns->NO3[layer] = 0.0;
		}

	}


	/***********************************************************************************************************************************/
	/* IV: WATER STATE VARIABLES */

	/* multilayer soil */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (ws->soilw[layer] < 0 && fabs(ws->soilw[layer]) < CRIT_PREC)
		{
			ws->Wprec_snk += ws->soilw[layer];
			ws->soilw[layer] = 0.0;
		}
	}

	if (ws->snoww < 0 && fabs(ws->snoww) < CRIT_PREC)
	{
		ws->Wprec_snk += ws->snoww;
		ws->snoww = 0.0;
	}


	if (ws->canopyw < 0 && fabs(ws->canopyw) < CRIT_PREC)
	{
		ws->Wprec_snk += ws->canopyw;
		ws->canopyw = 0.0;
	}

	if (ws->pondw < 0 && fabs(ws->pondw) < CRIT_PREC)
	{
		ws->Wprec_snk += ws->pondw;
		ws->pondw = 0.0;
	}

	

	return(errorCode);
}	


int precision_reset(int dm, int layer, control_struct* ctrl, soilprop_struct* sprop, cstate_struct* cs, nstate_struct* ns, soilInfo_struct* soilInfo, double* cpool, double* npool, int* errorCode)
{

	double litrClayer;



	/* examination of ratios */
	if (fabs(*cpool) / fabs(*npool) < 1 && *npool > CRIT_PREC_lenient && dm > 1) ctrl->CNratio_flag = 1;

	/* examination of cpool and npool values */
	if ((fabs(*cpool) < CRIT_PREC_RIG && *cpool != 0) || (fabs(*npool) < CRIT_PREC_RIG && *npool != 0))
	{


		/* extra control fo litter */
		if (dm == 0)
		{
			litrClayer = cs->litrCabove[layer] + cs->litrCbelow[layer];
			if (litrClayer)
			{
				cs->litrCabove[layer] -= *cpool * (cs->litrCabove[layer] / litrClayer);
				cs->litrCbelow[layer] -= *cpool * (cs->litrCbelow[layer] / litrClayer);
			}

		}

		/* extra control fo cwd */
		if (dm == 1)
		{
			cs->cwdCabove[layer] = 0.0;
			cs->cwdCbelow[layer] = 0.0;
		}

		/* in GWlayer no precision control */

		if (layer != (int)sprop->GWlayer)
		{
			cs->Cprec_snk += *cpool;
			ns->Nprec_snk += *npool;
			*cpool = 0.0;
			*npool = 0.0;

			/* updating of soilIfno struct with change of soil pools (dm=2 to 9*/
			if (dm > 1)
			{
				soilInfo->content_soil[dm][layer] = 0;
				soilInfo->contentBOUND_soil[dm][layer] = 0;
				soilInfo->contentDISSOLV_soil[dm][layer] = 0;

				soilInfo->content_soil[dm + 4][layer] = 0;
				soilInfo->contentBOUND_soil[dm + 4][layer] = 0;
				soilInfo->contentDISSOLV_soil[dm + 4][layer] = 0;


				if (layer == (int)sprop->CFlayer && sprop->CFlayer != sprop->GWlayer)
				{
					soilInfo->contentBOUND_CAPILcf[dm] = 0;
					soilInfo->contentDISSOLV_CAPILcf[dm] = 0;
					soilInfo->content_CAPILcf[dm] = 0;
					soilInfo->contentBOUND_NORMcf[dm] = 0;
					soilInfo->contentDISSOLV_NORMcf[dm] = 0;
					soilInfo->content_NORMcf[dm] = 0;
					soilInfo->contentBOUND_CAPILcf[dm + 4] = 0;
					soilInfo->contentDISSOLV_CAPILcf[dm + 4] = 0;
					soilInfo->content_CAPILcf[dm + 4] = 0;
					soilInfo->contentBOUND_NORMcf[dm + 4] = 0;
					soilInfo->contentDISSOLV_NORMcf[dm + 4] = 0;
					soilInfo->content_NORMcf[dm + 4] = 0;
				}
			}
		}

	}


	return(*errorCode);
}