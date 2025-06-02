/* 
precision_control.c
Detects very low values in state variable structures, and forces them to
0.0, in order to avoid rounding and overflow errors.

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2022, D. Hidy [dori.hidy@gmail.com]
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

int precision_control(wstate_struct* ws, cstate_struct* cs, nstate_struct* ns, soilprop_struct* sprop, soilInfo_struct* soilInfo)
{
	int errorCode = 0;
	int layer, CFlayer, GWlayer, dm;
	double litrClayer;
	CFlayer = (int)sprop->CFlayer;
	GWlayer = (int)sprop->GWlayer;

	/* CARBON AND NITROGEN STATE VARIABLES */
	/* force very low leaf C to 0.0, to avoid roundoff
	error in canopy radiation routines. Send excess to litter 1.
	Fine root C and N follow leaf C and N. This control is triggered
	at a higher value than the others because leafc is used in exp.c
	in radtrans, and so can cause rounding error at larger values. */

	/************************/
	if ((fabs(cs->leafc) < 1e-8 && fabs(cs->leafc) != 0) || (fabs(ns->leafn) < 1e-9 && fabs(ns->leafn) != 0))
	{
		cs->litr1c[0] += cs->leafc;
		ns->litr1n[0] += ns->leafn;
		cs->litrCabove[0] += cs->leafc;
		cs->leafc = 0.0;
		ns->leafn = 0.0;


	}

	if ((cs->leafc_storage != 0 && fabs(cs->leafc_storage) < CRIT_PREC_RIG) || (ns->leafn_storage != 0 && fabs(ns->leafn_storage) < CRIT_PREC_RIG))
	{
		cs->litr1c[0] += cs->leafc_storage;
		ns->litr1n[0] += ns->leafn_storage;
		cs->litrCabove[0] += cs->leafc_storage;
		cs->leafc_storage = 0.0;
		ns->leafn_storage = 0.0;
	}

	if ((cs->leafc_transfer != 0 && fabs(cs->leafc_transfer) < CRIT_PREC_RIG) || (ns->leafn_transfer != 0 && fabs(ns->leafn_transfer) < CRIT_PREC_RIG))
	{
		cs->litr1c[0] += cs->leafc_transfer;
		ns->litr1n[0] += ns->leafn_transfer;
		cs->litrCabove[0] += cs->leafc_transfer;
		cs->leafc_transfer = 0.0;
		ns->leafn_transfer = 0.0;
	}

	/************************/
	if ((fabs(cs->frootc) < CRIT_PREC && fabs(cs->frootc) != 0) || (fabs(ns->frootn) < CRIT_PREC && fabs(ns->frootn) != 0))
	{
		cs->litr1c[0] += cs->frootc;
		ns->litr1n[0] += ns->frootn;
		cs->litrCbelow[0] += cs->frootc;
		cs->frootc = 0.0;
		ns->frootn = 0.0;
	}

	if ((cs->frootc_storage != 0 && fabs(cs->frootc_storage) < CRIT_PREC_RIG) || (ns->frootn_storage != 0 && fabs(ns->frootn_storage) < CRIT_PREC_RIG))
	{
		cs->litr1c[0] += cs->frootc_storage;
		ns->litr1n[0] += ns->frootn_storage;
		cs->litrCbelow[0] += cs->frootc_storage;
		cs->frootc_storage = 0.0;
		ns->frootn_storage = 0.0;
	}

	if ((cs->frootc_transfer != 0 && fabs(cs->frootc_transfer) < CRIT_PREC_RIG) || (ns->frootn_transfer != 0 && fabs(ns->frootn_transfer) < CRIT_PREC_RIG))
	{
		cs->litr1c[0] += cs->frootc_transfer;
		ns->litr1n[0] += ns->frootn_transfer;
		cs->litrCbelow[0] += cs->frootc_transfer;
		cs->frootc_transfer = 0.0;
		ns->frootn_transfer = 0.0;
	}

	/************************/
	if ((fabs(cs->yieldc) < CRIT_PREC && fabs(cs->yieldc) != 0) || (fabs(ns->yieldn) < CRIT_PREC && fabs(ns->yieldn) != 0))
	{
		cs->litr1c[0] += cs->yieldc;
		ns->litr1n[0] += ns->yieldn;
		cs->litrCabove[0] += cs->yieldc;
		cs->yieldc = 0.0;
		ns->yieldn = 0.0;
	}

	if ((cs->yieldc_storage != 0 && fabs(cs->yieldc_storage) < CRIT_PREC_RIG) || (ns->yieldn_storage != 0 && fabs(ns->yieldn_storage) < CRIT_PREC_RIG))
	{
		cs->litr1c[0] += cs->yieldc_storage;
		ns->litr1n[0] += ns->yieldn_storage;
		cs->litrCabove[0] += cs->yieldc_storage;
		cs->yieldc_storage = 0.0;
		ns->yieldn_storage = 0.0;
	}

	if ((cs->yieldc_transfer != 0 && fabs(cs->yieldc_transfer) < CRIT_PREC_RIG) || (ns->yieldn_transfer < 0 && fabs(ns->yieldn_transfer) < CRIT_PREC_RIG))
	{
		cs->litr1c[0] += cs->yieldc_transfer;
		ns->litr1n[0] += ns->yieldn_transfer;
		cs->litrCabove[0] += cs->yieldc_transfer;
		cs->yieldc_transfer = 0.0;
		ns->yieldn_transfer = 0.0;
	}

	/************************/
	if ((fabs(cs->softstemc) != 0 && fabs(cs->softstemc) < CRIT_PREC) || (fabs(ns->softstemn) != 0 && fabs(ns->softstemn) < CRIT_PREC))
	{
		cs->litr1c[0] += cs->softstemc;
		ns->litr1n[0] += ns->softstemn;
		cs->litrCabove[0] += cs->softstemc;
		cs->softstemc = 0.0;
		ns->softstemn = 0.0;
	}


	if ((cs->softstemc_storage != 0 && fabs(cs->softstemc_storage) < CRIT_PREC_RIG) || (ns->softstemn_storage != 0 && fabs(ns->softstemn_storage) < CRIT_PREC_RIG))
	{
		cs->litr1c[0] += cs->softstemc_storage;
		ns->litr1n[0] += ns->softstemn_storage;
		cs->litrCabove[0] += cs->softstemc_storage;
		cs->softstemc_storage = 0.0;
		ns->softstemn_storage = 0.0;
	}

	if ((cs->softstemc_transfer != 0 && fabs(cs->softstemc_transfer) < CRIT_PREC_RIG) || (ns->softstemn_transfer != 0 && fabs(ns->softstemn_transfer) < CRIT_PREC_RIG))
	{
		cs->litr1c[0] += cs->softstemc_transfer;
		ns->litr1n[0] += ns->softstemn_transfer;
		cs->litrCabove[0] += cs->softstemc_transfer;
		cs->softstemc_transfer = 0.0;
		ns->softstemn_transfer = 0.0;
	}

	/************************/
	if ((fabs(cs->livestemc) < CRIT_PREC && fabs(cs->livestemc) != 0) || (fabs(ns->livestemn) < CRIT_PREC && fabs(ns->livestemn) != 0))
	{
		cs->litr1c[0] += cs->livestemc;
		ns->litr1n[0] += ns->livestemn;
		cs->litrCabove[0] += cs->livestemc;
		cs->livestemc = 0.0;
		ns->livestemn = 0.0;
	}

	if ((cs->livestemc_storage != 0 && fabs(cs->livestemc_storage) < CRIT_PREC_RIG) || (ns->livestemn_storage != 0 && fabs(ns->livestemn_storage) < CRIT_PREC_RIG))
	{
		cs->litr1c[0] += cs->livestemc_storage;
		ns->litr1n[0] += ns->livestemn_storage;
		cs->litrCabove[0] += cs->livestemc_storage;
		cs->livestemc_storage = 0.0;
		ns->livestemn_storage = 0.0;
	}

	if ((cs->livestemc_transfer != 0 && fabs(cs->livestemc_transfer) < CRIT_PREC_RIG) || (ns->livestemn_transfer != 0 && fabs(ns->livestemn_transfer) < CRIT_PREC_RIG))
	{
		cs->litr1c[0] += cs->livestemc_transfer;
		ns->litr1n[0] += ns->livestemn_transfer;
		cs->litrCabove[0] += cs->livestemc_transfer;
		cs->livestemc_transfer = 0.0;
		ns->livestemn_transfer = 0.0;
	}

	/************************/
	if ((fabs(cs->deadstemc) < CRIT_PREC && fabs(cs->deadstemc) != 0) || (fabs(ns->deadstemn) < CRIT_PREC && fabs(ns->deadstemn) != 0))
	{
		cs->litr1c[0] += cs->deadstemc;
		ns->litr1n[0] += ns->deadstemn;
		cs->litrCabove[0] += cs->deadstemc;
		cs->deadstemc = 0.0;
		ns->deadstemn = 0.0;
	}

	if ((cs->deadstemc_storage != 0 && fabs(cs->deadstemc_storage) < CRIT_PREC_RIG) || (ns->deadstemn_storage != 0 && fabs(ns->deadstemn_storage) < CRIT_PREC_RIG))
	{
		cs->litr1c[0] += cs->deadstemc_storage;
		ns->litr1n[0] += ns->deadstemn_storage;
		cs->litrCabove[0] += cs->deadstemc_storage;
		cs->deadstemc_storage = 0.0;
		ns->deadstemn_storage = 0.0;
	}

	if ((cs->deadstemc_transfer != 0 && fabs(cs->deadstemc_transfer) < CRIT_PREC_RIG) || (ns->deadstemn_transfer != 0 && fabs(ns->deadstemn_transfer) < CRIT_PREC_RIG))
	{
		cs->litr1c[0] += cs->deadstemc_transfer;
		ns->litr1n[0] += ns->deadstemn_transfer;
		cs->litrCabove[0] += cs->deadstemc_transfer;
		cs->deadstemc_transfer = 0.0;
		ns->deadstemn_transfer = 0.0;
	}

	/************************/
	if ((fabs(cs->livecrootc) < CRIT_PREC && fabs(cs->livecrootc) != 0) || (fabs(ns->livecrootn) < CRIT_PREC && fabs(ns->livecrootn) != 0))
	{
		cs->litr1c[0] += cs->livecrootc;
		ns->litr1n[0] += ns->livecrootn;
		cs->litrCbelow[0] += cs->livecrootc;
		cs->livecrootc = 0.0;
		ns->livecrootn = 0.0;
	}

	if ((cs->livecrootc_storage < 0 && fabs(cs->livecrootc_storage) < CRIT_PREC) || (ns->livecrootn_storage < 0 && fabs(ns->livecrootn_storage) < CRIT_PREC))
	{
		cs->litr1c[0] += cs->livecrootc_storage;
		ns->litr1n[0] += ns->livecrootn_storage;
		cs->litrCbelow[0] += cs->livecrootc_storage;
		cs->livecrootc_storage = 0.0;
		ns->livecrootn_storage = 0.0;
	}

	if ((cs->livecrootc_transfer < 0 && fabs(cs->livecrootc_transfer) < CRIT_PREC) || (ns->livecrootn_transfer < 0 && fabs(ns->livecrootn_transfer) < CRIT_PREC))
	{
		cs->litr1c[0] += cs->livecrootc_transfer;
		ns->litr1n[0] += ns->livecrootn_transfer;
		cs->litrCbelow[0] += cs->livecrootc_transfer;
		cs->livecrootc_transfer = 0.0;
		ns->livecrootn_transfer = 0.0;
	}

	/************************/
	if ((fabs(cs->deadcrootc) < CRIT_PREC && fabs(cs->deadcrootc) != 0) || (fabs(ns->deadcrootn) < CRIT_PREC && fabs(ns->deadcrootn) != 0))
	{
		cs->litr1c[0] += cs->deadcrootc;
		ns->litr1n[0] += ns->deadcrootn;
		cs->litrCbelow[0] += cs->deadcrootc;
		cs->deadcrootc = 0.0;
		ns->deadcrootn = 0.0;
	}

	if ((cs->deadcrootc_storage < 0 && fabs(cs->deadcrootc_storage) < CRIT_PREC) || (ns->deadcrootn_storage < 0 && fabs(ns->deadcrootn_storage) < CRIT_PREC))
	{
		cs->litr1c[0] += cs->deadcrootc_storage;
		ns->litr1n[0] += ns->deadcrootn_storage;
		cs->litrCbelow[0] += cs->deadcrootc_storage;
		cs->deadcrootc_storage = 0.0;
		ns->deadcrootn_storage = 0.0;
	}

	if ((cs->deadcrootc_transfer < 0 && fabs(cs->deadcrootc_transfer) < CRIT_PREC) || (ns->deadcrootn_transfer < 0 && fabs(ns->deadcrootn_transfer) < CRIT_PREC))
	{
		cs->litr1c[0] += cs->deadcrootc_transfer;
		ns->litr1n[0] += ns->deadcrootn_transfer;
		cs->litrCbelow[0] += cs->deadcrootc_transfer;
		cs->deadcrootc_transfer = 0.0;
		ns->deadcrootn_transfer = 0.0;
	}

	/************************/
	if (fabs(cs->gresp_transfer) < CRIT_PREC && cs->gresp_transfer != 0)
	{
		cs->litr1c[0] += cs->gresp_transfer;
		cs->litrCabove[0] += cs->gresp_transfer;
		cs->gresp_transfer = 0.0;

	}

	if (fabs(cs->cpool) < CRIT_PREC && cs->cpool != 0)
	{
		cs->litr1c[0] += cs->cpool;
		cs->litrCabove[0] += cs->cpool;
		cs->cpool = 0.0;

	}

	if (fabs(cs->NSCnw) < CRIT_PREC && cs->NSCnw != 0)
	{
		cs->litr1c[0] += cs->NSCnw;
		cs->litrCabove[0] += cs->NSCnw;
		cs->NSCnw = 0.0;

	}

	if (fabs(cs->NSCw) < CRIT_PREC && cs->NSCw != 0)
	{
		cs->litr1c[0] += cs->NSCw;
		cs->litrCabove[0] += cs->NSCw;
		cs->NSCw = 0.0;

	}

	if (fabs(cs->SCnw) < CRIT_PREC && cs->SCnw != 0)
	{
		cs->litr1c[0] += cs->SCnw;
		cs->litrCabove[0] += cs->SCnw;
		cs->SCnw = 0.0;

	}

	if (fabs(cs->SCw) < CRIT_PREC && cs->SCw != 0)
	{
		cs->litr1c[0] += cs->SCw;
		cs->litrCabove[0] += cs->SCw;
		cs->SCw = 0.0;

	}

	/************************/
	if ((fabs(cs->STDBc_leaf) < CRIT_PREC && fabs(cs->STDBc_leaf) != 0) || (fabs(ns->STDBn_leaf) < CRIT_PREC && fabs(ns->STDBn_leaf) != 0))
	{
		cs->litr1c[0] += cs->STDBc_leaf;
		ns->litr1n[0] += ns->STDBn_leaf;
		cs->litrCabove[0] += cs->STDBc_leaf;
		cs->STDBc_leaf = 0.0;
		ns->STDBn_leaf = 0.0;
	}

	if ((fabs(cs->STDBc_froot) < CRIT_PREC && fabs(cs->STDBc_froot) != 0) || (fabs(ns->STDBn_froot) < CRIT_PREC && fabs(ns->STDBn_froot) != 0))
	{
		cs->litr1c[0] += cs->STDBc_froot;
		ns->litr1n[0] += ns->STDBn_froot;
		cs->litrCbelow[0] += cs->STDBc_froot;
		cs->STDBc_froot = 0.0;
		ns->STDBn_froot = 0.0;
	}

	if ((fabs(cs->STDBc_softstem) < CRIT_PREC && fabs(cs->STDBc_softstem) != 0) || (fabs(ns->STDBn_softstem) < CRIT_PREC && fabs(ns->STDBn_softstem) != 0))
	{
		cs->litr1c[0] += cs->STDBc_softstem;
		ns->litr1n[0] += ns->STDBn_softstem;
		cs->litrCabove[0] += cs->STDBc_softstem;
		cs->STDBc_softstem = 0.0;
		ns->STDBn_softstem = 0.0;
	}

	if ((fabs(cs->STDBc_yield) < CRIT_PREC && fabs(cs->STDBc_yield) != 0) || (fabs(ns->STDBn_yield) < CRIT_PREC && fabs(ns->STDBn_yield) != 0))
	{
		cs->litr1c[0] += cs->STDBc_yield;
		ns->litr1n[0] += ns->STDBn_yield;
		cs->litrCabove[0] += cs->STDBc_yield;
		cs->STDBc_yield = 0.0;
		ns->STDBn_yield = 0.0;
	}



	/************************/
	if ((fabs(cs->CTDBc_leaf) < CRIT_PREC && fabs(cs->CTDBc_leaf) != 0) || (fabs(ns->CTDBn_leaf) < CRIT_PREC && fabs(ns->CTDBn_leaf) != 0))
	{
		cs->litr1c[0] += cs->CTDBc_leaf;
		ns->litr1n[0] += ns->CTDBn_leaf;
		cs->litrCabove[0] += cs->CTDBc_leaf;
		cs->CTDBc_leaf = 0.0;
		ns->CTDBn_leaf = 0.0;
	}

	if ((fabs(cs->CTDBc_froot) < CRIT_PREC && fabs(cs->CTDBc_froot) != 0) || (fabs(ns->CTDBn_froot) < CRIT_PREC && fabs(ns->CTDBn_froot) != 0))
	{
		cs->litr1c[0] += cs->CTDBc_froot;
		ns->litr1n[0] += ns->CTDBn_froot;
		cs->litrCbelow[0] += cs->CTDBc_froot;
		cs->CTDBc_froot = 0.0;
		ns->CTDBn_froot = 0.0;
	}

	if ((fabs(cs->CTDBc_softstem) < CRIT_PREC && fabs(cs->CTDBc_softstem) != 0) || (fabs(ns->CTDBn_softstem) < CRIT_PREC && fabs(ns->CTDBn_softstem) != 0))
	{
		cs->litr1c[0] += cs->CTDBc_softstem;
		ns->litr1n[0] += ns->CTDBn_softstem;
		cs->litrCabove[0] += cs->CTDBc_softstem;
		cs->CTDBc_softstem = 0.0;
		ns->CTDBn_softstem = 0.0;
	}

	if ((fabs(cs->CTDBc_yield) < CRIT_PREC && fabs(cs->CTDBc_yield) != 0) || (fabs(ns->CTDBn_yield) < CRIT_PREC && fabs(ns->CTDBn_yield) != 0))
	{
		cs->litr1c[0] += cs->CTDBc_yield;
		ns->litr1n[0] += ns->CTDBn_yield;
		cs->litrCabove[0] += cs->CTDBc_yield;
		cs->CTDBc_yield = 0.0;
		ns->CTDBn_yield = 0.0;
	}



	/************************/
	/* test for litter and soil poils in multilayer soil. Excess goes to HR sink (C) or volatilized sink (N) */


	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if ((cs->soil1c[layer] != 0 && fabs(cs->soil1c[layer]) < CRIT_PREC) || (ns->soil1n[layer] != 0 && fabs(ns->soil1n[layer]) < CRIT_PREC))
		{
			cs->HRsoil1_snk += cs->soil1c[layer];
			ns->Nprec_snk += ns->soil1n[layer];
			cs->soil1c[layer] = 0.0;
			ns->soil1n[layer] = 0.0;
		}
		if ((cs->soil2c[layer] != 0 && fabs(cs->soil2c[layer]) < CRIT_PREC) || (ns->soil2n[layer] != 0 && fabs(ns->soil2n[layer]) < CRIT_PREC))
		{
			cs->HRsoil2_snk += cs->soil2c[layer];
			ns->Nprec_snk += ns->soil2n[layer];
			cs->soil2c[layer] = 0.0;
			ns->soil2n[layer] = 0.0;
		}
		if ((cs->soil3c[layer] != 0 && fabs(cs->soil3c[layer]) < CRIT_PREC) || (ns->soil3n[layer] != 0 && fabs(ns->soil3n[layer]) < CRIT_PREC))
		{
			cs->HRsoil3_snk += cs->soil3c[layer];
			ns->Nprec_snk += ns->soil3n[layer];
			cs->soil3c[layer] = 0.0;
			ns->soil3n[layer] = 0.0;
		}
		if ((cs->soil4c[layer] != 0 && fabs(cs->soil4c[layer]) < CRIT_PREC) || (ns->soil4n[layer] != 0 && fabs(ns->soil4n[layer]) < CRIT_PREC))
		{
			cs->HRsoil4_snk += cs->soil4c[layer];
			ns->Nprec_snk += ns->soil4n[layer];
			cs->soil4c[layer] = 0.0;
			ns->soil4n[layer] = 0.0;
		}


		if ((cs->litr1c[layer] != 0 && fabs(cs->litr1c[layer]) < CRIT_PREC) || (ns->litr1n[layer] != 0 && fabs(ns->litr1n[layer]) < CRIT_PREC))
		{
			cs->HRlitr1_snk += cs->litr1c[layer];
			ns->Nprec_snk += ns->litr1n[layer];

			litrClayer = cs->litrCabove[layer] + cs->litrCbelow[layer];
			cs->litrCabove[layer] -= cs->litr1c[layer] * (cs->litrCabove[layer] / litrClayer);
			cs->litrCbelow[layer] -= cs->litr1c[layer] * (cs->litrCabove[layer] / litrClayer);

			if (cs->litrCabove[layer] < 0 || cs->litrCbelow[layer] < 0 || (cs->litrCabove[layer] && cs->litrCabove[layer] / cs->litrCabove[layer] != 1) || (cs->litrCbelow[layer] && cs->litrCbelow[layer] / cs->litrCbelow[layer] != 1))
			{
				if (fabs(cs->litrCabove[layer]) < CRIT_PREC) cs->litrCabove[layer] = 0;
				if (fabs(cs->litrCbelow[layer]) < CRIT_PREC) cs->litrCbelow[layer] = 0;

				if (cs->litrCabove[layer] < 0 || cs->litrCbelow[layer] < 0 || (cs->litrCabove[layer] && cs->litrCabove[layer] / cs->litrCabove[layer] != 1) || (cs->litrCbelow[layer] && cs->litrCbelow[layer] / cs->litrCbelow[layer] != 1))
				{
					printf("\n");
					printf("ERROR: negative above/below litr in precision_control.c\n");
					errorCode = 1;
				}
			}

			cs->litr1c[layer] = 0.0;
			ns->litr1n[layer] = 0.0;
		}
		if ((cs->litr2c[layer] != 0 && fabs(cs->litr2c[layer]) < CRIT_PREC) || (ns->litr2n[layer] != 0 && fabs(ns->litr2n[layer]) < CRIT_PREC))
		{
			cs->HRlitr2_snk += cs->litr2c[layer];
			ns->Nprec_snk += ns->litr2n[layer];

			litrClayer = cs->litrCabove[layer] + cs->litrCbelow[layer];
			cs->litrCabove[layer] -= cs->litr2c[layer] * (cs->litrCabove[layer] / litrClayer);
			cs->litrCbelow[layer] -= cs->litr2c[layer] * (cs->litrCbelow[layer] / litrClayer);

			if (cs->litrCabove[layer] < 0 || cs->litrCbelow[layer] < 0 || (cs->litrCabove[layer] && cs->litrCabove[layer] / cs->litrCabove[layer] != 1) || (cs->litrCbelow[layer] && cs->litrCbelow[layer] / cs->litrCbelow[layer] != 1))
			{
				if (fabs(cs->litrCabove[layer]) < CRIT_PREC) cs->litrCabove[layer] = 0;
				if (fabs(cs->litrCbelow[layer]) < CRIT_PREC) cs->litrCbelow[layer] = 0;

				if (cs->litrCabove[layer] < 0 || cs->litrCbelow[layer] < 0 || (cs->litrCabove[layer] && cs->litrCabove[layer] / cs->litrCabove[layer] != 1) || (cs->litrCbelow[layer] && cs->litrCbelow[layer] / cs->litrCbelow[layer] != 1))

				{
					printf("\n");
					printf("ERROR: negative above/below litr in precision_control.c\n");
					errorCode = 1;
				}
			}
			cs->litr2c[layer] = 0.0;
			ns->litr2n[layer] = 0.0;
		}
		if ((cs->litr3c[layer] != 0 && fabs(cs->litr3c[layer]) < CRIT_PREC) || (ns->litr3n[layer] != 0 && ns->litr3n[layer] < 0 && fabs(ns->litr3n[layer]) < CRIT_PREC))
		{
			cs->HRlitr4_snk += cs->litr3c[layer]; /* NO LITR3C HR SINK */
			ns->Nprec_snk += ns->litr3n[layer];

			litrClayer = cs->litrCabove[layer] + cs->litrCbelow[layer];
			cs->litrCabove[layer] -= cs->litr3c[layer] * (cs->litrCabove[layer] / litrClayer);
			cs->litrCbelow[layer] -= cs->litr3c[layer] * (cs->litrCbelow[layer] / litrClayer);

			if (cs->litrCabove[layer] < 0 || cs->litrCbelow[layer] < 0 || (cs->litrCabove[layer] && cs->litrCabove[layer] / cs->litrCabove[layer] != 1) || (cs->litrCbelow[layer] && cs->litrCbelow[layer] / cs->litrCbelow[layer] != 1))
			{
				if (fabs(cs->litrCabove[layer]) < CRIT_PREC) cs->litrCabove[layer] = 0;
				if (fabs(cs->litrCbelow[layer]) < CRIT_PREC) cs->litrCbelow[layer] = 0;

				if (cs->litrCabove[layer] < 0 || cs->litrCbelow[layer] < 0 || (cs->litrCabove[layer] && cs->litrCabove[layer] / cs->litrCabove[layer] != 1) || (cs->litrCbelow[layer] && cs->litrCbelow[layer] / cs->litrCbelow[layer] != 1))
				{
					printf("\n");
					printf("ERROR: negative above/below litr in precision_control.c\n");
					errorCode = 1;
				}
			}

			cs->litr3c[layer] = 0.0;
			ns->litr3n[layer] = 0.0;
		}
		if ((cs->litr4c[layer] != 0 && fabs(cs->litr4c[layer]) < CRIT_PREC) || (ns->litr4n[layer] != 0 && fabs(ns->litr4n[layer]) < CRIT_PREC))
		{
			cs->HRlitr4_snk += cs->litr4c[layer];
			ns->Nprec_snk += ns->litr4n[layer];

			litrClayer = cs->litrCabove[layer] + cs->litrCbelow[layer];
			cs->litrCabove[layer] -= cs->litr4c[layer] * (cs->litrCabove[layer] / litrClayer);
			cs->litrCbelow[layer] -= cs->litr4c[layer] * (cs->litrCbelow[layer] / litrClayer);

			if (cs->litrCabove[layer] < 0 || cs->litrCbelow[layer] < 0 || (cs->litrCabove[layer] && cs->litrCabove[layer] / cs->litrCabove[layer] != 1) || (cs->litrCbelow[layer] && cs->litrCbelow[layer] / cs->litrCbelow[layer] != 1))
			{
				if (fabs(cs->litrCabove[layer]) < CRIT_PREC) cs->litrCabove[layer] = 0;
				if (fabs(cs->litrCbelow[layer]) < CRIT_PREC) cs->litrCbelow[layer] = 0;

				if (cs->litrCabove[layer] < 0 || cs->litrCbelow[layer] < 0 || (cs->litrCabove[layer] && cs->litrCabove[layer] / cs->litrCabove[layer] != 1) || (cs->litrCbelow[layer] && cs->litrCbelow[layer] / cs->litrCbelow[layer] != 1))
				{
					printf("\n");
					printf("ERROR: negative above/below litr in precision_control.c\n");
					errorCode = 1;
				}
			}

			cs->litr4c[layer] = 0.0;
			ns->litr4n[layer] = 0.0;
		}

		if ((cs->cwdc[layer] != 0 && fabs(cs->cwdc[layer]) < CRIT_PREC) || (ns->cwdn[layer] != 0 && fabs(ns->cwdn[layer]) < CRIT_PREC))
		{
			cs->HRlitr4_snk += cs->cwdc[layer];
			ns->Nprec_snk += ns->cwdn[layer];
			cs->cwdc[layer] = 0.0;
			ns->cwdn[layer] = 0.0;
			cs->cwdCabove[layer] = 0.0;
			cs->cwdCbelow[layer] = 0.0;
		}

	}


	if (fabs(ns->retransn) < CRIT_PREC && ns->retransn != 0)
	{
		ns->Nprec_snk += ns->retransn;
		ns->retransn = 0.0;
	}
	/* additional tests for soil mineral N */

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




	/* WATER STATE VARIABLES */

	/* multilayer soil */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (ws->soilw[layer] < 0 && fabs(ws->soilw[layer]) < CRIT_PREC)
		{
			ws->soilEVP_snk += ws->soilw[layer];
			ws->soilw[layer] = 0.0;
		}

	}

	if (ws->snoww < 0 && fabs(ws->snoww) < CRIT_PREC)
	{
		ws->snowSUBL_snk += ws->snoww;
		ws->snoww = 0.0;
	}


	if (ws->canopyw < 0 && fabs(ws->canopyw) < CRIT_PREC)
	{
		ws->EVPcanopyw_snk += ws->canopyw;
		ws->canopyw = 0.0;
	}

	if (ws->pondw < 0 && fabs(ws->pondw) < CRIT_PREC)
	{
		ws->pondEVP_snk += ws->pondw;
		ws->pondw = 0.0;
	}

	/* transfer value: NH4, NO3, DOC, DON - > content_array  etc.*/
	if (!errorCode && check_soilcontent(-1, 0, sprop, cs, ns, soilInfo))
	{
		printf("ERROR in check_soilcontent.c for groundwater_preproc.c\n");
		errorCode = 1;
	}

	if (GWlayer > DATA_GAP && GWlayer < N_SOILLAYERS)
	{ 
		for (dm = 0; dm < N_DISSOLVMATER; dm++)
		{
			if (soilInfo->content_soil[dm][CFlayer] == 0)
			{
				soilInfo->content_CAPILcf[dm] = 0;
				soilInfo->content_NORMcf[dm] = 0;
			}
			if (soilInfo->content_soil[dm][GWlayer] == 0)
			{
				soilInfo->content_SATgw[dm] = 0;
				soilInfo->content_CAPILgw[dm] = 0;
				soilInfo->content_NORMgw[dm] = 0;

			}
		}
	}
	return(errorCode);
}	
