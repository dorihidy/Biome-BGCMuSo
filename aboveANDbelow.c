 /* 
aboveANDbelow.c
CALCULATING aboveground and belowground biomass 

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

int aboveANDbelow(soilprop_struct* sprop, epvar_struct* epv, cstate_struct* cs, cflux_struct* cf)
{

	/* soilCovering parameters */
	
	int layer, errorCode; 
	double litrINabove, litrINbelow, litrOUTabove, cwdINabove, cwdINbelow, cwdOUTabove, litrC, cwdC;


	errorCode = 0;
	litrINabove=litrINbelow= litrOUTabove = cwdINabove = cwdINbelow = cwdOUTabove = litrC = cwdC = 0;



	/* Summaraizing above- and belowground biomass */


	cs->litrCabove_total = 0;
	cs->litrCbelow_total = 0;
	cs->cwdCabove_total = 0;
	cs->cwdCbelow_total = 0;

	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{ 
	

		litrINabove = (cf->leafc_to_litr1c + cf->leafc_to_litr2c + cf->leafc_to_litr3c + cf->leafc_to_litr4c +
			cf->softstemc_to_litr1c + cf->softstemc_to_litr2c + cf->softstemc_to_litr3c + cf->softstemc_to_litr4c +
			cf->yieldc_to_litr1c + cf->yieldc_to_litr2c + cf->yieldc_to_litr3c + cf->yieldc_to_litr4c +
			cf->m_leafc_to_litr1c + cf->m_leafc_to_litr2c + cf->m_leafc_to_litr3c + cf->m_leafc_to_litr4c +
			cf->m_softstemc_to_litr1c + cf->m_softstemc_to_litr2c + cf->m_softstemc_to_litr3c + cf->m_softstemc_to_litr4c +
			cf->m_yieldc_to_litr1c + cf->m_yieldc_to_litr2c + cf->m_yieldc_to_litr3c + cf->m_yieldc_to_litr4c +
			cf->m_leafc_storage_to_litr1c + cf->m_yieldc_storage_to_litr1c + cf->m_softstemc_storage_to_litr1c + cf->m_livestemc_storage_to_litr1c + cf->m_deadstemc_storage_to_litr1c +
			cf->m_leafc_transfer_to_litr1c + cf->m_yieldc_transfer_to_litr1c + cf->m_softstemc_transfer_to_litr1c + cf->m_livestemc_transfer_to_litr1c + cf->m_deadstemc_transfer_to_litr1c +
			cf->m_gresp_storage_to_litr1c + cf->m_gresp_transfer_to_litr1c +
			cf->GRZ_to_litr1c + cf->GRZ_to_litr2c + cf->GRZ_to_litr3c + cf->GRZ_to_litr4c +
			cf->FRZ_to_litr1c + cf->FRZ_to_litr2c + cf->FRZ_to_litr3c + cf->FRZ_to_litr4c +
			cf->STDBc_leaf_to_litr + cf->STDBc_softstem_to_litr + cf->STDBc_yield_to_litr + cf->CTDBc_leaf_to_litr + cf->CTDBc_softstem_to_litr + cf->CTDBc_yield_to_litr) * sprop->PROPlayerDC[layer];

		litrINbelow = (cf->m_frootc_storage_to_litr1c + cf->m_frootc_transfer_to_litr1c) * sprop->PROPlayerDC[layer] +
		              (cf->frootc_to_litr1c + cf->frootc_to_litr2c + cf->frootc_to_litr3c + cf->frootc_to_litr4c +
			           cf->m_frootc_to_litr1c + cf->m_frootc_to_litr2c + cf->m_frootc_to_litr3c + cf->m_frootc_to_litr4c +
					   cf->STDBc_froot_to_litr + cf->CTDBc_froot_to_litr) * epv->rootlengthProp[layer];

		cwdINabove = (cf->m_livestemc_to_cwdc + cf->m_deadstemc_to_cwdc + cf->CTDBc_cstem_to_cwd) * sprop->PROPlayerDC[layer];

		cwdINbelow = (cf->m_livecrootc_to_cwdc + cf->m_deadcrootc_to_cwdc + cf->CTDBc_croot_to_cwd) * epv->rootlengthProp[layer];


		/* if proportion of decompositon is greater than 0 - also abovegrdound biomass (else: onlya belowground) */
		if (sprop->PROPlayerDC[layer] > 0)
		{
			if (cs->litrCabove[layer] + cs->litrCbelow[layer])
			{
				litrOUTabove = (cf->litr1c_to_soil1c[layer] + cf->litr1_hr[layer] + cf->litr1c_to_release[layer] + cf->m_litr1c_to_fire[layer] +
					cf->litr2c_to_soil2c[layer] + cf->litr2_hr[layer] + cf->litr2c_to_release[layer] + cf->m_litr2c_to_fire[layer] +
					cf->m_litr3c_to_fire[layer] +
					cf->litr4c_to_soil3c[layer] + cf->litr4_hr[layer] + cf->litr4c_to_release[layer] + cf->m_litr4c_to_fire[layer]) * (cs->litrCabove[layer] / (cs->litrCabove[layer] + cs->litrCbelow[layer]));
			}
			else
			{
				if (litrINabove + litrINbelow)
				{
					litrOUTabove = (cf->litr1c_to_soil1c[layer] + cf->litr1_hr[layer] + cf->litr1c_to_release[layer] + cf->m_litr1c_to_fire[layer] +
						cf->litr2c_to_soil2c[layer] + cf->litr2_hr[layer] + cf->litr2c_to_release[layer] + cf->m_litr2c_to_fire[layer] +
						cf->m_litr3c_to_fire[layer] +
						cf->litr4c_to_soil3c[layer] + cf->litr4_hr[layer] + cf->litr4c_to_release[layer] + cf->m_litr4c_to_fire[layer]) * (litrINabove / (litrINabove + litrINbelow));
				}
			}

			if (cs->cwdCabove[layer] + cs->cwdCbelow[layer])
			{
				cwdOUTabove = (cf->cwdc_to_litr2c[layer] + cf->cwdc_to_litr3c[layer] + cf->cwdc_to_litr4c[layer] + cf->m_cwdc_to_fire[layer]) * (cs->cwdCabove[layer] / (cs->cwdCabove[layer] + cs->cwdCbelow[layer]));
			}
			else
			{
				if (cwdINabove + cwdINbelow) cwdOUTabove = (cf->cwdc_to_litr2c[layer] + cf->cwdc_to_litr3c[layer] + cf->cwdc_to_litr4c[layer] + cf->m_cwdc_to_fire[layer]) * (cwdINabove / (cwdINabove + cwdINbelow));
			}

		}
		else
		{
			litrINabove = 0;
			litrOUTabove = 0;
			cwdINabove = 0;
			cwdOUTabove = 0;
		}




		cs->litrCabove[layer] += litrINabove - litrOUTabove;
		cs->litrCbelow[layer] = cs->litr1c[layer] + cs->litr2c[layer] + cs->litr3c[layer] + cs->litr4c[layer] - cs->litrCabove[layer];


		cs->cwdCabove[layer] += cwdINabove - cwdOUTabove;
		cs->cwdCbelow[layer] = cs->cwdc[layer] - cs->cwdCabove[layer];

		if (fabs(cs->litrCabove[layer]) < CRIT_PREC) cs->litrCabove[layer] = 0;
		if (fabs(cs->litrCbelow[layer]) < CRIT_PREC) cs->litrCbelow[layer] = 0;
		if (fabs(cs->cwdCabove[layer]) < CRIT_PREC) cs->cwdCabove[layer] = 0;
		if (fabs(cs->cwdCbelow[layer]) < CRIT_PREC) cs->cwdCbelow[layer] = 0;

		/* conrtol */
		if (cs->litrCabove[layer] < 0 || cs->litrCbelow[layer] < 0 || cs->cwdCabove[layer] < 0 || cs->cwdCbelow[layer] < 0)
		{
			printf("\n");
			printf("ERROR: negative above/below litr in aboveANDbelow.c\n");
			errorCode = 1;
		}

		if (cs->litrCabove[layer] - (cs->litr1c[layer] + cs->litr2c[layer] + cs->litr3c[layer] + cs->litr4c[layer]) > 0)
		{
			if (cs->litrCabove[layer] - (cs->litr1c[layer] + cs->litr2c[layer] + cs->litr3c[layer] + cs->litr4c[layer]) > CRIT_PREC)
			{
				printf("\n");
				printf("ERROR: negative above/below litr in aboveANDbelow.c\n");
				errorCode = 1;
			}
			else
				cs->litrCabove[layer] = cs->litr1c[layer] + cs->litr2c[layer] + cs->litr3c[layer] + cs->litr4c[layer];

		}
		if (cs->cwdCabove[layer] - cs->cwdc[layer] > 0)
		{
			if (cs->cwdCabove[layer] - cs->cwdc[layer] > CRIT_PREC)
			{
				printf("\n");
				printf("ERROR: negative above/below litr in aboveANDbelow.c\n");
				errorCode = 1;
			}
			else
				cs->cwdCabove[layer] = cs->cwdc[layer];

		}

		cs->litrCabove_total += cs->litrCabove[layer];
		cs->litrCbelow_total += cs->litrCbelow[layer];
        cs->cwdCabove_total  += cs->cwdCabove[layer];
		cs->cwdCbelow_total += cs->cwdCbelow[layer];

		litrC += cs->litr1c[layer] +  cs->litr2c[layer] + cs->litr3c[layer] + cs->litr4c[layer];
		cwdC  += cs->cwdc[layer];

	}

	

	/* control */
	if (fabs(cs->litrCabove_total + cs->litrCbelow_total - litrC) > CRIT_PREC_lenient)
	{
		printf("ERROR: in above/below litter calculation in aboveANDbelow.c\n");
		errorCode = 1;
	}

	if (fabs(cs->cwdCabove_total + cs->cwdCbelow_total - cwdC) > CRIT_PREC_lenient)
	{
		printf("ERROR: in above/below cwdc calculation in aboveANDbelow.c\n");
		errorCode = 1;
	}



   return (errorCode);
}
	