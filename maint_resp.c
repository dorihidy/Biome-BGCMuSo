/* 
maint_resp.c
daily maintenance respiration

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2022, D. Hidy [dori.hidy@gmail.com]
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
#include "bgc_func.h"
#include "bgc_constants.h"

int maint_resp(const control_struct* ctrl, const planting_struct* PLT, const cstate_struct* cs, const nstate_struct* ns, const epconst_struct* epc, const metvar_struct* metv,
	          epvar_struct* epv, cflux_struct* cf)
{
	/*
	maintenance respiration routine
	Uses reference values at 20 deg C and an empirical relationship between
	tissue N content and respiration rate given in:

	Ryan, M.G., 1991. Effects of climate change on plant respiration.
	Ecological Applications, 1(2):157-167.
	
	Uses the same value of Q_10 (2.0) for all compartments, leaf, stem, 
	coarse and fine roots.
	
	From Ryan's figures and regressions equations, the maintenance respiration
	in kgC/day per kg of tissue N is:
	MRperN = 0.218 (kgC/kgN/d)
	
	Leaf maintenance respiration is calculated separately for day and
	night, since the PSN routine needs the daylight value.
	
	Leaf and fine root respiration are dependent on phenology.*/
	
	int errorCode=0;
	double t1;
	double exponent;
	double n_area_sun, n_area_shade, dlmr_area_sun, dlmr_area_shade;
	int layer;
	double Tsoil, frootn_layer, livecrootn_layer, q10, froot_MR,livecroot_MR,softstem_MR,livestem_MR;


	double MRperN = epc->MRperN;
	double acclim_const = -0.00794; 	/* 1. calculate the Tsoil exponenet regarding to froot_resp in multilayer soil */

	froot_MR = livecroot_MR = softstem_MR = livestem_MR = 0;


	/* ********************************************************* */
	/* 1. possibility to use temperature dependent Q10 */ 

	if (ctrl->q10depend_flag)
		q10= 3.22 - 0.046 * metv->Tavg;
	else
		q10 = Q10_VALUE;

	/* ********************************************************* */
	/* 2. Leaf day and night maintenance respiration when leaves on */

	if (cs->leafc && epv->projSLA_sun && epv->projSLA_shade)
	{
		t1 = ns->leafn * MRperN;
		
		/* leaf, day */
		exponent = (metv->Tday - 20.0) / 10.0;
		cf->leaf_day_MR = t1 * pow(q10, exponent) * metv->dayl / nSEC_IN_DAY;

		/* for day respiration, also determine rates of maintenance respiration
		per unit of projected leaf area in the sunlit and shaded portions of
		the canopy, for use in the photosynthesis routine */
		/* first, calculate the mass of N per unit of projected leaf area
		in each canopy fraction (kg N/m2 projected area) */
		n_area_sun   = 1.0/(epv->projSLA_sun * epc->leaf_cn);
		n_area_shade = 1.0/(epv->projSLA_shade * epc->leaf_cn);
		/* convert to respiration flux in kg C/m2 projected area/day, and
		correct for temperature */
		dlmr_area_sun   = n_area_sun * MRperN * pow(q10, exponent);
		dlmr_area_shade = n_area_shade * MRperN * pow(q10, exponent);
		/* finally, convert from mass to molar units, and from a daily rate to 
		a rate per second */
		epv->dlmr_area_sun = dlmr_area_sun/(nSEC_IN_DAY * 12.011e-9);
		epv->dlmr_area_shade = dlmr_area_shade/(nSEC_IN_DAY * 12.011e-9);
		
		/* leaf, night */
		exponent = (metv->Tnight - 20.0) / 10.0;
		cf->leaf_night_MR = t1 * pow(q10, exponent) * (nSEC_IN_DAY - metv->dayl) / nSEC_IN_DAY;
	}
	else /* no leaves on */
	{
		cf->leaf_day_MR = 0.0;
		epv->dlmr_area_sun = 0.0;
		epv->dlmr_area_shade = 0.0;
		cf->leaf_night_MR = 0.0;
	}

	/* ********************************************************* */
	/* 2. fine root maintenance respiration when fine roots on */
	
	if (cs->frootc)
	{
		for (layer = epv->germ_layer; layer < epv->n_rootlayers; layer++)
		{
			Tsoil = metv->Tsoil[layer];
			
			frootn_layer = ns->frootn * epv->rootlengthProp[layer];
			exponent = (Tsoil - 20.0) / 10.0;
			t1 = pow(q10, exponent);
			froot_MR += frootn_layer * MRperN * t1;
		}
		
	}
	else froot_MR = 0;/* no fine roots on */

	cf->froot_MR = froot_MR;

	/* ********************************************************* */
	/* 3. yield simulation maintenance respiration when yields on */
   
	if (cs->yieldc && PLT->PLT_num == 0)
	{
		exponent = (metv->Tavg - 20.0) / 10.0;
		t1 = pow(q10, exponent);
		cf->yield_MR = ns->yieldn * MRperN * t1;
	}
	/* no yields on */
	else cf->yield_MR = 0.0;
		

	/* ********************************************************* */
	/* 4. softstem maintenance respiration */
	
	if (cs->softstemc)
	{
		exponent = (metv->Tavg - 20.0) / 10.0;
		t1 = pow(q10, exponent);
		cf->softstem_MR = ns->softstemn * MRperN * t1;
	}

	/* ********************************************************* */
	/* 5. live stem maintenance respiration */
	
	if (cs->livestemc)
	{
		exponent = (metv->Tavg - 20.0) / 10.0;
		t1 = pow(q10, exponent);
		cf->livestem_MR = ns->livestemn * MRperN * t1;
	}
	
    /* ********************************************************* */
	/* 6. live coarse root maintenance respiration */
	
	if (cs->livecrootc)
	{
		for (layer = epv->germ_layer; layer < epv->n_rootlayers; layer++)
		{
			Tsoil = metv->Tsoil[layer];

			livecrootn_layer = ns->livecrootn * epv->rootlengthProp[layer];
			
			exponent = (Tsoil - 20.0) / 10.0;
			t1 = pow(q10, exponent);
			livecroot_MR += livecrootn_layer * MRperN * t1;
		}
		
	}
	
	cf->livecroot_MR = livecroot_MR;
	
	
	/* ********************************************************* */
	/* 7. acclimation (acc_flag=1 - only respiration is acclimated, acc_flag=3 - respiration and photosynt. are acclimated) */
	
	
	if (ctrl->resp_acclim_flag)
	{
		cf->leaf_day_MR   = cf->leaf_day_MR   * pow(10,(acclim_const * (metv->TavgRA10-20.0)));
		cf->leaf_night_MR = cf->leaf_night_MR * pow(10,(acclim_const * (metv->TavgRA10-20.0)));
		cf->froot_MR      = cf->froot_MR      * pow(10,(acclim_const * (metv->TavgRA10-20.0)));
		cf->yield_MR      = cf->yield_MR      * pow(10,(acclim_const * (metv->TavgRA10-20.0)));
		cf->livestem_MR   = cf->livestem_MR   * pow(10,(acclim_const * (metv->TavgRA10-20.0)));
		cf->livecroot_MR  = cf->livecroot_MR  * pow(10,(acclim_const * (metv->TavgRA10-20.0)));
		cf->softstem_MR   = cf->softstem_MR   * pow(10,(acclim_const * (metv->TavgRA10-20.0)));
	}


	
	
	return (errorCode);
}
