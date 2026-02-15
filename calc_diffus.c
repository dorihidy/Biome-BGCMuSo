 /*
calc_diffus.c
Calculation of diffusion flux between two soil layers

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

int calc_diffus(int layer, const soilprop_struct* sprop, double dz0, double VWC0, double VWC0_sat, double VWC0_EqFC, double VWC0_wp, double VWC0_hw,
	double dz1, double VWC1, double VWC1_sat, double VWC1_EqFC, double VWC1_wp, double VWC1_hw,
	double fluxLimit, double dLk, double* DBAR, double* soilwDiffus)
{

	int errorCode = 0;
	double ESWi0, ESWi1, THETi0, THETi1, innerTHET, innerESW, DBARlimited, GRAD, FLOW;
	double dz0_cm, dz1_cm, soilw0, soilw1, soilwDiffus_act;
	double soilw_hw0, soilw_hw1, soilw_sat0, soilw_sat1;

	dz0_cm = dz0 * 100;
	dz1_cm = dz1 * 100;


	soilw0 = VWC0 * dz0 * water_density;
	soilw1 = VWC1 * dz1 * water_density;

	soilw_sat0 = VWC0_sat * dz0 * water_density;
	soilw_sat1 = VWC1_sat * dz1 * water_density;

	soilw_hw0 = VWC0_hw * dz0 * water_density;
	soilw_hw1 = VWC1_hw * dz1 * water_density;


	/* the plant-extractable soil water */
	ESWi0 = (VWC0_EqFC - VWC0_wp);
	ESWi1 = (VWC1_EqFC - VWC1_wp);


	/* interation for determine diffusion 	*/
	THETi0 = MIN(VWC0_EqFC, VWC0) - VWC0_wp;
	THETi1 = MIN(VWC1_EqFC, VWC1) - VWC1_wp;


	innerTHET = (THETi0 * 0.5 + THETi1 * 0.5);
	//innerTHET = (VWC0 * 0.5 + VWC1 * 0.5);
	*DBAR = sprop->p1diffus[layer] * exp(sprop->p2diffus[layer] * innerTHET);


	DBARlimited = MIN(*DBAR, sprop->p3diffus[layer]);


	innerESW = (ESWi0 * 0.5 + ESWi1 * 0.5);
	GRAD = innerESW * (THETi1 / ESWi1 - THETi0 / ESWi0);

	if (dLk == DATA_GAP)
		FLOW = DBARlimited * GRAD / ((dz0_cm + dz1_cm) / 2.);
	else
		FLOW = DBARlimited * GRAD / dLk;

	if (fabs(FLOW) > CRIT_PREC)
		soilwDiffus_act = -1 * (FLOW / mm_to_cm);
	else
		soilwDiffus_act = 0;

 	/*  diffusion limitation: fluxLimit, hw and sat */
	if (soilwDiffus_act > 0)
	{
		if (fabs(soilwDiffus_act) > fabs(fluxLimit)) soilwDiffus_act = fabs(fluxLimit);
		if (soilw0 - soilw_hw0 - soilwDiffus_act < 0) soilwDiffus_act = soilw0 - soilw_hw0;
		if ((soilw_sat1 - soilw_hw1 > CRIT_PREC) && (soilw_sat1 - (soilw1 + soilwDiffus_act) < 0)) soilwDiffus_act = soilw_sat1 - soilw1;
	}
	else
	{
		if (fabs(soilwDiffus_act) > fabs(fluxLimit)) soilwDiffus_act = -1 * fabs(fluxLimit);
		if ((soilw_sat1 - soilw_hw1 > CRIT_PREC) && (soilw1 - soilw_hw1 + soilwDiffus_act < 0)) soilwDiffus_act = -1 * (soilw1 - soilw_hw1);
		if (soilw_sat0 - (soilw0 - soilwDiffus_act) < 0) soilwDiffus_act = -1 * (soilw_sat0 - soilw0);
	}



	/* update water content */
	soilw0 -= soilwDiffus_act;

	/* control to avoid oversaturation*/
	if (!errorCode && (soilw0 - soilw_sat0 > CRIT_PREC_lenient || soilw0 <= 0))
	{
		printf("\n");
		printf("ERROR in calc_diffus.c (soilw0 - soilw_sat0 > 0) %12.2f%12.2f%12.2f\n", soilw0, soilw_sat0, soilw0 - soilw_sat0);
		errorCode = 1;
	}

	soilw1 += soilwDiffus_act;

	/* control to avoid oversaturation (except of GWlayer: VWCsat=VWChw */
	if (!errorCode && (soilw_sat1 - soilw_hw1 > CRIT_PREC) && (soilw1 - soilw_sat1 > CRIT_PREC_lenient || soilw1 <= 0))
	{
		printf("\n");
		printf("ERROR in calc_diffus.c (soilw1 - soilw_sat1) %12.2f%12.2f\n", soilw1,soilw_sat1);
		errorCode = 1;
	}

	
			
	VWC0  =  soilw0 / dz0 / water_density;
	VWC1  =  soilw1 / dz1 / water_density;
		
	*soilwDiffus = soilwDiffus_act;
	




	return (errorCode);

}
