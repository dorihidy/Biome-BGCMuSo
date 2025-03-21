/* 
multilayer_rootDepth.c
calculation of changing rooting depth based on empirical function 
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
#include "bgc_func.h"
#include "bgc_constants.h"

int multilayer_rootDepth(const epconst_struct* epc, const soilprop_struct* sprop, const cstate_struct* cs, siteconst_struct* sitec,  epvar_struct* epv)
{

	int errorCode=0;
	int layer;

	double RLprop_sum1, RLprop_sum2, RLprop_sum0, frootc;
	double maxRD;


	/* initalizing internal variables */

	RLprop_sum1=RLprop_sum2=RLprop_sum0=0.0;

	if (sprop->soildepth < epc->rootzoneDepth_max)
		maxRD = sprop->soildepth;
	else
		maxRD = epc->rootzoneDepth_max;

	frootc = cs->frootc;



	/* ***************************************************************************************************** */	
	/* 1. Calculating rooting depth in case of non-wwody ecosystems (based on Campbell and Diaz, 1988) 
	      actual rooting depth determines the rootzone depth (epv->n_rootlayers) */
	

	if (frootc > CRIT_PREC) 
	{
		if (frootc < epc->rootlength_par1)
		{
			/* par1: root weight corresponding to max root depth, par2: root depth function shape parameter */
			epv->rootDepth = epv->germDepth + (maxRD-epv->germDepth) * pow(frootc / epc->rootlength_par1, epc->rootlength_par2);
		}
		else
			epv->rootDepth = maxRD;

	}
	else
		epv->rootDepth = 0;

	if (epc->woody) epv->rootDepth = maxRD;

	if (epv->germ_layer && epv->rootDepth)
		epv->rootlength = epv->rootDepth - sitec->soillayer_depth[epv->germ_layer-1];
	else
		epv->rootlength = epv->rootDepth;

	if (epv->rootlength > 0 && epv->rootlength < 0.01)
	{
		epv->rootlength = 0.01;
		epv->rootDepth = epv->germDepth + epv->rootlength;
	}

	/* ***************************************************************************************************** */	
	/* 3. Calculating the number of the soil layers in which root can be found. It determines the rootzone depth (epv->n_rootlayers) */
	
	if (!errorCode && calc_nrootlayers(1, maxRD, frootc, sitec, epv))
	{
		printf("\n");
		printf("ERROR in calc_nrootlayers.c for multilayer_rootDepth.c\n");
		errorCode=1;
	}
	
	/* ***************************************************************************************************** */	
	/* 4. Calculating the distribution of the root in the soil layers based on empirical function (Jarvis, 1989)*/
	

	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		/* live root distribution for soil water calculation */
		if (layer < epv->n_rootlayers && layer >= epv->germ_layer)
		{
			if (epv->n_rootlayers > 1)
				epv->rootlengthProp[layer] = epc->rootdistrib_param * (sitec->soillayer_thickness[layer] / epv->rootlength) * exp(-epc->rootdistrib_param * (sitec->soillayer_midpoint[layer] / epv->rootlength));
			else
				epv->rootlengthProp[layer] = 1;

			RLprop_sum1 += epv->rootlengthProp[layer];

		}
		else
			epv->rootlengthProp[layer]   = 0;

		/* dead root distribution for soil water calculation */
		if (epv->n_maxrootlayers)
		{
			if (layer < epv->n_maxrootlayers-1)
				epv->rootlengthLandD_prop[layer] = sitec->soillayer_thickness[layer] / maxRD;
			else
			{
				if (layer == epv->n_maxrootlayers-1 && layer > 1) 
					epv->rootlengthLandD_prop[layer] = (maxRD-sitec->soillayer_depth[layer-1]) / maxRD;
				else
					epv->rootlengthLandD_prop[layer] = 0;
			}
			 
		}
		else
			epv->rootlengthLandD_prop[layer]   = 0;
		
		RLprop_sum0 += epv->rootlengthLandD_prop[layer];
	}

	if ((epv->n_maxrootlayers && fabs(1-RLprop_sum0) > CRIT_PREC) || (!epv->n_maxrootlayers && RLprop_sum0 != 0)) 
	{
		printf("\n");
		printf("ERROR in multilayer_rootDepth: sum of soillayer_RZportion is 0.0\n");
		errorCode=1;
	}

	if (RLprop_sum1 == 0) 
	{
		if (epv->rootDepth)
		{
			printf("\n");
			printf("ERROR in multilayer_rootDepth: sum of soillayer_RZportion is 0.0\n");
			errorCode=1;
		}
	}

	/* correction */
	if (RLprop_sum1)
	{
		for (layer =0; layer < N_SOILLAYERS; layer++)
		{	
			epv->rootlengthProp[layer] = epv->rootlengthProp[layer] / RLprop_sum1;
			RLprop_sum2                 += epv->rootlengthProp[layer];
		}

		if (fabs(1. - RLprop_sum2) > CRIT_PREC)
		{
			printf("\n");
			printf("ERROR in multilayer_rootDepth: sum of soillayer_RZportion is not equal to 1.0\n");
			errorCode=1;
		}
	}
	

	/* ***************************************************************************************************** */	
	/* 5. calculation of plant height (based on 4M)*/

	/* par1: stem weight corresponding to max plant height, par2: plant height function shape parameter */
	if (epc->woody)
	{
		epv->plantHeight = epc->plantHeight_max*(1-exp((-5/epc->plantHeight_par1)*(cs->livestemc+cs->deadstemc)));
	}
	else
	{
		
		epv->plantHeight = epc->plantHeight_max * pow(((cs->softstemc+cs->STDBc_softstem)/epc->plantHeight_par1),epc->plantHeight_par2);
		
	}
	if (epv->plantHeight > epc->plantHeight_max) epv->plantHeight = epc->plantHeight_max;

	


	return(errorCode);
}

