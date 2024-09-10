 /*
diffusCalc.c
Calculation of diffusion flux between two soil layers

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2022, D. Hidy [dori.hidy@gmail.com]
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

int diffusCalc(const soilprop_struct* sprop, double dz0, double VWC0, double VWC0_sat, double VWC0_fc, double VWC0_wp, double VWC0_limit,
	                                         double dz1, double VWC1, double VWC1_sat, double VWC1_fc, double VWC1_wp, double VWC1_limit, double* soilwDiffus)
{

	int errorCode = 0;
	int noVWC1limit = 0;
	double ESWi0, ESWi1, THETi0, THETi1, innerTHET, innerESW, DBAR, GRAD, FLOW;
	double dz0_cm, dz1_cm, soilw0, soilw1, soilw_sat0, soilw_sat1, soilw_fc0, soilw_fc1, soilwDiffus_act;
	double soilw0_limit, soilw1_limit, diff1, diff2, diff;

	dz0_cm = dz0 * 100;
	dz1_cm = dz1 * 100;




	soilw0_limit = VWC0_limit * dz0 * water_density;
	soilw1_limit = VWC1_limit * dz1 * water_density;

	soilw0 = VWC0 * dz0 * water_density;
	soilw1 = VWC1 * dz1 * water_density;

	/* the plant-extractable soil water */
	ESWi0 = (VWC0_fc - VWC0_wp);
	ESWi1 = (VWC1_fc - VWC1_wp);


	/* interation for determine diffusion 	*/
	THETi0 = MIN(VWC0 - VWC0_wp, ESWi0);
	THETi1 = MIN(VWC1 - VWC1_wp, ESWi1);


	THETi0 = MAX(THETi0, 0);
	THETi1 = MAX(THETi1, 0);

	innerTHET = (THETi0 * 0.5 + THETi1 * 0.5);
	DBAR = sprop->p1diffus_tipping * exp(sprop->p2diffus_tipping * innerTHET);

	DBAR = MIN(DBAR, sprop->p3diffus_tipping);



	innerESW = (ESWi0 * 0.5 + ESWi1 * 0.5);
	GRAD     = innerESW * (THETi1/ESWi1 - THETi0/ESWi0);
	
	FLOW   = DBAR * GRAD / ((dz0_cm + dz1_cm) / 2.);

	if (fabs(FLOW) > CRIT_PREC)
		soilwDiffus_act = -1*(FLOW / m_to_cm) * water_density;
	else
		soilwDiffus_act = 0;
	
	/* tipping diffusion limitation */
	soilw_fc0  = VWC0_fc  * dz0 * water_density;
	soilw_fc1  = VWC1_fc  * dz1 * water_density;
	soilw_sat0 = VWC0_sat * dz0 * water_density;
	soilw_sat1 = VWC1_sat * dz1 * water_density;

	/* control to remain balance and avoid increasing VWC above equalization value or field caparity */
	if (soilwDiffus_act != 0)
	{
		/* downward movement */
		if (soilwDiffus_act > 0)
		{
			/* equalization control */
			diff1 = soilw1 + soilwDiffus_act - soilw1_limit;
			diff2 = soilw0 - soilwDiffus_act - soilw0_limit;

			if (diff1 > 0 || diff2 < 0)
			{
				if (diff1 > 0 && diff2 < 0)
				{
					if (fabs(diff1) > fabs(diff2))
						diff = fabs(diff1);
					else
						diff = fabs(diff2);
				}
				else
				{
					if (diff1 > 0)
						diff = diff1;
					else
						diff = fabs(diff2);
				}
				soilwDiffus_act -= diff;
				
				/* FC-control */
				diff = soilw1 + soilwDiffus_act - soilw_fc1;
				if (diff > 0)
				{
					soilwDiffus_act -= diff;
				}

			}
			if (soilwDiffus_act < 0)
			{
				if (fabs(soilwDiffus_act) > CRIT_PRECwater && !errorCode)
				{
					printf("\n");
					printf("ERROR in diffusCalc.c from tipping.c\n");
					errorCode = 1;
				}
				else
					soilwDiffus_act = 0;
			}
			

		}
		/* upward movement */
		else
		{

			diff1 = soilw0 - soilwDiffus_act - soilw0_limit;
			diff2 = soilw1 + soilwDiffus_act - soilw1_limit;

			if (noVWC1limit == 1) diff2 = 0;

			if (diff1 > 0 || diff2 < 0)
			{
				if (diff1 > 0 && diff2 < 0)
				{
					if (fabs(diff1) > fabs(diff2))
						diff = fabs(diff1);
					else
						diff = fabs(diff2);
				}
				else
				{
					if (diff1 > 0)
						diff = diff1;
					else
						diff = fabs(diff2);
				}

				soilwDiffus_act += diff;
			}

			/* FC-control */
			diff = soilw0 - soilwDiffus_act - soilw_fc0;
			if (diff > 0)
			{
				soilwDiffus_act += diff;
			}

			if (soilwDiffus_act > 0)
			{
				if (soilwDiffus_act > CRIT_PRECwater && !errorCode)
				{
					printf("\n");
					printf("ERROR in diffusCalc.c from tipping.c\n");
					errorCode = 1;
				}
				else
					soilwDiffus_act = 0;

			}

		}
	}

	soilw0  -= soilwDiffus_act;
			
	/* control to avoid oversaturation*/
	if (soilw0 - soilw_sat0 > CRIT_PRECwater && !errorCode)
	{
		printf("\n");
		printf("ERROR in diffusCalc.c from tipping.c\n");
		errorCode = 1;
	}
			
	soilw1 += soilwDiffus_act;
			
	/* control to avoid oversaturation*/ 
	if (soilw1 - soilw_sat1 > CRIT_PRECwater && !errorCode)
	{
		printf("\n");
		printf("ERROR in diffusCalc.c from tipping.c\n");
		errorCode = 1;
	}

	
			
	VWC0  =  soilw0 / dz0 / water_density;
	VWC1  =  soilw1 / dz1 / water_density;
		
	*soilwDiffus = soilwDiffus_act;
	




	return (errorCode);

}
