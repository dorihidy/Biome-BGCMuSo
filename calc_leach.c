 /*
calc_leach.c
Calculation of leaching flux between two soil layers (layer0 and layer1; materials: NO3, NH4, DON and DOC)

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

int calc_leach(control_struct* ctrl, soilInfo_struct* soilInfo, double wflux, double wstate0, double wstate1, double state0[N_DISSOLVMATER], double state1[N_DISSOLVMATER],
	                                                             double conc0[N_DISSOLVMATER], double conc1[N_DISSOLVMATER], double leachFlux[N_DISSOLVMATER])
{
	int dm, nc;
	int errorCode = 0;
	double conc0_act[N_DISSOLVMATER], conc1_act[N_DISSOLVMATER], leachFlux_act[N_DISSOLVMATER], state0_act[N_DISSOLVMATER], state1_act[N_DISSOLVMATER];
	double pool0, pool1;


	/* calculation of concentration */

	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		state0_act[dm] = state0[dm];
		state1_act[dm] = state1[dm];
		conc0_act[dm] = soilInfo->dissolv_prop[dm] * (state0_act[dm] / wstate0);
		conc1_act[dm] = soilInfo->dissolv_prop[dm] * (state1_act[dm] / wstate1);

		if (wflux > 0)
			leachFlux_act[dm] = conc0_act[dm] * wflux;
		else
			leachFlux_act[dm] = conc1_act[dm] * wflux;

		pool0 = state0_act[dm] - leachFlux_act[dm];
		pool1 = state1_act[dm] + leachFlux_act[dm];

		if (pool0 < 0)
		{
			leachFlux_act[dm] += pool0;
			/* limitleach_flag: flag of WARNING writing in log file (only at first time) */
			if (fabs(pool0) > CRIT_PREC && !ctrl->limitleach_flag) ctrl->limitleach_flag = 1;
		}

		if (pool1 < 0)
		{
			leachFlux_act[dm] -= pool1;
			/* limitdiffus_flag: flag of WARNING writing in log file (only at first time) */
			if (fabs(pool1) > CRIT_PREC && !ctrl->limitdiffus_flag) ctrl->limitdiffus_flag = 1;
		}

		if (fabs(leachFlux_act[dm]) < CRIT_PREC) leachFlux_act[dm] = 0;
	
		/* control CNratio of fluxes */
		for (nc = N_DISSOLVN; nc<N_DISSOLVMATER; nc++)
		{ 
			if (leachFlux_act[nc - 4] == 0) leachFlux_act[nc] = 0;
		}
	

		state0[dm] = state0_act[dm] - leachFlux_act[dm];
		state1[dm] = state1_act[dm] + leachFlux_act[dm];

		conc0[dm] = soilInfo->dissolv_prop[dm] * ((state0_act[dm] - leachFlux_act[dm]) / wstate0);
		conc1[dm] = soilInfo->dissolv_prop[dm] * ((state1_act[dm] + leachFlux_act[dm]) / wstate1);

		leachFlux[dm] = leachFlux_act[dm];

		if (state0[dm] < 0)
		{
			if (fabs(state0[dm]) > CRIT_PREC)
			{
				printf("\n");
				printf("ERROR calc_leach.c negative storage values.c\n");
				errorCode = 1;
			}
			else
			{
				leachFlux_act[dm] += state0[dm];
				state0[dm] = 0;
			}
		}

		if (state1[dm] < 0)
		{
			if (fabs(state1[dm]) > CRIT_PREC)
			{
				printf("\n");
				printf("ERROR calc_leach.c negative storage values.c\n");
				errorCode = 1;
			}
			else
			{
				leachFlux_act[dm] -= state1[dm];
				state1[dm] = 0;
			}
		}

			
	}




	return (errorCode);

}
