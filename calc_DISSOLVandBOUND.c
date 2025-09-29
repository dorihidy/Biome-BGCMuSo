/* 
calc_DISSOLVandBOUND.c
transformation and checking o fbounded and dissolved content of soil material
firsttime_flag: before calculation - 0: note pre values for calculating difference, 1: calculating difference
partlyORtotal_flag: 0 - DISSOLV or BOUND is affected, 1 - both are affected


*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2025, D. Hidy [dori.hidy@gmail.com]
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

int calc_DISSOLVandBOUND(int firsttime_flag, int partlyORtotal_flag, const soilprop_struct* sprop, soilInfo_struct* soilInfo)
{
	/*  intitialiazion*/

	int errorCode, layer, dm, CFlayer, GWlayer, N_layer;
	static double content_soil_pre[N_DISSOLVMATER][N_SOILLAYERS];
	static double content_NORMcf_pre[N_DISSOLVMATER], content_CAPILcf_pre[N_DISSOLVMATER];
	errorCode = 0;

	CFlayer = (int)sprop->CFlayer;
	GWlayer = (int)sprop->GWlayer;

	if (GWlayer == DATA_GAP)
		N_layer = N_SOILLAYERS;
	else
		N_layer = GWlayer;

	/* A. in case firsttimeflag=0 Note the initial value for calculating the difference.  */
	if (firsttime_flag == 0)
	{ 	
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				for (layer = 0; layer < N_SOILLAYERS; layer++) content_soil_pre[dm][layer] = soilInfo->content_soil[dm][layer];
				content_NORMcf_pre[dm] = soilInfo->content_NORMcf[dm];
				content_CAPILcf_pre[dm] = soilInfo->content_CAPILcf[dm];

			}

	}
	/* B. in case firsttimeflag=1 calculate changes */
	else
	{ 
		/* B1. if change of BOUND and DISSOLV is calculated outside, the content is the sum of BOUND and DISSOLV */
		if (partlyORtotal_flag == 0)
		{
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				for (layer = 0; layer < N_SOILLAYERS; layer++) soilInfo->content_soil[dm][layer] = soilInfo->contentBOUND_soil[dm][layer] + soilInfo->contentDISSOLV_soil[dm][layer];

				if (GWlayer != DATA_GAP)
				{
					soilInfo->content_NORMcf[dm] = soilInfo->contentBOUND_NORMcf[dm] + soilInfo->contentDISSOLV_NORMcf[dm];
					soilInfo->content_CAPILcf[dm] = soilInfo->contentBOUND_CAPILcf[dm] + soilInfo->contentDISSOLV_CAPILcf[dm];
					soilInfo->content_NORMgw[dm] = soilInfo->contentBOUND_NORMgw[dm] + soilInfo->contentDISSOLV_NORMgw[dm];
					soilInfo->content_CAPILgw[dm] = soilInfo->contentBOUND_CAPILgw[dm] + soilInfo->contentDISSOLV_CAPILgw[dm];
					soilInfo->content_SATgw[dm] = soilInfo->contentBOUND_SATgw[dm] + soilInfo->contentDISSOLV_SATgw[dm];
				}
			}
		}
		/* B2. if change of total content is calculated outsuide, changes of BOUND and DISSOLV are calculated proportionally. Exception: GWlayer and below */
		else
		{
			/* B2.1 above GW layer: proportionally change (exception: CFlayer, must be calculated separately CAPILcf and NORMcf */
			for (layer = 0; layer < N_layer; layer++)
			{
				if (layer != CFlayer)
				{
					for (dm = 0; dm < N_DISSOLVMATER; dm++)
					{
						soilInfo->contentBOUND_soil[dm][layer] += (soilInfo->content_soil[dm][layer] - content_soil_pre[dm][layer]) * (1 - soilInfo->dissolv_prop[dm]);
						soilInfo->contentDISSOLV_soil[dm][layer] += (soilInfo->content_soil[dm][layer] - content_soil_pre[dm][layer]) * soilInfo->dissolv_prop[dm];

						/* avoid negative pools */
						if (soilInfo->contentBOUND_soil[dm][layer] < 0)
						{
							soilInfo->contentDISSOLV_soil[dm][layer] += soilInfo->contentBOUND_soil[dm][layer];
							soilInfo->contentBOUND_soil[dm][layer] = 0;
						}
						if (soilInfo->contentDISSOLV_soil[dm][layer] < 0)
						{
							soilInfo->contentBOUND_soil[dm][layer] += soilInfo->contentDISSOLV_soil[dm][layer];
							soilInfo->contentDISSOLV_soil[dm][layer] = 0;
						}
						if (soilInfo->contentBOUND_soil[dm][layer] < 0)
						{
							if (fabs(soilInfo->contentBOUND_soil[dm][layer]) > CRIT_PREC)
							{
								printf("\n");
								printf("ERROR: negative BOUND soil content problem (calc_DISSOLVandBOUND.c)\n");
								errorCode = 1;
							}
							else
								soilInfo->contentBOUND_soil[dm][layer] = 0;

						}

						if (soilInfo->contentDISSOLV_soil[dm][layer] < 0)
						{
							if (fabs(soilInfo->contentDISSOLV_soil[dm][layer]) > CRIT_PREC)
							{
								printf("\n");
								printf("ERROR: negative DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
								errorCode = 1;
							}
							else
								soilInfo->contentDISSOLV_soil[dm][layer] = 0;

						}
					}
				}
				else
				{
					for (dm = 0; dm < N_DISSOLVMATER; dm++)
					{
						/*---------------*/
						/* NORMcf */
						if (sprop->dz_NORMcf)
						{
							soilInfo->contentBOUND_NORMcf[dm] += (soilInfo->content_NORMcf[dm] - content_NORMcf_pre[dm]) * (1 - soilInfo->dissolv_prop[dm]);
							soilInfo->contentDISSOLV_NORMcf[dm] += (soilInfo->content_NORMcf[dm] - content_NORMcf_pre[dm]) * soilInfo->dissolv_prop[dm];

							/* avoid negative pools */
							if (soilInfo->contentBOUND_NORMcf[dm] < 0)
							{
								soilInfo->contentDISSOLV_NORMcf[dm] += soilInfo->contentBOUND_NORMcf[dm];
								soilInfo->contentDISSOLV_soil[dm][CFlayer] += soilInfo->contentBOUND_NORMcf[dm];
								soilInfo->contentBOUND_soil[dm][CFlayer] -= soilInfo->contentBOUND_NORMcf[dm];
								soilInfo->contentBOUND_NORMcf[dm] = 0;
							}
							if (soilInfo->contentDISSOLV_NORMcf[dm] < 0)
							{
								soilInfo->contentBOUND_NORMcf[dm] += soilInfo->contentDISSOLV_NORMcf[dm];
								soilInfo->contentBOUND_soil[dm][CFlayer] += soilInfo->contentDISSOLV_NORMcf[dm];
								soilInfo->contentDISSOLV_soil[dm][CFlayer] -= soilInfo->contentDISSOLV_NORMcf[dm];
								soilInfo->contentDISSOLV_NORMcf[dm] = 0;
							}
							if (soilInfo->contentBOUND_NORMcf[dm] < 0 || soilInfo->contentDISSOLV_NORMcf[dm] < 0)
							{
								printf("\n");
								printf("ERROR: negative BOUND or DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
								errorCode = 1;
							}
						}

						/*---------------*/
						/* CAPILcf */
						if (sprop->dz_CAPILcf)
						{
							soilInfo->contentBOUND_CAPILcf[dm] += (soilInfo->content_CAPILcf[dm] - content_CAPILcf_pre[dm]) * (1 - soilInfo->dissolv_prop[dm]);
							soilInfo->contentDISSOLV_CAPILcf[dm] += (soilInfo->content_CAPILcf[dm] - content_CAPILcf_pre[dm]) * soilInfo->dissolv_prop[dm];

							/* avoid negative pools */
							if (soilInfo->contentBOUND_CAPILcf[dm] < 0)
							{
								soilInfo->contentDISSOLV_CAPILcf[dm] += soilInfo->contentBOUND_CAPILcf[dm];
								soilInfo->contentDISSOLV_soil[dm][CFlayer] += soilInfo->contentBOUND_CAPILcf[dm];
								soilInfo->contentBOUND_soil[dm][CFlayer] -= soilInfo->contentBOUND_CAPILcf[dm];
								soilInfo->contentBOUND_CAPILcf[dm] = 0;
							}
							if (soilInfo->contentDISSOLV_CAPILcf[dm] < 0)
							{
								soilInfo->contentBOUND_CAPILcf[dm] += soilInfo->contentDISSOLV_CAPILcf[dm];
								soilInfo->contentBOUND_soil[dm][CFlayer] += soilInfo->contentDISSOLV_CAPILcf[dm];
								soilInfo->contentDISSOLV_soil[dm][CFlayer] -= soilInfo->contentDISSOLV_CAPILcf[dm];
								soilInfo->contentDISSOLV_CAPILcf[dm] = 0;
							}
							if (soilInfo->contentBOUND_CAPILcf[dm] < 0 || soilInfo->contentDISSOLV_CAPILcf[dm] < 0)
							{
								printf("\n");
								printf("ERROR: negative BOUND or DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
								errorCode = 1;
							}
						}

						soilInfo->contentBOUND_soil[dm][layer] = soilInfo->contentBOUND_CAPILcf[dm] + soilInfo->contentBOUND_NORMcf[dm];
						soilInfo->contentDISSOLV_soil[dm][layer] = soilInfo->contentDISSOLV_CAPILcf[dm] + soilInfo->contentDISSOLV_NORMcf[dm];

						soilInfo->content_soil[dm][layer] = soilInfo->contentBOUND_soil[dm][layer] + soilInfo->contentDISSOLV_soil[dm][layer];

					}
				} /* endelse CFlayer*/
			} /* endfor above GWLayers */

			/* B2.2 in GW layer: calculation from DISSOLV and BOUND*/	
			if (GWlayer != DATA_GAP)
			{ 
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
						
					soilInfo->content_NORMgw[dm] = soilInfo->contentBOUND_NORMgw[dm] + soilInfo->contentDISSOLV_NORMgw[dm];
					soilInfo->content_CAPILgw[dm] = soilInfo->contentBOUND_CAPILgw[dm] + soilInfo->contentDISSOLV_CAPILgw[dm];
					soilInfo->content_SATgw[dm] = soilInfo->contentBOUND_SATgw[dm] + soilInfo->contentDISSOLV_SATgw[dm];
					
					soilInfo->contentBOUND_soil[dm][GWlayer] = soilInfo->contentBOUND_SATgw[dm] + soilInfo->contentBOUND_CAPILgw[dm] + soilInfo->contentBOUND_NORMgw[dm];
					soilInfo->contentDISSOLV_soil[dm][GWlayer] = soilInfo->contentDISSOLV_SATgw[dm] + soilInfo->contentDISSOLV_CAPILgw[dm] + soilInfo->contentDISSOLV_NORMgw[dm];

					soilInfo->content_soil[dm][GWlayer] = soilInfo->contentBOUND_soil[dm][GWlayer] + soilInfo->contentDISSOLV_soil[dm][GWlayer];
				
				}
			}

			/* B2.3 below GWlayer: calculation from DISSOLV and BOUND*/
			for (layer = N_layer+1; layer < N_SOILLAYERS; layer++)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++) soilInfo->content_soil[dm][layer] = soilInfo->contentBOUND_soil[dm][layer] + soilInfo->contentDISSOLV_soil[dm][layer];
			
			}
		}
	}

	
		
	/* control for negative pool */
	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{
			if (!errorCode && fabs(soilInfo->content_soil[dm][layer] - soilInfo->contentBOUND_soil[dm][layer] - soilInfo->contentDISSOLV_soil[dm][layer]) > CRIT_PREC)
			{
 				printf("\n");
				printf("ERROR: BOUND+DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
				errorCode = 1;
			}

			if (!errorCode && soilInfo->contentDISSOLV_soil[dm][layer] < 0)
			{
				printf("\n");
				printf("ERROR: BOUND+DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
				errorCode = 1;
			}
			if (!errorCode && soilInfo->contentBOUND_soil[dm][layer] < 0)
			{
				printf("\n");
				printf("ERROR: BOUND+BOUND soil content problem (calc_DISSOLVandBOUND.c)\n");
				errorCode = 1;
			}

	
			if (sprop->CFlayer == layer && sprop->CFlayer != sprop->GWlayer)
			{
				if (!errorCode && fabs(soilInfo->contentDISSOLV_soil[dm][layer] - soilInfo->contentDISSOLV_CAPILcf[dm] - soilInfo->contentDISSOLV_NORMcf[dm]) > CRIT_PREC)
				{
					printf("\n");
					printf("ERROR: BOUND+DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
					errorCode = 1;
				}

				if (!errorCode && (soilInfo->contentDISSOLV_CAPILcf[dm] < 0 || soilInfo->contentDISSOLV_NORMcf[dm] < 0))
				{
					printf("\n");
					printf("ERROR: BOUND+DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
					errorCode = 1;
				}

				if (!errorCode && (soilInfo->contentBOUND_CAPILcf[dm] < 0 || soilInfo->contentBOUND_NORMcf[dm] < 0))
				{
					printf("\n");
					printf("ERROR: BOUND+BOUND soil content problem (calc_DISSOLVandBOUND.c)\n");
					errorCode = 1;
				}

			}

			if (sprop->GWlayer == layer)
			{
				if (!errorCode && fabs(soilInfo->contentDISSOLV_soil[dm][layer] - soilInfo->contentDISSOLV_SATgw[dm] - soilInfo->contentDISSOLV_CAPILgw[dm] - soilInfo->contentDISSOLV_NORMgw[dm]) > CRIT_PREC)
				{
					printf("\n");
					printf("ERROR: BOUND+DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
					errorCode = 1;
				}

				if (!errorCode && (soilInfo->contentDISSOLV_SATgw[dm] < 0 || soilInfo->contentDISSOLV_CAPILgw[dm] < 0 || soilInfo->contentDISSOLV_NORMgw[dm] < 0))
				{
					printf("\n");
					printf("ERROR: BOUND+DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
					errorCode = 1;
				}

				if (!errorCode && (soilInfo->contentBOUND_SATgw[dm] < 0 || soilInfo->contentBOUND_CAPILgw[dm] < 0 || soilInfo->contentBOUND_NORMgw[dm] < 0))
				{
					printf("\n");
					printf("ERROR: BOUND+BOUND soil content problem (calc_DISSOLVandBOUND.c)\n");
					errorCode = 1;
				}
			}
		}
		if (!errorCode && fabs(soilInfo->content_NORMcf[dm] - soilInfo->contentBOUND_NORMcf[dm] - soilInfo->contentDISSOLV_NORMcf[dm]) > CRIT_PREC)
		{
			printf("\n");
			printf("ERROR: BOUND+DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
			errorCode = 1;
		}
		if (!errorCode && fabs(soilInfo->content_CAPILcf[dm] - soilInfo->contentBOUND_CAPILcf[dm] - soilInfo->contentDISSOLV_CAPILcf[dm]) > CRIT_PREC)
		{
			printf("\n");
			printf("ERROR: BOUND+DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
			errorCode = 1;
		}
		if (!errorCode && fabs(soilInfo->content_NORMgw[dm] - soilInfo->contentBOUND_NORMgw[dm] - soilInfo->contentDISSOLV_NORMgw[dm]) > CRIT_PREC)
		{
			printf("\n");
			printf("ERROR: BOUND+DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
			errorCode = 1;
		}
		if (!errorCode && fabs(soilInfo->content_CAPILgw[dm] - soilInfo->contentBOUND_CAPILgw[dm] - soilInfo->contentDISSOLV_CAPILgw[dm]) > CRIT_PREC)
		{
			printf("\n");
			printf("ERROR: BOUND+DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
			errorCode = 1;
		}
		if (!errorCode && fabs(soilInfo->content_SATgw[dm] - soilInfo->contentBOUND_SATgw[dm] - soilInfo->contentDISSOLV_SATgw[dm]) > CRIT_PREC)
		{
			printf("\n");
			printf("ERROR: BOUND+DISSOLV soil content problem (calc_DISSOLVandBOUND.c)\n");
			errorCode = 1;
		}



	}

	return(errorCode);
}

