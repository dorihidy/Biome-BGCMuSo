/* 
check_soilcontent.c
transformation and checking of soil content data from normal pools (ns, cs) to multi array
content2pool_flag = 0: NH4, NO3, DOC, DON - > content_array
content2pool_flag = 1: content_array -> NH4, NO3, DOC, DON

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

int check_soilcontent(int layerFlag, int content2pool_flag, const soilprop_struct* sprop, cstate_struct* cs, nstate_struct* ns, soilInfo_struct* soilInfo)
{
	int errorCode, layer, layerS, layerE, dm;
	errorCode = 0;

	/* calculation for all layers (layerFlag == -1) or only the given layer (layerFlag > -1) */
	if (layerFlag == -1)
	{
		layerS = 0;
		layerE = N_SOILLAYERS;
	}
	else
	{
		layerS = layerFlag;
		layerE = layerFlag + 1;
	}

	for (layer = layerS; layer < layerE; layer++)
	{
		if (content2pool_flag == 0)
		{
			soilInfo->content_soil[0][layer] = ns->NH4[layer];
			soilInfo->content_soil[1][layer] = ns->NO3[layer];
			soilInfo->content_soil[2][layer] = ns->soil1n[layer];
			soilInfo->content_soil[3][layer] = ns->soil2n[layer];
			soilInfo->content_soil[4][layer] = ns->soil3n[layer];
			soilInfo->content_soil[5][layer] = ns->soil4n[layer];
			soilInfo->content_soil[6][layer] = cs->soil1c[layer];
			soilInfo->content_soil[7][layer] = cs->soil2c[layer];
			soilInfo->content_soil[8][layer] = cs->soil3c[layer];
			soilInfo->content_soil[9][layer] = cs->soil4c[layer];
		}
		else
		{
			ns->NH4[layer] = soilInfo->content_soil[0][layer];
			ns->NO3[layer] = soilInfo->content_soil[1][layer];
			ns->soil1n[layer] = soilInfo->content_soil[2][layer];
			ns->soil2n[layer] = soilInfo->content_soil[3][layer];
			ns->soil3n[layer] = soilInfo->content_soil[4][layer];
			ns->soil4n[layer] = soilInfo->content_soil[5][layer];
			cs->soil1c[layer] = soilInfo->content_soil[6][layer];
			cs->soil2c[layer] = soilInfo->content_soil[7][layer];
			cs->soil3c[layer] = soilInfo->content_soil[8][layer];
			cs->soil4c[layer] = soilInfo->content_soil[9][layer];
		}

		
		/* control for negative pool */
		for (dm = 0; dm < N_DISSOLVMATER; dm++)
		{
			if (!errorCode && soilInfo->content_soil[dm][layer] < 0.0)
			{
				if (fabs(soilInfo->content_soil[dm][layer]) > CRIT_PREC && !errorCode)
				{
					printf("\n");
					printf("ERROR: negative nitrogen pool (check_soilcontent.c)\n");
					errorCode = 1;
				}
				else
					soilInfo->content_soil[dm][layer] = 0;
			}

			if (!errorCode && soilInfo->content_soil[dm][layer] && soilInfo->content_soil[dm][layer] / soilInfo->content_soil[dm][layer] != 1)
			{
				printf("\n");
				printf("ERROR: invalid nitrogen pool (check_soilcontent.c)\n");
				errorCode = 1;
			}

		}


		/* control for unsat CF-layer is it is avaialbe  */
		if (sprop->GWlayer != DATA_GAP && layer == (int)sprop->CFlayer && sprop->dz_NORMcf && content2pool_flag == 1)
		{
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
			
				if (soilInfo->content_NORMcf[dm] < 0)
				{
					if (fabs(soilInfo->content_NORMcf[dm]) < CRIT_PREC)
					{
						soilInfo->content_soil[dm][layer] -= soilInfo->content_NORMcf[dm];
						soilInfo->content_NORMcf[dm] = 0;
					}
					else
					{
						printf("\n");
						printf("ERROR: negative pool in zone_NORM (check_soilcontent.c)\n");
						errorCode = 1;
					}

				}
			
				
				if (soilInfo->content_CAPILcf[dm] < 0)
				{
					if (fabs(soilInfo->content_CAPILcf[dm]) < CRIT_PREC)
					{
						soilInfo->content_soil[dm][layer] -= soilInfo->content_CAPILcf[dm];
						soilInfo->content_CAPILcf[dm] = 0;
					}
					else
					{ 
						printf("\n");
						printf("ERROR: negative pool in zone_CAPIL (check_soilcontent.c)\n");
						errorCode = 1;
					}
				}

		
				if (!errorCode && soilInfo->content_soil[dm][layer] - (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm]) > CRIT_PREC * 10)
				{
					printf("\n");
					printf("ERROR: sum of zone_NORM, zone_CAPIL  is not equal to content of unsat CFlayer (check_soilcontent.c)\n");
					errorCode = 1;
				}

			}
		}

		/* control for GWlayer */
		if (sprop->GWlayer != DATA_GAP  && layer == (int)sprop->GWlayer && sprop->dz_CAPILgw && content2pool_flag == 1)
		{
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				
				if (soilInfo->content_NORMgw[dm] < 0)
				{
					if (fabs(soilInfo->content_NORMgw[dm]) < CRIT_PREC)
					{
						soilInfo->content_soil[dm][layer] -= soilInfo->content_NORMgw[dm];
						soilInfo->content_NORMgw[dm] = 0;
					}
					else
					{ 
						printf("\n");
						printf("ERROR: negative pool in zone_NORM (check_soilcontent.c)\n");
						errorCode = 1;
					}
				}

			
			
				if (soilInfo->content_CAPILgw[dm] < 0)
				{
					if (fabs(soilInfo->content_CAPILgw[dm]) < CRIT_PREC)
					{
						soilInfo->content_soil[dm][layer] -= soilInfo->content_CAPILgw[dm];
						soilInfo->content_CAPILgw[dm] = 0;
					}
					else
					{ 
						printf("\n");
						printf("ERROR: negative pool in zone_CAPIL (check_soilcontent.c)\n");
						errorCode = 1;
					}
				}

	
				
				if (soilInfo->content_SATgw[dm] < 0)
				{
					printf("\n");
					printf("ERROR: negative pool in zone_SAT (check_soilcontent.c)\n");
					errorCode = 1;
				}

				if (!errorCode && soilInfo->content_soil[dm][layer] - (soilInfo->content_NORMgw[dm] + soilInfo->content_CAPILgw[dm] + soilInfo->content_SATgw[dm]) > CRIT_PREC*100000)
				{
					printf("\n");
					printf("ERROR: sum of zone_NORM, zone_CAPIL and zone_SAT is not equal to content of GWlayer (check_soilcontent.c)\n");
					errorCode = 1;
				}

		
			

			}
		}

		

	
	}



	return(errorCode);
}

