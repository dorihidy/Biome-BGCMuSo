/* 
scc_init.c
Initialize the scalar climate change parameters for bgc simulation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group (NTSG)
School of Forestry, University of Montana
Missoula, MT 59812
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ini.h"
#include "bgc_struct.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"

int scc_init(file init, climchange_struct* scc)
{
	int errorCode=0;
	char key1[] = "CLIM_CHANGE";
	char keyword[STRINGSIZE];

	/********************************************************************
	**                                                                 **
	** Begin reading initialization file block starting with keyword:  **
	** CLIM_CHANGE                                                     ** 
	**                                                                 **
	********************************************************************/

	/* scan for the climate change block keyword, exit if not next */
	if (!errorCode && scan_value(init, keyword, 's'))
	{
		printf("ERROR reading keyword, scc_init.c\n");
		errorCode=214;
	}
	if (!errorCode && strcmp(keyword,key1))
	{
		printf("Expecting keyword --> %s in %s\n",key1,init.name);
		errorCode=214;
	}

	/* begin reading climate change data */
	if (!errorCode && scan_value(init, &scc->s_Tmax, 'd'))
	{
		printf("ERROR reading scalar for Tmax, scc_init.c\n");
		errorCode=21401;
	}
	if (!errorCode && scan_value(init, &scc->s_Tmin, 'd'))
	{
		printf("ERROR reading scalar for Tmin, scc_init.c\n");
		errorCode=21402;
	}
	if (!errorCode && scan_value(init, &scc->s_prcp, 'd'))
	{
		printf("ERROR reading scalar for prcp, scc_init.c\n");
		errorCode=21403;
	}
	if (!errorCode && scan_value(init, &scc->s_vpd, 'd'))
	{
		printf("ERROR reading scalar for vpd, scc_init.c\n");
		errorCode=21404;
	}
	if (!errorCode && scan_value(init, &scc->s_swavgfd, 'd'))
	{
		printf("ERROR reading scalar for swavgfd, scc_init.c\n");
		errorCode=21405;
	}
	
	/* some error checking on scalar climate change values */
	if (scc->s_prcp < 0.0)
	{
		printf("ERROR in scc_init.c: prcp scalar must be positive\n");
		errorCode=21403;
	}
	if (scc->s_vpd < 0.0)
	{
		printf("ERROR in scc_init.c: vpd scalar must be positive\n");
		errorCode=21404;
	}
	if (scc->s_swavgfd < 0.0)
	{
		printf("ERROR in scc_init.c: swavgfd scalar must be positive\n");
		errorCode=21405;
	}

	return (errorCode);
}
