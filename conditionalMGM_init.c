/* 
conditionalMGM_init.c
read conditional management information for pointbgc simulation

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
#include "bgc_struct.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"


int conditionalMGM_init(file init, control_struct* ctrl, irrigating_struct* IRG, mowing_struct* MOW)
{
	
	char key[] = "CONDITIONAL_MANAGEMENT_STRATEGIES";
	char keyword[STRINGSIZE];

	int errorCode=0;

	/********************************************************************
	**                                                                 **
	** Begin reading initialization file block starting with keyword:  **
	** irrigating                                                       ** 
	**                                                                 **
	********************************************************************/
	
	/* scan for the conditional_management keyword, exit if not next */
	if (!errorCode && scan_value(init, keyword, 's'))
	{
		printf("ERROR reading keyword for control data\n");
		errorCode=215;
	}
	if (!errorCode && strcmp(keyword, key))
	{
		printf("Expecting keyword --> %s in file %s\n",key,init.name);
		errorCode=215;
	}
	
	/* conditional mowing parameters */	
	if (!errorCode && scan_value(init, &MOW->condMOW_flag, 'i'))
	{
		printf("ERROR reading condMOW_flag\n");
		errorCode=21501;
	}

	/* control of condMOW and MOW_flag */
	if (!errorCode && MOW->condMOW_flag && MOW->MOW_num)
	{
		MOW->MOW_num = 0;
		ctrl->condMOWerr_flag = 1;

	}

	if (!errorCode && scan_value(init, &MOW->fixLAIbef_condMOW, 'd'))
	{
		printf("ERROR reading fixLAIbef_condMOW\n");
		errorCode=21502;
	}

	if (!errorCode && scan_value(init, &MOW->fixLAIaft_condMOW, 'd'))
	{
		printf("ERROR reading fixLAIaft_condMOW\n");
		errorCode=21503;
	}

	if (!errorCode && scan_value(init, &MOW->transpCOEFF_condMOW, 'd'))
	{
		printf("ERROR reading transpCOEFF_condMOW\n");
		errorCode=21504;
	}



	/* conditional irrigating parameters */	
	if (!errorCode && scan_value(init, &IRG->condIRG_flag, 'i'))
	{
		printf("ERROR reading condIRG_flag\n");
		errorCode=21505;
	}

	/* control of condIRG and IRG_flag */
	if (!errorCode && IRG->condIRG_flag && IRG->IRG_num)
	{
		IRG->IRG_num = 0;
		ctrl->condIRGerr_flag = 1;
	}


	if (!errorCode && scan_value(init, &IRG->nLayer_condIRG, 'd'))
	{
		printf("ERROR reading nLayer_condIRG\n");
		errorCode=21506;
	}

	
	if (!errorCode && scan_value(init, &IRG->startPoint_condIRG, 'd'))
	{
		printf("ERROR reading startPoint_condIRG\n");
		errorCode=21507;
	}

	
	if (!errorCode && scan_value(init, &IRG->aftVWCratio_condIRG, 'd'))
	{
		printf("ERROR reading aftVWCratio_condIRG\n");
		errorCode=21508;
	}

	if (!errorCode && scan_value(init, &IRG->maxAMOUNT_condIRG, 'd'))
	{
		printf("ERROR reading maxAMOUNT_condIRG\n");
		errorCode=21509;
	}

	/* control */
	if (!errorCode && IRG->condIRG_flag)
	{

		if (IRG->condIRG_flag &&  IRG->condIRG_flag > 0 && IRG->maxAMOUNT_condIRG == 0)
		{
			printf("ERROR in conditional irrigating: maxAMOUNT_condIRG must greater than 0\n");
			errorCode=2150501;
		}

		if (IRG->condIRG_flag == 2 && IRG->startPoint_condIRG > IRG->aftVWCratio_condIRG)
		{
			printf("ERROR in VWCratio values of conditional irrigating: befVWC_ratio must less then aftVWC_ratio\n");
			errorCode=2150502;
		}

	}


	
	return (errorCode);
}
