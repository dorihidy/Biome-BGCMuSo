/* 
mgm_init.c
read mgm file for pointbgc simulation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2022, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
20.03.01 Galina Churkina added variable "sum" substituting  t1+t2+t3 in IF statement,
which gave an error.
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
#include "bgc_func.h"
#include "bgc_constants.h"



int mgm_init(file init, control_struct* ctrl,  
	         fertilizing_struct* FRZ, grazing_struct* GRZ, harvesting_struct* HRV, mowing_struct* MOW, planting_struct* PLT, ploughing_struct* PLG, 
			 thinning_struct* THN, irrigating_struct* IRG, mulching_struct* MUL, CWDextract_struct* CWE)
{


	int errorCode=0;
	int dofilecloseMANAGEMENT = 1;
	file mgm_file;
	char key1[] = "MANAGEMENT_FILE";
	char key2[] = "none";
	char keyword[STRINGSIZE];
	char header[STRINGSIZE];

	int mgm, PLTyday, HRVyday, GRZstart_yday,GRZend_yday,doy,leap;
	int* mondays=0;
	int* enddays=0;


	/********************************************************************
	**                                                                 **
	** Begin reading initialization file block starting with keyword:  **
	** MANAGEMENT_FILE                                                       ** 
	**                                                                 **
	********************************************************************/
	
	/* initialization */
	PLT->mgmdPLT = -1;
	THN->THN_num = 0;
	MOW->MOW_num = 0;
	HRV->HRV_num = 0;
	PLG->PLG_num = 0;
	GRZ->GRZ_num = 0;
	FRZ->FRZ_num = 0;
	IRG->IRG_num = 0;
	MUL->MUL_num = 0;
	CWE->CWE_num = 0;
	GRZ->trampleff_act	= DATA_GAP;

	
		/* scan for the MANAGEMENT file keyword, exit if not next */
		if (!errorCode && scan_value(init, keyword, 's'))
		{
			printf("ERROR reading keyword for control data\n");
			errorCode=210;
			dofilecloseMANAGEMENT = 0;
		}
		if (!errorCode && strcmp(keyword, key1))
		{
			printf("Expecting keyword --> %s in file %s\n",key1,init.name);
			errorCode=210;
			dofilecloseMANAGEMENT = 0;
		}
	/* open simple MANAGEMENT file  */
	if (!errorCode && scan_open(init,&mgm_file,'r',0)) 
	{
		
			if (!errorCode && strcmp(mgm_file.name, key2))
			{
				printf("ERROR opening mgm_file file from INI file, mgm_init.c\n");
				errorCode=210;
			}
			dofilecloseMANAGEMENT = 0;

	}
	
	/* if management data file is open: read management data */
	if (dofilecloseMANAGEMENT)
	{
		/* first scan sprop header to ensure proper *.sprop format */
		if (!errorCode && scan_value(mgm_file, header, 's'))
		{
			printf("ERROR reading header, mgm_init.c\n");
			errorCode=210;
		}
		
		/* allocate space for enddays and mondays */	
		if (!errorCode) 
		{
			enddays = (int*) malloc(nMONTHS_OF_YEAR * sizeof(int));
			if (!enddays)
			{
				printf("ERROR allocating for enddays in bgc.c\n");
				errorCode=2100;
			}
		}

		if (!errorCode) 
		{
			mondays = (int*) malloc(nMONTHS_OF_YEAR * sizeof(int));
			if (!mondays)
			{
				printf("ERROR allocating for enddays in bgc.c\n");
				errorCode=2100;
			}
		}

		/* -------------------------------------------------------------------------*/
		/* MANAGEMENT SECTION  */
	
		/* read the planting information */
		if (!errorCode && planting_init(mgm_file, ctrl, PLT))
		{
			printf("ERROR in call to planting_init.c from mgm_init.c... Exiting\n");
			errorCode=2101;
		}

		/* read the thinning information */
		if (!errorCode && thinning_init(mgm_file, ctrl, THN))
		{
			printf("ERROR in call to thinning_init.c from mgm_init.c... Exiting\n");
			errorCode=2102;
		}

		/* read the mowing  information */
		if (!errorCode && mowing_init(mgm_file, ctrl, MOW))
		{
			printf("ERROR in call to mowing_init.c from mgm_init.c... Exiting\n");
			errorCode=2103;
		}

		/* read the grazing information */
		if (!errorCode && grazing_init(mgm_file, ctrl, GRZ))
		{
			printf("ERROR in call to grazing_init.c from mgm_init.c... Exiting\n");
			errorCode=2104;
		}

		/* read the harvesting information */
		if (!errorCode && harvesting_init(mgm_file, ctrl, HRV))
		{
			printf("ERROR in call to harvesting_init.c from mgm_init.c... Exiting\n");
			errorCode=2105;
		}

		/* read the harvesting information */
		if (!errorCode && ploughing_init(mgm_file, ctrl, PLG))
		{
			printf("ERROR in call to ploughing_init.c from mgm_init.c... Exiting\n");
			errorCode=2106;
		}

		/* read the fertilizing  information */
		if (!errorCode && fertilizing_init(mgm_file, ctrl, FRZ))
		{
			printf("ERROR in call to fertilizing_init.c from mgm_init.c... Exiting\n");
			errorCode=2107;
		}

		/* read the irrigating information */
		if (!errorCode && irrigating_init(mgm_file, ctrl, IRG))
		{
			printf("ERROR in call to irrigating_init.c from mgm_init.c... Exiting\n");
			errorCode=2108;
		}

		/* read the mulching information */
		if (!errorCode && mulching_init(mgm_file, ctrl, MUL))
		{
			printf("ERROR in call to mulching_init.c from mgm_init.c... Exiting\n");
			errorCode=2109;
		}

		/* read the CWDextract information */
		if (!errorCode && CWDextract_init(mgm_file, ctrl, CWE))
		{
			printf("ERROR in call to CWDextract_init.c from mgm_init.c... Exiting\n");
			errorCode=2110;
		}

		




		fclose(mgm_file.ptr);

		/* CONTROL OF MANAGEMENT DATA */
		
		/* planting and harvest setting*/
		if (!errorCode && (PLT->PLT_num != HRV->HRV_num))
		{
			if (PLT->PLT_num == HRV->HRV_num + 1)
			{
				if (ctrl->onscreen) printf("WARNING: last planting is inactive (because of the lack of harvest)\n");
				PLT->PLT_num -= 1;
			}
			else
			{
				printf("ERROR in management data: number planting and harvest action should be equal\n");
				errorCode=2100001;
			}
		}

		/* planting and harvest date */
		if (!errorCode && PLT->PLT_num)
		{
			for (mgm = 0; mgm < PLT->PLT_num; mgm++)
			{
				if (!errorCode && leapControl(PLT->PLTyear_array[mgm], enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from mgm_init.c\n");
					errorCode=2100002;
				}
				doy = date_to_doy(mondays, PLT->PLTmonth_array[mgm], PLT->PLTday_array[mgm]);
				PLTyday = PLT->PLTyear_array[mgm] * nDAYS_OF_YEAR + doy;

				if (leap == 1 && PLT->PLTmonth_array[mgm] == 12 && PLT->PLTday_array[mgm] == 31)
				{
					printf("ERROR in planting date in mgm_init.c: data from 31 December in a leap year is found in planting file\n");
					printf("Please read the manual and modify the input data\n");
					errorCode=2100003;
				}

				if (!errorCode && leapControl(HRV->HRVyear_array[mgm], enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from mgm_init.c\n");
					errorCode=2100004;
				}
				doy = date_to_doy(mondays, HRV->HRVmonth_array[mgm], HRV->HRVday_array[mgm]);
				HRVyday = HRV->HRVyear_array[mgm] * nDAYS_OF_YEAR + doy;
				
				if (leap == 1 && HRV->HRVmonth_array[mgm] == 12 && HRV->HRVday_array[mgm] == 31)
				{
					printf("ERROR in harvesting date in mgm_init.c: data from 31 December in a leap year is found in harvesting file\n");
					printf("Please read the manual and modify the input data\n");
					errorCode=2100005;
				}

				if (HRVyday <= PLTyday)
				{
					printf("ERROR in management data: PLANTING must be before HARVESTING date\n");
					errorCode=2100006;
				}


			}
		}


		/* grazing date */
		if (!errorCode && GRZ->GRZ_num)
		{
			for (mgm = 0; mgm < GRZ->GRZ_num; mgm++)
			{
				if (!errorCode && leapControl(GRZ->GRZstart_year_array[mgm], enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from mgm_init.c\n");
					errorCode=2100007;
				}
				GRZstart_yday = GRZ->GRZstart_year_array[mgm] * nDAYS_OF_YEAR + date_to_doy(mondays, GRZ->GRZstart_month_array[mgm], GRZ->GRZstart_day_array[mgm]);
				
				if (leap == 1 && GRZ->GRZstart_month_array[mgm] == 12 && GRZ->GRZstart_day_array[mgm] == 31)
				{
					printf("ERROR in grazing start date in mgm_init.c: data from 31 December in a leap year is found in grazing file\n");
					printf("Please read the manual and modify the input data\n");
					errorCode=2100008;
				}

				if (!errorCode && leapControl(GRZ->GRZend_year_array[mgm], enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from mgm_init.c\n");
					errorCode=2100009;
				}
				GRZend_yday   = GRZ->GRZend_year_array[mgm] * nDAYS_OF_YEAR + date_to_doy(mondays, GRZ->GRZend_month_array[mgm], GRZ->GRZend_day_array[mgm]);
				
				if (leap == 1 && GRZ->GRZend_month_array[mgm] == 12 && GRZ->GRZend_day_array[mgm] == 31)
				{
					printf("ERROR in grazing end date in mgm_init.c: data from 31 December in a leap year is found in grazing file\n");
					printf("Please read the manual and modify the input data\n");
					errorCode=2100010;
				}

				if (GRZend_yday < GRZstart_yday)
				{
					printf("ERROR in management date: start of grazing must be before end of grazing\n");
					errorCode=2100011;
				}
			}
		}

		/* leap year */
		if (!errorCode && THN->THN_num)
		{
			for (mgm = 0; mgm < THN->THN_num; mgm++)
			{
				if (!errorCode && leapControl(THN->THNyear_array[mgm], enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from mgm_init.c\n");
					errorCode=2100012;
				}
				if (leap == 1 && THN->THNmonth_array[mgm] == 12 && THN->THNday_array[mgm] == 31)
				{
					printf("ERROR in thinning date in mgm_init.c: data from 31 December in a leap year is found in thinning file\n");
					printf("Please read the manual and modify the input data\n");
					errorCode=2100012;
				}
			}
		}

		if (!errorCode && MOW->MOW_num)
		{
			for (mgm = 0; mgm < MOW->MOW_num; mgm++)
			{
				if (!errorCode && leapControl(MOW->MOWyear_array[mgm], enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from mgm_init.c\n");
					errorCode=2100013;
				}
				if (leap == 1 && MOW->MOWmonth_array[mgm] == 12 && MOW->MOWday_array[mgm] == 31)
				{
					printf("ERROR in mowing date in mgm_init.c: data from 31 December in a leap year is found in mowing file\n");
					printf("Please read the manual and modify the input data\n");
					errorCode=2100013;
				}
			}
		}

		if (!errorCode && PLG->PLG_num)
		{
			for (mgm = 0; mgm < PLG->PLG_num; mgm++)
			{
				if (!errorCode && leapControl(PLG->PLGyear_array[mgm], enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from mgm_init.c\n");
					errorCode=2100013;
				}
				if (leap == 1 && PLG->PLGmonth_array[mgm] == 12 && PLG->PLGday_array[mgm] == 31)
				{
					printf("ERROR in plouging date in mgm_init.c: data from 31 December in a leap year is found in plouging file\n");
					printf("Please read the manual and modify the input data\n");
					errorCode=2100013;
				}
			}
		}

		if (!errorCode && FRZ->FRZ_num)
		{
			for (mgm = 0; mgm < FRZ->FRZ_num; mgm++)
			{
				if (!errorCode && leapControl(FRZ->FRZyear_array[mgm], enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from mgm_init.c\n");
					errorCode=2100014;
				}
				if (leap == 1 && FRZ->FRZmonth_array[mgm] == 12 && FRZ->FRZday_array[mgm] == 31)
				{
					printf("ERROR in fertilizing date in mgm_init.c: data from 31 December in a leap year is found in fertilizing file\n");
					printf("Please read the manual and modify the input data\n");
					errorCode=2100014;
				}
			}
		}

		if (!errorCode && IRG->IRG_num)
		{
			for (mgm = 0; mgm < IRG->IRG_num; mgm++)
			{
				if (!errorCode && leapControl(IRG->IRGyear_array[mgm], enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from mgm_init.c\n");
					errorCode=2100015;
				}
				if (leap == 1 && IRG->IRGmonth_array[mgm] == 12 && IRG->IRGday_array[mgm] == 31)
				{
					printf("ERROR in irrigation date in mgm_init.c: data from 31 December in a leap year is found in irrigation file\n");
					printf("Please read the manual and modify the input data\n");
					errorCode=2100015;
				}
			}
		}

	

	}
	/* if no management */
	else
	{
		PLT->PLT_num = 0;
		THN->THN_num = 0;
		MOW->MOW_num = 0;
		GRZ->GRZ_num = 0;
		HRV->HRV_num = 0;
		PLG->PLG_num = 0;
		FRZ->FRZ_num = 0;
		IRG->IRG_num = 0;

	}
	
	if (!errorCode || errorCode > 2100) free(enddays);
	if (!errorCode || errorCode > 2100) free(mondays);

	return (errorCode);

}
