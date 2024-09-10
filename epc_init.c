/* 
epc_init.c
read epc file for pointbgc simulation

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
#include "bgc_constants.h"



int epc_init(file init, epconst_struct* epc, control_struct* ctrl, int EPCfromINI)
{
	int errorCode=0;
	int dofilecloseEPC = 1;
	double t1 = 0;
	double t2 = 0;
	double t3 = 0;
	double t4,r1;
	double sum, maxIND, diff;
	int phenphase, scanflag;
	file epc_file, FM_file, WPM_file, MSC_file, SGS_file, EGS_file; 	
	char key[] = "EPC_FILE";
	char keyword[STRINGSIZE];
	char header[STRINGSIZE];

	int p1;
	double p2;

	int* OPTINyear_array;							
	int dataread;
	int ndata = 0;

	sum=maxIND=diff=0;

	ctrl->planttypeName = (char*) malloc(STRINGSIZE * sizeof(char));
			
	/* allocate space for the optional input temporal variables */
	OPTINyear_array    = (int*) malloc(ctrl->simyears*sizeof(int));  
	
	
	/********************************************************************
	**                                                                 **
	** Begin reading initialization file block starting with keyword:  **
	** EPC_FILES                                                       ** 
	**                                                                 **
	********************************************************************/
	
	/* if EPC file reading occurs from INI  */
	
	if (EPCfromINI)
	{
		/* scan for the EPC file keyword, exit if not next */
		if (!errorCode && scan_value(init, keyword, 's'))
		{
			printf("ERROR reading keyword for control data\n");
			errorCode=209;
			dofilecloseEPC = 0;
		}
		if (!errorCode && strcmp(keyword, key))
		{
			printf("Expecting keyword --> %s in file %s\n",key,init.name);
			errorCode=209;
		}
		/* open simple EPC file  */
		if (!errorCode && scan_open(init,&epc_file,'r',1)) 
		{
			printf("ERROR opening epconst file, epc_init.c\n");
			dofilecloseEPC = 0;
			errorCode=20900;
		}
	}
	else
	{
		epc_file = init;
		/* open simple EPC file */
		if (!errorCode && file_open(&epc_file,'r',1)) 
		{
			printf("ERROR opening epconst file, epc_init.c\n");
			dofilecloseEPC = 0;
			errorCode=20900;
		}
	}


	/* first scan epc keyword to ensure proper *.epc format */
	if (!errorCode && scan_value(epc_file, ctrl->planttypeName, 's'))
	{
		printf("ERROR reading header (planttypeName), epc_init.c\n");
		errorCode=209001;
	}


	/****************************************************************************************************************/
	/* dividing line from file */ 
	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading 1. dividing line, epc_init.c\n");
		errorCode=209002;
	}

	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading 1. block title, epc_init.c\n");
		errorCode=209003;
	}
	/****************************************************************************************************************/
	/* FLAGS */
	/****************************************************************************************************************/
	if (!errorCode && scan_value(epc_file, &epc->woody, 'i'))
	{
		printf("ERROR reading woody/non-woody flag, epc_init.c\n");
		errorCode=209004;
	}

	if (!errorCode && scan_value(epc_file, &epc->evergreen, 'i'))
	{
		printf("ERROR reading evergreen/deciduous flag, epc_init.c\n");
		errorCode=209005;
	}

	if (!errorCode && scan_value(epc_file, &epc->c3_flag, 'i'))
	{
		printf("ERROR reading C3/C4 flag, epc_init.c\n");
		errorCode=209006;
	}

	/****************************************************************************************************************/
	/* dividing line from file */ 
	
	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading 2. dividing line, epc_init.c\n");
		errorCode=209007;
	}

	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading 2. block title, epc_init.c\n");
		errorCode=209008;
	}

	/****************************************************************************************************************/
	/* PLANT FUNCTIONING PARAMETERS */
	/****************************************************************************************************************/
	
	/*-------------------------------------------------------------*/
	/* using varying onday values (in transient or in normal run)*/
	if (!errorCode && scan_value(epc_file, &epc->ondayUSER, 'i'))
	{
		printf("ERROR reading onday, epc_init.c\n");
		errorCode=209009;
	}

	if (ctrl->spinup == 0) 
		strcpy(SGS_file.name,  "onday_normal.txt");
	else
	{
		if (ctrl->spinup == 1) 
			strcpy(SGS_file.name,  "onday_spinup.txt");
		else
			strcpy(SGS_file.name,  "onday_transient.txt");
	}
	
	/* SGS flag: constans or varying SGS from file */
	if (!errorCode && !file_open(&SGS_file,'j',1)) 
		ctrl->varSGS_flag = 1;
	else
		ctrl->varSGS_flag = 0;

	if (ctrl->varSGS_flag == 1 && epc->ondayUSER ==  DATA_GAP)
	{
		printf("ERROR in using annual varying onday data, if user-defined bareground simulation is defined (onday=-9999), epc_init.c\n");
		errorCode=2090091;
	}

	if (!errorCode && ctrl->varSGS_flag) 
	{
		/* allocate space for the annual SGS array */
		epc->SGS_array = (double*) malloc(ctrl->simyears * sizeof(double));
		if (!epc->SGS_array)
		{
			printf("ERROR allocating for annual SGS array, epc_init.c\n");
			errorCode=2090091;
		}
	
		/* initialization with default value */
		for (ndata = 0; ndata < ctrl->simyears; ndata++) epc->SGS_array[ndata] = epc->ondayUSER;
		
		ndata=0;
		while (!errorCode && !(dataread = scan_array (SGS_file, &p1, 'i', 0, 0)))
		{
			dataread = fscanf(SGS_file.ptr, "%lf%*[^\n]", &p2);
				
			if (p1 >= ctrl->simstartyear && p1 < ctrl->simstartyear + ctrl->simyears)
			{
				OPTINyear_array[ndata] = p1;
				epc->SGS_array[ndata]  = p2;

				if (epc->SGS_array[ndata] < 0.0 && epc->SGS_array[ndata] != DATA_GAP)
				{
					printf("ERROR in using annual varying onday data, if user-defined bareground simulation is defined (offday=-9999), epc_init.c\n");
					errorCode=2090091;
				}

				if (epc->SGS_array[ndata] > 364)
				{
					printf("ERROR in epc_init.c: SGS must less than 365\n");
					errorCode=2090091;
				}
			
				ndata += 1;
			}
		}
		/* read year and SGS for each simday in each simyear */
		
		fclose(SGS_file.ptr);
	}	
	else /* if no changing data constant EPC parameter are used */
	{
		epc->SGS_array = 0;
	}	

	/*-------------------------------------------------------------*/
	/* using varying offday values (in transient or in normal run)*/

	if (!errorCode && scan_value(epc_file, &epc->offdayUSER, 'i'))
	{
		printf("ERROR reading offday, epc_init.c\n");
		errorCode=209010;
	}

	if (ctrl->spinup == 0) 
		strcpy(EGS_file.name, "offday_normal.txt");
	else
	{
		if (ctrl->spinup == 1) 
			strcpy(EGS_file.name, "offday_spinup.txt");
		else
			strcpy(EGS_file.name, "offday_transient.txt");
	}
	
	/* EGS flag: constans or varying EGS from file */
	if (!errorCode && !file_open(&EGS_file,'j',1)) 
		ctrl->varEGS_flag = 1;
	else
		ctrl->varEGS_flag = 0;


	if (ctrl->varEGS_flag == 1 && epc->offdayUSER ==  DATA_GAP)
	{
		printf("ERROR in using annual varying offday data, if user-defined bareground simulation is defined (offday=-9999), epc_init.c\n");
		errorCode=2090101;
	}

	if (ctrl->varEGS_flag > 364)
	{
		printf("ERROR in using annual varying offday data: EGS must less than 365\n");
		errorCode=2090101;
	}

	if (!errorCode && ctrl->varEGS_flag) 
	{
		/* allocate space for the annual EGS array */
		epc->EGS_array = (double*) malloc(ctrl->simyears * sizeof(double));
		if (!epc->EGS_array)
		{
			printf("ERROR allocating for annual EGS array, epc_init.c\n");
			errorCode=2090101;
		}
	
		/* initialization with default value */
		for (ndata = 0; ndata < ctrl->simyears; ndata++) epc->EGS_array[ndata] = epc->offdayUSER;
		
		ndata=0;
		while (!errorCode && !(dataread = scan_array (EGS_file, &p1, 'i', 0, 0)))
		{
			dataread = fscanf(EGS_file.ptr, "%lf%*[^\n]", &p2);
				
			if (p1 >= ctrl->simstartyear && p1 < ctrl->simstartyear + ctrl->simyears)
			{
				OPTINyear_array[ndata] = p1;
				epc->EGS_array[ndata]  = p2;

				if (epc->EGS_array[ndata] < 0.0 && epc->EGS_array[ndata] != DATA_GAP)
				{
					printf("ERROR in epc_init.c: EGS must be positive\n");
					errorCode=2090101;
					
				}
			
				ndata += 1;
			}
		}
		/* read year and EGS for each simday in each simyear */
		
		fclose(EGS_file.ptr);
	}	
	else /* if no changing data constant EPC parameter are used */
	{
		epc->EGS_array = 0;
	}	

	/*-------------------------------------------------------------*/
	/* transfer growth and litterfall period */
	if (!errorCode && scan_value(epc_file, &epc->transfer_pdays, 'd'))
	{
		printf("ERROR reading transfer_pdays, epc_init.c\n");
		errorCode=209011;
		
	}


	if (!errorCode && scan_value(epc_file, &epc->litfall_pdays, 'd'))
	{
		printf("ERROR reading litfall_pdays, epc_init.c\n");
		errorCode=209012;
		
	}

	/*-------------------------------------------------------------*/
	/* base temperature for calculation GDD / heatsum  */
	if (!errorCode && scan_value(epc_file, &epc->base_temp, 'd'))
	{
		printf("ERROR reading base_temp parameter: epc_init.c\n");
		errorCode=209013;
		
	}

	/*-------------------------------------------------------------*/
	/* minimum/optimal/maximum temperature for growth displayed on current day (-9999: no T-dependence of allocation) */
	if (!errorCode && scan_value(epc_file, &epc->pnow_minT, 'd'))
	{
		printf("ERROR reading pnow_minT parameter: epc_init.c\n");
		errorCode=209014;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->pnow_opt1T, 'd'))
	{
		printf("ERROR reading pnow_opt1T parameter: epc_init.c\n");
		errorCode=209015;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->pnow_opt2T, 'd'))
	{
		printf("ERROR reading pnow_opt1T parameter: epc_init.c\n");
		errorCode=209016;
		
	}

	/* temperature for growth data with control */
	if (!errorCode && scan_value(epc_file, &epc->pnow_maxT, 'd'))
	{
		printf("ERROR reading pnow_maxT parameter: epc_init.c\n");
		errorCode=209017;
		
	}
	
	if (!errorCode)
	{
		if (epc->pnow_minT == DATA_GAP ||epc->pnow_opt1T == DATA_GAP || epc->pnow_opt2T == DATA_GAP || epc->pnow_maxT == DATA_GAP)
		{
			if (epc->pnow_minT != DATA_GAP ||epc->pnow_opt1T != DATA_GAP || epc->pnow_opt2T != DATA_GAP || epc->pnow_maxT != DATA_GAP) 
			{
				printf("ERROR in minimum/optimal/maximum temperature for growth data in EPC file\n");
				printf("All or none temperature data should to be set by the user\n");
				errorCode=2090171;
				
			}
		}
		if (epc->pnow_minT > epc->pnow_opt1T ||epc->pnow_opt1T > epc->pnow_opt2T || epc->pnow_opt2T > epc->pnow_maxT) 
		{
				printf("ERROR in minimum/optimal/maximum temperature for pnow data in EPC file\n");
				printf("Correct temperature data: minT <= opt1T <= opt2T <= maxT\n");
				errorCode=2090171;
				
		}
	}


	/*-------------------------------------------------------------*/
	/* minimum/optimal/maximum temperature for C-assimilation displayed on current day (-9999: no limitation) */
	
	if (!errorCode && scan_value(epc_file, &epc->assim_minT, 'd'))
	{
		printf("ERROR reading assim_minT parameter: epc_init.c\n");
		errorCode=209018;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->assim_opt1T, 'd'))
	{
		printf("ERROR reading assim_opt1T parameter: epc_init.c\n");
		errorCode=209019;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->assim_opt2T, 'd'))
	{
		printf("ERROR reading assim_opt1T parameter: epc_init.c\n");
		errorCode=209020;
		
	}

	/* temperature data for C-assimilation with control */
	if (!errorCode && scan_value(epc_file, &epc->assim_maxT, 'd'))
	{
		printf("ERROR reading assim_maxT parameter: epc_init.c\n");
		errorCode=209021;
		
	}

	if (!errorCode)
	{
		if (epc->assim_minT == DATA_GAP || epc->assim_opt1T == DATA_GAP || epc->assim_opt2T == DATA_GAP || epc->assim_maxT == DATA_GAP)
		{
			if (epc->assim_minT != DATA_GAP ||epc->assim_opt1T != DATA_GAP || epc->assim_opt2T != DATA_GAP || epc->assim_maxT != DATA_GAP) 
			{
				printf("ERROR in minimum/optimal/maximum temperature for assim data in EPC file\n");
				printf("All or none temperature data should to be set by the user\n");
				errorCode=2090211;
				
			}
		}
		if (epc->assim_minT > epc->assim_opt1T ||epc->assim_opt1T > epc->assim_opt2T || epc->assim_opt2T > epc->assim_maxT) 
		{
				printf("ERROR in minimum/optimal/maximum temperature for assim data in EPC file\n");
				printf("Correct temperature data: minT <= opt1T <= opt2T <= maxT\n");
				errorCode=2090211;
				
		}
	}

	/*-------------------------------------------------------------*/

	if (!errorCode && scan_value(epc_file, &epc->PT_ETcritT, 'd'))
	{
		printf("ERROR reading ET_critT parameter: epc_init.c\n");
		errorCode=209022;
		
	}


	/*-------------------------------------------------------------*/
	/*  leaf turnover fraction with force leaf turnover fraction to 1.0 if deciduous */
	if (!errorCode && scan_value(epc_file, &epc->nonwoody_turnover, 'd'))
	{
		printf("ERROR reading non-woody biomass turnover, epc_init.c\n");
		errorCode=209023;
		
	}
	
	if (!epc->evergreen)
	{
		epc->nonwoody_turnover = 1.0;
	}
	
	if (!errorCode && scan_value(epc_file, &epc->woody_turnover, 'd'))
	{
		printf("ERROR reading woody biomass turnover, epc_init.c\n");
		errorCode=209024;
		
	}
	/* ------------------------------------------------------ */
	/* fire mortality */

	if (!errorCode && scan_value(epc_file, &epc->FMyr, 'd'))
	{
		printf("ERROR reading fire mortality, epc_init.c\n");
		errorCode=209025;
		
	}

	/* using varying fire mortality values */
	if (ctrl->spinup == 0)
		strcpy(FM_file.name, "FM_normal.txt");
	else
	{
		if (ctrl->spinup == 1) 
			strcpy(FM_file.name, "FM_spinup.txt");
		else
			strcpy(FM_file.name, "FM_transient.txt");
	}

	/* FM flag: constans or varying EGS from file */
	if (!errorCode && !file_open(&FM_file,'j',1)) 
		ctrl->varFM_flag = 1;
	else
		ctrl->varFM_flag = 0;
	
	if (!errorCode && ctrl->varFM_flag) 
	{
		/* allocate space for the annual FM array */
		epc->FMyr_array = (double*) malloc(ctrl->simyears * sizeof(double));
		if (!epc->FMyr_array)
		{
			printf("ERROR allocating for annual FM array, epc_init.c\n");
			errorCode=2090251;
			
		}
	
		/* initialization with default value */
		for (ndata = 0; ndata < ctrl->simyears; ndata++) epc->FMyr_array[ndata] = epc->FMyr;
		
		ndata=0;
		while (!errorCode && !(dataread = scan_array (FM_file, &p1, 'i', 0, 0)))
		{
			dataread = fscanf(FM_file.ptr, "%lf%*[^\n]", &p2);
				
			if (p1 >= ctrl->simstartyear && p1 < ctrl->simstartyear + ctrl->simyears)
			{
				OPTINyear_array[ndata] = p1;
				epc->FMyr_array[ndata]  = p2;

				if (epc->FMyr_array[ndata] < 0.0 && epc->FMyr_array[ndata] != DATA_GAP)
				{
					printf("ERROR in epc_init.c: FM must be positive\n");
					errorCode=2090251;
					
				}
			
				ndata += 1;
			}
		}
		/* read year and FM for each simday in each simyear */
		
		fclose(FM_file.ptr);
	}	
	else /* if no changing data constant EPC parameter are used */
	{
		epc->FMyr_array = 0;
	}	
	
	/*-------------------------------------------------------------*/
	/* whole plant mortality */

	if (!errorCode && scan_value(epc_file, &epc->WPMyr, 'd'))
	{
		printf("ERROR reading whole-plant mortality, epc_init.c\n");
		errorCode=209026;
		
	}

	/* using varying whole plant mortality values */
	if (ctrl->spinup == 0)
		strcpy(WPM_file.name, "WPM_normal.txt");
	else
	{
		if (ctrl->spinup == 1) 
			strcpy(WPM_file.name, "WPM_spinup.txt");
		else
			strcpy(WPM_file.name, "WPM_transient.txt");
	}

	/* WPM flag: constans or varying EGS from file */
	if (!errorCode && !file_open(&WPM_file,'j',1)) 
		ctrl->varWPM_flag = 1;
	else
		ctrl->varWPM_flag = 0;
	
	if (!errorCode && ctrl->varWPM_flag) 
	{
		/* allocate space for the annual WPM array */
		epc->WPMyr_array = (double*) malloc(ctrl->simyears * sizeof(double));
		if (!epc->WPMyr_array)
		{
			printf("ERROR allocating for annual WPM array, epc_init.c\n");
			errorCode=2090261;
			
		}
	
		/* initialization with default value */
		for (ndata = 0; ndata < ctrl->simyears; ndata++) epc->WPMyr_array[ndata] = epc->WPMyr;
		
		ndata=0;
		while (!errorCode && !(dataread = scan_array (WPM_file, &p1, 'i', 0, 0)))
		{
			dataread = fscanf(WPM_file.ptr, "%lf%*[^\n]", &p2);
				
			if (p1 >= ctrl->simstartyear && p1 < ctrl->simstartyear + ctrl->simyears)
			{
				OPTINyear_array[ndata] = p1;
				epc->WPMyr_array[ndata]  = p2;

				if (epc->WPMyr_array[ndata] < 0.0 && epc->WPMyr_array[ndata] != DATA_GAP)
				{
					printf("ERROR in epc_init.c: WPM must be positive\n");
					errorCode=2090262;
					
				}
			
				ndata += 1;
			}
		}
		/* read year and WPM for each simday in each simyear */
		
		fclose(WPM_file.ptr);
	}	
	else /* if no changing data constant EPC parameter are used */
	{
		epc->WPMyr_array = 0;
	}	


	/* ------------------------------------------------------ */
	/* dead stem biomass combustion proportion */
	if (!errorCode && scan_value(epc_file, &epc->dscp, 'd'))
	{
		printf("ERROR reading dead stem biomass combustion proportion, epc_init.c\n");
		errorCode=209027;
		
	}

	/* ------------------------------------------------------ */
	/* dead stem biomass combustion proportion */
	if (!errorCode && scan_value(epc_file, &epc->cwcp, 'd'))
	{
		printf("ERROR reading coarse woody biomass combustion proportion, epc_init.c\n");
		errorCode=209028;
		
	}

	/*--------------------------------------------------------------------*/
	/* C:N ratios */

	/* leaf C:N > 0 */
	if (!errorCode && scan_value(epc_file, &epc->leaf_cn, 'd'))
	{
		printf("ERROR reading average leaf C:N, epc_init.c\n");
		errorCode=209029;
		
	}

	if (!errorCode && epc->leaf_cn <= 0	)
	{
		printf("ERROR in epc_init.c:leaf  C:N must be > 0\n");
		printf("change the values in EPC file\n");
		errorCode=2090291;
		
	}

	/* leaflitter C:N > leaf C:N  */
	if (!errorCode && scan_value(epc_file, &epc->leaflitr_cn, 'd'))
	{
		printf("ERROR reading leaf litter C:N, epc_init.c\n");
		errorCode=209030;
		
	}
	
	if (!errorCode && epc->leaflitr_cn < epc->leaf_cn)
	{
		printf("ERROR in epc_init.c:leaf litter C:N must be >= leaf C:N\n");
		printf("change the values in EPC file\n");
		errorCode=2090301;
		
	}

	/* froot C:N > leaf C:N  */
	if (!errorCode && scan_value(epc_file, &epc->froot_cn, 'd'))
	{
		printf("ERROR reading initial fine root C:N, epc_init.c\n");
		errorCode=209031;
		
	}
	if (!errorCode && epc->froot_cn < epc->leaf_cn)
	{
		printf("ERROR in epc_init.c:fine root C:N must be >= leaf C:N\n");
		printf("change the values in EPC file\n");
		errorCode=209031;
		
	}

	/* yield C:N > leaf C:N  */
	if (!errorCode && scan_value(epc_file, &epc->yield_cn, 'd'))
	{
		printf("ERROR reading initial yield C:N, epc_init.c\n");
		errorCode=209032;
		
	}
	if (!errorCode && epc->yield_cn > 0 && epc->yield_cn < epc->leaf_cn)
	{
		printf("ERROR in epc_init.c:yield C:N must be >= leaf C:N\n");
		printf("change the values in EPC file\n");
		errorCode=2090321;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->softstem_cn, 'd'))
	{
		printf("ERROR reading initial softstem C:N, epc_init.c\n");
		errorCode=209033;
		
	}
	if (!errorCode && epc->softstem_cn > 0 && epc->softstem_cn < epc->leaf_cn)
	{
		printf("ERROR in epc_init.c:softstem C:N must be >= leaf C:N\n");
		printf("change the values in EPC file\n");
		errorCode=2090331;
		
	}


	/* livewood C:N setting for woody biomes */
	if (!errorCode && scan_value(epc_file, &epc->livewood_cn, 'd'))
	{
		printf("ERROR reading initial livewood C:N, epc_init.c\n");
		errorCode=209034;
		
	}
	if (!errorCode && epc->woody && epc->livewood_cn == 0)
	{
		printf("ERROR in epc_init.c:livewood C:N must be > 0 in woody simulation \n");
		printf("change the values in EPC file\n");
		errorCode=2090341;
		
	}

	/* deadwood C:N > livewood C:N */
	if (!errorCode && scan_value(epc_file, &epc->deadwood_cn, 'd'))
	{
		printf("ERROR reading initial deadwood C:N, epc_init.c\n");
		errorCode=2090351;
		
	}
	if (!errorCode && epc->deadwood_cn < epc->livewood_cn)
	{
		printf("ERROR in epc_init.c:deadwood C:N must be >= livewood C:N\n");
		printf("change the values in EPC file\n");
		errorCode=2090351;
		
	}


	/*--------------------------------------------------------------------*/
	/* C:DM ratios */
	if (!errorCode && scan_value(epc_file, &epc->leafC_DM, 'd'))
	{
		printf("ERROR reading dry matter carbon content of leaves, epc_init.c\n");
		errorCode=209036;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->leaflitrC_DM, 'd'))
	{
		printf("ERROR reading dry matter carbon content of leaf litter, epc_init.c\n");
		errorCode=209037;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->frootC_DM, 'd'))
	{
		printf("ERROR reading dry matter carbon content of fine root, epc_init.c\n");
		errorCode=209038;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->yield_DM, 'd'))
	{
		printf("ERROR reading dry matter carbon content of yield, epc_init.c\n");
		errorCode=209039;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->softstemC_DM, 'd'))
	{
		printf("ERROR reading dry matter carbon content of softstem, epc_init.c\n");
		errorCode=209040;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->livewoodC_DM, 'd'))
	{
		printf("ERROR reading dry matter carbon content of live wood, epc_init.c\n");
		errorCode=209041;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->deadwoodC_DM, 'd'))
	{
		printf("ERROR reading dry matter carbon content of dead wood, epc_init.c\n");
		errorCode=209042;
		
	}

	/*--------------------------------------------------------------------*/

	/* LEAF LITTER PROPORTION */
	if (!errorCode && scan_value(epc_file, &t1, 'd'))
	{
		printf("ERROR reading leaf litter labile proportion, epc_init.c\n");
		errorCode=209043;
		
	}

	if (!errorCode) epc->leaflitr_flab = t1;

	if (!errorCode && scan_value(epc_file, &t2, 'd'))
	{
		printf("ERROR reading leaf litter cellulose proportion, epc_init.c\n");
		errorCode=209044;
		
	}

	/* livecroot lignin proportion is calculated from labile and cellulose */
	t3 = 1 - t2 - t1 ;

    /* test for fractions sum to 1.0 */
	if (!errorCode && ((t1+t2) > 1))
	{
		printf("ERROR in fractions: %.2f %.2f %.2f %.3f\n",t1, t2, t3, t1+t2+t3);
		printf("leaf litter proportions of labile and cellulose \n");
		printf("must less than 1.0. Check initialization file and try again.\n");
		errorCode=2090441;
		
	}

	if (!errorCode) epc->leaflitr_flig = t3;

	/* calculate shielded and unshielded cellulose fraction */
	if (!errorCode)
	{
		r1 = t3/t2;
		if (r1 <= 0.45)
		{
			epc->leaflitr_fscel = 0.0;
			epc->leaflitr_fucel = t2;
		}
		else if (r1 > 0.45 && r1 < 0.7)
		{
			t4 = (r1 - 0.45)*3.2;
			epc->leaflitr_fscel = t4*t2;
			epc->leaflitr_fucel = (1.0 - t4)*t2;
		}
		else
		{
			epc->leaflitr_fscel = 0.8*t2;
			epc->leaflitr_fucel = 0.2*t2;
		}
	}

	/*--------------------------------------------------------------------*/
	/* FROOT LITTER PROPORTION */
	if (!errorCode && scan_value(epc_file, &t1, 'd'))
	{
		printf("ERROR reading froot litter labile proportion, epc_init.c\n");
		errorCode=209045;
		
	}

	epc->frootlitr_flab = t1;
	if (!errorCode && scan_value(epc_file, &t2, 'd'))
	{
		printf("ERROR reading froot litter cellulose proportion, epc_init.c\n");
		errorCode=209046;
		
	}
		
	/* livecroot lignin proportion is calculated from labile and cellulose */
	t3 = 1 - t2 - t1;
	
     
	/* test for fractions sum to 1.0 */
	if (!errorCode && ((t1+t2) > 1))
	{
		printf("ERROR in fractions: %.2f %.2f %.2f %.3f\n",t1, t2, t3, t1+t2+t3);
		printf("froot proportions of labile and cellulose \n");
		printf("must less than 1.0. Check initialization file and try again.\n");
		errorCode=2090461;
		
	}

	if (!errorCode) epc->frootlitr_flig = t3;
	

	/* calculate shielded and unshielded cellulose fraction */
	if (!errorCode)
	{
		r1 = t3/t2;
		if (r1 <= 0.45)
		{
			epc->frootlitr_fscel = 0.0;
			epc->frootlitr_fucel = t2;
		}
		else if (r1 > 0.45 && r1 < 0.7)
		{
			t4 = (r1 - 0.45)*3.2;
			epc->frootlitr_fscel = t4*t2;
			epc->frootlitr_fucel = (1.0 - t4)*t2;
		}
		else
		{
			epc->frootlitr_fscel = 0.8*t2;
			epc->frootlitr_fucel = 0.2*t2;
		}
	}

	/*--------------------------------------------------------------------*/
	/* yield LITTER PROPORTION */

	if (!errorCode && scan_value(epc_file, &t1, 'd'))
	{
		printf("ERROR reading yield litter labile proportion, epc_init.c\n");
		errorCode=209047;
		
	}

	epc->yieldlitr_flab = t1;
	if (!errorCode && scan_value(epc_file, &t2, 'd'))
	{
		printf("ERROR reading yield litter cellulose proportion, epc_init.c\n");
		errorCode=209048;
		
	}
			
	/* livecroot lignin proportion is calculated from labile and cellulose */
	t3 = 1 - t2 - t1;

     
	/* test for fractions sum to 1.0 */
	if (!errorCode && ((t1+t2) > 1))
	{
		printf("ERROR in fractions: %.2f %.2f %.2f %.3f\n",t1, t2, t3, t1+t2+t3);
		printf("froot proportions of labile and cellulose \n");
		printf("must less than 1.0. Check initialization file and try again.\n");
		errorCode=2090481;
		
	}

	if (!errorCode) epc->yieldlitr_flig = t3;
		
	/* calculate shielded and unshielded cellulose fraction */
	if (!errorCode)
	{
		r1 = t3/t2;
		if (r1 <= 0.45)
		{
			epc->yieldlitr_fscel = 0.0;
			epc->yieldlitr_fucel = t2;
		}
		else if (r1 > 0.45 && r1 < 0.7)
		{
			t4 = (r1 - 0.45)*3.2;
			epc->yieldlitr_fscel = t4*t2;
			epc->yieldlitr_fucel = (1.0 - t4)*t2;
		}
		else
		{
			epc->yieldlitr_fscel = 0.8*t2;
			epc->yieldlitr_fucel = 0.2*t2;
		}
	}

	/*--------------------------------------------------------------------*/
	/* SOFT STEM LITTER PROPORTION */
	
	if (!errorCode && scan_value(epc_file, &t1, 'd'))
	{
		printf("ERROR reading soft stem litter labile proportion, epc_init.c\n");
		errorCode=209049;
		
	}
	epc->softstemlitr_flab = t1;


	if (!errorCode && scan_value(epc_file, &t2, 'd'))
	{
		printf("ERROR reading soft stem  litter cellulose proportion, epc_init.c\n");
		errorCode=209050;
		
	}
	
	/* livecroot lignin proportion is calculated from labile and cellulose */
	t3 = 1 - t2 - t1;

	 /* test for fractions sum to 1.0 */
	if (!errorCode && ((t1+t2) > 1))
	{
		printf("ERROR in fractions: %.2f %.2f %.2f %.3f\n",t1, t2, t3, t1+t2+t3);
		printf("softstem proportions of labile and cellulose \n");
		printf("must less than 1.0. Check initialization file and try again.\n");
		errorCode=2090501;
		
	}


	if (!errorCode) epc->softstemlitr_flig = t3;

	/* calculate shielded and unshielded cellulose fraction */
	if (!errorCode)
	{
		r1 = t3/t2;
		if (r1 <= 0.45)
		{
			epc->softstemlitr_fscel = 0.0;
			epc->softstemlitr_fucel = t2;
		}
		else if (r1 > 0.45 && r1 < 0.7)
		{
			t4 = (r1 - 0.45)*3.2;
			epc->softstemlitr_fscel = t4*t2;
			epc->softstemlitr_fucel = (1.0 - t4)*t2;
		}
		else
		{
			epc->softstemlitr_fscel = 0.8*t2;
			epc->softstemlitr_fucel = 0.2*t2;
		}
	}

	/*--------------------------------------------------------------------*/
	/* DEAD WOOD LITTER PROPORTION */

	if (!errorCode && scan_value(epc_file, &t1, 'd'))
	{
		printf("ERROR reading dead wood %% cellulose, epc_init.c\n");
		errorCode=209051;
		
	}

     /* lignin proportion is calculated from cellulose */
	t2 = 1 - t1;


      /* test for fractions sum to 1.0 */
	if (!errorCode && t1 > 1)
	{
		printf("ERROR in fractions: %.2f %.2f %.3f\n",t1, t2, t1+t2);
		printf("dead wood cellulose, and lignin\n");
		printf("must less than 1.0. Check initialization file and try again.\n");
		errorCode=2090511;
		
	}

	if (!errorCode) epc->deadwood_flig = t2;

	if (!errorCode && epc->woody && epc->deadwood_flig == 0)
	{
		printf("ERROR in epc_init.c:dead wood cellulose proportion must be > 0 in woody simulation \n");
		printf("change the values in EPC file\n");
		errorCode=2090511;
		
	}
	
	/* calculate shielded and unshielded cellulose fraction */
	if (!errorCode)
	{
		r1 = t2/t1;
		if (r1 <= 0.45)
		{
			epc->deadwood_fscel = 0.0;
			epc->deadwood_fucel = t1;
		}
		else if (r1 > 0.45 && r1 < 0.7)
		{
			t4 = (r1 - 0.45)*3.2;
			epc->deadwood_fscel = t4*t1;
			epc->deadwood_fucel = (1.0 - t4)*t1;
		}
		else
		{
			epc->deadwood_fscel = 0.8*t1;
			epc->deadwood_fucel = 0.2*t1;
		}
	}

	/* ---------------------------------------------- */
	/* radiation and LAI-parameters */

	if (!errorCode && scan_value(epc_file, &epc->int_coef, 'd'))
	{
		printf("ERROR reading canopy water int coef, epc_init.c\n");
		errorCode=209052;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->ext_coef, 'd'))
	{
		printf("ERROR reading canopy light ext.coef, epc_init.c\n");
		errorCode=209053;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->potRUE, 'd'))
	{
		printf("ERROR reading potential radiation use efficiency, epc_init.c\n");
		errorCode=209054;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->rad_param1, 'd'))
	{
		printf("ERROR reading radiation parameter1, epc_init.c\n");
		errorCode=209055;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->rad_param2, 'd'))
	{
		printf("ERROR reading radiation parameter2, epc_init.c\n");
		errorCode=209056;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->lai_ratio, 'd'))
	{
		printf("ERROR reading all to projected LAI ratio, epc_init.c\n");
		errorCode=209057;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->sla_ratio, 'd'))
	{
		printf("ERROR reading shaded to sunlit SLA ratio, epc_init.c\n");
		errorCode=209058;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->flnr, 'd'))
	{
		printf("ERROR reading Rubisco N fraction, epc_init.c\n");
		errorCode=209059;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->flnp, 'd'))
	{
		printf("ERROR reading PeP N fraction, epc_init.c\n");
		errorCode=209060;
		
	}

	/* ------------------------------------------------------ */
	/* maximum stomatal conductance (optional annual varying value) */

	if (!errorCode && scan_value(epc_file, &epc->gl_sMAX, 'd'))
	{
		printf("ERROR reading gl_sMAX, epc_init.c\n");
		errorCode=209061;
		
	}

	if (ctrl->spinup == 0)
		strcpy(MSC_file.name, "conductance_normal.txt");
	else
	{
		if (ctrl->spinup == 1) 
			strcpy(MSC_file.name, "conductance_spinup.txt");
		else
			strcpy(MSC_file.name, "conductance_transient.txt");
	}
	
	/* MSC flag: constans or varying MSC from file */
	if (!errorCode && !file_open(&MSC_file,'j',1)) 
		ctrl->varMSC_flag = 1;
	else
		ctrl->varMSC_flag = 0;

	if (!errorCode && ctrl->varMSC_flag) 
	{
		/* allocate space for the annual MSC array */
		epc->MSC_array = (double*) malloc(ctrl->simyears * sizeof(double));
		if (!epc->MSC_array)
		{
			printf("ERROR allocating for annual MSC array, epc_init.c\n");
			errorCode=2090611;
			
		}
	
		/* initialization with default value */
		for (ndata = 0; ndata < ctrl->simyears; ndata++) epc->MSC_array[ndata] = epc->gl_sMAX;
		
		ndata=0;
		while (!errorCode && !(dataread = scan_array (MSC_file, &p1, 'i', 0, 0)))
		{
			dataread = fscanf(MSC_file.ptr, "%lf%*[^\n]", &p2);
				
			if (p1 >= ctrl->simstartyear && p1 < ctrl->simstartyear + ctrl->simyears)
			{
				OPTINyear_array[ndata] = p1;
				epc->MSC_array[ndata]  = p2;

				if (epc->MSC_array[ndata] < 0.0 && epc->MSC_array[ndata] != DATA_GAP)
				{
					printf("ERROR in epc_init.c: MSC must be positive\n");
					errorCode=2090611;
					
				}
			
				ndata += 1;
			}
		}
		/* read year and MSC for each simday in each simyear */
		
		fclose(MSC_file.ptr);
	}	
	else /* if no changing data constant EPC parameter are used */
	{
		epc->MSC_array = 0;
	}	
	
	
	/* ------------------------------------------------------ */
	if (!errorCode && scan_value(epc_file, &epc->gl_c, 'd'))
	{
		printf("ERROR reading gl_c, epc_init.c\n");
		errorCode=209062;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->gl_bl, 'd'))
	{
		printf("ERROR reading gl_bl, epc_init.c\n");
		errorCode=209063;
		
	}

	/* plantHeight_max */
	if (!errorCode && scan_value(epc_file, &epc->plantHeight_max, 'd'))
	{
		printf("ERROR reading number of plantHeight_max\n");
		errorCode=209064;
		
	}

	/* plantHeight parameter 1 (weight of stem at maximum height */  
	if (!errorCode && scan_value(epc_file, &epc->plantHeight_par1, 'd'))
	{
		printf("ERROR reading plantHeight_par1\n");
		errorCode=209065;
		
	}

	/* plantHeight parameter 2 (slope) */  
	if (!errorCode && scan_value(epc_file, &epc->plantHeight_par2, 'd'))
	{
		printf("ERROR reading plantHeight_par2\n");
		errorCode=209066;
		
	}

	/* rootzoneDepth_max  */
	if (!errorCode && scan_value(epc_file, &epc->rootzoneDepth_max, 'd'))
	{
		printf("ERROR reading number of rootzoneDepth_max\n");
		errorCode=209067;
		
	}
	/* CONTROL to avoid negative data */
 	if (!errorCode && epc->rootzoneDepth_max <= 0)
	{
		printf("ERROR in epc data in EPC file: negative or zero rootzone_depth, epc_init.c\n");
		errorCode=209067;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->rootdistrib_param, 'd'))
	{
		printf("ERROR reading rootdistrib_param: epc_init.c\n");
		errorCode=209068;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->rootlength_par1, 'd'))
	{
		printf("ERROR reading rootlength_par1: epc_init.c\n");
		errorCode=209069;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->rootlength_par2, 'd'))
	{
		printf("ERROR reading rootlength_par2: epc_init.c\n");
		errorCode=209070;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->rootweight2length, 'd'))
	{
		printf("ERROR reading rootweight2length: epc_init.c\n");
		errorCode=209071;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->GR_ratio, 'd'))
	{
		printf("ERROR reading growth resp.ratio: epc_init.c\n");
		errorCode=209072;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->MRperN, 'd'))
	{
		printf("ERROR reading MRperN: epc_init.c\n");
		errorCode=209073;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->prop_NSCvsSC, 'd'))
	{
		printf("ERROR reading prop_NSCvsSC: epc_init.c\n");
		errorCode=209074;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->prop_NSCcrit, 'd'))
	{
		printf("ERROR reading prop_NSCcrit: epc_init.c\n");
		errorCode=209075;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->nfix, 'd'))
	{
		printf("ERROR reading N fixation, epc_init.c\n");
		errorCode=209076;
		
	}
	
	/* CONTROL to avoid irrealistic temperature data */
	if (!errorCode && scan_value(epc_file, &epc->tau, 'd'))
	{
		printf("ERROR reading tau (time delay photosynthesis acclim.), epc_init.c\n");
		errorCode=209077;
		
	}
	
 	if (!errorCode && epc->tau < 1 && epc->tau > 0)
	{
		if (epc->tau > 0.5) 
			epc->tau = 1;
		else
			epc->tau = 0;
	}

	/*********************************************************************************************************************/
	/* dividing line from file */ 
	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading 3. dividing line, epc_init.c\n");
		errorCode=209078;
		
	}

	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading 3. block title, epc_init.c\n");
		errorCode=209079;
		
	}

	/*********************************************************************************************************************/
	/* CROP SPECIFIC PARAMETERS */
	/*********************************************************************************************************************/

	/* number of phenophase of germination */
	if (!errorCode && scan_value(epc_file, &epc->n_germ_phenophase, 'i'))
	{
		printf("ERROR reading number of phenophase of germination \n");
		errorCode=209080;
		
	}

	/* number of phenophase of emergence */
	if (!errorCode && scan_value(epc_file, &epc->n_emerg_phenophase, 'i'))
	{
		printf("ERROR reading number of phenophase of emergence\n");
		errorCode=209081;
		
	}
	
	/* critical relative VWC (prop. to FC) in germination */
	if (!errorCode && scan_value(epc_file, &epc->grmn_paramVWC, 'd'))
	{
		printf("ERROR reading critical VWC ratio in germination, epc_init.c\n");
		errorCode=209082;
		
	}

	/* number of phenophase of photoperiodic slowing effect */
	if (!errorCode && scan_value(epc_file, &epc->n_phpsl_phenophase, 'i'))
	{
		printf("ERROR reading number of phenophase of photoperiodic slowing effect \n");
		errorCode=209083;
		
	}

	/* critical photoslow daylength */
	if (!errorCode && scan_value(epc_file, &epc->phpsl_parDL, 'd'))
	{
		printf("ERROR reading critical photoslow daylength, epc_init.c\n");
		errorCode=209084;
		
	}
	
	/* slope of relative photoslow development rate */
	if (!errorCode && scan_value(epc_file, &epc->phpsl_parDR, 'd'))
	{
		printf("ERROR reading slope of relative photoslow development rate, epc_init.c\n");
		errorCode=209085;
		
	}

	/* number of phenophase of vernalization */
	if (!errorCode && scan_value(epc_file, &epc->n_vern_phenophase, 'i'))
	{
		printf("ERROR reading number of vernalization phenophase\n");
		errorCode=209086;
		
	}

	/* critical vernalization temperature data */
	if (!errorCode && scan_value(epc_file, &epc->vern_parT1, 'd'))
	{
		printf("ERROR reading critical vernalization temperature 1, epc_init.c\n");
		errorCode=209087;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->vern_parT2, 'd'))
	{
		printf("ERROR reading critical vernalization temperature 2, epc_init.c\n");
		errorCode=209088;
		
	}

	if (!errorCode && scan_value(epc_file, &epc->vern_parT3, 'd'))
	{
		printf("ERROR reading critical vernalization temperature 3, epc_init.c\n");
		errorCode=209089;
	}

	if (!errorCode && scan_value(epc_file, &epc->vern_parT4, 'd'))
	{
		printf("ERROR reading critical vernalization temperature 4, epc_init.c\n");
		errorCode=209090;
	}
	
	/* control of vernalization temperature parameters */
	if (!errorCode)
	{
		if (epc->vern_parT1 > epc->vern_parT2 ||epc->vern_parT2 > epc->vern_parT3 || epc->vern_parT3 > epc->vern_parT4) 
		{
				printf("ERROR in vernalization temperature data in EPC file\n");
				printf("Correct temperature data: vern_parT1 <= vern_parT2 <= vern_parT3 <= vern_parT4\n");
				errorCode=2090901;
		}
	}

	/* slope of relative vernalization development rate  */
	if (!errorCode && scan_value(epc_file, &epc->vern_parDR1, 'd'))
	{
		printf("ERROR reading slope of relative vernalization development rate, epc_init.c\n");
		errorCode=209091;
	}

	/* required vernalization days (in vernalization development rate  */
	if (!errorCode && scan_value(epc_file, &epc->vern_parDR2, 'd'))
	{
		printf("ERROR reading required vernalization days (in vernalization development rate, epc_init.c\n");
		errorCode=209092;
	}

	/* number of phenophase of flowering heat stress */
	if (!errorCode && scan_value(epc_file, &epc->n_flowHS_phenophase, 'i'))
	{
		printf("ERROR reading number of phenophase of flowering heat stress\n");
		errorCode=209093;
	}

	/* critical flowering heat stress temperatures */
	if (!errorCode && scan_value(epc_file, &epc->flowHS_parT1, 'd'))
	{
		printf("ERROR reading critical flowering heat stress temperature 1, epc_init.c\n");
		errorCode=209094;
	}
	
	/* flowering heat stress temperature parameter with control */
	if (!errorCode && scan_value(epc_file, &epc->flowHS_parT2, 'd'))
	{
		printf("ERROR reading critical flowering heat stress temperature 2, epc_init.c\n");
		errorCode=209095;
	}
	if (epc->flowHS_parT1 > epc->flowHS_parT2) 
	{
		printf("ERROR in flowering heat stress temperature data in EPC file\n");
		printf("Correct temperature data: flowHS_parT1 <= flowHS_parT2 \n");
		errorCode=2090951;
	}

	/* mortality parameter of flowering heat stress */
	if (!errorCode && scan_value(epc_file, &epc->flowHS_parMORT, 'd'))
	{
		printf("ERROR reading mortality parameter of flowering heat stress, epc_init.c\n");
		errorCode=209096;
	}

	/*********************************************************************************************************************/
	/* dividing line from file */ 
	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading 4. dividing line, epc_init.c\n");
		errorCode=209097;
	}

	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading 4. block titke, epc_init.c\n");
		errorCode=209098;
	}

	/*********************************************************************************************************************/
	/* STRESS AND SENESCENCE PARAMETERS */
	/*********************************************************************************************************************/

	if (!errorCode && scan_value(epc_file, &epc->VWCratio_WScrit1, 'd'))
	{
		printf("ERROR reading VWCratio_WScrit1, epc_init.c\n");
		errorCode=209099;
	}
	
	if (!errorCode && scan_value(epc_file, &epc->VWCratio_WScrit2, 'd'))
	{
		printf("ERROR reading VWCratio_WScrit2, epc_init.c\n");
		errorCode=209100;
	}

	if (!errorCode && scan_value(epc_file, &epc->nday_fullWS, 'd'))
	{
		printf("ERROR reading nday_fullWS: epc_init.c\n");
		errorCode = 209101;
	}

	/* control: nday_fullWS must greater than 1 -> number of days */
	if (!errorCode && epc->nday_fullWS < 1)
	{
		printf("ERROR in critical value of total water stress length parameter in EPC file\n");
		printf("Correct parameter value: must greater than 1 \n");
		errorCode = 2091011;
	}
	if (!errorCode && scan_value(epc_file, &epc->vpd_open, 'd')) 
	{
		printf("ERROR reading vpd_max, epc_init.c\n");
		errorCode=209102;
	}

	if (!errorCode && scan_value(epc_file, &epc->vpd_close, 'd'))
	{
		printf("ERROR reading vpd_min, epc_init.c\n");
		errorCode=209103;
	}

	if (!errorCode && scan_value(epc_file, &epc->maxSNSCmort_leaf, 'd'))
	{
		printf("ERROR reading maximal senescence mortality parameter of leaf: epc_init.c\n");
		errorCode=209104;
	}

	if (!errorCode && scan_value(epc_file, &epc->maxSNSCmort_other, 'd'))
	{
		printf("ERROR reading maximal senescence mortality of softstem and froot: epc_init.c\n");
		errorCode=209105;
	}

	if (!errorCode && scan_value(epc_file, &epc->m_nscSNSCmort, 'd'))
	{
		printf("ERROR reading multiplier of senescence mortality calculation of non-stuctured biomass: epc_init.c\n");
		errorCode=209106;
	}


	if (!errorCode && scan_value(epc_file, &epc->SNSC_extremT1, 'd'))
	{
		printf("ERROR reading lower limit extreme high temperature effect on senesncene mortality: epc_init.c\n");
		errorCode=209107;
	}

	/* flowering heat stress temperature parameter with control */
	if (!errorCode && scan_value(epc_file, &epc->SNSC_extremT2, 'd'))
	{
		printf("ERROR reading lower limit extreme high temperature effect on senesncene mortality: epc_init.c\n");
		errorCode=209108;
	}
	if (!errorCode && epc->SNSC_extremT1 >= epc->SNSC_extremT2) 
	{
		printf("ERROR in extreme high temperature effect data in EPC file\n");
		printf("Correct temperature data: SNSC_extremT1 < SNSC_extremT2 \n");
		errorCode=2091081;
	}

	if (!errorCode && scan_value(epc_file, &epc->mort_SNSC_to_litter, 'd'))
	{
		printf("ERROR reading turnover rate of wilted standing biomass to litter parameter: epc_init.c\n");
		errorCode=209109;
	}

	if (!errorCode && scan_value(epc_file, &epc->mort_CnW_to_litter, 'd'))
	{
		printf("ERROR reading turnover rate of cut-down non-woody biomass to litter parameter: epc_init.c\n");
		errorCode=209110;
	}

	if (!errorCode && scan_value(epc_file, &epc->mort_CW_to_litter, 'd'))
	{
		printf("ERROR reading turnover rate of cut-down woody biomass to litter parameter: epc_init.c\n");
		errorCode=209111;
	}

	if (!errorCode && scan_value(epc_file, &epc->WSlenght_crit, 'd'))
	{
		printf("ERROR reading WSlenght_crit: epc_init.c\n");
		errorCode=209112;
	}

	if (!errorCode && scan_value(epc_file, &epc->photoSTRESSeffect, 'd'))
	{
		printf("ERROR reading photoSTRESSeffect: epc_init.c\n");
		errorCode=209113;
	}

	/*********************************************************************************************************************/
	/* dividing line from file */ 
	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading 5. dividing line, epc_init.c\n");
		errorCode=209114;
	}

	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading 5. block title, epc_init.c\n");
		errorCode=209115;
	}

	/*********************************************************************************************************************/
	/* GROWING SEASON PARAMETERS */
	/*********************************************************************************************************************/
	
	/* snowcover_limit for calculation heatsum */
	if (!errorCode && scan_value(epc_file, &epc->snowcover_limit, 'd'))
	{
		printf("ERROR reading snowcover_limit parameter: epc_init.c\n");
		errorCode=209116;
	}

	/* heatsum_limit for calculation heatsum index */
	if (!errorCode && scan_value(epc_file, &epc->heatsum_limit1, 'd'))
	{
		printf("ERROR reading heatsum_limit1 parameter: epc_init.c\n");
		errorCode=209117;
	}

	/* heatsum_limit2 for calculation heatsum index with control */
	if (!errorCode && scan_value(epc_file, &epc->heatsum_limit2, 'd'))
	{
		printf("ERROR reading heatsum_limit2 parameter: epc_init.c\n");
		errorCode=209118;
	}

	if (!errorCode && epc->heatsum_limit1 > epc->heatsum_limit2) 
	{
		printf("ERROR in heatsum_limit data in EPC file\n");
		printf("Correct data: limit1 < limit2 \n");
		errorCode=2091181;
	}

	/* Tmin_limit1 for calculation Tmin index */
	if (!errorCode && scan_value(epc_file, &epc->Tmin_limit1, 'd'))
	{
		printf("ERROR reading Tmin_limit1 parameter: epc_init.c\n");
		errorCode=209119;
	}

	/* Tmin_limit1 for calculation Tmin index */
	if (!errorCode && scan_value(epc_file, &epc->Tmin_limit2, 'd'))
	{
		printf("ERROR reading Tmin_limit2 parameter: epc_init.c\n");
		errorCode=209120;
	}

	if (!errorCode && epc->Tmin_limit1 > epc->Tmin_limit2) 
	{
		printf("ERROR in Tmin_limit data in EPC file\n");
		printf("Correct data: limit1 < limit2 \n");
		errorCode=2091201;
	}

	/* vpd_limit for calculation vpd index */
	if (!errorCode && scan_value(epc_file, &epc->vpd_limit1, 'd'))
	{
		printf("ERROR reading vpd_limit1 parameter: epc_init.c\n");
		errorCode=209121;
	}
	
	if (!errorCode && scan_value(epc_file, &epc->vpd_limit2, 'd'))
	{
		printf("ERROR reading vpd_limit2 parameter: epc_init.c\n");
		errorCode=209122;
	}
	
	if (!errorCode && epc->vpd_limit1 < epc->vpd_limit2) 
	{
		printf("ERROR in vpd_limit data in EPC file\n");
		printf("Correct data: limit1 > limit2 \n");
		errorCode=2091221;
	}

	/* dayl_limit for calculation dayl index */
	if (!errorCode && scan_value(epc_file, &epc->dayl_limit1, 'd'))
	{
		printf("ERROR reading dayl_limit1 parameter: epc_init.c\n");
		errorCode=209123;
	}
	

	if (!errorCode && scan_value(epc_file, &epc->dayl_limit2, 'd'))
	{
		printf("ERROR reading dayl_limit2 parameter: epc_init.c\n");
		errorCode=209124;
	}

	if (!errorCode && epc->dayl_limit1 > epc->dayl_limit2) 
	{
		printf("ERROR in vpd_limit data in EPC file\n");
		printf("Correct data: limit1 > limit2 \n");
		errorCode=2091241;
	}

	/* n_moving_avg for calculation moving average from indexes */
	if (!errorCode && scan_value(epc_file, &epc->n_moving_avg, 'i'))
	{
		printf("ERROR reading n_moving_avg parameter: epc_init.c\n");
		errorCode=209125;
	}
	if (!errorCode && (epc->n_moving_avg < 2))
	{
		printf("ERROR in moving average parameter of epc calculation:\n");
		printf("must greater then 2. Check EPC file and try again.\n");
		errorCode=2091251;
	}

	/* epc_limit for calculation yday of start growing season  */
	if (!errorCode && scan_value(epc_file, &epc->GSI_limit_SGS, 'd'))
	{
		printf("ERROR reading epc_limit_SGS parameter: epc_init.c\n");
		errorCode=209126;
	}

	/* epc_limit for calculation yday of end growing season  */
	if (!errorCode && scan_value(epc_file, &epc->GSI_limit_EGS, 'd'))
	{
		printf("ERROR reading epc_limit_EGS parameter: epc_init.c\n");
		errorCode=209127;
	}

	/*********************************************************************************************************************/
	/* dividing line from file */ 
	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading 6. dividing line from, epc_init.c\n");
		errorCode=209128;
	}

	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading header, epc_init.c\n");
		errorCode=209129;
	}

	/*********************************************************************************************************************/
	/* PHENOLOGICAL (ALLOCATION) PARAMETERS (7 phenological phases) */
    /*********************************************************************************************************************/

	if (!errorCode && scan_value(epc_file, header, 's'))
	{
		printf("ERROR reading name of the phenophases\n");
		errorCode=209130;
	}

	scanflag=0; 
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		if (phenphase==N_PHENPHASES-1) scanflag=1;
		if (!errorCode && scan_array(epc_file, &(epc->phenophase_length[phenphase]), 'd', scanflag, 1))
		{
			printf("ERROR reading phenophase_length in phenophase %i, epc_init.c\n", phenphase+1);
			errorCode=209131;
		}
	}

	scanflag=0; 
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		if (phenphase==N_PHENPHASES-1) scanflag=1;
		if (!errorCode && scan_array(epc_file, &(epc->alloc_leafc[phenphase]), 'd', scanflag, 1))
		{
			printf("ERROR reading alloc_leafc in phenophase %i, epc_init.c\n", phenphase+1);
			errorCode=209132;
		}
	}

	scanflag=0; 
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		if (phenphase==N_PHENPHASES-1) scanflag=1;
		if (!errorCode && scan_array(epc_file, &(epc->alloc_frootc[phenphase]), 'd', scanflag, 1))
		{
			printf("ERROR reading alloc_frootc in phenophase %i, epc_init.c\n", phenphase+1);
			errorCode=209133;
		}
	}

	scanflag=0; 
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		if (phenphase==N_PHENPHASES-1) scanflag=1;
		if (!errorCode && scan_array(epc_file, &(epc->alloc_yield[phenphase]), 'd', scanflag, 1))
		{
			printf("ERROR reading alloc_yield in phenophase %i, epc_init.c\n", phenphase+1);
			errorCode=209134;
		}
	}

	scanflag=0; 
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		if (phenphase==N_PHENPHASES-1) scanflag=1;
		if (!errorCode && scan_array(epc_file, &(epc->alloc_softstemc[phenphase]), 'd', scanflag, 1))
		{
			printf("ERROR reading alloc_softstemc in phenophase %i, epc_init.c\n", phenphase+1);
			errorCode=209135;
		}
	}

	scanflag=0; 
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		if (phenphase==N_PHENPHASES-1) scanflag=1;
		if (!errorCode && scan_array(epc_file, &(epc->alloc_livestemc[phenphase]), 'd', scanflag, 1))
		{
			printf("ERROR reading alloc_woodystemc_live in phenophase %i, epc_init.c\n", phenphase+1);
			errorCode=209136;
		}
	}

	scanflag=0; 
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		if (phenphase==N_PHENPHASES-1) scanflag=1;
		if (!errorCode && scan_array(epc_file, &(epc->alloc_deadstemc[phenphase]), 'd', scanflag, 1))
		{
			printf("ERROR reading alloc_woodystemc_dead in phenophase %i, epc_init.c\n", phenphase+1);
			errorCode=209137;
		}
	}

	scanflag=0; 
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		if (phenphase==N_PHENPHASES-1) scanflag=1;
		if (!errorCode && scan_array(epc_file, &(epc->alloc_livecrootc[phenphase]), 'd', scanflag, 1))
		{
			printf("ERROR reading alloc_crootc_live in phenophase %i, epc_init.c\n", phenphase+1);
			errorCode=209138;
		}
	}

	scanflag=0; 
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		if (phenphase==N_PHENPHASES-1) scanflag=1;
		if (!errorCode && scan_array(epc_file, &(epc->alloc_deadcrootc[phenphase]), 'd', scanflag, 1))
		{
			printf("ERROR reading alloc_crootc_dead in phenophase %i, epc_init.c\n", phenphase+1);
			errorCode=209139;
		}
	}


	/* control of allocation parameter: test for alloc.fractions sum to 1.0  in every phenophases */
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		sum = epc->alloc_leafc[phenphase] + epc->alloc_frootc[phenphase] + epc->alloc_yield[phenphase] + epc->alloc_softstemc[phenphase] +
			  epc->alloc_livestemc[phenphase]   + epc->alloc_deadstemc[phenphase] + 
			  epc->alloc_livecrootc[phenphase]  + epc->alloc_deadcrootc[phenphase];

		/* control of allocation parameter settings */
		if (epc->alloc_leafc[phenphase]      < 1e-04)  epc->alloc_leafc[phenphase]      = 0;
		if (epc->alloc_frootc[phenphase]     < 1e-04)  epc->alloc_frootc[phenphase]     = 0;
		if (epc->alloc_yield[phenphase]     < 1e-04)   epc->alloc_yield[phenphase]     = 0;
		if (epc->alloc_softstemc[phenphase]  < 1e-04)  epc->alloc_softstemc[phenphase]  = 0;
		if (epc->alloc_livestemc[phenphase]  < 1e-04)  epc->alloc_livestemc[phenphase]  = 0;
		if (epc->alloc_deadstemc[phenphase]  < 1e-04)  epc->alloc_deadstemc[phenphase]  = 0;
		if (epc->alloc_livecrootc[phenphase] < 1e-04)  epc->alloc_livecrootc[phenphase] = 0;
		if (epc->alloc_deadcrootc[phenphase] < 1e-04)  epc->alloc_deadcrootc[phenphase] = 0;

		/*  searching of the maximal allocation */
		if (epc->alloc_leafc[phenphase]      > maxIND) maxIND=1;
		if (epc->alloc_frootc[phenphase]     > maxIND) maxIND=2;
		if (epc->alloc_yield[phenphase]     > maxIND) maxIND=3;
		if (epc->alloc_softstemc[phenphase]  > maxIND) maxIND=4;
		if (epc->alloc_livestemc[phenphase]  > maxIND) maxIND=5;
		if (epc->alloc_deadstemc[phenphase]  > maxIND) maxIND=6;
		if (epc->alloc_livecrootc[phenphase] > maxIND) maxIND=7;
		if (epc->alloc_deadcrootc[phenphase] > maxIND) maxIND=8;

		diff=fabs(sum - 1.);
		if (!errorCode && sum != 0 && diff > CRIT_PREC)
		{
			if (diff < 1e-04)
			{
				if (maxIND == 1) epc->alloc_leafc[phenphase]      -= diff;
				if (maxIND == 2) epc->alloc_frootc[phenphase]     -= diff;
				if (maxIND == 3) epc->alloc_yield[phenphase]     -= diff;
				if (maxIND == 4) epc->alloc_softstemc[phenphase]  -= diff;
				if (maxIND == 5) epc->alloc_livestemc[phenphase]  -= diff;
				if (maxIND == 6) epc->alloc_deadstemc[phenphase]  -= diff;
				if (maxIND == 7) epc->alloc_livecrootc[phenphase] -= diff;
				if (maxIND == 8) epc->alloc_deadcrootc[phenphase] -= diff;
				ctrl->allocControl_flag = 1;
			}
			else
			{
				printf("ERROR in allocation parameters in phenophase %i, epc_init.c\n", phenphase+1);
				printf("Allocation parameters must sum to 1.0 in every phenophase. Check EPC file and try again.\n");
				errorCode=2091391;
			}
		
		}

		/* control: after germination allocation parameter setting is necessary */
		if (!errorCode && (phenphase == 0)
			   && (epc->alloc_leafc[phenphase] == 0 || epc->alloc_frootc[phenphase] == 0))
		{
			printf("ERROR in allocation parameters in phenophase %i, epc_init.c\n", phenphase+1);
			printf("In the first phenophase leaf and froot allocation parameter setting is necessary. Check EPC file and try again.\n");
			errorCode=2091321;
		}

		/* control: woody biomes - no softstem */
		if (!errorCode && epc->woody && epc->alloc_softstemc[phenphase] > 0)
		{
			printf("ERROR in allocation parameters in phenophase %i, epc_init.c\n", phenphase+1);
			printf("No softstem allocation in case of woody biomass. Check EPC file and try again.\n");
			errorCode=2091361;
		}

		/* control: non-woody biomes - no woody allocation */
		if (!errorCode && !epc->woody && 
		   (epc->alloc_livecrootc[phenphase] > 0 || epc->alloc_livestemc[phenphase] > 0 || epc->alloc_deadcrootc[phenphase] > 0 || epc->alloc_deadstemc[phenphase] > 0))
		{
			printf("ERROR in allocation parameters in phenophase %i, epc_init.c\n", phenphase+1);
			printf("No woody allocation in case of non-woody biomass. Check EPC file and try again.\n");
			errorCode=2091381;
		}

	}

	scanflag=0; 
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		if (phenphase==N_PHENPHASES-1) scanflag=1;
		if (!errorCode && scan_array(epc_file, &(epc->avg_proj_sla[phenphase]), 'd', scanflag, 1))
		{
			printf("ERROR reading avg_proj_sla in phenophase %i, epc_init.c\n", phenphase+1);
			errorCode=209140;
		}

		if (!errorCode && epc->avg_proj_sla[phenphase] < 0)
		{
			printf("ERROR in specific leaf area in phenophase %i, epc_init.c\n", phenphase+1);
			printf("Specific leaf area must be positive in every phenophase. Check EPC file and try again.\n");
			errorCode=2091401;
		
		}

		if (!errorCode && epc->avg_proj_sla[phenphase] == 0 && epc->alloc_leafc[phenphase] > 0)
		{
			printf("ERROR in specific leaf area in phenophase %i, epc_init.c\n", phenphase+1);
			printf("Specific leaf area must greater than to 0.0 in case of leaf allocation. Check EPC file and try again.\n");
			errorCode=2091401;
		
		}
	}

	scanflag=0; 
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		if (phenphase==N_PHENPHASES-1) scanflag=1;
		if (!errorCode && scan_array(epc_file, &(epc->curgrowth_prop[phenphase]), 'd', scanflag, 1))
		{
			printf("ERROR reading curgrowth_prop in phenophase %i, epc_init.c\n", phenphase+1);
			errorCode=209141;
		}
	}

	scanflag=0; 
	for (phenphase=0; phenphase<N_PHENPHASES; phenphase++)
	{
		if (phenphase==N_PHENPHASES-1) scanflag=1;
		if (!errorCode && scan_array(epc_file, &(epc->max_plantlifetime[phenphase]), 'd', scanflag, 1))
		{
			printf("ERROR reading max_plantlifetime in phenophase %i, epc_init.c\n", phenphase+1);
			errorCode=209142;
		}


		if (!errorCode && ((epc->max_plantlifetime[0] == DATA_GAP && epc->max_plantlifetime[phenphase] != DATA_GAP) ||
			(epc->max_plantlifetime[0] != DATA_GAP && epc->max_plantlifetime[phenphase] == DATA_GAP)))
		{
			printf("ERROR in max_plantlifetime data: All or none maximal lifetime of plant tissue variables should be set at same time. Check EPC file and try again. \n");
			errorCode=2091421;
		}
	}

	

	free(OPTINyear_array);		

	/* -------------------------------------------*/
	if (dofilecloseEPC) fclose(epc_file.ptr);

		
	return (errorCode);

}
