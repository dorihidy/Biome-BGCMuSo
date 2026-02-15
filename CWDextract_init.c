/* 
CWDextract_init.c
read CWDextract information for pointbgc simulation

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


int CWDextract_init(file init, control_struct* ctrl, CWDextract_struct* CWE)
{
	
	char header[STRINGSIZE];
	char CWE_filename[STRINGSIZE];
	file CWE_file;

	int errorCode=0;
	int okFILE = 1;


	int mgmread;
	int nmgm = 0;

	int p1,p2,p3;
	double p4;
	char tempvar;

	int n_CWEparam, maxCWE_num;

	int* CWEyear_array;					
	int* CWEmonth_array;				
    int* CWEday_array;							
	double* removePROP_CWE;			      
		


	maxCWE_num=1000;


	/********************************************************************
	**                                                                 **
	** Begin reading initialization file block starting with keyword:  **
	** CWDextract                                                       ** 
	**                                                                 **
	********************************************************************/
	
	/* header reading */
	if (!errorCode && scan_value(init, header, 's'))
	{
		printf("ERROR reading header, CWDextract_init.c\n");
		errorCode=1;
	}
	

	/* keyword control */
	if (!errorCode && scan_value(init, header, 's'))
	{
		printf("ERROR reading header, CWDextract_init.c\n");
		errorCode=1;
	}
	
	/* number of management action */
	if (!errorCode && scan_value(init, &CWE->CWE_num, 'i'))
	{
		printf("ERROR reading number of CWDextract in CWDextract section\n");
		errorCode=1;
	}


	/* if CWE_num = 1 -> CWDextract */
	if (!errorCode && CWE->CWE_num)
	{
		/* allocate space for the MGM array */
		CWEyear_array         = (int*) malloc(maxCWE_num*sizeof(int));  
		CWEmonth_array        = (int*) malloc(maxCWE_num*sizeof(int)); 
		CWEday_array          = (int*) malloc(maxCWE_num*sizeof(int)); 
		removePROP_CWE        = (double*) malloc(maxCWE_num*sizeof(double)); 
		
		if (!errorCode && scan_value(init, CWE_filename, 's'))
		{
			printf("ERROR reading CWDextract calculating file\n");
			errorCode=1;
		}
		
		strcpy(CWE_file.name, CWE_filename);
		
		/* open the main init file for ascii read and check for errors */
		if (file_open(&CWE_file,'i',1))
		{
			printf("ERROR opening CWE_file, CWDextract_int.c\n");
			errorCode=1;
			okFILE=0;
		}

		if (!errorCode && scan_value(CWE_file, header, 's'))
		{
			printf("ERROR reading header for CWDextract section in MANAGMENET file\n");
			errorCode=1;
		}

	
		while (!errorCode && !(mgmread = scan_array (CWE_file, &p1, 'i', 0, 0)))
		{
			n_CWEparam = 5;

			mgmread = fscanf(CWE_file.ptr, "%c%d%c%d%lf%*[^\n]",&tempvar,&p2,&tempvar,&p3,&p4);

			if (mgmread != n_CWEparam)
			{
				printf("ERROR reading CWDextract parameters from CWDextract file\n");
				errorCode=1;
			}
			

			if (p1 >= ctrl->simstartyear && p1 < ctrl->simstartyear + ctrl->simyears)
			{
				CWEyear_array[nmgm]         = p1;
				CWEmonth_array[nmgm]        = p2;
				CWEday_array[nmgm]          = p3;
				removePROP_CWE[nmgm]        = p4;

				nmgm += 1;

			}
		}

		CWE->CWE_num = nmgm;
		nmgm = 0;
			

		/* allocate space for the temporary MGM array */
		CWE->CWEyear_array         = (int*) malloc(CWE->CWE_num*sizeof(int));  
		CWE->CWEmonth_array        = (int*) malloc(CWE->CWE_num*sizeof(int)); 
		CWE->CWEday_array          = (int*) malloc(CWE->CWE_num*sizeof(int)); 
		CWE->removePROP_CWE        = (double*) malloc(CWE->CWE_num*sizeof(double)); 
	
		for (nmgm = 0; nmgm < CWE->CWE_num; nmgm++)
		{

			CWE->CWEyear_array[nmgm]         = CWEyear_array[nmgm];
			CWE->CWEmonth_array[nmgm]        = CWEmonth_array[nmgm];
			CWE->CWEday_array[nmgm]          = CWEday_array[nmgm];
			CWE->removePROP_CWE[nmgm]		 = removePROP_CWE[nmgm];
		}

		/* close CWDextract file  and free temporary memory*/
		if (okFILE) fclose (CWE_file.ptr);

		free(CWEyear_array);
		free(CWEmonth_array);
		free(CWEday_array);
		free(removePROP_CWE);
	}
	else
	{
		/* reading the line of management file into a temporary variable */
		if (!errorCode && scan_value(init, header, 's'))
		{
			printf("ERROR reading line of management file (in case of no management)\n");
			errorCode=1;
		}
	}


	
	
	CWE->mgmdCWE = 0;

	
	return (errorCode);
 }
