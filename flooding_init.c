/* 
flooding_init.c
read flooding depth information if it is available

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
#include <time.h>
#include <math.h>
#include "ini.h"
#include "bgc_struct.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"
#include "bgc_func.h"


int flooding_init(flooding_struct* FLS, control_struct* ctrl)
{
	int errorCode=0;
	file FL_file;	
	char header[STRINGSIZE];

	int dataread, leap, n_FLparam;
	int ndata = 0;

	int p1,p2,p3,p4,p5,p6, maxFLnum, nmgm;
	double p7, p8, p9, p10, p11, p12, p13, p14, p15, p16,p17;
	char tempvar;

	int* FLstart_year_array;			
	int* FLstart_month_array;						
	int* FLstart_day_array;
	int* FLend_year_array;			
	int* FLend_month_array;						
	int* FLend_day_array;
	double* FLheight_array;
	double* FL_NH4ppm_array;
	double* FL_NO3ppm_array;
	double* FL_DON1ppm_array;
	double* FL_DON2ppm_array;
	double* FL_DON3ppm_array;
	double* FL_DON4ppm_array;
	double* FL_DOC1ppm_array;
	double* FL_DOC2ppm_array;
	double* FL_DOC3ppm_array;
	double* FL_DOC4ppm_array;


	int* mondays=0;
	int* enddays=0;

	nmgm=0;
	maxFLnum=ctrl->simyears* nDAYS_OF_YEAR;
	FLS->FLnum = 0;

	/********************************************************************
	**                                                                 **
	** Reading FL data if available                                   ** 
	**                                                                 **
	********************************************************************/

	if (!errorCode)
	{
		if (ctrl->spinup == 0)   /* normal run */
		{
			strcpy(FL_file.name, "flooding_normal.txt");
			if (!file_open(&FL_file,'j',1)) FLS->FLnum = 1;
		}
		else                     /* spinup and transient run */        
		{ 	
			strcpy(FL_file.name, "flooding_transient.txt");
			if (!file_open(&FL_file,'j',1)) FLS->FLnum = 1;	
		}
	}


	if (!errorCode && FLS->FLnum > 0) 
	{		

		if (!errorCode && scan_value(FL_file, header, 's'))
		{
			printf("ERROR reading header for FLOODING file\n");
			errorCode = 1;
		}
		
		if (!errorCode) 
		{
			enddays = (int*) malloc(nMONTHS_OF_YEAR * sizeof(int));
			if (!enddays)
			{
				printf("ERROR allocating for enddays in bgc.c\n");
				errorCode=1;
			}
		}

		if (!errorCode) 
		{
			mondays = (int*) malloc(nMONTHS_OF_YEAR * sizeof(int));
			if (!mondays)
			{
				printf("ERROR allocating for enddays in bgc.c\n");
				errorCode=1;
			}
		}

		/* allocate space for the temporary MGM array */
		FLstart_year_array  = (int*) malloc(maxFLnum*sizeof(int));  
		FLstart_month_array = (int*) malloc(maxFLnum*sizeof(int)); 
		FLstart_day_array   = (int*) malloc(maxFLnum*sizeof(int)); 
		FLend_year_array    = (int*) malloc(maxFLnum*sizeof(int));  
		FLend_month_array   = (int*) malloc(maxFLnum*sizeof(int)); 
		FLend_day_array     = (int*) malloc(maxFLnum*sizeof(int)); 
        FLheight_array      = (double*) malloc(maxFLnum*sizeof(double)); 
		FL_NH4ppm_array     = (double*) malloc(maxFLnum * sizeof(double));
		FL_NO3ppm_array     = (double*) malloc(maxFLnum * sizeof(double));
		FL_DON1ppm_array = (double*)malloc(maxFLnum * sizeof(double));
		FL_DON2ppm_array = (double*)malloc(maxFLnum * sizeof(double));
		FL_DON3ppm_array = (double*)malloc(maxFLnum * sizeof(double));
		FL_DON4ppm_array = (double*)malloc(maxFLnum * sizeof(double));
		FL_DOC1ppm_array = (double*)malloc(maxFLnum * sizeof(double));
		FL_DOC2ppm_array = (double*)malloc(maxFLnum * sizeof(double));
		FL_DOC3ppm_array = (double*)malloc(maxFLnum * sizeof(double));
		FL_DOC4ppm_array = (double*)malloc(maxFLnum * sizeof(double));


		
		ndata=0;
		while (!errorCode && !(dataread = scan_array (FL_file, &p1, 'i', 0, 0)))
		{
			n_FLparam = 20;

			dataread = fscanf(FL_file.ptr, "%c %d %c %d %d %c %d %c %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf[^\n]", &tempvar,&p2,&tempvar,&p3,&p4,&tempvar,&p5,&tempvar,&p6,&p7,&p8,&p9,&p10,&p11,&p12,&p13,&p14,&p15,&p16,&p17);
				
			if (dataread != n_FLparam)
			{
				printf("ERROR reading FLOODING data from flooding file  file\n");
				errorCode = 1;
			}

			if (p1 >= ctrl->simstartyear && p1 < ctrl->simstartyear + ctrl->simyears)
			{
				FLstart_year_array[ndata]     = p1;
				FLstart_month_array[ndata]    = p2;
				FLstart_day_array[ndata]      = p3;
				FLend_year_array[ndata]       = p4;
				FLend_month_array[ndata]      = p5;
				FLend_day_array[ndata]        = p6;
				FLheight_array[ndata]         = p7;
				FL_NH4ppm_array[ndata]        = p8;
				FL_NO3ppm_array[ndata]        = p9;
				FL_DON1ppm_array[ndata] = p10;
				FL_DON2ppm_array[ndata] = p11;
				FL_DON3ppm_array[ndata] = p12;
				FL_DON4ppm_array[ndata] = p13;
				FL_DOC1ppm_array[ndata] = p14;
				FL_DOC2ppm_array[ndata] = p15;
				FL_DOC3ppm_array[ndata] = p16;
				FL_DOC4ppm_array[ndata] = p17;

				if (!errorCode && leapControl(FLstart_year_array[ndata], enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from flooding_init.c\n");
					errorCode=1;
				}
				if (leap == 1 && FLstart_month_array[ndata] == 12 && FLstart_day_array[ndata] == 31)
				{
					printf("ERROR in flooding date in flooding_init.c: data from 31 December in a leap year is found in flooding file\n");
					printf("Please read the manual and modify the input data\n");
					errorCode=1;
				}

				if (!errorCode && leapControl(FLend_year_array[ndata], enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from flooding_init.c\n");
					errorCode=1;
				}
				if (leap == 1 && FLend_month_array[ndata] == 12 && FLend_day_array[ndata] == 31)
				{
					printf("ERROR in flooding date in flooding_init.c: data from 31 December in a leap year is found in flooding file\n");
					printf("Please read the manual and modify the input data\n");
					errorCode=1;
				}

                nmgm += 1;
				ndata += 1;
			}
		}

		FLS->FLnum = nmgm;
		nmgm = 0;
	
		FLS->FLstart_year_array      = (int*) malloc(FLS->FLnum*sizeof(double));  
		FLS->FLstart_month_array     = (int*) malloc(FLS->FLnum*sizeof(double)); 
		FLS->FLstart_day_array       = (int*) malloc(FLS->FLnum*sizeof(double)); 
		FLS->FLend_year_array        = (int*) malloc(FLS->FLnum*sizeof(double));  
		FLS->FLend_month_array       = (int*) malloc(FLS->FLnum*sizeof(double)); 
		FLS->FLend_day_array         = (int*) malloc(FLS->FLnum*sizeof(double)); 
		FLS->FLheight_array          = (double*) malloc(FLS->FLnum*sizeof(double)); 
		FLS->FL_NH4ppm_array         = (double*) malloc(FLS->FLnum * sizeof(double));
		FLS->FL_NO3ppm_array         = (double*) malloc(FLS->FLnum * sizeof(double));
		FLS->FL_DON1ppm_array = (double*)malloc(FLS->FLnum * sizeof(double));
		FLS->FL_DON2ppm_array = (double*)malloc(FLS->FLnum * sizeof(double));
		FLS->FL_DON3ppm_array = (double*)malloc(FLS->FLnum * sizeof(double));
		FLS->FL_DON4ppm_array = (double*)malloc(FLS->FLnum * sizeof(double));
		FLS->FL_DOC1ppm_array = (double*)malloc(FLS->FLnum * sizeof(double));
		FLS->FL_DOC2ppm_array = (double*)malloc(FLS->FLnum * sizeof(double));
		FLS->FL_DOC3ppm_array = (double*)malloc(FLS->FLnum * sizeof(double));
		FLS->FL_DOC4ppm_array = (double*)malloc(FLS->FLnum * sizeof(double));

		for (nmgm = 0; nmgm < FLS->FLnum; nmgm++)
		{
			FLS->FLstart_year_array[nmgm]    = FLstart_year_array[nmgm];
			FLS->FLstart_month_array[nmgm]   = FLstart_month_array[nmgm];
			FLS->FLstart_day_array[nmgm]     = FLstart_day_array[nmgm];

			FLS->FLend_year_array[nmgm]      = FLend_year_array[nmgm];
			FLS->FLend_month_array[nmgm]     = FLend_month_array[nmgm];
			FLS->FLend_day_array[nmgm]       = FLend_day_array[nmgm];

			FLS->FLheight_array[nmgm]        = FLheight_array[nmgm];
			FLS->FL_NH4ppm_array[nmgm]       = FL_NH4ppm_array[nmgm];
			FLS->FL_NO3ppm_array[nmgm]       = FL_NO3ppm_array[nmgm];
			FLS->FL_DON1ppm_array[nmgm] = FL_DON1ppm_array[nmgm];
			FLS->FL_DON2ppm_array[nmgm] = FL_DON2ppm_array[nmgm];
			FLS->FL_DON3ppm_array[nmgm] = FL_DON3ppm_array[nmgm];
			FLS->FL_DON4ppm_array[nmgm] = FL_DON4ppm_array[nmgm];
			FLS->FL_DOC1ppm_array[nmgm] = FL_DOC1ppm_array[nmgm];
			FLS->FL_DOC2ppm_array[nmgm] = FL_DOC2ppm_array[nmgm];
			FLS->FL_DOC3ppm_array[nmgm] = FL_DOC3ppm_array[nmgm];
			FLS->FL_DOC4ppm_array[nmgm] = FL_DOC4ppm_array[nmgm];
		}

		if (nmgm > maxFLnum)
		{
			printf("ERROR in flooding data reading flooding_init.c\n");
			errorCode=1;
		}

		/* read year and FLD for each simday in each simyear */
		
		free(FLstart_year_array);	
		free(FLstart_month_array);	
		free(FLstart_day_array);	
		free(FLend_year_array);	
		free(FLend_month_array);	
		free(FLend_day_array);	
        free(FLheight_array);	
		free(FL_NH4ppm_array);
		free(FL_NO3ppm_array);
		free(FL_DON1ppm_array);
		free(FL_DON2ppm_array);
		free(FL_DON3ppm_array);
		free(FL_DON4ppm_array);
		free(FL_DOC1ppm_array);
		free(FL_DOC2ppm_array);
		free(FL_DOC3ppm_array);
		free(FL_DOC4ppm_array);

		fclose(FL_file.ptr);
	}	

	

	FLS->mgmdFL = 0;

	if (!errorCode && FLS->FLnum > 0) 
	{
		free(enddays);
		free(mondays);
	}
		
	return (errorCode);
}
