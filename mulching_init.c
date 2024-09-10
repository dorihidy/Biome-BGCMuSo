/* 
MULCHING_init.c
read MULCHING information for pointbgc simulation

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


int mulching_init(file init, control_struct* ctrl, mulching_struct* MUL)
{
	
	char header[STRINGSIZE];
	char MUL_filename[STRINGSIZE];
	file MUL_file;



	int errorCode=0;
	int okFILE = 1;

	

	int mgmread;
	int nmgm = 0;

	int p1,p2,p3;
	double p4,p5,p6,p7;
	char tempvar;

	int n_MULparam, maxMUL_num;

	int* MULyear_array;					
	int* MULmonth_array;				
    int* MULday_array;							
	double* litrCabove_MUL;			      
	double* litrCNabove_MUL;				    
	double* cwdCabove_MUL;				
	double* cwdCNabove_MUL;				

	

	maxMUL_num=1000;


	/********************************************************************
	**                                                                 **
	** Begin reading initialization file block starting with keyword:  **
	** MULCHING                                                       ** 
	**                                                                 **
	********************************************************************/
	
	/* header reading */
	if (!errorCode && scan_value(init, header, 's'))
	{
		printf("ERROR reading header, MULCHING_init.c\n");
		errorCode=1;
	}
	

	/* keyword control */
	if (!errorCode && scan_value(init, header, 's'))
	{
		printf("ERROR reading header, MULCHING_init.c\n");
		errorCode=1;
	}
	
	/* number of management action */
	if (!errorCode && scan_value(init, &MUL->MUL_num, 'i'))
	{
		printf("ERROR reading number of MULCHING in MULCHING section\n");
		errorCode=1;
	}


	/* if MUL_num = 1 -> MULCHING */
	if (!errorCode && MUL->MUL_num)
	{
		/* allocate space for the MGM array */
		MULyear_array         = (int*) malloc(maxMUL_num*sizeof(int));  
		MULmonth_array        = (int*) malloc(maxMUL_num*sizeof(int)); 
		MULday_array          = (int*) malloc(maxMUL_num*sizeof(int)); 
		litrCabove_MUL        = (double*) malloc(maxMUL_num*sizeof(double)); 
		litrCNabove_MUL       = (double*) malloc(maxMUL_num*sizeof(double)); 
		cwdCabove_MUL         = (double*) malloc(maxMUL_num*sizeof(double)); 
		cwdCNabove_MUL        = (double*) malloc(maxMUL_num*sizeof(double)); 
		
		if (!errorCode && scan_value(init, MUL_filename, 's'))
		{
			printf("ERROR reading MULCHING calculating file\n");
			errorCode=1;
		}
		
		strcpy(MUL_file.name, MUL_filename);
		
		/* open the main init file for ascii read and check for errors */
		if (file_open(&MUL_file,'i',1))
		{
			printf("ERROR opening MUL_file, MULCHING_int.c\n");
			errorCode=1;
			okFILE=0;
		}

		if (!errorCode && scan_value(MUL_file, header, 's'))
		{
			printf("ERROR reading header for MULCHING section in MANAGMENET file\n");
			errorCode=1;
		}

	
		while (!errorCode && !(mgmread = scan_array (MUL_file, &p1, 'i', 0, 0)))
		{
			n_MULparam = 8;

			mgmread = fscanf(MUL_file.ptr, "%c%d%c%d%lf%lf%lf%lf%*[^\n]",&tempvar,&p2,&tempvar,&p3,&p4,&p5,&p6,&p7);

			if (mgmread != n_MULparam)
			{
				printf("ERROR reading MULCHING parameters from MULCHING file\n");
				errorCode=1;
			}
			

			if (p1 >= ctrl->simstartyear && p1 < ctrl->simstartyear + ctrl->simyears)
			{
				MULyear_array[nmgm]         = p1;
				MULmonth_array[nmgm]        = p2;
				MULday_array[nmgm]          = p3;
				litrCabove_MUL[nmgm]        = p4;
				litrCNabove_MUL[nmgm]       = p5;
				cwdCabove_MUL[nmgm]         = p6;
				cwdCNabove_MUL[nmgm]        = p7;

				nmgm += 1;

			}
		}

		MUL->MUL_num = nmgm;
		nmgm = 0;
			

		/* allocate space for the temporary MGM array */
		MUL->MULyear_array         = (int*) malloc(MUL->MUL_num*sizeof(int));  
		MUL->MULmonth_array        = (int*) malloc(MUL->MUL_num*sizeof(int)); 
		MUL->MULday_array          = (int*) malloc(MUL->MUL_num*sizeof(int)); 
		MUL->litrCabove_MUL        = (double*) malloc(MUL->MUL_num*sizeof(double)); 
		MUL->litrCNabove_MUL       = (double*) malloc(MUL->MUL_num*sizeof(double)); 
		MUL->cwdCabove_MUL         = (double*) malloc(MUL->MUL_num*sizeof(double)); 
		MUL->cwdCNabove_MUL        = (double*) malloc(MUL->MUL_num*sizeof(double)); 

		for (nmgm = 0; nmgm < MUL->MUL_num; nmgm++)
		{

			MUL->MULyear_array[nmgm]         = MULyear_array[nmgm];
			MUL->MULmonth_array[nmgm]        = MULmonth_array[nmgm];
			MUL->MULday_array[nmgm]          = MULday_array[nmgm];
			MUL->litrCabove_MUL[nmgm]		 = litrCabove_MUL[nmgm];
			MUL->litrCNabove_MUL[nmgm]       = litrCNabove_MUL[nmgm] ;
			MUL->cwdCabove_MUL[nmgm]         = cwdCabove_MUL[nmgm];
			MUL->cwdCNabove_MUL[nmgm]        = cwdCNabove_MUL[nmgm];	

			if (MUL->litrCabove_MUL[nmgm] > 0 && MUL->litrCNabove_MUL[nmgm] == 0)
			{
				printf("ERROR in mulching parameters: if non-woody carbon content is set, non-woody CNratio must be greater than 0\n");
				errorCode=1;
			}

			if (MUL->cwdCabove_MUL[nmgm] > 0 && MUL->cwdCNabove_MUL[nmgm] == 0)
			{
				printf("ERROR in mulching parameters: if non-woody carbon content is set, non-woody CNratio must be greater than 0\n");
				errorCode=1;
			}
		}

		/* close MULCHING file  and free temporary memory*/
		if (okFILE) fclose (MUL_file.ptr);

		free(MULyear_array);
		free(MULmonth_array);
		free(MULday_array);
		free(litrCabove_MUL);
		free(litrCNabove_MUL);
		free(cwdCabove_MUL);
		free(cwdCNabove_MUL);
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


	
	
	MUL->mgmdMUL = 0;

	
	return (errorCode);
 }
