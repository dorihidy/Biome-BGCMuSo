/*
prephenology.c
Initialize phenology arrays, called prior to annual loop in bgc.c

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2022, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"    
#include "bgc_func.h" 
#include "bgc_constants.h"
#include "misc_func.h"

int prephenology(file logfile, const epconst_struct* epc, const metarr_struct* metarr, 
	             planting_struct* PLT, harvesting_struct* HRV, control_struct* ctrl, phenarray_struct* phenarr)
{
	int errorCode=0;
	
	int model,woody,evergreen;
	int countONOFFDAY, leap;
	double t1;
	int i,pday,ndays,py;
	int nyears,phenyears;
	int onday,offday;
	/* phenology model variables */
	int *onday_arr = 0;
	int *offday_arr = 0;
	int fall_Tavg_count;
	int onset_day, offset_day;
	double mean_Tavg,fall_Tavg;
	double phensoilt,phendayl;
	double onset_critsum, sum_soilt;
	double criTdayl = 39300.0; /* seconds */
	/* grass model parameters */
	double ann_prcp;
	double sum_prcp, phenprcp;
	double grass_stsumcrit;
	double grass_prcpcrit;
	double grass_stsummax = 1380.0;
	double grass_stsummid = 900.0;
	double grass_stsummin = 418.0;
	double grass_a = 32.9;
	double grass_k = 0.15;
	double grass_tmid = 9.0;
	double grass_prcpyear[nDAYS_OF_YEAR];
	double grass_prcpprevcrit = 1.14;
	double grass_prcpprev;
	double grass_prcpnextcrit = 0.97;
	double grass_prcpnext;
	double grass_Tmaxyear[nDAYS_OF_YEAR];
	double grass_Tminyear[nDAYS_OF_YEAR];
	double grass_3dayTmin[nDAYS_OF_YEAR];
	int psum_starTday, psum_stopday;
	double Tmax_ann, Tmax, new_Tmax;
	double Tmin_annavg;

	int onday_min, onday_max, offday_min, offday_max;
	int PLTyday, HRVyday, PLTyear, HRVyear;
	int* mondays=0;
	int* enddays=0;


	onday = offday = PLTyday = HRVyday = PLTyear = HRVyear = 0;

	/* in case of natural ecosystem, nyears = number of meteorological year, in case agricultural system: nyears = number of plantings (plus fallow years) */
	if (!PLT->PLT_num)
		nyears = ctrl->simyears;
	else
		nyears = PLT->PLT_num + ctrl->simyears;

	ndays = nDAYS_OF_YEAR * nyears;

	onday_min = offday_min = 364;
	onday_max = offday_max = 0;
	

	if (!errorCode)
	{
		onday_arr = (int*) malloc((nyears+1) * sizeof(int));
		if (!onday_arr)
		{
			printf("ERROR allocating for onday_arr, prephenology.c\n");
			errorCode=1;
		}
	}

	if (!errorCode)
	{
		offday_arr = (int*) malloc((nyears+1) * sizeof(int));
		if (!offday_arr)
		{
			printf("ERROR allocating for offday_arr, prephenology.c\n");
			errorCode=1;
		}
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

	/* set some local flags to control the phenology model behavior */
	/* model=1 --> use phenology model   model=0 --> user specified phenology */
	/* woody=1 --> woody veg type        woody=0 --> non-woody veg type */
	/* evergreen=1 --> evergreen type    evergreen=0 --> deciduous type */
	/* south=1 --> southern hemisphere   south=0 --> northern hemisphere */

	model = ctrl->phenology_flag;
	
	woody = epc->woody;
	evergreen = epc->evergreen;
	
	
	/* before MuSo 6.0.4: for southern hemisphere sites, use an extra phenology year
	if (south) phenyears = nyears+1;
	else phenyears = nyears;  */
	phenyears = nyears;

	fprintf(logfile.ptr, "INFORMATION ABOUT SGS AND EGS VALUES (yday of onday and offday)\n");



	/************************************************************************************************************/
	/* I. CONSTANT onday and offday
	      define the phenology signals for cases in which the phenology signals are constant between years */
	/************************************************************************************************************/

	
	if (evergreen || !model || PLT->PLT_num)
	{
		
		if (PLT->PLT_num)
		{
			/* onday and offday values from yearly variable from  management data */
			onday_min = 364;
			onday_max = 0;
			offday_min = 364;
			offday_max = 0;
			countONOFFDAY = 0;
			
			for (py=0 ; py<PLT->PLT_num; py++)
			{
				PLTyear = PLT->PLTyear_array[countONOFFDAY];
				if (!errorCode && leapControl(PLTyear, enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from prephenology_init.c\n");
					errorCode=1;
				}
				PLTyday = date_to_doy(mondays, PLT->PLTmonth_array[countONOFFDAY], PLT->PLTday_array[countONOFFDAY]);
				
				
				HRVyear = HRV->HRVyear_array[countONOFFDAY];
				if (!errorCode && leapControl(HRVyear, enddays, mondays, &leap))
				{
					printf("ERROR in call to leapControl.c from prephenology_init.c\n");
					errorCode=1;
				}
				HRVyday = date_to_doy(mondays, HRV->HRVmonth_array[countONOFFDAY], HRV->HRVday_array[countONOFFDAY]);
					
				
				onday_arr[countONOFFDAY] = PLTyday;
				offday_arr[countONOFFDAY] = HRVyday;
						
				phenarr->onday_arr[countONOFFDAY][0] = PLTyear;
				phenarr->onday_arr[countONOFFDAY][1] = onday_arr[countONOFFDAY];

				phenarr->offday_arr[countONOFFDAY][0] = HRVyear;
				phenarr->offday_arr[countONOFFDAY][1] = offday_arr[countONOFFDAY];

				if (onday_arr[countONOFFDAY]  < onday_min)  onday_min  = onday_arr[countONOFFDAY];
				if (onday_arr[countONOFFDAY]  > onday_max)  onday_max  = onday_arr[countONOFFDAY];
						

				if (offday_arr[countONOFFDAY]  < offday_min)  offday_min  = offday_arr[countONOFFDAY];
				if (offday_arr[countONOFFDAY]  > offday_max)  offday_max  = offday_arr[countONOFFDAY];
				
				countONOFFDAY += 1;

			} 
			

		}/* endfor mgm */
		else
		{
		
			/*  2. define USER-DEFINED onday and offdays (from EPC file or SGS/EGS file) */
			for (py=0 ; py<phenyears ; py++)
			{
	
				if (ctrl->varSGS_flag)
					onday = (int) epc->SGS_array[py];
				else
					onday = epc->ondayUSER;
			
				if (ctrl->varEGS_flag)
					offday = (int) epc->EGS_array[py];
				else
					offday = epc->offdayUSER;

				if (evergreen)
				{
					onday=0;
					offday=364;
				}

				if (onday == DATA_GAP && offday == DATA_GAP) ctrl->bareground_flag = 1;


				phenarr->onday_arr[py][0]  = ctrl->simstartyear+py;
				phenarr->offday_arr[py][0] = ctrl->simstartyear+py;
				
				phenarr->onday_arr[py][1]  = onday;
				phenarr->offday_arr[py][1] = offday;
			}
			/* for wrinting out log file */
			if (onday  < onday_min  && onday != DATA_GAP)  onday_min  = onday;
			if (onday  > onday_max  && onday != DATA_GAP)  onday_max  = onday;
			if (offday < offday_min && offday != DATA_GAP) offday_min = offday;
			if (offday > offday_max && offday != DATA_GAP) offday_max = offday;

		
		} /* end of user defined on and off days  */

		
		
		/* 3. wrinting log file */
		fprintf(logfile.ptr, "SGS value (min and max): %6i %6i\n", onday_min, onday_max);
		fprintf(logfile.ptr, "EGS value (min and max): %6i %6i\n", offday_min, offday_max);
		fprintf(logfile.ptr, " \n");

	} /* end if constant phenological signals */

	/************************************************************************************************************/
	/* II. MODEL-CALCULATED onday and offday */
	/************************************************************************************************************/
	else
	{
		/*  1. define MODEL-DEFINED onday and offdays (original or GSI-method */
		/**************************************************************/
		/* A: ORIGINAL METHOD for calculating onset and offset day */
		/**************************************************************/
		if (!ctrl->GSI_flag)
		{
			/* Cases that have variable phenological signals between years */
			/* Use the phenology model described in White et al., 1997 */
			/* the two cases that make it to this block are:
			model, deciduous, woody   and
			model, deciduous, non-woody (grass), which are the two cases
			handled by the White et al. paper */
			if (woody)
			{
				/* use DECIDUOUS TREE PHENOLOGY MODEL */
				/* loop through the entire Tavg timeseries to calculate long-term
				average Tavg */
				mean_Tavg = 0.0;
				for (i=0 ; i<ndays ; i++)
				{
					mean_Tavg += (metarr->Tmax_array[i] + metarr->Tmin_array[i])/2.;
				}
				mean_Tavg /= (double)ndays;
				/* tree onset equation from Mike White, Aug. 1997 */
				onset_critsum = exp(4.795 + 0.129*mean_Tavg);
				
				/* now go through the phenological years and generate expansionand litterfall arrays. Some complications for Southern hemisphere sites... */
				/* calculate fall_Tavg, the mean Tavg from phenyday 244-304 */
				fall_Tavg = 0.0;
				fall_Tavg_count = 0;
				for (py=0 ; py<phenyears ; py++)
				{
					for (pday=244 ; pday<305 ; pday++)
					{
				
						phensoilt = metarr->TavgRA11_array[py*nDAYS_OF_YEAR+pday];
						
						fall_Tavg += phensoilt;
						fall_Tavg_count++;
						
					} /* end pday loop */
				} /* end py loop */
				fall_Tavg /= (double)fall_Tavg_count;
				
				/* loop through phenyears again, fill onset and offset arrays */
				for (py=0 ; py<phenyears ; py++)
				{
					sum_soilt = 0.0;
					onset_day = offset_day = -1;
					for (pday=0 ; pday<nDAYS_OF_YEAR ; pday++)
					{
					
						phensoilt = metarr->TavgRA11_array[py*nDAYS_OF_YEAR+pday];
						phendayl = metarr->dayl_array[py*nDAYS_OF_YEAR+pday];
						
						/* tree onset test */
						if (onset_day == -1)
						{
							if (phensoilt > 0.0) sum_soilt += phensoilt;
							if (sum_soilt >= onset_critsum) onset_day = pday;
						}
						
						/* tree offset test */
						if (onset_day != -1 && offset_day == -1)
						{
							if ((pday > 182) && (((phendayl <= criTdayl) && (phensoilt <= fall_Tavg)) || (phensoilt <= 2.0)))
							{
								offset_day = pday;
							}
						}
						
					} /* end pday loop */
					
					/* now do some exception handling for this year's phenology */
					if (onset_day != -1)
					{
						/* leaves are turned on sometime this year */
						/* subtract 15 days from onset day to approximate the
						start of the new growth period, instead of the middle of
						the new growth period, as is used in the White et al. ms. */
						if (onset_day >= 15)
						{
							onset_day -= 15;
							/*onset_day-=1;*/
						}
						else onset_day = 0;

						/* if leaves never got turned off, force off on last day */
						if (offset_day == -1) offset_day = 364;
						/* add 15 days to offset day to approximate the
						end of the litterfall period, instead of the middle
						as in the White et al. ms. */
						if (offset_day <= 349)
						{
							offset_day += 15;
							/*offset_day+=1;*/
						}
						else offset_day = 364;
						
						/* force onset and offset to be at least one day apart */
						if (onset_day == offset_day)
						{
							if (onset_day > 0) onset_day--;
							else offset_day++;
						}
					}
					else
					{
						/* leaves never got turned on, this is a non-growth
						year.  This probably indicates a problem with the
						phenology model */
						onset_day = -1;
						offset_day = -1;
					}

					/* save these onset and offset days and go to the next
					phenological year */
					onday_arr[py] = onset_day;
					offday_arr[py] = offset_day;

			
					
				} /* end phenyears loop for filling onset and offset arrays */
			} /* end if woody (tree phenology model) */
			else
			{
			

				/* non-woody, use the grass PHENOLOGY MODEL to calculate the
				array of onset and offset days */
				/* loop through the entire Tavg timeseries to calculate long-term
				average Tavg and long-term average annual total precip */
				mean_Tavg = 0.0;
				ann_prcp = 0.0;
				for (i=0 ; i<ndays ; i++)
				{
					mean_Tavg += metarr->Tday_array[i];
					ann_prcp += metarr->prcp_array[i];
				}
				mean_Tavg /= (double)ndays;
				ann_prcp /= (double)ndays / nDAYS_OF_YEAR;
				
				/* grass onset equation from White et al., 1997, with parameter
				values specified by Mike White, Aug. 1997 */
				t1 = exp(grass_a * (mean_Tavg - grass_tmid));
				grass_stsumcrit = ((grass_stsummax - grass_stsummin)* 0.5 *
					((t1-1)/(t1+1))) + grass_stsummid;
				grass_prcpcrit = ann_prcp * grass_k;
				
				/* now go through the phenological years and generate onset
				and offset days */
				
				/* calculate the long-term average annual high temperature
				for use in grass offset prediction */
				Tmax_ann = 0.0;
				Tmin_annavg = 0.0;
				for (py=0 ; py<phenyears ; py++)
				{
					new_Tmax = -1000.0;
					for (pday=0 ; pday<nDAYS_OF_YEAR ; pday++)
					{
						Tmax = metarr->Tmax_array[py*nDAYS_OF_YEAR+pday];
						Tmin_annavg += metarr->Tmin_array[py*nDAYS_OF_YEAR+pday];
						
						if (Tmax > new_Tmax) new_Tmax = Tmax;
						
					} /* end pday loop */
					
					Tmax_ann += new_Tmax;
				} /* end py loop */
				Tmax_ann /= (double) phenyears;
				/* 92% of Tmax_ann is the threshold used in grass offset below */
				Tmax_ann *= 0.92;
				Tmin_annavg /= (double) phenyears * nDAYS_OF_YEAR;
				
				/* loop through phenyears again, fill onset and offset arrays */
				for (py=0 ; py<phenyears ; py++)
				{
					sum_soilt = 0.0;
					sum_prcp = 0.0;
					onset_day = offset_day = -1;
					for (pday=0 ; pday<nDAYS_OF_YEAR ; pday++)
					{
					
						phensoilt = metarr->TavgRA11_array[py*nDAYS_OF_YEAR+pday];
						phenprcp = metarr->prcp_array[py*nDAYS_OF_YEAR+pday];
						grass_prcpyear[pday] = phenprcp;
						grass_Tminyear[pday] = metarr->Tmin_array[py*nDAYS_OF_YEAR+pday];
						grass_Tmaxyear[pday] = metarr->Tmax_array[py*nDAYS_OF_YEAR+pday];
						
						/* grass onset test */
						if (onset_day == -1)
						{
							if (phensoilt > 0.0) sum_soilt += phensoilt;
							sum_prcp += phenprcp;
							if (sum_soilt >= grass_stsumcrit &&
								sum_prcp >= grass_prcpcrit) onset_day = pday;
						}

						
					} /* end pday loop */
					
					/* do averaging operations on grass_prcpyear and grass_Tminyear,
					and do tests for offset day. Offset due to hot & dry can't
					happen within one month after the onset day, and offset due
					to cold can't happen before midyear (yearday 182) */
					if (onset_day != -1)
					{
						/* calculate three-day boxcar average of Tmin */
						if (boxcar_smooth(grass_Tminyear, grass_3dayTmin, nDAYS_OF_YEAR,3,0))
						{
							printf("ERROR in prephenology.c call to boxcar.c\n");
							errorCode=1;
						}
						
						for (pday=onset_day+30 ; pday<nDAYS_OF_YEAR ; pday++)
						{
							/* calculate the previous 31-day prcp total */
							psum_starTday = pday - 30;
							grass_prcpprev = 0.0;
							for (i=psum_starTday ; i<=pday ; i++)
							{
								grass_prcpprev += grass_prcpyear[i];
							}

							/* calculate the next 7-day prcp total */
							if (pday > 358) psum_stopday = 364;
							else psum_stopday = pday + 6;
							grass_prcpnext = 0.0;
							for (i=pday ; i<=psum_stopday ; i++)
							{
								grass_prcpnext += grass_prcpyear[i];
							}
							
							/* test for hot and dry conditions */
							if (offset_day == -1)
							{
								if (grass_prcpprev < grass_prcpprevcrit && 
									grass_prcpnext < grass_prcpnextcrit &&
									grass_Tmaxyear[pday] > Tmax_ann)
									offset_day = pday;
							}
							
							/* test for cold offset condition */
							if (offset_day == -1)
							{
								if (pday > 182 &&
									grass_3dayTmin[pday] <= Tmin_annavg)
									offset_day = pday;
							}
							
						} /* end of pdays loop for grass offset testing */
					} /* end of if onset_day != -1 block */
						
					/* now do some exception handling for this year's phenology */
					if (onset_day != -1)
					{
						/* leaves are turned on sometime this year */
						/* subtract 15 days from onset day to approximate the
						start of the new growth period, instead of the middle of
						the new growth period, as is used in the White et al. ms. */
						if (onset_day >= 15)
						{
							onset_day -= 15;
						}
						else onset_day = 0;

						/* if leaves never got turned off, force off on last day */
						if (offset_day == -1) offset_day = 364;

						/* force onset and offset to be at least one day apart */
						if (onset_day == offset_day)
						{
							if (onset_day > 0) onset_day--;
							else offset_day++;
						}
					}
					else
					{
						/* leaves never got turned on, this is a non-growth
						year.  This probably indicates a problem with the
						phenology model */
						onset_day = -1;
						offset_day = -1;
					}
					
					/* save these onset and offset days and go to the next
					phenological year */
					onday_arr[py] = onset_day;
					offday_arr[py] = offset_day;
					
				} 	/* end phenyears loop for filling onset and offset arrays */
			} /* end else !woody (grass phenology model) */

			/* fix the onday and offday value */
			for (py=0 ; py<phenyears ; py++)
			{
				phenarr->onday_arr[py][0]  = ctrl->simstartyear+py;
				phenarr->offday_arr[py][0] = ctrl->simstartyear+py;

				phenarr->onday_arr[py][1]  = onday_arr[py];
				phenarr->offday_arr[py][1] = offday_arr[py];

				/* wrinting out log file */
				if (onday_arr[py]  < onday_min)  onday_min  = onday_arr[py];
				if (onday_arr[py]  > onday_max)  onday_max  = onday_arr[py];
				if (offday_arr[py] < offday_min) offday_min = offday_arr[py];
				if (offday_arr[py] > offday_max) offday_max = offday_arr[py];

			}

		}
		/**************************************************************/
		/* B. GSI-METHOD for calculating onset and offset day */
		/**************************************************************/
		else
		{
			for (py=0 ; py<phenyears ; py++)
			{
				onday_arr[py] = phenarr->onday_arr[py][1];
				offday_arr[py] = phenarr->offday_arr[py][1];

				/* wrinting out log file */
				if (onday_arr[py]  < onday_min)  onday_min  = onday_arr[py];
				if (onday_arr[py]  > onday_max)  onday_max  = onday_arr[py];
				if (offday_arr[py] < offday_min) offday_min = offday_arr[py];
				if (offday_arr[py] > offday_max) offday_max = offday_arr[py];

				
			}
		}

	
		/* 3. wrinting log file */
		fprintf(logfile.ptr, "SGS value (min and max): %6i %6i\n", onday_min, onday_max);
		fprintf(logfile.ptr, "EGS value (min and max): %6i %6i\n", offday_min, offday_max);
		fprintf(logfile.ptr, " \n");

		
		/* 4. now the onset and offset days are established for each phenyear,either by the deciduous tree or the grass model.  
		Next loop through phenyears filling the phenology signal arrays and copying them to the permanent phen struct arrays */
		for (py=0 ; py<phenyears ; py++)
		{
			onday = onday_arr[py];
			offday = offday_arr[py];
		
		} /* end phenyears loop for filling permanent arrays */
	} /* end else phenology model block */
	
	
	/* writing in logfile */
	if (PLT->PLT_num) 
		fprintf(logfile.ptr, "START OF PHENPHASES (MM/DD)\n");
	else
		fprintf(logfile.ptr, " \n");
	
	/* free the local array memory */
	free(onday_arr); 
	free(offday_arr);
	free(enddays);
	free(mondays);

	return (errorCode);
}

