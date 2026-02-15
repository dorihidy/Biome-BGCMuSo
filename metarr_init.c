/*
metarr_init.c
Initialize meteorological data arrays for pointbgc simulation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2025, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_constants.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"
#include "misc_func.h"
  

/* It is assumed here that the meteorological datafile contains the following 
list of variables, with the indicated units.  All other variables in the metv
arrays are derived from this basic set:

VARIABLE    UNITS
yday        (none) (yearday)
prcp        cm     (daily total precipitation, water equivalent)  
Tmax        deg C  (daily maximum temperature) 
Tmin        deg C  (daily minimum temperature)
VPD         Pa     (daylight average VPD)
swavgfd     W/m2   (daylight average shortwave flux density)
daylength   s      (daylight duration)

*/

int metarr_init(point_struct* point, metarr_struct* metarr, const climchange_struct* scc,const siteconst_struct* sitec, const control_struct* ctrl) 
{
	int errorCode=0;
	int sd,j, metread, shift,y,m,d,nm;
	int ndays, nsimdays, dd, n_METvar;
	int nyears, year, day,gapday;
	double Tmax, Tmin, Tday, Tavg, prcp, vpd, swavgfd, dayl;
	double sw_MJ;
	double annTavgSUM, annTrangeSUM,monthTavgSUM, minTmonth, maxTmonth;

	int nMONTHday_array[]={31,28,31,30,31,30,31,31,30,31,30,31};

	n_METvar = 9;
	Tmax=Tmin=Tday=Tavg=prcp=vpd=swavgfd=dayl=0;
	nyears = ctrl->simyears;
	ndays    = nDAYS_OF_YEAR * nyears;
	nsimdays = nDAYS_OF_YEAR * (nyears-1) + point->nday_lastsimyear;

	/* allocate space for the metv arrays */
	if (!errorCode)
	{
		metarr->Tmax_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->Tmax_array)
		{
			printf("ERROR allocating for Tmax array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->Tmin_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->Tmin_array)
		{
			printf("ERROR allocating for Tmin array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->prcp_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->prcp_array)
		{
			printf("ERROR allocating for prcp array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->vpd_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->vpd_array)
		{
			printf("ERROR allocating for vpd array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->Tday_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->Tday_array)
		{
			printf("ERROR allocating for Tday array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->Tavg_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->Tavg_array)
		{
			printf("ERROR allocating for Tavg11 array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->TavgRA11_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->TavgRA11_array)
		{
			printf("ERROR allocating for TavgRA11 array\n");
			errorCode=218;
		}
	}

	/* new arrays */
	
	if (!errorCode)
	{
		metarr->TavgRA30_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->TavgRA30_array)
		{
			printf("ERROR allocating for TavgRA30 array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->TavgRA10_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->TavgRA10_array)
		{
			printf("ERROR allocating for TavgRA10 array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->tempradF_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->tempradF_array)
		{
			printf("ERROR allocating for tempradF array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->tempradFra_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->tempradFra_array)
		{
			printf("ERROR allocating for tempradFra array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->swavgfd_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->swavgfd_array)
		{
			printf("ERROR allocating for swavgfd array\n");
			errorCode=218;
		}
	}
	
	if (!errorCode)
	{
		metarr->par_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->par_array)
		{
			printf("ERROR allocating for par array\n");
			errorCode=218;
		}
	}
	
	if (!errorCode)
	{
		metarr->dayl_array = (double*) malloc(ndays * sizeof(double));
		if (!metarr->dayl_array)
		{
			printf("ERROR allocating for dayl_array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->annTavg_array = (double*) malloc(nyears * sizeof(double));
		if (!metarr->annTavg_array)
		{
			printf("ERROR allocating for annTavg_array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->monTavg_array = (double*) malloc(nyears * 12 * sizeof(double));
		if (!metarr->monTavg_array)
		{
			printf("ERROR allocating for monTavg_array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->annTrange_array= (double*) malloc(nyears * sizeof(double));
		if (!metarr->annTrange_array)
		{
			printf("ERROR allocating for annTrange_array array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->annTavgRA_array = (double*) malloc(nyears * sizeof(double));
		if (!metarr->annTavgRA_array)
		{
			printf("ERROR allocating for annTavgRA_array\n");
			errorCode=218;
		}
	}

	if (!errorCode)
	{
		metarr->annTrangeRA_array= (double*) malloc(nyears * sizeof(double));
		if (!metarr->annTrangeRA_array)
		{
			printf("ERROR allocating for annTrangeRA_array array\n");
			errorCode=218;
		}
	}


	/* southen hemisphere: year from 1th of July - last simulation year: truncated year */
	shift = ctrl->south_shift;
	point->nday_lastsimyear -= shift;
	if (point->nday_lastsimyear < 0) 
	{
		printf("ERROR in meteorological file: in southen hemisphere the number of days in the last simulation year must greater than 182\n");
		errorCode=218;
	}
	
	sd=0;
	gapday=0;
	/* begin daily loop: read input file, generate array values */
	while (!errorCode && sd < ndays)
	{
		metread = fscanf(point->metf.ptr,"%i%i%lf%lf%lf%lf%lf%lf%lf%*[^\n]",&year,&day,&Tmax,&Tmin,&Tday,&prcp,&vpd,&swavgfd,&dayl);
		
		if (sd+shift < nsimdays && metread != n_METvar)
		{
			printf("ERROR reading MET variables (must be 9 columns: year,day,Tmax,Tmin,Tday,prcp,vpd,swavgfd,dayl)\n");
			errorCode=218;
		}

		if (sd == 0 && year > ctrl->simstartyear)
		{
			printf("ERROR reading MET variables: file must contain meteorological data for each simulation day (from the beginning)\n");
			errorCode=218;
		}
			
		/* shifting meteorological data in southen hemisphere */
		if (point->nday_lastsimyear < nDAYS_OF_YEAR && sd+shift >= nsimdays)
		{
			Tmax=Tmin=Tday=prcp=vpd=swavgfd=dayl=0;
			for (j=0; j<nyears-1;j++)
			{
				dd = (j * nDAYS_OF_YEAR) + (point->nday_lastsimyear + gapday - 1);
				Tmax     += metarr->Tmax_array[dd];
				Tmin     += metarr->Tmin_array[dd];
				Tday     += metarr->Tday_array[dd];
				prcp     += metarr->prcp_array[dd];
				vpd      += metarr->vpd_array[dd];
				swavgfd  += metarr->swavgfd_array[dd];
				dayl     += metarr->dayl_array[dd];

			}
			Tmax    = Tmax    / (nyears-1);
			Tmin    = Tmin    / (nyears-1);
			Tday    = Tday    / (nyears-1);
			prcp    = prcp    / (nyears-1);
			vpd     = vpd     / (nyears-1);
			swavgfd = swavgfd / (nyears-1);
			dayl    = dayl    / (nyears-1);
			gapday      += 1;
		}
	

		/* read Tmax, Tmin, prcp, vpd, and srad */
		/* the following scan statement discards the Tday field in the standard MTCLIM version 3.1 input file */
		if (year >= ctrl->simstartyear && year < ctrl->simstartyear + ctrl->simyears)
		{
		
			/* control to avoid negative meteorological data */
 			if (prcp < 0 || vpd < 0 || swavgfd <= 0 || dayl <= 0)
			{
				printf("ERROR in met file: negative prcp/vpd/swavgfd/dayl or zero swavgfd/dayl data, metv_init.c\n");
				errorCode=218;
			}
		}
		else
		{
			if (sd == 0 && year > ctrl->simstartyear)
			{
				printf("ERROR reading met array, metarr_init.c\n");
				printf("Note: file must contain meteorological data for each simulation day\n");
				errorCode=218;
			}
		}

		
		/* apply the climate change scenario and store */
		if (!errorCode && (year > ctrl->simstartyear || (year == ctrl->simstartyear && day > shift)))
		{
			metarr->Tmax_array[sd] = Tmax + scc->s_Tmax;
			metarr->Tmin_array[sd] = Tmin + scc->s_Tmin;
			Tavg = (Tmax + Tmin) / 2.0;
			
			metarr->Tavg_array[sd] = (metarr->Tmax_array[sd] + metarr->Tmin_array[sd]) / 2.0;
			metarr->Tday_array[sd] = Tday + (metarr->Tavg_array[sd] - Tavg); 

			metarr->prcp_array[sd] = prcp * scc->s_prcp;
			metarr->vpd_array[sd] = vpd * scc->s_vpd;
			metarr->swavgfd_array[sd] = swavgfd * scc->s_swavgfd;
			metarr->par_array[sd] = swavgfd * 0.45 * scc->s_swavgfd;
			metarr->dayl_array[sd] = dayl;

			/* factor for soil temp. calculation */
			sw_MJ = swavgfd * dayl / 1000000;
			if (sd == 0)
				metarr->tempradF_array[sd] = metarr->Tavg_array[sd] +  (Tmax - metarr->Tavg_array[sd]) * sqrt(sw_MJ*0.03);
			else
				metarr->tempradF_array[sd] = (1 - sitec->albedo_sw) * (metarr->Tavg_array[sd] +  (Tmax - metarr->Tavg_array[sd]) * sqrt(sw_MJ*0.03))
										+ sitec->albedo_sw * metarr->tempradF_array[sd-1];

						
			/* counter of simdays*/			
			sd += 1;
		}
		

	}

	/* calculating annual air temperature and annual air temperature range */

	for (y=0; y<nyears;y++)
	{
		annTavgSUM = 0;
		monthTavgSUM = 0;
		minTmonth = 0;
		maxTmonth = 0;
		m=0;
		nm=0;

		for (d=0; d<nDAYS_OF_YEAR-1;d++)
		{
			annTavgSUM += metarr->Tavg_array[y*nDAYS_OF_YEAR+d];

			if (nm < nMONTHday_array[m])
			{
				monthTavgSUM += metarr->Tavg_array[y*nDAYS_OF_YEAR+d];
				nm=nm+1;
			}
			else
			{
				metarr->monTavg_array[y*nMONTHS_OF_YEAR+m] = monthTavgSUM/nMONTHday_array[m];

				if (monthTavgSUM/nMONTHday_array[m] < minTmonth) minTmonth = monthTavgSUM/nMONTHday_array[m];
				if (monthTavgSUM/nMONTHday_array[m] > maxTmonth) maxTmonth = monthTavgSUM/nMONTHday_array[m];

				m=m+1;
				nm=0;
				monthTavgSUM = 0;
			}
				
		}
		metarr->annTavg_array[y]    = annTavgSUM/nDAYS_OF_YEAR;
		metarr->annTrange_array[y] = maxTmonth - minTmonth;
	}

	/* perform running averages of daily average temperature for 
	use in soil temperature routine. 

	This implementation uses a linearly ramped 11-day running average 
	of daily mean air temperature, with days 1-10 based on a 1-10 day
	running average, respectively. 
	*/
	
	
	if (!errorCode && run_avg(metarr->Tavg_array, metarr->TavgRA11_array, ndays, 11, 1))
	{
		printf("ERROR in run_avg.c for metv_init.c \n");
		errorCode=218;
	}
	
	/*  new averages */
	if (!errorCode && run_avg(metarr->Tavg_array, metarr->TavgRA30_array, ndays, 30, 0))
	{
		printf("ERROR in run_avg.c for metv_init.c \n");
		errorCode=218;
	}

	if (!errorCode && run_avg(metarr->Tavg_array, metarr->TavgRA10_array, ndays, 10, 0))
	{
		printf("ERROR in run_avg.c for metv_init.c \n");
		errorCode=218;
	}

	if (!errorCode && run_avg(metarr->tempradF_array, metarr->tempradFra_array, ndays, 5, 0))
	{
		printf("ERROR in run_avg.c for metv_init.c \n");
		errorCode=218;
	}


	/* calculation of 10-year running average of mean annual air temperature and mean annual air tempearture range */
	if (nyears > 10)
	{
		if (!errorCode && run_avg(metarr->annTavg_array, metarr->annTavgRA_array, nyears, 10, 0))
		{
			printf("ERROR in run_avg.c for metv_init.c \n");
			errorCode=218;
		}
		if (!errorCode && run_avg(metarr->annTrange_array, metarr->annTrangeRA_array, nyears, 10, 0))
		{
			printf("ERROR in run_avg.c for metv_init.c \n");
			errorCode=218;
		}
	}
	else
	{

		annTavgSUM=0;
		annTrangeSUM=0;
		for (y=0; y<nyears; y++)
		{
			annTavgSUM   += metarr->annTavg_array[y];
			annTrangeSUM += metarr->annTrange_array[y];
		}
		for (y=0; y<nyears; y++)
		{
			metarr->annTavgRA_array[y]   = annTavgSUM/nyears;
			metarr->annTrangeRA_array[y] = annTrangeSUM/nyears;
		}
	}

	return (errorCode);
}
