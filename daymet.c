/*daymet.c
transfer one day of meteorological data from metarr struct to metv struct

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2025, D. Hidy [dori.hidy@gmail.com]
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

int daymet(const control_struct* ctrl,const metarr_struct* metarr, const epconst_struct* epc, metvar_struct* metv, double snoww)
{
	/* generates daily meteorological variables from the metarray struct */
	double Tmax,Tmin,Tavg,TavgRA11,TavgRA30,TavgRA10,Tday,tdiff, Tsoil_top;
	int errorCode=0;
	static double Tacclim_pre;


	/* convert prcp from cm --> kg/m2 */
	metv->prcp = metarr->prcp_array[ctrl->metday] * 10.0;
	

	/* air temperature calculations (all temperatures deg C) */
	metv->Tmax			= Tmax    = metarr->Tmax_array[ctrl->metday];
	metv->Tmin			= Tmin    = metarr->Tmin_array[ctrl->metday];
	metv->Tavg		    = Tavg    = metarr->Tavg_array[ctrl->metday];

	metv->Tday			= Tday	= metarr->Tday_array[ctrl->metday];

	metv->Tnight		= (Tday + Tmin) / 2.0;
	metv->TavgRA11	    = TavgRA11 = metarr->TavgRA11_array[ctrl->metday];
	metv->TavgRA30	    = TavgRA30 = metarr->TavgRA30_array[ctrl->metday];
	metv->TavgRA10	    = TavgRA10 = metarr->TavgRA10_array[ctrl->metday];
	
	metv->annTavg       = metarr->annTavg_array[ctrl->simyr];
	metv->annTavgRA     = metarr->annTavgRA_array[ctrl->simyr];
	metv->annTrange     = metarr->annTrange_array[ctrl->simyr];
	metv->annTrangeRA   = metarr->annTrangeRA_array[ctrl->simyr];
	
	metv->tempradF     = metarr->tempradF_array[ctrl->metday];
	metv->tempradFra   = metarr->tempradFra_array[ctrl->metday];

	if (!ctrl->metday)
	{
		metv->Tacclim            = metv->Tday;
		Tacclim_pre = metv->Tday;
	}
	else
	{
		if (epc->tau)
			metv->Tacclim = Tacclim_pre + ((metv->Tday - Tacclim_pre) / epc->tau);
		else
			metv->Tacclim = metv->Tday; 
		
		Tacclim_pre = metv->Tacclim;
	
	}

	
	/* **********************************************************************************/
	/* new estimation of Tsoil - on the first day original method is used */
	
	/* ORIGINAL: for this version, an 11-day running weighted average of daily average temperature is used as the soil temperature at 10 cm.
	For days 1-10, a 1-10 day running weighted average is used instead.The tail of the running average is weighted linearly from 1 to 11.
	There are no corrections for snowpack or vegetation cover.	*/

	if (ctrl->metday < 1)
	{
		Tsoil_top = metv->TavgRA11;
		/* soil temperature correction using difference from annual average tair */
		tdiff =  metv->annTavg - Tsoil_top;
		
		if (snoww)
		{
			Tsoil_top += 0.2 * tdiff;
		}
		else
		{
			Tsoil_top += 0.1 * tdiff;
		}

		metv->Tsoil_surface     = Tsoil_top;
	}

	/* 3 m below the ground surface (last layer) is specified by the annual mean surface air temperature */
	metv->Tsoil[N_SOILLAYERS-1] = metv->annTavgRA;
	

	
	/* **********************************************************************************/
	
	/* daylight average vapor pressure deficit (Pa) */
	metv->vpd = metarr->vpd_array[ctrl->metday];

	/* daylight average	shortwave flux density (W/m2) */
	metv->swavgfd =  metarr->swavgfd_array[ctrl->metday];
	
	/* PAR (W/m2) */
	metv->par = metarr->par_array[ctrl->metday];

	/* daylength (s) */
	metv->dayl = metarr->dayl_array[ctrl->metday];

	

	return (errorCode);
}
