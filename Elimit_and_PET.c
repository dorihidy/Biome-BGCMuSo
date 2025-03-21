/* 
Elimit_and_PET.c
daily bare soil evaporation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2022, D. Hidy [dori.hidy@gmail.com]
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

int Elimit_and_PET(const control_struct* ctrl, const epconst_struct* epc, const soilprop_struct* sprop, const metvar_struct* metv, epvar_struct* epv, wflux_struct* wf)
{
	int errorCode=0;
	double rbl;					                /* (m/s) boundary layer resistance */
	double potEVPandSUBLsurface, potETcanopy;	/* (kg/m2/s) potential evaporation (daytime) */
	double rcorr;				                /* correction factor for temp and pressure */
	double lhvap; 
	pmet_struct pmet_in;		                /* input structure for penmon function */

	double EEQ;			        /* internal variable of DSSAT model */
	double wisp = 3;			/* wind speed (no input data - constans value = 3m/s) */
	double potETS = 0; 

	/*---------------------------------------------------------------*/
	/* 0. ENERGETIC CONTROL - maximum energy */
	/*---------------------------------------------------------------*/
	
	lhvap = 2.5023e6 - 2430.54 * metv->Tday;
	wf->ET_Elimit = ((metv->swabs+metv->swtrans) * metv->dayl)/lhvap;
	if (wf->ET_Elimit < 0) wf->ET_Elimit = 0; 

	/*---------------------------------------------------------------*/
	/* I. Penman-Montieth */
	/*---------------------------------------------------------------*/
	if (ctrl->ET_flag == 0)
	{
		/* correct conductances for temperature and pressure based on Jones (1992)
		with standard conditions assumed to be 20 deg C, 101300 Pa */
		rcorr = 1.0/(pow((metv->Tday+273.15)/293.15, 1.75) * 101300/metv->pa);

		/* new bare-soil evaporation routine */
		/* first calculate potential evaporation, assuming the resistance for vapor transport is equal to the resistance for sensible heat transport.  
		That is, no additional resistance for vapor transport to the soil surface. This represents evaporation from a wet surface with
		a specified aerodynamic resistance (= boundary layer resistance). The aerodynamic resistance is for now set as a constant, and is
		taken from observations over bare soil in tiger-bush in south-west Niger: rbl = 107 s m-1 (Wallace and Holwill, 1997). */
		rbl = rcorr * sprop->aerodyn_resist;
	
		/* normal run: aerodyn_resist is set to a measured value (test run: -9999 -> no evaporation, transpiration is calculated using a default value: 107m/s) */
		if (sprop->aerodyn_resist > 0)
		{
			/*---------------------------------------------------------------*/
			/* 1. Limit of Surface evaporation with transmitted radiation */
			/*---------------------------------------------------------------*/

			/* fill the pmet_in structure */
			pmet_in.ta = metv->Tday;
			pmet_in.pa = metv->pa;
			pmet_in.vpd = metv->vpd;
			pmet_in.irad = metv->swtrans;
			pmet_in.rv = rbl;
			pmet_in.rh = rbl;

			/* calculate potEVP in kg/m2/s */
			penmon(&pmet_in, 0, &potEVPandSUBLsurface);
		

			/*---------------------------------------------------------------*/
			/* 2. Limit of canopy evaporation with absorbed radiation */
			/*---------------------------------------------------------------*/

			/* fill the pmet_in structure */
			pmet_in.ta = metv->Tday;
			pmet_in.pa = metv->pa;
			pmet_in.vpd = metv->vpd;
			pmet_in.irad = metv->swabs;
			pmet_in.rv = rbl;
			pmet_in.rh = rbl;

			/* calculate potEVP in kg/m2/s */
			penmon(&pmet_in, 0, &potETcanopy);

			if (epv->projLAI < LAIlimit_potET)
			{	
				potETcanopy *= epv->projLAI / LAIlimit_potET;
			}
		}
		else
		{
			potEVPandSUBLsurface = 0;
			potETcanopy          = 0;
		
		}

	
		/* convert to daily total kg/m2 */
		wf->potEVPandSUBLsurface = potEVPandSUBLsurface * metv->dayl;
		wf->potETcanopy          = potETcanopy * metv->dayl;

		
		/*---------------------------------------------------------------*/
		/* 3. maximum energy limits potential */
		/*---------------------------------------------------------------*/
		
		potETS        = wf->potEVPandSUBLsurface+wf->potETcanopy;


		if (potETS > wf->ET_Elimit)
		{
			wf->potEVPandSUBLsurface *= (wf->ET_Elimit/potETS);
			wf->potETcanopy          *= (wf->ET_Elimit/potETS);
		}

		wf->PET = wf->potEVPandSUBLsurface+wf->potETcanopy;


	}
	/*---------------------------------------------------------------*/
	/* II. Priestley-Taylor */
	/*---------------------------------------------------------------*/
	else
	{
		/*-----------------------*/
		/* 1. DSSAT calculation  */
		/*-----------------------*/
		EEQ = metv->swavgfd * (nSEC_IN_DAY / 1000000.)  * (4.88E-3 - 4.37E-3*epv->albedo_LAI)*(0.6*metv->Tmax+0.4*metv->Tmin+29); //unit change: W/m2 to MJ/m2/day

		wf->PET = 1.1*EEQ*(1+wisp/100);

		if (metv->Tmax > epc->PT_ETcritT)  wf->PET=EEQ*((metv->Tmax-epc->PT_ETcritT)*0.05+1.1);
		if (metv->Tmax < 5)                wf->PET=EEQ*0.01*exp(0.18*(metv->Tmax+20.0));

		/*-----------------------*/
		/* 2. energy limitation */
		/*-----------------------*/
		if (wf->PET > wf->ET_Elimit)
		{
			wf->PET *= (wf->ET_Elimit/wf->PET);
		}

		/*-----------------------*/
		/* 3. separation */
		/*-----------------------*/

		if (metv->swabs + metv->swtrans > 0)
		{
			wf->potEVPandSUBLsurface = wf->PET * (metv->swtrans / (metv->swabs + metv->swtrans));
			wf->potETcanopy          = wf->PET * (metv->swabs   / (metv->swabs + metv->swtrans));
		}
		else
		{
			wf->potEVPandSUBLsurface = 0;
			wf->potETcanopy          = 0;
		}
	}

		
	/*---------------------------------------------------------------*/
	/* Surface evaporation - Snow sublimation energy decreases energy to evaporation (canopyET and soilEVP is zero belo 0 Celsius) */
	/*---------------------------------------------------------------*/
		
	if (metv->Tday > 0)
	{
		if (wf->SUBLsnoww - wf->potEVPandSUBLsurface > CRIT_PREC_lenient)
		{
			wf->SUBLsnoww     = wf->potEVPandSUBLsurface;
			wf->potEVPsurface = 0;
		}
		else
			wf->potEVPsurface = wf->potEVPandSUBLsurface - wf->SUBLsnoww;
	}
	else
	{
		if (wf->SUBLsnoww - wf->potEVPandSUBLsurface > CRIT_PREC_lenient) wf->SUBLsnoww = wf->potEVPandSUBLsurface;
		wf->potEVPsurface = 0;
		wf->potETcanopy   = 0;
		wf->PET           = wf->potEVPandSUBLsurface + wf->potETcanopy;
	}

	
	

	
	return(errorCode);
}

