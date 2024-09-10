/* 
canopy_et.c
A single-function treatment of canopy evaporation and transpiration
fluxes.  

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group (NTSG)
School of Forestry, University of Montana
Missoula, MT 59812
*-*-*-*-*-*-*-*-*-*-*-*-*-*
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

int canopy_et(const control_struct* ctrl, const metvar_struct* metv, epvar_struct* epv, wflux_struct* wf)
{
	int errorCode=0;
	double e, cwe, t, TRP, TRPsun, TRPshade, TRP_pot, TRPsun_pot, TRPshade_pot,e_dayl,t_dayl,diff;
	
	pmet_struct pmet_in;
    
	e=cwe=t=TRP=TRPsun=TRPshade=e_dayl=t_dayl=TRP_pot=TRPsun_pot=TRPshade_pot=0;


	/* Assign values in pmet_in that don't change */
	pmet_in.ta = metv->Tday;
	pmet_in.pa = metv->pa;
	pmet_in.vpd = metv->vpd;


	/* Canopy evaporation, if any water was intercepted */
	/* Calculate Penman-Monteith evaporation, given the canopy conductances to
	evaporated water and sensible heat.  Calculate the time required to 
	evaporate all the canopy water at the daily average conditions, and 
	subtract that time from the daylength to get the effective daylength for
	transpiration. */
	if (wf->prcp_to_canopyw)
	{
		/* assign appropriate resistance and radiation for pmet_in */
		pmet_in.rv = 1.0/epv->gc_e_wv;
		pmet_in.rh = 1.0/epv->gc_sh;

		/* choose radiation calculation method */
		if (ctrl->radiation_flag == 0)
			pmet_in.irad = metv->swabs;
		else
			pmet_in.irad = metv->RADnet;
		
		/* call penman-monteith function, returns e in kg/m2/s */
		if (penmon(&pmet_in, 0, &e))
		{
			printf("ERROR in penmon.c for canopy EVP in canopy_et.c\n");
			errorCode=1;
		}
		
		/* calculate the time required to evaporate all the canopy water */
		e_dayl = wf->prcp_to_canopyw/e;
		
		if (e_dayl > metv->dayl)  
		{
			/* day not long enough to EVP. all int. water */
			TRP = 0.0;    /* no time left for transpiration */
			cwe = e * metv->dayl;   /* daylength limits canopy evaporation */
		}
		else                
		{
			/* all intercepted water evaporated */
			cwe = wf->prcp_to_canopyw;
			
			/* adjust daylength for transpiration */
			t_dayl = metv->dayl - e_dayl;
			 
			/**********ACTUAL******************/
			/* calculate  transpiration using adjusted daylength */
			/* first for sunlit canopy fraction */
			pmet_in.rv = 1.0/epv->gl_t_wv_sun;
			pmet_in.rh = 1.0/epv->gl_sh;

			/* choose radiation calculation method */
			if (ctrl->radiation_flag == 0)
				pmet_in.irad = metv->swabs_per_plaisun;
			else
				pmet_in.irad = metv->RADnet_per_plaisun;

			/* call Penman-Monthieth function */
			if (penmon(&pmet_in, 0, &t))
			{
				printf("ERROR in penmon.c for adjusted transpiration in canopy_et.c\n");
				errorCode=1;
			}
			TRPsun = t * t_dayl * epv->plaisun;
			
			/* next for shaded canopy fraction */
			pmet_in.rv = 1.0/epv->gl_t_wv_shade;
			pmet_in.rh = 1.0/epv->gl_sh;
			
			/* choose radiation calculation method */
			if (ctrl->radiation_flag == 0)
				pmet_in.irad = metv->swabs_per_plaishade;
			else
				pmet_in.irad = metv->RADnet_per_plaishade;
		
			/* call Penman-Monthieth function */
			if (penmon(&pmet_in, 0, &t))
			{
				printf("ERROR in penmon.c for adjusted transpiration in canopy_et.c\n");
				errorCode=1;
			}
			TRPshade = t * t_dayl * epv->plaishade;
			TRP = TRPsun + TRPshade;

			/**********POTENTIAL******************/
			/* calculate  transpiration using adjusted daylength */
			/* first for sunlit canopy fraction */
			pmet_in.rv = 1.0/epv->gl_t_wv_sunPOT;
			pmet_in.rh = 1.0/epv->gl_sh;

			/* choose radiation calculation method */
			if (ctrl->radiation_flag == 0)
				pmet_in.irad = metv->swabs_per_plaisun;
			else
				pmet_in.irad = metv->RADnet_per_plaisun;

			/* call Penman-Monthieth function */
			if (penmon(&pmet_in, 0, &t))
			{
				printf("ERROR in penmon.c for adjusted transpiration in canopy_et.c\n");
				errorCode=1;
			}
			TRPsun_pot = t * t_dayl * epv->plaisun;
			
			/* next for shaded canopy fraction */
			pmet_in.rv = 1.0/epv->gl_t_wv_shadePOT;
			pmet_in.rh = 1.0/epv->gl_sh;
			
			/* choose radiation calculation method */
			if (ctrl->radiation_flag == 0)
				pmet_in.irad = metv->swabs_per_plaishade;
			else
				pmet_in.irad = metv->RADnet_per_plaishade;
		
			/* call Penman-Monthieth function */
			if (penmon(&pmet_in, 0, &t))
			{
				printf("ERROR in penmon.c for adjusted transpiration in canopy_et.c\n");
				errorCode=1;
			}
			TRPshade_pot = t * t_dayl * epv->plaishade;
			TRP_pot = TRPsun_pot + TRPshade_pot;
		}
		wf->EVPcanopyw = cwe;
	}    /* end if canopy_water */
	
	else /* no canopy water, transpiration with unadjusted daylength */
	{
		wf->EVPcanopyw = 0;

		/**********ACTUAL******************/
		/* first for sunlit canopy fraction */
		pmet_in.rv = 1.0/epv->gl_t_wv_sun;
		pmet_in.rh = 1.0/epv->gl_sh;
	
		/* choose radiation calculation method */
		if (ctrl->radiation_flag == 0)
			pmet_in.irad = metv->swabs_per_plaisun;
		else
			pmet_in.irad = metv->RADnet_per_plaisun;
			

		/* call Penman-Monthieth function */
		if (penmon(&pmet_in, 0, &t))
		{
			printf("ERROR in penmon.c for adjusted transpiration in canopy_et.c\n");
			errorCode=1;
		}
		
		TRPsun = t * metv->dayl * epv->plaisun;
		
		/* next for shaded canopy fraction */
		pmet_in.rv = 1.0/epv->gl_t_wv_shade;
		pmet_in.rh = 1.0/epv->gl_sh;
		
		/* choose radiation calculation method */
		if (ctrl->radiation_flag == 0)
			pmet_in.irad = metv->swabs_per_plaishade;
		else
			pmet_in.irad = metv->RADnet_per_plaishade;
		

		/* call Penman-Monthieth function */
		if (penmon(&pmet_in, 0, &t))
		{
			printf("ERROR in penmon.c for adjusted transpiration in canopy_et.c\n");
			errorCode=1;
		}
		TRPshade = t * metv->dayl * epv->plaishade;
		TRP = TRPsun + TRPshade;

		/**********POTENTIAL******************/
		/* first for sunlit canopy fraction */
		pmet_in.rv = 1.0/epv->gl_t_wv_sunPOT;
		pmet_in.rh = 1.0/epv->gl_sh;
	
		/* choose radiation calculation method */
		if (ctrl->radiation_flag == 0)
			pmet_in.irad = metv->swabs_per_plaisun;
		else
			pmet_in.irad = metv->RADnet_per_plaisun;
			

		/* call Penman-Monthieth function */
		if (penmon(&pmet_in, 0, &t))
		{
			printf("ERROR in penmon.c for adjusted transpiration in canopy_et.c\n");
			errorCode=1;
		}
		
		TRPsun_pot = t * metv->dayl * epv->plaisun;
		
		/* next for shaded canopy fraction */
		pmet_in.rv = 1.0/epv->gl_t_wv_shadePOT;
		pmet_in.rh = 1.0/epv->gl_sh;
		
		/* choose radiation calculation method */
		if (ctrl->radiation_flag == 0)
			pmet_in.irad = metv->swabs_per_plaishade;
		else
			pmet_in.irad = metv->RADnet_per_plaishade;
		

		/* call Penman-Monthieth function */
		if (penmon(&pmet_in, 0, &t))
		{
			printf("ERROR in penmon.c for adjusted transpiration in canopy_et.c\n");
			errorCode=1;
		}
		TRPshade_pot = t * metv->dayl * epv->plaishade;
		TRP_pot = TRPsun_pot + TRPshade_pot;
	}
	wf->potTRPsoilw        = TRP_pot;

	/*---------------------------------------------------------------*/
	/* 5. Energy limit  */
	/*---------------------------------------------------------------*/
	diff = wf->EVPcanopyw + wf->potTRPsoilw - wf->potETcanopy;
	if (diff > CRIT_PRECwater)
	{
		if (wf->EVPcanopyw  > wf->potETcanopy)
		{
			wf->EVPcanopyw   = wf->potETcanopy;
			wf->potTRPsoilw = 0;
		}
		else
		{
			wf->potTRPsoilw -= diff;
		}
	}


	return (errorCode);
}