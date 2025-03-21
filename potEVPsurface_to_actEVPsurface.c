/* 
potEVPsurface_to_actEVPsurface.c
calculation of actual soil evaporation from potential soil evaporation regarding to dry day limit (Joe Rictchie, 1970)

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
#include "bgc_constants.h"
#include "bgc_func.h"    
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))   

int potEVPsurface_to_actEVPsurface(soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{

	/* internal variables */
	int errorCode=0;
	double infiltPOT;

	
	/* assign water fluxes, all excess not evaporated goes to soil water compartment */
    wf->canopyw_to_soilw = wf->prcp_to_canopyw - wf->EVPcanopyw;

	infiltPOT     = wf->prcp_to_soilSurface + wf->snoww_to_soilw + wf->canopyw_to_soilw + wf->IRG_to_prcp + wf->soilwFLuxFromBelow;

	
	/*-----------------------------------------------------------------------------*/
	/* 1. first evaporation phase - no limit  if CF layer is at soil surface */
	if (ws->EVPsurface1cum < sprop->soilEVPcrit || sprop->CFlayer == 0)
	{
		epv->DSR = 0;
	
		if (infiltPOT >= ws->EVPsurface1cum)
			ws->EVPsurface1cum = 0;
		else
			ws->EVPsurface1cum -= infiltPOT;
		
		if (!errorCode && EVPphase1TOphase2(sprop, epv, ws, wf))
		{
			printf("\n");
			printf("ERROR in EVPphase1TOphase2.c for potEVPsurface_to_actEVPsurface.c\n");
			errorCode=1; 
		} 

	}
	/*-----------------------------------------------------------------------------*/
	/* 2. second evaporation phase - limit of DSR: ws->EVPsurface1cum > soilEVPcrit  */
	else
	{
		if (infiltPOT >= ws->EVPsurface2cum) 
		{
			epv->DSR = 0;

			/* berszivárgó csapadék feltölt CUM2-t, a maradék meg a CUM1-bõl amennyit csak tud */
			infiltPOT             -= ws->EVPsurface2cum;
			ws->EVPsurface2cum   =0;

			if (infiltPOT > sprop->soilEVPcrit)
				ws->EVPsurface1cum  = 0;
			else
				ws->EVPsurface1cum  -= infiltPOT;

			if (!errorCode && EVPphase1TOphase2(sprop, epv, ws, wf))
			{
				printf("\n");
				printf("ERROR in EVPphase1TOphase2.c for potEVPsurface_to_actEVPsurface.c\n");
				errorCode=1; 
			} 
		}
		else
		{

			ws->EVPsurface2cum -= infiltPOT;
			epv->DSR = pow((ws->EVPsurface2cum/sprop->coeff_EVPcum),2); 
			epv->DSR += 1;
			wf->EVPsoilw = sprop->coeff_EVPcum * pow(epv->DSR, 0.5) - ws->EVPsurface2cum;
			
			if (wf->EVPsoilw > wf->potEVPsurface) wf->EVPsoilw = wf->potEVPsurface;
		
			ws->EVPsurface2cum += wf->EVPsoilw;
			epv->DSR = pow((ws->EVPsurface2cum/sprop->coeff_EVPcum),2); 
		}
	}


	/* control */
	if (wf->EVPsoilw < 0 ||  epv->DSR < 0) 
	{
		printf("ERROR in soil evaporation calculation (potEVPsurface_to_actEVPsurface)\n");
		errorCode=1;
	}

	
	wf->EVPsoilw *= epv->SC_EVPred;

	return (errorCode);
}

int EVPphase1TOphase2(const soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{
	/* internal variables */
	int errorCode=0;


	ws->EVPsurface1cum += wf->potEVPsurface;
	if (ws->EVPsurface1cum > sprop->soilEVPcrit)
	{
		wf->EVPsoilw = wf->potEVPsurface - sprop->coeff_EVPlim*(ws->EVPsurface1cum - sprop->soilEVPcrit);
		ws->EVPsurface2cum = (1-sprop->coeff_EVPlim)*(ws->EVPsurface1cum - sprop->soilEVPcrit);
		epv->DSR = pow(ws->EVPsurface2cum/sprop->coeff_EVPcum,2);
		ws->EVPsurface1cum = sprop->soilEVPcrit;
	}
	else
	{
		wf->EVPsoilw = wf->potEVPsurface; 
	}





return (errorCode);
}