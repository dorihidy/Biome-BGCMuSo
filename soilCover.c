 /* 
soilCover.c
CALCULATING the effect of soilCover on the top soil layers: evaporation limitation effect + UNIFORM DISTRIBUTION OF VWC and tsoil 

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
#include "pointbgc_struct.h"
#include "bgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"

int soilCover(siteconst_struct* sitec, soilprop_struct* sprop, metvar_struct* metv,  epvar_struct* epv, cstate_struct* cs)
{

	/* soilCovering parameters */
	int layer_SCeff, layer, errorCode;
	double soilw_SUM, tsoilSUM, soilCover;	 


	errorCode=0;
	tsoilSUM=soilw_SUM=0;	 
	

	layer_SCeff = ( int) sprop->pLAY_soilCover;

	/**********************************************************************************************/
	/* I. CALCULATING the effect of soil cover on the top soil layers: UNIFORM DISTRIBUTION OF tsoil */

	soilCover = cs->litrCabove_total + cs->cwdCabove_total;
	epv->SCpercent = sprop->p1_soilCover * pow((soilCover)/sprop->p2_soilCover, sprop->p3_soilCover);
	epv->SC_EVPred   = pow(sprop->pRED_soilCover, epv->SCpercent/100);

	if (soilCover > sprop->pCRIT_soilCover)
	{
 		for (layer = 0; layer<layer_SCeff; layer++)
		{
			tsoilSUM += metv->tsoil[layer] * sitec->soillayer_thickness[layer] / (sitec->soillayer_depth[layer_SCeff-1]);
		}


		for (layer = 0; layer<layer_SCeff; layer++)
		{
			metv->tsoil[layer] = tsoilSUM;
		}

		metv->tsoil_surface     = metv->tsoil[0];
	}

   return (errorCode);
}
	