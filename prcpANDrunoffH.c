/* 
prcpANDrunoffH.c
routing of daily precipitation to canopy, soil, snowpack and Hortonian runoff

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group (NTSG)
School of Forestry, University of Montana
Missoula, MT 59812
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

int prcpANDrunoffH(const control_struct* ctrl, const wstate_struct* ws, const metvar_struct* metv, const soilprop_struct* sprop,  const epconst_struct* epc,
	               epvar_struct* epv, wflux_struct* wf)
{
	/* Precipitation routing, as either rain or snow.  Rain can be intercepted on the canopy, and amount in excess of interception is
	routed to the soil.  The intercepted volume is saved in a temporary variable that is later passed to the canopy evaporation routine,
	but there is no day-to-day storage of water on the canopy: any that is not evaporated in one day is routed to the soil in canopy_et.c.
	
	There is no canopy interception of snow.

	There is no runoff in absence of snow
	*/
	
	int errorCode=0;
	double max_int;
	double prcp, through;
	double coeff_soiltype, coeff_soilmoist, coeff_runoff;

	prcp = metv->prcp + wf->IRG_to_prcp;

	/* 1. maximum daily canopy interception  */
	/* Original method: with allLAI (kg intercepted/kg rain/unit all-sided LAI/day)  */
	/* New method: with single LAI and (1/LAI/d) canopy water interception coefficient*/
	
	if (ctrl->interception_flag == 1)
		max_int = epc->int_coef * prcp                   * epv->allLAI;
	else
		max_int = epc->int_coef * (1-exp(-0.1*(prcp+1))) * epv->projLAI;

	
	/* 2. rain vs. snow, and canopy interception */
	if (metv->Tavg > 0.0 && ws->snoww == 0)             /* rain */
	{
		if (prcp <= max_int)          /* all intercepted */
		{
			wf->prcp_to_canopyw = prcp; 
			through = 0.0;
		}
		else                          /* canopy limits interception */
		{
			wf->prcp_to_canopyw = max_int;
			through = prcp - max_int;
		}
		
		/* 3. throughfall to soil water and Hortonian runoff */
		/* when the precipitation at the surface exceeds the max. infiltration rate, the excess water is put into surface runoff (Balsamo et al. 20008; Eq.(7)) */
		if (sprop->RCN > 0 && ws->snoww == 0 && ws->pondw == 0)
		{
			coeff_soiltype  = 254*(100 / sprop->RCN - 1);

			coeff_soilmoist = 0.15 * ((sprop->VWCsat[0] - epv->VWC[0]) / (sprop->VWCsat[0]  - sprop->VWChw[0]));

			coeff_runoff = coeff_soiltype * coeff_soilmoist;

			if (through > coeff_runoff)
				wf->prcp_to_runoff = pow(through - coeff_runoff, 2) / (through + (1 - coeff_soilmoist)*coeff_soiltype);
			else
				wf->prcp_to_runoff = 0;
		}
		else	
			wf->prcp_to_runoff = 0;

		if (through > wf->prcp_to_runoff)
			wf->prcp_to_soilSurface = through - wf->prcp_to_runoff;  
		else
			wf->prcp_to_soilSurface = 0;  
	}
	else                              /* snow */
	{
		wf->prcp_to_snoww = prcp;     /* no interception */
	}
	
	return(errorCode);
}
