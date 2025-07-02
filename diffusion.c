 /*
diffusion.c
Calculation of diffusion flux between soil layerts between the difference of their equilibrium temperature

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2025, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_func.h"
#include "bgc_constants.h"
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

int diffusion(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{

	int errorCode = 0;
	int layer, N_NAGlayers;

	double soilwDiffus_act, DBAR;
	double VWC0, VWC0_sat, VWC0_eq, VWC0_wp, VWC0_EqFC, VWC1, VWC1_sat, VWC1_eq, VWC1_wp, VWC1_EqFC;
	double rVWC0, rVWC1, rVWC_limit, VWC0_limit, VWC1_limit;

	double dz0, dz1;
	double dLk;


	soilwDiffus_act = 0;

	/* diffusion.c is used only for layers without GW */
	if (sprop->GWlayer == DATA_GAP)
		N_NAGlayers = N_SOILLAYERS-1;
	else
		N_NAGlayers = (int)sprop->CFlayer-1;

	/* --------------------------------------------------------------------------------------------------------*/
	/* diffusion fluxes between noGWlayers: last noGWlayer and first GWlayers is handled in groundwater_diffusion.c  */

	for (layer = N_NAGlayers - 1; layer >= 0; layer--)
	{
		dz0 = sitec->soillayer_thickness[layer];
		VWC0 = epv->VWC[layer];
		VWC0_sat = sprop->VWCsat[layer];
		VWC0_eq = sprop->VWCeq[layer];
		VWC0_wp = sprop->VWCwp[layer];
		VWC0_EqFC = MAX(VWC0_eq, sprop->VWCfc[layer]);

		dz1 = sitec->soillayer_thickness[layer + 1];
		VWC1 = epv->VWC[layer + 1];
		VWC1_sat = sprop->VWCsat[layer + 1];
		VWC1_eq = sprop->VWCeq[layer + 1];
		VWC1_wp = sprop->VWCwp[layer + 1];
		VWC1_EqFC = MAX(VWC1_eq, sprop->VWCfc[layer+1]);
			
		rVWC0 = (VWC0 - VWC0_wp) / (VWC0_EqFC - VWC0_wp);
		rVWC1 = (VWC1 - VWC1_wp) / (VWC1_EqFC - VWC1_wp);
		rVWC_limit = (rVWC0 + rVWC1) / 2;
		VWC0_limit = rVWC_limit * (VWC0_EqFC - VWC0_wp) + VWC0_wp;
		VWC1_limit = rVWC_limit * (VWC1_EqFC - VWC1_wp) + VWC1_wp;

		if (VWC0_limit > VWC0_sat) VWC0_limit = VWC0_sat;
		if (VWC1_limit > VWC1_sat) VWC1_limit = VWC1_sat;

		dLk = DATA_GAP; // only interpreted in GW and CF layers

	

		if (!errorCode && calc_diffus(layer, sprop, dz0, VWC0, VWC0_sat, VWC0_EqFC, VWC0_wp, VWC0_limit, dz1, VWC1, VWC1_sat, VWC1_EqFC, VWC1_wp, VWC1_limit, dLk, &DBAR, &soilwDiffus_act))
		{
			printf("\n");
			printf("ERROR in calc_diffus.c for diffusion.c\n");
			errorCode = 1;
		}

		/* element layer */
		sprop->DBARarray[layer] = DBAR;

		if (fabs(soilwDiffus_act) < CRIT_PREC_lenient) soilwDiffus_act = 0;

		wf->soilwDiffus[layer] = soilwDiffus_act;

		/* state update */
		ws->soilw[layer] -= wf->soilwDiffus[layer];
		ws->soilw[layer + 1] += wf->soilwDiffus[layer];


		epv->VWC[layer] = ws->soilw[layer] / dz0 / water_density;
		epv->VWC[layer + 1] = ws->soilw[layer + 1] / dz1 / water_density;

	}


	

	/* --------------------------------------------------------------------------------------------------------*/
	/* BOTTOM LAYER IS SPECIAL 	*/

	if (sprop->GWD == DATA_GAP)
	{
		wf->soilwDiffus[N_SOILLAYERS - 1] = wf->soilwDiffus[N_SOILLAYERS - 2];
		ws->soilw[N_SOILLAYERS - 1] -= wf->soilwDiffus[N_SOILLAYERS - 1];
		epv->VWC[N_SOILLAYERS - 1] = ws->soilw[N_SOILLAYERS - 1] / sitec->soillayer_thickness[N_SOILLAYERS - 1] / water_density;
	}








	return (errorCode);

}
