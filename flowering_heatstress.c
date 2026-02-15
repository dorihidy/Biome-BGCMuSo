/* 
flowering_heatstress.c
Calculate daily crop-specific phenology fluxes

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2025, D. Hidy [dori.hidy@gmail.com]
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

int flowering_heatstress(const epconst_struct* epc, const metvar_struct* metv, cstate_struct* cs, epvar_struct* epv, cflux_struct* cf, nflux_struct* nf)

{
	double mort_act;


	int errorCode=0;
	
	/* phenological stages */
	int n_actphen          = (int) epv->n_actphen;
	int n_flowHS_phenophase = (int) epc->n_flowHS_phenophase;
	
	/* in flowering phase heat stress can affect the seed fertilization */
	if (n_actphen == n_flowHS_phenophase)
	{
		if(metv->Tmax < epc->flowHS_parT1)
			mort_act = 0;
		else
		{
			if (metv->Tmax < epc->flowHS_parT2) 
				mort_act = (metv->Tmax - epc->flowHS_parT1) / (epc->flowHS_parT2 - epc->flowHS_parT1);
			else
				mort_act = 1;
		}

		epv->flowHS_mort = mort_act;
		


	}
	else
		epv->flowHS_mort = 0;
		/* in phenophase of yield filling the heat stress affects the material flow to yield (parameter: flowHS_mort)
		if (n_actphen == n_flowHS_phenophase+1 && epc->yield_cn) 
		{
			cf->yieldc_to_flowHS = epv->flowHS_mort * epc->flowHS_parMORT * cs->yieldc;
			nf->yieldn_to_flowHS = cf->yieldc_to_flowHS / epc->yield_cn;
		}
		else */
			

	cf->yieldc_to_flowHS = epv->flowHS_mort * epc->flowHS_parMORT * cs->yieldc;
	nf->yieldn_to_flowHS = cf->yieldc_to_flowHS / epc->yield_cn;
	cs->calc_flowHS += cf->yieldc_to_flowHS;
	cs->flowHSsnk_C += cf->yieldc_to_flowHS;



	return (errorCode);
}

