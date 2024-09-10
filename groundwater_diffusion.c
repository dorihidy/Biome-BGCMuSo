/*
groundwater_diffusion.c
UPWARD WATER MOVEMENT in groundwater layers

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2024, D. Hidy [dori.hidy@gmail.com]
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

int groundwater_diffusion(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)
{

	int errorCode = 0;
	int GWlayer;

	double soilwDiffus_act;
	double VWC0, VWC0_sat, VWC0_fc, VWC0_wp, VWC1, VWC1_sat, VWC1_fc, VWC1_wp;
	double rVWC0, rVWC1, rVWC_limit, VWC0_limit, VWC1_limit;

	double dz0, dz1;


	soilwDiffus_act = 0;
	GWlayer = (int)sprop->GWlayer;

	

	/* --------------------------------------------------------------------------------------------------------*/
	/* 1. between zoneSAT vs. zoneCAP -  */

	dz0 = sprop->dz_zoneCAPIL;
	VWC0 = sprop->VWC_zoneCAPIL;;
	VWC0_sat = sprop->VWCsat[GWlayer];
	VWC0_fc = sprop->VWCsat[GWlayer];
	VWC0_wp = sprop->VWCwp[GWlayer];

	
	dz1 = sprop->dz_zoneSAT;
	VWC1 = sprop->VWCsat[GWlayer];
	VWC1_sat = sprop->VWCsat[GWlayer];
	VWC1_fc = sprop->VWCsat[GWlayer];
	VWC1_wp = sprop->VWCwp[GWlayer];



	rVWC0 = (VWC0 - VWC0_wp) / (VWC0_fc - VWC0_wp);
	rVWC1 = (VWC1 - VWC1_wp) / (VWC1_fc - VWC1_wp);
	rVWC_limit = (rVWC0 + rVWC1) / 2;
	VWC0_limit = rVWC_limit * (VWC0_fc - VWC0_wp) + VWC0_wp;
	VWC1_limit = rVWC_limit * (VWC1_fc - VWC1_wp) + VWC1_wp;


	if (!errorCode && diffusCalc(sprop, dz0, VWC0, VWC0_sat, VWC0_fc, VWC0_wp, VWC0_limit, dz1, VWC1, VWC1_sat, VWC1_fc, VWC1_wp, VWC1_limit, &soilwDiffus_act))
	{
		printf("\n");
		printf("ERROR in diffusCalc.c for tipping.c\n");
		errorCode = 1;
	}

	if (fabs(soilwDiffus_act) < CRIT_PRECwater) soilwDiffus_act = 0;


	if (soilwDiffus_act > 0)
	{
		printf("\n");
		printf("ERROR in groundwater_diffusion.c: no diffusion from zoneCAPIL to zoneSAT.c\n");
		errorCode = 1;
	}

	wf->GWdischarge = -1 * soilwDiffus_act;


	sprop->soilw_zoneCAPIL -= soilwDiffus_act;
	sprop->VWC_zoneCAPIL = sprop->soilw_zoneCAPIL / dz0 / water_density;

	/* --------------------------------------------------------------------------------------------------------*/
	/* 2. if zoneNORM is present:  -  zoneCAPIL vs. zoneNORM, zoneNORM vs unsat; if zoneNORM is not present:  -  zoneCAPIL vs unsat  */

	if (sprop->dz_zoneNORM)
	{
		/* ------------------------------------ */
		/* 2.1.1 between zoneCAPIL vs. zoneNORM  */
		/* ------------------------------------ */
		dz0 = sprop->dz_zoneNORM;
		VWC0 = sprop->VWC_zoneNORM;;
		VWC0_sat = sprop->VWCsat[GWlayer];
		VWC0_fc = sprop->VWCfc_base[GWlayer];
		VWC0_wp = sprop->VWCwp[GWlayer];

		dz1 = sprop->dz_zoneCAPIL;
		VWC1 = sprop->VWC_zoneCAPIL;;
		VWC1_sat = sprop->VWCsat[GWlayer];
		VWC1_fc = sprop->VWCfc_base[GWlayer];
		VWC1_wp = sprop->VWCwp[GWlayer];

		rVWC0 = (VWC0 - VWC0_wp) / (VWC0_fc - VWC0_wp);
		rVWC1 = (VWC1 - VWC1_wp) / (VWC1_fc - VWC1_wp);
		rVWC_limit = (rVWC0 + rVWC1) / 2;
		VWC0_limit = rVWC_limit * (VWC0_fc - VWC0_wp) + VWC0_wp;
		VWC1_limit = rVWC_limit * (VWC1_fc - VWC1_wp) + VWC1_wp;


		if (!errorCode && diffusCalc(sprop, dz0, VWC0, VWC0_sat, VWC0_fc, VWC0_wp, VWC0_limit, dz1, VWC1, VWC1_sat, VWC1_fc, VWC1_wp, VWC1_limit, &soilwDiffus_act))
		{
			printf("\n");
			printf("ERROR in diffusCalc.c for tipping.c\n");
			errorCode = 1;
		}

		if (fabs(soilwDiffus_act) < CRIT_PRECwater) soilwDiffus_act = 0;


		wf->soilwDiffus_NORMvsCAPIL = soilwDiffus_act;


		sprop->soilw_zoneNORM  -= soilwDiffus_act;
		sprop->soilw_zoneCAPIL += soilwDiffus_act;

		sprop->VWC_zoneNORM  = sprop->soilw_zoneNORM / dz0 / water_density;
		sprop->VWC_zoneCAPIL = sprop->soilw_zoneCAPIL / dz1 / water_density;

		/* ------------------------------------ */
		/* 2.1.2 between zoneNORM vs. noGWlayer -  */
		/* ------------------------------------ */
		dz0 = sitec->soillayer_thickness[GWlayer-1];
		VWC0 = epv->VWC[GWlayer - 1];
		VWC0_sat = sprop->VWCsat[GWlayer - 1];
		VWC0_fc = sprop->VWCfc_base[GWlayer - 1];
		VWC0_wp = sprop->VWCwp[GWlayer - 1];

		dz1 = sprop->dz_zoneNORM;
		VWC1 = sprop->VWC_zoneNORM;;
		VWC1_sat = sprop->VWCsat[GWlayer];
		VWC1_fc = sprop->VWCfc_base[GWlayer];
		VWC1_wp = sprop->VWCwp[GWlayer];

		rVWC0 = (VWC0 - VWC0_wp) / (VWC0_fc - VWC0_wp);
		rVWC1 = (VWC1 - VWC1_wp) / (VWC1_fc - VWC1_wp);
		rVWC_limit = (rVWC0 + rVWC1) / 2;
		VWC0_limit = rVWC_limit * (VWC0_fc - VWC0_wp) + VWC0_wp;
		VWC1_limit = rVWC_limit * (VWC1_fc - VWC1_wp) + VWC1_wp;


		if (!errorCode && diffusCalc(sprop, dz0, VWC0, VWC0_sat, VWC0_fc, VWC0_wp, VWC0_limit, dz1, VWC1, VWC1_sat, VWC1_fc, VWC1_wp, VWC1_limit, &soilwDiffus_act))
		{
			printf("\n");
			printf("ERROR in diffusCalc.c for tipping.c\n");
			errorCode = 1;
		}

		if (fabs(soilwDiffus_act) < CRIT_PRECwater) soilwDiffus_act = 0;

		wf->soilwDiffus[GWlayer-1] = soilwDiffus_act;
		ws->soilw[GWlayer-1]      -= soilwDiffus_act;
		sprop->soilw_zoneNORM     += soilwDiffus_act;

		epv->VWC[GWlayer - 1] = ws->soilw[GWlayer - 1] / dz0 / water_density;
		sprop->VWC_zoneNORM = sprop->soilw_zoneNORM / dz1 / water_density;

	}
	/* ------------------------------------ */
	/* 2.2. between zoneCF (no zoneNORM) vs. noGWlayer -  */
	/* ------------------------------------ */
	else
	{
		dz0 = sitec->soillayer_thickness[GWlayer - 1];
		VWC0 = epv->VWC[GWlayer - 1];
		VWC0_sat = sprop->VWCsat[GWlayer - 1];
		VWC0_fc = sprop->VWCfc_base[GWlayer - 1];
		VWC0_wp = sprop->VWCwp[GWlayer - 1];

		dz1 = sprop->dz_zoneCAPIL;
		VWC1 = sprop->VWC_zoneCAPIL;;
		VWC1_sat = sprop->VWCsat[GWlayer];
		VWC1_fc = sprop->VWCfc_base[GWlayer];
		VWC1_wp = sprop->VWCwp[GWlayer];

		rVWC0 = (VWC0 - VWC0_wp) / (VWC0_fc - VWC0_wp);
		rVWC1 = (VWC1 - VWC1_wp) / (VWC1_fc - VWC1_wp);
		rVWC_limit = (rVWC0 + rVWC1) / 2;
		VWC0_limit = rVWC_limit * (VWC0_fc - VWC0_wp) + VWC0_wp;
		VWC1_limit = rVWC_limit * (VWC1_fc - VWC1_wp) + VWC1_wp;


		if (!errorCode && diffusCalc(sprop, dz0, VWC0, VWC0_sat, VWC0_fc, VWC0_wp, VWC0_limit, dz1, VWC1, VWC1_sat, VWC1_fc, VWC1_wp, VWC1_limit, &soilwDiffus_act))
		{
			printf("\n");
			printf("ERROR in diffusCalc.c for tipping.c\n");
			errorCode = 1;
		}

		if (fabs(soilwDiffus_act) < CRIT_PRECwater) soilwDiffus_act = 0;


		wf->soilwDiffus[GWlayer-1] = soilwDiffus_act;


		ws->soilw[GWlayer - 1] -= soilwDiffus_act;
		sprop->soilw_zoneCAPIL += soilwDiffus_act;

		epv->VWC[GWlayer - 1] = ws->soilw[GWlayer - 1] / dz0 / water_density;
		sprop->VWC_zoneCAPIL = sprop->soilw_zoneCAPIL / dz1 / water_density;
	}

	/* calculation net water transport */
	wf->soilwFlux[GWlayer]    = wf->soilwPercol[GWlayer] + wf->soilwDiffus[GWlayer];
	wf->soilwFlux_NORMvsCAPIL = wf->soilwPercol_NORMvsCAPIL + wf->soilwDiffus_NORMvsCAPIL;


	ws->soilw[GWlayer] = sprop->soilw_zoneNORM + sprop->soilw_zoneCAPIL + sprop->soilw_zoneSAT;
	epv->VWC[GWlayer] = ws->soilw[GWlayer] / sitec->soillayer_thickness[GWlayer] / water_density;



	return (errorCode);

}
