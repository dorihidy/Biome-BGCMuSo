/* 
groundwater_movement.c
calculate the deoth of the virtual groundsater layers and VWC and nutrient change due to GWD-moving

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
#include <time.h>
#include <math.h>
#include "ini.h"
#include "bgc_struct.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"
#include "bgc_func.h"


int groundwater_movement(const siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, soilInfo_struct* soilInfo,
	wstate_struct* ws, wflux_struct* wf, cstate_struct* cs, nstate_struct* ns)
{
	/* ------------------------------------------------------------------------------------------------------ */
	/* initialization */

	int errorCode, vlayer, layer, CFlayer, GWlayer, vlt, vlb, ll;
	int Itop, Ibot, FLAGtop, FLAGbot, dm;
	errorCode = vlayer = FLAGtop = FLAGbot = Itop = Ibot = 0;
	double tmp1, tmp2;
	double preDENS[N_DISSOLVMATER+1][N_SOILvirtLAYERS], actDENS[N_DISSOLVMATER+1][N_SOILvirtLAYERS]; // N_DISSOLVEDMATER+ VWC in all  layers (also virtual) - only dissolved matter!
	double VWCsat;
	double content_dissolved, content_bounded;
	double ratioNORM, ratioCAPIL, ratioSAT;

	static double Ztop_pre[N_SOILLAYERS + 2];
	static double Zbot_pre[N_SOILLAYERS + 2];
	static int virtLayer_pre[N_SOILvirtLAYERS];
	static int Nvl_act;
	static int Nvl_pre;


	GWlayer = (int)sprop->GWlayer;
	CFlayer = (int)sprop->CFlayer;

	for (vlayer = 0; vlayer < N_SOILvirtLAYERS; vlayer++)
	{
		sprop->Zbot[vlayer] = DATA_GAP;
		sprop->Ztop[vlayer] = DATA_GAP;
		sprop->virtLayer[vlayer] = DATA_GAP;
		for (dm = 0; dm < N_DISSOLVMATER + 1; dm++) preDENS[dm][vlayer] = DATA_GAP;
		for (dm = 0; dm < N_DISSOLVMATER + 1; dm++) actDENS[dm][vlayer] = DATA_GAP;
	}

	/* Save soilw data for leaching calculation  */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		ws->soilw_pre[layer] = ws->soilw[layer];
		for (dm = 0; dm < N_DISSOLVMATER; dm++) soilInfo->content_soil_pre[dm][layer] = soilInfo->content_soil[dm][layer];
	}


	/*----------------------------------------------------------------------------------------------*/
	/*  I. layer numbering regarding pre arrays (on first day, simple VWC-change is calculated): 0, 1, 2, 3 - normal layers, -1: NORMcf, -2 CAPILcf, -3: NORMgw, -4: CAPILgw */
	/*----------------------------------------------------------------------------------------------*/

	vlayer = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		/* if GW layer  */
		if (layer == GWlayer)
		{

			/* if NORMgw exists */
			if (sprop->dz_NORMgw)
			{
				/* NORMgw */
				sprop->Zbot[vlayer] = sprop->CFD;
				sprop->Ztop[vlayer] = sitec->soillayer_depth[layer] - sitec->soillayer_thickness[layer];
				sprop->virtLayer[vlayer] = -3;
				vlayer += 1;
			}

			/* if CAPILgw exists */
			if (sprop->dz_CAPILgw)
			{
				/* CAPILgw */
				sprop->Zbot[vlayer] = sprop->GWD;
				if (sprop->dz_NORMgw)
					sprop->Ztop[vlayer] = sprop->CFD;
				else
					sprop->Ztop[vlayer] = sitec->soillayer_depth[layer] - sitec->soillayer_thickness[layer];

				sprop->virtLayer[vlayer] = -4;
				vlayer += 1;
			}

			/* SATgw */
			sprop->Zbot[vlayer] = sitec->soillayer_depth[layer];
			sprop->Ztop[vlayer] = sprop->GWD;
			sprop->virtLayer[vlayer] = -5;
			vlayer += 1;

		}
		else
		{
			/* if CF layer  */
			if (layer == CFlayer && CFlayer != GWlayer)
			{
				/* if NORMcf exists */
				if (sprop->dz_NORMcf)
				{
					/* NORMcf */
					sprop->Zbot[vlayer] = sprop->CFD;
					sprop->Ztop[vlayer] = sitec->soillayer_depth[layer] - sitec->soillayer_thickness[layer];
					sprop->virtLayer[vlayer] = -1;
					vlayer += 1;
				}

				/* CAPILcf */
				sprop->Zbot[vlayer] = sitec->soillayer_depth[layer];
				sprop->Ztop[vlayer] = sprop->CFD;
				sprop->virtLayer[vlayer] = -2;
				vlayer += 1;
			}
			/* normal layers   */
			else
			{
				sprop->Zbot[vlayer] = sitec->soillayer_depth[layer];
				sprop->Ztop[vlayer] = sitec->soillayer_depth[layer] - sitec->soillayer_thickness[layer];
				sprop->virtLayer[vlayer] = layer;
				vlayer += 1;
			}

		}

		if (layer == 9) Nvl_act = vlayer;

	}

	/*----------------------------------------------------------------------------------------------*/
	/* II. on  first GW-simulation day: initialization: : saturated layer below GWD (DENSentration from GW), mixesd layer in GW-layer, and unchanged normal layer */
	/*----------------------------------------------------------------------------------------------*/
	if (sprop->GWD_pre == DATA_GAP || sprop->GWD_pre >= 10)
	{
		Nvl_pre = N_SOILLAYERS;
		for (vlayer = 0; vlayer < N_SOILLAYERS; vlayer++)
		{
			virtLayer_pre[vlayer] = vlayer;
			Zbot_pre[vlayer] = sitec->soillayer_depth[vlayer];
			Ztop_pre[vlayer] = sitec->soillayer_depth[vlayer] - sitec->soillayer_thickness[vlayer];

			if (vlayer == CFlayer && vlayer != GWlayer)
			{
				ratioNORM = sprop->dz_NORMcf / (sprop->dz_CAPILcf + sprop->dz_NORMcf);
				ratioCAPIL = sprop->dz_CAPILcf / (sprop->dz_CAPILcf + sprop->dz_NORMcf);
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					soilInfo->content_NORMcf[dm] = soilInfo->content_soil[dm][vlayer] * ratioNORM;
					soilInfo->content_CAPILcf[dm] = soilInfo->content_soil[dm][vlayer] * ratioCAPIL;
				}
			}

			if (vlayer == GWlayer)
			{
				ratioNORM = sprop->dz_NORMgw / (sprop->dz_CAPILgw + sprop->dz_NORMgw + sprop->dz_SATgw);
				ratioCAPIL = sprop->dz_CAPILgw / (sprop->dz_CAPILgw + sprop->dz_NORMgw + sprop->dz_SATgw);
				ratioSAT = sprop->dz_SATgw / (sprop->dz_CAPILgw + sprop->dz_NORMgw + sprop->dz_SATgw);
				for (dm=0; dm< N_DISSOLVMATER; dm++)
				{ 
					soilInfo->content_NORMgw[dm] = soilInfo->content_soil[dm][vlayer] * ratioNORM;
					soilInfo->content_CAPILgw[dm] = soilInfo->content_soil[dm][vlayer] * ratioCAPIL;
					soilInfo->content_SATgw[dm] = soilInfo->content_soil[dm][vlayer] * ratioSAT;
				}
			}
		}
		for (vlayer = N_SOILLAYERS; vlayer < N_SOILvirtLAYERS; vlayer++)
		{
			virtLayer_pre[vlayer] = DATA_GAP;
			Zbot_pre[vlayer] = DATA_GAP;
			Ztop_pre[vlayer] = DATA_GAP;
		}

	}
	/*----------------------------------------------------------------------------------------------*/
	/* III. on  non-first GW-simulation day: calculation of mixture: pre and act values */
	/*----------------------------------------------------------------------------------------------*/
	/* -------------------------------- */
	/* III./1. conservation of pre values  */
	layer = 0;
	for (vlayer = 0; vlayer < Nvl_pre; vlayer++)
	{

		/* dissolved materials (DENSentration from content) and VWC: pre values  */
		if (virtLayer_pre[vlayer] >= 0)
		{
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				preDENS[dm][vlayer] = (soilInfo->content_soil[dm][layer] * soilInfo->dissolv_prop[dm]) / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
			}
			preDENS[dm][vlayer] = epv->VWC[layer] * water_density;
			layer += 1;
		}
		else
		{
			/* NORMcf */
			if (virtLayer_pre[vlayer] == -1)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++) preDENS[dm][vlayer] = (soilInfo->content_NORMcf_pre[dm] * soilInfo->dissolv_prop[dm]) / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				preDENS[dm][vlayer] = sprop->soilw_NORMcf_pre / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
			}

			/* CAPILcf */
			if (virtLayer_pre[vlayer] == -2)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++) preDENS[dm][vlayer] = (soilInfo->content_CAPILcf_pre[dm] * soilInfo->dissolv_prop[dm]) / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				preDENS[dm][vlayer] = sprop->soilw_CAPILcf_pre / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				layer = layer + 1;

			}
			if (virtLayer_pre[vlayer] == -3)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++) preDENS[dm][vlayer] = (soilInfo->content_NORMgw_pre[dm] * soilInfo->dissolv_prop[dm]) / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				preDENS[dm][vlayer] = sprop->soilw_NORMgw_pre / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
			}
			if (virtLayer_pre[vlayer] == -4)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++) preDENS[dm][vlayer] = (soilInfo->content_CAPILgw_pre[dm] * soilInfo->dissolv_prop[dm]) / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				preDENS[dm][vlayer] = sprop->soilw_CAPILgw_pre / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
			}
			if (virtLayer_pre[vlayer] == -5)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++) preDENS[dm][vlayer] = (soilInfo->content_SATgw_pre[dm] * soilInfo->dissolv_prop[dm]) / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				preDENS[dm][vlayer] = sprop->soilw_SATgw_pre / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				layer += 1;
			}
		}

	}


	/* -------------------------------- */
	/* III/2. calculation of new dissolves material and VWC values */
	for (vlayer = 0; vlayer < Nvl_act; vlayer++)
	{
		/* new layer in pre layer */
		FLAGtop = 0;
		FLAGbot = 0;
		vlt = 0;
	
		while (FLAGtop == 0)
		{
			if (sprop->Ztop[vlayer] >= Ztop_pre[vlt] && sprop->Ztop[vlayer] < Zbot_pre[vlt])
			{
				Itop = vlt;
				FLAGtop = 1;
				vlb = vlt;
				while (FLAGbot == 0)
				{
					if ((sprop->Zbot[vlayer] >= Ztop_pre[vlb] || Ztop_pre[vlb] == DATA_GAP) && (sprop->Zbot[vlayer] <= Zbot_pre[vlb] || Zbot_pre[vlb] == DATA_GAP))
					{
						Ibot = vlb;
						FLAGbot = 1;
					}
					else
						vlb += 1;
				}
			}
			else
				vlt += 1;
		}
		
		/* actual values */
		/* in saturated zone of GW-layer and below GWlayer: saturation values, otherwise: combination of layers */
		if (sprop->virtLayer[vlayer] == -5 || (sprop->virtLayer[vlayer] > sprop->GWlayer))
		{
			if (sprop->virtLayer[vlayer] == -5)
				VWCsat = sprop->VWCsat[GWlayer];
			else
				VWCsat = sprop->VWCsat[sprop->virtLayer[vlayer]];

			for (dm = 0; dm < N_DISSOLVMATER + 1; dm++)
			{
				if (dm < N_DISSOLVMATER)
					actDENS[dm][vlayer] = (soilInfo->GWconc[dm] * VWCsat * water_density);
				else
					actDENS[dm][vlayer] = VWCsat * water_density;
			}
		}
		else
		{
			if (Ibot == Itop)
			{
				for (dm = 0; dm < N_DISSOLVMATER + 1; dm++) actDENS[dm][vlayer] = preDENS[dm][Ibot];

			}
			else
			{
				for (dm = 0; dm < N_DISSOLVMATER + 1; dm++)
				{

					tmp1 = ((Zbot_pre[Itop] - sprop->Ztop[vlayer]) * preDENS[dm][Itop] + (sprop->Zbot[vlayer] - Ztop_pre[Ibot]) * preDENS[dm][Ibot]);
					tmp2 = 0;
					for (ll = Itop + 1; ll < Ibot; ll++)
					{
						tmp2 += (sprop->Zbot[ll] - sprop->Ztop[ll]) * preDENS[dm][ll];
					}
					actDENS[dm][vlayer] = (tmp1 + tmp2) / (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
				}

			}
	
		}

	}

	/* -------------------------------- */
	/* III/3. calculation of changes due to GWmoving */

		
	for (dm = 0; dm < N_DISSOLVMATER + 1; dm++)
	{
		layer = 0;
		for (vlayer = 0; vlayer < Nvl_act; vlayer++)
		{
			/* dissolved materials and VWC: pre values  */
			if (sprop->virtLayer[vlayer] == layer)
			{
				if (dm < N_DISSOLVMATER)
				{ 
					content_bounded = soilInfo->content_soil[dm][layer] * (1 - soilInfo->dissolv_prop[dm]);
					content_dissolved = actDENS[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
					soilInfo->content_soil[dm][layer] = content_bounded + content_dissolved;
				}
				else
				{
					ws->soilw[layer] = actDENS[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
					epv->VWC[layer] = actDENS[dm][vlayer] / water_density;
				}
				layer += 1;
			}
			else
			{
				if (sprop->virtLayer[vlayer] == -1)
				{
					if (dm < N_DISSOLVMATER)
					{
						content_bounded = soilInfo->content_NORMcf[dm] * (1 - soilInfo->dissolv_prop[dm]);
						content_dissolved = actDENS[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
						soilInfo->content_NORMcf[dm] = content_bounded + content_dissolved;
					}			
					else
					{				
						sprop->soilw_NORMcf = actDENS[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]) ;
						sprop->VWC_NORMcf = actDENS[dm][vlayer] / water_density;

						if (fabs((sprop->Zbot[vlayer] - sprop->Ztop[vlayer] - sprop->dz_NORMcf)) > CRIT_PREC)
						{
							printf("\n");
							printf("ERROR in virtual layer calculation in groundwater_movement.c\n");
							errorCode = 1;
						}
					}

				}
				if (sprop->virtLayer[vlayer] == -2)
				{
					if (dm < N_DISSOLVMATER)
					{
						content_bounded = soilInfo->content_CAPILcf[dm] * (1 - soilInfo->dissolv_prop[dm]);
						content_dissolved = actDENS[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
						soilInfo->content_CAPILcf[dm] = content_bounded + content_dissolved;
						soilInfo->content_soil[dm][layer] = soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm];
					}
					else
					{		
						sprop->soilw_CAPILcf = actDENS[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
						sprop->VWC_CAPILcf = actDENS[dm][vlayer] / water_density;
						if (fabs((sprop->Zbot[vlayer] - sprop->Ztop[vlayer] - sprop->dz_CAPILcf)) > CRIT_PREC)
						{
							printf("\n");
							printf("ERROR in virtual layer calculation in groundwater_movement.c\n");
							errorCode = 1;
						}
						ws->soilw[layer] = sprop->soilw_NORMcf + sprop->soilw_CAPILcf;
						epv->VWC[layer] = ws->soilw[layer] / (water_density * sitec->soillayer_thickness[layer]);
							
					}
					layer += 1;
				}

				if (sprop->virtLayer[vlayer] == -3)
				{
					if (dm < N_DISSOLVMATER)
					{
						content_bounded = soilInfo->content_NORMgw[dm] * (1 - soilInfo->dissolv_prop[dm]);
						content_dissolved = actDENS[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
						soilInfo->content_NORMgw[dm] = content_bounded + content_dissolved;
					}
					else
					{
							
						sprop->soilw_NORMgw = actDENS[dm][vlayer] * sprop->dz_NORMgw;
						sprop->VWC_NORMgw = actDENS[dm][vlayer] / water_density;

						if (fabs((sprop->Zbot[vlayer] - sprop->Ztop[vlayer] - sprop->dz_NORMgw)) > CRIT_PREC)
						{
							printf("\n");
							printf("ERROR in virtual layer calculation in groundwater_movement.c\n");
							errorCode = 1;
						}
					}
				}
				if (sprop->virtLayer[vlayer] == -4)
				{
					if (dm < N_DISSOLVMATER)
					{
						content_bounded = soilInfo->content_CAPILgw[dm] * (1 - soilInfo->dissolv_prop[dm]);
						content_dissolved = actDENS[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
						soilInfo->content_CAPILgw[dm] = content_bounded + content_dissolved;
					}
					else
					{
							
						sprop->soilw_CAPILgw = actDENS[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
						sprop->VWC_CAPILgw = actDENS[dm][vlayer] / water_density;

						if (fabs((sprop->Zbot[vlayer] - sprop->Ztop[vlayer] - sprop->dz_CAPILgw)) > CRIT_PREC)
						{
							printf("\n");
							printf("ERROR in virtual layer calculation in groundwater_movement.c\n");
							errorCode = 1;
						}
					}
				}
				if (sprop->virtLayer[vlayer] == -5)
				{
					if (dm < N_DISSOLVMATER)
					{
						content_bounded = soilInfo->content_SATgw[dm] * (1 - soilInfo->dissolv_prop[dm]);
						content_dissolved = actDENS[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
						soilInfo->content_SATgw[dm] = content_bounded + content_dissolved;
						soilInfo->content_soil[dm][layer] = soilInfo->content_NORMgw[dm] + soilInfo->content_CAPILgw[dm] + soilInfo->content_SATgw[dm];
					}
					else
					{
						sprop->soilw_SATgw = sprop->VWCsat[layer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]) * water_density;
						if (fabs((sprop->Zbot[vlayer] - sprop->Ztop[vlayer] - sprop->dz_SATgw)) > CRIT_PREC)
						{
							printf("\n");
							printf("ERROR in virtual layer calculation in groundwater_movement.c\n");
							errorCode = 1;
						}
						ws->soilw[layer] = sprop->soilw_NORMgw + sprop->soilw_CAPILgw + sprop->soilw_SATgw;
						epv->VWC[layer] = ws->soilw[layer] / (water_density * sitec->soillayer_thickness[layer]);
					}
					layer += 1;
				}
			}
	
		} /* endfor: layer */

	} /* endfor: dm */


	/*----------------------------------------------------------------------------------------------*/
	/* IV. calculation of change due to GW-moving  */
	/*----------------------------------------------------------------------------------------------*/

	/* transfer value: content_array -> NH4, NO3, DOC, DON */
	if (!errorCode && check_soilcontent(-1, 1, sprop, cs, ns, soilInfo))
	{
		printf("ERROR in check_soilcontent.c for groundwater_movement.c\n");
		errorCode = 1;
	}

	for (layer = 0; layer < N_SOILLAYERS; layer++) wf->GWmovchange += ws->soilw[layer] - ws->soilw_pre[layer];
	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		for (layer = 0; layer < N_SOILLAYERS; layer++) soilInfo->dismatGWmovchange[dm] += soilInfo->content_soil[dm][layer] - soilInfo->content_soil_pre[dm][layer];
	}


	/*----------------------------------------------------------------------------------------------*/
	/* IV. src/snk variables  */
	/*----------------------------------------------------------------------------------------------*/
	
	/* src/snk variables*/
	for (dm = 0; dm < N_DISSOLVN; dm++) soilInfo->dismatGWmovchangeN_total += soilInfo->dismatGWmovchange[dm];
	for (dm = N_DISSOLVN; dm < N_DISSOLVMATER; dm++) soilInfo->dismatGWmovchangeC_total += soilInfo->dismatGWmovchange[dm];
	

	if (soilInfo->dismatGWmovchangeN_total > 0)
		ns->GWsrc_N += soilInfo->dismatGWmovchangeN_total;
	else
		ns->GWsnk_N += -1 * soilInfo->dismatGWmovchangeN_total;

	
	if (soilInfo->dismatGWmovchangeC_total > 0)
		cs->GWsrc_C += soilInfo->dismatGWmovchangeC_total;
	else
		cs->GWsnk_C += -1 * soilInfo->dismatGWmovchangeC_total;


	if (wf->GWmovchange > 0)
		ws->GWsrc_W += wf->GWmovchange;
	else
		ws->GWsnk_W += -1 * wf->GWmovchange;



	Nvl_pre = Nvl_act;

	for (vlayer = 0; vlayer < N_SOILvirtLAYERS; vlayer++)
	{
		Zbot_pre[vlayer] = sprop->Zbot[vlayer];
		Ztop_pre[vlayer] = sprop->Ztop[vlayer];
		virtLayer_pre[vlayer] = sprop->virtLayer[vlayer];


	}

	/* reset variables that are not current */
	if (sprop->dz_NORMcf == 0)
	{
		sprop->soilw_NORMcf = 0;
		sprop->VWC_NORMcf = 0;
		for (dm = 0; dm < N_DISSOLVMATER; dm++) soilInfo->content_NORMcf[dm] = 0;
	}

	if (sprop->dz_CAPILcf == 0)
	{
		sprop->soilw_CAPILcf = 0;
		sprop->VWC_CAPILcf = 0;
		for (dm = 0; dm < N_DISSOLVMATER; dm++) soilInfo->content_CAPILcf[dm] = 0;
	}

	if (sprop->dz_NORMgw == 0)
	{
		sprop->soilw_NORMgw = 0;
		sprop->VWC_NORMgw = 0;
		for (dm = 0; dm < N_DISSOLVMATER; dm++) soilInfo->content_NORMgw[dm] = 0;
	}

	if (sprop->dz_CAPILgw == 0)
	{
		sprop->soilw_CAPILgw = 0;
		sprop->VWC_CAPILgw = 0;
		for (dm = 0; dm < N_DISSOLVMATER; dm++) soilInfo->content_CAPILgw[dm] = 0;
	}




	return (errorCode);
}
