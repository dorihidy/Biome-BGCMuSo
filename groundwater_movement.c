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


int groundwater_movement(const siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, soilInfo_struct* soilInfo, wstate_struct* ws, wflux_struct* wf, cstate_struct* cs, nstate_struct* ns)
{
	/* ------------------------------------------------------------------------------------------------------ */
	/* initialization */

	int errorCode, vlayer, layer, CFlayer, GWlayer, vlt, vlb, ll;
	int Itop, Ibot, FLAGtop, FLAGbot, dm;
	errorCode = vlayer = FLAGtop = FLAGbot = Itop = Ibot = 0;
	double tmp1, tmp2;
	double preDENS_BOUND[N_DISSOLVMATER][N_SOILvirtLAYERS], actDENS_BOUND[N_DISSOLVMATER+1][N_SOILvirtLAYERS]; /* N_DISSOLVEDMATER in all  layers(also virtual) - only bound matter! */ 
	double preDENS_DISSOLV[N_DISSOLVMATER + 1][N_SOILvirtLAYERS], actDENS_DISSOLV[N_DISSOLVMATER + 1][N_SOILvirtLAYERS]; /*  N_DISSOLVEDMATER+ VWC in all  layers (also virtual) - only dissolved matter! */
	double VWCsat;
	double ratioNORM, ratioCAPIL, ratioSAT;

	static double Ztop_pre[N_SOILLAYERS + 2];
	static double Zbot_pre[N_SOILLAYERS + 2];
	static double soilw_NORMcf_pre;        
	static double  soilw_CAPILcf_pre;       
	static double  soilw_NORMgw_pre;          
	static double  soilw_CAPILgw_pre;         
	static double  soilw_SATgw_pre;        
	static double  contentBOUND_NORMgw_pre[N_DISSOLVMATER];                  
	static double  contentBOUND_CAPILgw_pre[N_DISSOLVMATER];
	static double  contentBOUND_SATgw_pre[N_DISSOLVMATER];
	static double  contentBOUND_NORMcf_pre[N_DISSOLVMATER];
	static double  contentBOUND_CAPILcf_pre[N_DISSOLVMATER];
	static double  contentDISSOLV_NORMgw_pre[N_DISSOLVMATER];
	static double  contentDISSOLV_CAPILgw_pre[N_DISSOLVMATER];
	static double  contentDISSOLV_SATgw_pre[N_DISSOLVMATER];
	static double  contentDISSOLV_NORMcf_pre[N_DISSOLVMATER];
	static double  contentDISSOLV_CAPILcf_pre[N_DISSOLVMATER];
	static int virtLayer_pre[N_SOILvirtLAYERS];
	static int Nvl_act;
	static int Nvl_pre;


	GWlayer = (int)sprop->GWlayer;
	CFlayer = (int)sprop->CFlayer;

	/* initalizing */
	for (vlayer = 0; vlayer < N_SOILvirtLAYERS; vlayer++)
	{
		sprop->Zbot[vlayer] = DATA_GAP;
		sprop->Ztop[vlayer] = DATA_GAP;
		sprop->virtLayer[vlayer] = DATA_GAP;
		for (dm = 0; dm < N_DISSOLVMATER + 1; dm++)
		{
			if (dm < N_DISSOLVMATER)
			{
				preDENS_BOUND[dm][vlayer] = DATA_GAP;
				actDENS_BOUND[dm][vlayer] = DATA_GAP;
			}
			preDENS_DISSOLV[dm][vlayer] = DATA_GAP;
			actDENS_DISSOLV[dm][vlayer] = DATA_GAP;
		}
	}

	/*----------------------------------------------------------------------------------------------*/
	/*  0. Initializing pre values and actual values */
	/*----------------------------------------------------------------------------------------------*/
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		ws->soilw_pre[layer] = ws->soilw[layer];
		for (dm = 0; dm < N_DISSOLVMATER; dm++) soilInfo->content_soil_pre[dm][layer] = soilInfo->content_soil[dm][layer];
	}

	/* initializing pre values and avtualvalues  */
	soilw_NORMcf_pre = sprop->soilw_NORMcf;
	soilw_CAPILcf_pre = sprop->soilw_CAPILcf;
	soilw_NORMgw_pre = sprop->soilw_NORMgw;
	soilw_CAPILgw_pre = sprop->soilw_CAPILgw;
	soilw_SATgw_pre = sprop->soilw_SATgw;

	sprop->soilw_NORMcf = 0;
	sprop->soilw_CAPILcf = 0;
	sprop->soilw_NORMgw = 0;
	sprop->soilw_CAPILgw = 0;
	sprop->soilw_SATgw = 0;

	sprop->VWC_NORMcf = 0;
	sprop->VWC_CAPILcf = 0;
	sprop->VWC_NORMgw = 0;
	sprop->VWC_CAPILgw = 0;
	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		contentBOUND_NORMcf_pre[dm] = soilInfo->contentBOUND_NORMcf[dm];
		contentBOUND_CAPILcf_pre[dm] = soilInfo->contentBOUND_CAPILcf[dm];
		contentBOUND_NORMgw_pre[dm] = soilInfo->contentBOUND_NORMgw[dm];
		contentBOUND_CAPILgw_pre[dm] = soilInfo->contentBOUND_CAPILgw[dm];
		contentBOUND_SATgw_pre[dm] = soilInfo->contentBOUND_SATgw[dm];

		contentDISSOLV_NORMcf_pre[dm] = soilInfo->contentDISSOLV_NORMcf[dm];
		contentDISSOLV_CAPILcf_pre[dm] = soilInfo->contentDISSOLV_CAPILcf[dm];
		contentDISSOLV_NORMgw_pre[dm] = soilInfo->contentDISSOLV_NORMgw[dm];
		contentDISSOLV_CAPILgw_pre[dm] = soilInfo->contentDISSOLV_CAPILgw[dm];
		contentDISSOLV_SATgw_pre[dm] = soilInfo->contentDISSOLV_SATgw[dm];

		soilInfo->contentBOUND_NORMcf[dm] = 0;
		soilInfo->contentBOUND_CAPILcf[dm] = 0;
		soilInfo->contentBOUND_NORMgw[dm] = 0;
		soilInfo->contentBOUND_CAPILgw[dm] = 0;
		soilInfo->contentBOUND_SATgw[dm] = 0;

		soilInfo->contentDISSOLV_NORMcf[dm] = 0;
		soilInfo->contentDISSOLV_CAPILcf[dm] = 0;
		soilInfo->contentDISSOLV_NORMgw[dm] = 0;
		soilInfo->contentDISSOLV_CAPILgw[dm] = 0;
		soilInfo->contentDISSOLV_SATgw[dm] = 0;

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


				/* if CAPILcf exists */
				if (sprop->dz_CAPILcf)
				{
					sprop->Zbot[vlayer] = sitec->soillayer_depth[layer];
					sprop->Ztop[vlayer] = sprop->CFD;
					sprop->virtLayer[vlayer] = -2;
					vlayer += 1;
				}

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
/* II. on  first GW-simulation day: initialization: : saturated layer below GWD (concentration from GW), mixesd layer in GW-layer, and unchanged normal layer */
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
					contentBOUND_NORMcf_pre[dm] = soilInfo->contentBOUND_soil[dm][vlayer] * ratioNORM;
					contentBOUND_CAPILcf_pre[dm] = soilInfo->contentBOUND_soil[dm][vlayer] * ratioCAPIL;

					contentDISSOLV_NORMcf_pre[dm] = soilInfo->contentDISSOLV_soil[dm][vlayer] * ratioNORM;
					contentDISSOLV_CAPILcf_pre[dm] = soilInfo->contentDISSOLV_soil[dm][vlayer] * ratioCAPIL;
				}
				soilw_NORMcf_pre = ws->soilw[vlayer] * ratioNORM;
				soilw_CAPILcf_pre = ws->soilw[vlayer] * ratioCAPIL;
			}

			if (vlayer == GWlayer)
			{
				ratioNORM = sprop->dz_NORMgw / (sprop->dz_CAPILgw + sprop->dz_NORMgw + sprop->dz_SATgw);
				ratioCAPIL = sprop->dz_CAPILgw / (sprop->dz_CAPILgw + sprop->dz_NORMgw + sprop->dz_SATgw);
				ratioSAT = sprop->dz_SATgw / (sprop->dz_CAPILgw + sprop->dz_NORMgw + sprop->dz_SATgw);
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					contentBOUND_NORMgw_pre[dm] = soilInfo->contentBOUND_soil[dm][vlayer] * ratioNORM;
					contentBOUND_CAPILgw_pre[dm] = soilInfo->contentBOUND_soil[dm][vlayer] * ratioCAPIL;
					contentBOUND_SATgw_pre[dm] = soilInfo->contentBOUND_soil[dm][vlayer] * ratioSAT;

					contentDISSOLV_NORMgw_pre[dm] = soilInfo->contentDISSOLV_soil[dm][vlayer] * ratioNORM;
					contentDISSOLV_CAPILgw_pre[dm] = soilInfo->contentDISSOLV_soil[dm][vlayer] * ratioCAPIL;
					contentDISSOLV_SATgw_pre[dm] = soilInfo->contentDISSOLV_soil[dm][vlayer] * ratioSAT;
				}
				soilw_NORMgw_pre = ws->soilw[vlayer] * ratioNORM;
				soilw_CAPILgw_pre = ws->soilw[vlayer] * ratioCAPIL;
				soilw_SATgw_pre = ws->soilw[vlayer] * ratioSAT;
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
	/* -------------------------------- */
	layer = 0;
	for (vlayer = 0; vlayer < Nvl_pre; vlayer++)
	{

		/* dissolved materials (concentration from content) and VWC: pre values  */
		if (virtLayer_pre[vlayer] >= 0)
		{
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				preDENS_BOUND[dm][vlayer] = soilInfo->contentBOUND_soil[dm][layer] / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				preDENS_DISSOLV[dm][vlayer] = soilInfo->contentDISSOLV_soil[dm][layer] / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
			}
			preDENS_DISSOLV[dm][vlayer] = epv->VWC[layer] * water_density;
			layer += 1;
		}
		else
		{
			/* NORMcf */
			if (virtLayer_pre[vlayer] == -1)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					preDENS_BOUND[dm][vlayer] = contentBOUND_NORMcf_pre[dm] / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
					preDENS_DISSOLV[dm][vlayer] = contentDISSOLV_NORMcf_pre[dm] / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				}
				preDENS_DISSOLV[dm][vlayer] = soilw_NORMcf_pre / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				if (soilw_CAPILcf_pre == 0) layer = layer + 1;
			}

			/* CAPILcf */
			if (virtLayer_pre[vlayer] == -2)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					preDENS_BOUND[dm][vlayer] = contentBOUND_CAPILcf_pre[dm] / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
					preDENS_DISSOLV[dm][vlayer] = contentDISSOLV_CAPILcf_pre[dm] / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				}
				preDENS_DISSOLV[dm][vlayer] = soilw_CAPILcf_pre / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				layer = layer + 1;

			}
			/* NORMgw */
			if (virtLayer_pre[vlayer] == -3)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					preDENS_BOUND[dm][vlayer] = contentBOUND_NORMgw_pre[dm] / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
					preDENS_DISSOLV[dm][vlayer] = contentDISSOLV_NORMgw_pre[dm] / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				}
				preDENS_DISSOLV[dm][vlayer] = soilw_NORMgw_pre / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
			}

			/* CAPILgw */
			if (virtLayer_pre[vlayer] == -4)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					preDENS_BOUND[dm][vlayer] = contentBOUND_CAPILgw_pre[dm] / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
					preDENS_DISSOLV[dm][vlayer] = contentDISSOLV_CAPILgw_pre[dm] / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				}
				preDENS_DISSOLV[dm][vlayer] = soilw_CAPILgw_pre / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
			}

			/* SATgw */
			if (virtLayer_pre[vlayer] == -5)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					preDENS_BOUND[dm][vlayer] = contentBOUND_SATgw_pre[dm] / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
					preDENS_DISSOLV[dm][vlayer] = contentDISSOLV_SATgw_pre[dm] / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				}
				preDENS_DISSOLV[dm][vlayer] = soilw_SATgw_pre / (Zbot_pre[vlayer] - Ztop_pre[vlayer]);
				layer += 1;
			}
		}

	}


	/* -------------------------------- */
	/* III/2. calculation of actDENS from preDENS values */
	/* -------------------------------- */
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
		if (Ibot == Itop)
		{
			/* BOUND material */
			for (dm = 0; dm < N_DISSOLVMATER; dm++) actDENS_BOUND[dm][vlayer] = preDENS_BOUND[dm][Ibot];
			
		
			/* DISSOLV material: in saturated zone of GW-layer and below GWlayer: saturation values, otherwise: preDENS */
			if (sprop->virtLayer[vlayer] == -5 || (sprop->virtLayer[vlayer] > sprop->GWlayer))
			{
				if (sprop->virtLayer[vlayer] == -5)
					VWCsat = sprop->VWCsat[GWlayer];
				else
					VWCsat = sprop->VWCsat[sprop->virtLayer[vlayer]];

				for (dm = 0; dm < N_DISSOLVMATER+1; dm++)
				{
					if (dm < N_DISSOLVMATER)
						actDENS_DISSOLV[dm][vlayer] = (soilInfo->GWconc[dm] * VWCsat * water_density);
					else
						actDENS_DISSOLV[dm][vlayer] = VWCsat * water_density;
				}
			}
			else
				for (dm = 0; dm < N_DISSOLVMATER+1; dm++) actDENS_DISSOLV[dm][vlayer] = preDENS_DISSOLV[dm][Ibot];
			

		}
		else
		{
			
			/* BOUND material */
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				tmp1 = ((Zbot_pre[Itop] - sprop->Ztop[vlayer]) * preDENS_BOUND[dm][Itop] + (sprop->Zbot[vlayer] - Ztop_pre[Ibot]) * preDENS_BOUND[dm][Ibot]);
				tmp2 = 0;
				for (ll = Itop + 1; ll < Ibot; ll++)
				{
					tmp2 += (sprop->Zbot[ll] - sprop->Ztop[ll]) * preDENS_BOUND[dm][ll];
				}
				actDENS_BOUND[dm][vlayer] = (tmp1 + tmp2) / (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
			}

			/* DISSOLV material: in saturated zone of GW-layer and below GWlayer: saturation values, otherwise: combination of layers */
			if (sprop->virtLayer[vlayer] == -5 || (sprop->virtLayer[vlayer] > sprop->GWlayer))
			{
				if (sprop->virtLayer[vlayer] == -5)
					VWCsat = sprop->VWCsat[GWlayer];
				else
					VWCsat = sprop->VWCsat[sprop->virtLayer[vlayer]];

				for (dm = 0; dm < N_DISSOLVMATER + 1; dm++)
				{
					if (dm < N_DISSOLVMATER)
						actDENS_DISSOLV[dm][vlayer] = (soilInfo->GWconc[dm] * VWCsat * water_density);
					else
						actDENS_DISSOLV[dm][vlayer] = VWCsat * water_density;
				}
			}
			else
			{
				for (dm = 0; dm < N_DISSOLVMATER + 1; dm++)
				{
					tmp1 = ((Zbot_pre[Itop] - sprop->Ztop[vlayer]) * preDENS_DISSOLV[dm][Itop] + (sprop->Zbot[vlayer] - Ztop_pre[Ibot]) * preDENS_DISSOLV[dm][Ibot]);
					tmp2 = 0;
					for (ll = Itop + 1; ll < Ibot; ll++)
					{
						tmp2 += (sprop->Zbot[ll] - sprop->Ztop[ll]) * preDENS_DISSOLV[dm][ll];
					}
					actDENS_DISSOLV[dm][vlayer] = (tmp1 + tmp2) / (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
				}
			}

		} /* endelse Ibot == Itop */
	}/* endfor layer */


		/* -------------------------------- */
	/* III/3. calculation of new dissolves material and VWC values */
	/* -------------------------------- */
	layer = 0;
	for (vlayer = 0; vlayer < Nvl_act; vlayer++)
	{

		/* in normal layers  */
		if (sprop->virtLayer[vlayer] == layer)
		{
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				soilInfo->contentBOUND_soil[dm][layer] = actDENS_BOUND[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
				soilInfo->contentDISSOLV_soil[dm][layer] = actDENS_DISSOLV[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
			}
			ws->soilw[layer] = actDENS_DISSOLV[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
			epv->VWC[layer] = ws->soilw[layer] / (water_density * sitec->soillayer_thickness[layer]);
			layer += 1;
		}
		else
		{

			/* calculation of special layers */
			/* NORMcf */
			if (sprop->virtLayer[vlayer] == -1)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					soilInfo->contentBOUND_NORMcf[dm] = actDENS_BOUND[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
					soilInfo->contentDISSOLV_NORMcf[dm] = actDENS_DISSOLV[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
				}
				sprop->soilw_NORMcf = actDENS_DISSOLV[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
				sprop->VWC_NORMcf = sprop->soilw_NORMcf / (water_density * sprop->dz_NORMcf);

				if (sprop->dz_CAPILcf == 0) layer = layer + 1;
			}

			/* CAPILcf */
			if (sprop->virtLayer[vlayer] == -2)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					soilInfo->contentBOUND_CAPILcf[dm] = actDENS_BOUND[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
					soilInfo->contentDISSOLV_CAPILcf[dm] = actDENS_DISSOLV[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
					soilInfo->contentBOUND_soil[dm][layer] = soilInfo->contentBOUND_NORMcf[dm] + soilInfo->contentBOUND_CAPILcf[dm];
					soilInfo->contentDISSOLV_soil[dm][layer] = soilInfo->contentDISSOLV_NORMcf[dm] + soilInfo->contentDISSOLV_CAPILcf[dm];
				}


				sprop->soilw_CAPILcf = actDENS_DISSOLV[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
				sprop->VWC_CAPILcf = sprop->soilw_CAPILcf / (water_density * sprop->dz_CAPILcf);

				ws->soilw[layer] = sprop->soilw_CAPILcf + sprop->soilw_NORMcf;
				epv->VWC[layer] = ws->soilw[layer] / (water_density * sitec->soillayer_thickness[layer]);

				layer = layer + 1;
			}
			/* NORMgw */
			if (sprop->virtLayer[vlayer] == -3)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					soilInfo->contentBOUND_NORMgw[dm] = actDENS_BOUND[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
					soilInfo->contentDISSOLV_NORMgw[dm] = actDENS_DISSOLV[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
				}
				sprop->soilw_NORMgw = actDENS_DISSOLV[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
				sprop->VWC_NORMgw = sprop->soilw_NORMgw / (water_density * sprop->dz_NORMgw);
			}

			/* CAPILgw */
			if (sprop->virtLayer[vlayer] == -4)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					soilInfo->contentBOUND_CAPILgw[dm] = actDENS_BOUND[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
					soilInfo->contentDISSOLV_CAPILgw[dm] = actDENS_DISSOLV[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
				}
				sprop->soilw_CAPILgw = actDENS_DISSOLV[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
				sprop->VWC_CAPILgw = sprop->soilw_CAPILgw / (water_density * sprop->dz_CAPILgw);
			}

			/* SATgw */
			if (sprop->virtLayer[vlayer] == -5)
			{
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					soilInfo->contentBOUND_SATgw[dm] = actDENS_BOUND[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
					soilInfo->contentDISSOLV_SATgw[dm] = actDENS_DISSOLV[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
					soilInfo->contentBOUND_soil[dm][layer] = soilInfo->contentBOUND_NORMgw[dm] + soilInfo->contentBOUND_CAPILgw[dm] + soilInfo->contentBOUND_SATgw[dm];
					soilInfo->contentDISSOLV_soil[dm][layer] = soilInfo->contentDISSOLV_NORMgw[dm] + soilInfo->contentDISSOLV_CAPILgw[dm] + soilInfo->contentDISSOLV_SATgw[dm];
				}


				sprop->soilw_SATgw = actDENS_DISSOLV[dm][vlayer] * (sprop->Zbot[vlayer] - sprop->Ztop[vlayer]);
				ws->soilw[layer] = sprop->soilw_SATgw + sprop->soilw_CAPILgw + sprop->soilw_NORMgw;
				epv->VWC[layer] = ws->soilw[layer] / (water_density * sitec->soillayer_thickness[layer]);

				layer = layer + 1;
			}
		} /* endelse: special layers */
	} /* endfor vlayers */

	/*----------------------------------------------------------------------------------------------*/
	/* IV. calculation of change due to GW-moving  */
	/*----------------------------------------------------------------------------------------------*/

	/*---------------------------------------------------------------------------------*/
	/* calculating the difference */
	if (!errorCode && calc_DISSOLVandBOUND(1, 0, sprop, soilInfo))
	{
		printf("ERROR in calc_DISSOLVandBOUND.c for multilayer_leaching.c\n");
		errorCode = 1;
	}

	/* transfer value: content_array -> NH4, NO3, DOC, DON */
	if (!errorCode && check_soilcontent(-1, 1, sprop, cs, ns, soilInfo))
	{
		printf("ERROR in check_soilcontent.c for groundwater_movement.c\n");
		errorCode = 1;
	}

	for (layer = 0; layer < N_SOILLAYERS; layer++) wf->GWmovchange += ws->soilw[layer] - ws->soilw_pre[layer];
	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		for (layer = 0; layer < N_SOILLAYERS; layer++) soilInfo->dismatGWmovchange[dm][layer] = soilInfo->content_soil[dm][layer] - soilInfo->content_soil_pre[dm][layer];
	}


	/*----------------------------------------------------------------------------------------------*/
	/* V. src/snk variables  */
	/*----------------------------------------------------------------------------------------------*/

	/* src/snk variables*/
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		for (dm = 0; dm < N_DISSOLVN; dm++) soilInfo->dismatGWmovchangeN_total += soilInfo->dismatGWmovchange[dm][layer];
		for (dm = N_DISSOLVN; dm < N_DISSOLVMATER; dm++) soilInfo->dismatGWmovchangeC_total += soilInfo->dismatGWmovchange[dm][layer];
	}


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



	/* pre values */
	Nvl_pre = Nvl_act;

	for (vlayer = 0; vlayer < N_SOILvirtLAYERS; vlayer++)
	{
		Zbot_pre[vlayer] = sprop->Zbot[vlayer];
		Ztop_pre[vlayer] = sprop->Ztop[vlayer];
		virtLayer_pre[vlayer] = sprop->virtLayer[vlayer];


	}




	return (errorCode);
}
