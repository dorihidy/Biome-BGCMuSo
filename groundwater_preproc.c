/* 
groundwater_preproc.c
calculate GWD, CFD, GWeff, CFeff and the GWdicharge/GWrecharge/GWmovchange in function of GWD

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
#include "bgc_struct.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"
#include "bgc_func.h"


int groundwater_preproc(control_struct* ctrl, const siteconst_struct* sitec, const groundwaterINIT_struct* GWS, soilprop_struct* sprop, epvar_struct* epv, soilInfo_struct* soilInfo, 
	                     wstate_struct* ws, wflux_struct* wf, cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf)
{
	int layerSAT,errorCode,layer,md, year, GWlayer, CFlayer;
	double soilw_sat, CFact, VWCfc_CF;
	double GWboundL, GWboundU, ratio, GWflux;;
	double VWCnew, GWeff, CFeff;

	/*------------------------------------*/
	/* Initialization */
	
	layer=layerSAT=errorCode=0;
	GWlayer=CFlayer = DATA_GAP;

	md = GWS->mgmdGW-1;
	year = ctrl->simstartyear + ctrl->simyr;


	
	/*------------------------------------*/
	/* II. Groundwater calculation only if GWD data is available from GW-file */

	if (GWS->GWnum != 0)	
	{	
		/*-----------------------------------------------------------------------------*/
		/* 1. Determination of groundwater days */
	
		if (md >= 0 && year == GWS->GWyear_array[md] && ctrl->month == GWS->GWmonth_array[md] && ctrl->day == GWS->GWday_array[md])
		{
			sprop->GWD = GWS->GWdepth_array[md];

			/* temporary solution in order to avoid the coincidence of GWD and soil layer */
			for (layer = 0; layer < N_SOILLAYERS; layer++)
			{
				if (fabs(sprop->GWD - sitec->soillayer_depth[layer]) < 0.0001)
				{
					sprop->GWD = sitec->soillayer_depth[layer] + 0.0001;
				}
			}
	
			/* GW above the surface: pond water formation (limitation: maximal pond water height */
			if (sprop->GWD < 0)
			{	
				ws->GW_waterlogging = -1*sprop->GWD*1000;
				if (ws->pondw < ws->GW_waterlogging)
					wf->GW_to_pondw = -1 * sprop->GWD * 1000 - ws->pondw;
				else
					wf->GW_to_pondw = 0;

				sprop->GWD          = 0;
			}
		}
		else
		{
			sprop->GWD = 100;
			for (layer = 0; layer < N_SOILLAYERS; layer++) sprop->VWCfc[layer] = sprop->VWCfc_base[layer];
			sprop->CFlayer = DATA_GAP;
			sprop->preGWD = DATA_GAP;
			sprop->preGWlayer = DATA_GAP;
		}


		/*-----------------------------------------------------------------------------*/
		/* 2. Calculation of actual value of the depth of groundwater and capillary zone and GWmovchange - necessary only if GWD has changed */

		if (sprop->GWD != sprop->preGWD)
		{

			/* lower and upper boundary of the laye */
			layer = N_SOILLAYERS-1;
			while (GWlayer == DATA_GAP && layer >=0)
			{
				GWboundL = sitec->soillayer_depth[layer];
				if (layer == 0)
					GWboundU  = -1*CRIT_PREC;
				else
					GWboundU  = sitec->soillayer_depth[layer-1];
				
				/* if groundwater table is in actual layer (above lower boundary): GWlayer   */
				if (sprop->GWD >= sitec->soillayer_depth[N_SOILLAYERS-1] || (sprop->GWD <= GWboundL && sprop->GWD > GWboundU)) 
				{
					/* if GW-table is under 10 meter */
					if (sprop->GWD > sitec->soillayer_depth[N_SOILLAYERS-1])
					{
						GWlayer=N_SOILLAYERS;
						CFact=sprop->CapillFringe[N_SOILLAYERS-1];
					}
					else
					{
						GWlayer=layer;
						CFact=sprop->CapillFringe[layer];
					}
			
					/* searching of upper boundary of capillary zone */
					if (sprop->GWD - CFact < GWboundU)
					{			
						layer -= 1;
						CFact  = sprop->CapillFringe[layer] - (sprop->GWD - GWboundU);
						while (CFlayer == DATA_GAP && layer >=0)
						{
							GWboundL = sitec->soillayer_depth[layer];
							if (layer == 0)
								GWboundU  = -1*CRIT_PREC;
							else
								GWboundU  = sitec->soillayer_depth[layer-1];

							if (GWboundL - CFact < GWboundU)
							{
								layer -= 1;
								CFact  = sprop->CapillFringe[layer] - (sprop->GWD - GWboundU);
								if (CFact < 0)
								{
									CFlayer = layer+1;
									sprop->CFD = GWboundU; 
								}
							}
							else
							{
								CFlayer = layer;
								if (CFact > 0)
									sprop->CFD = GWboundL - CFact;
								else
									sprop->CFD = GWboundL;
							}
						}
					}
					else
					{
						CFlayer = layer;
						sprop->CFD = sprop->GWD - CFact;
						if (sprop->CFD >= sitec->soillayer_depth[N_SOILLAYERS-1]) CFlayer = N_SOILLAYERS;

					}

				}
				else
					layer -= 1;

			}
			
			/* 3. define GWlayer and actual CF (maximum of capillary fringe upper layer is the top soil zone */
			if (CFlayer == DATA_GAP) 
			{
				CFlayer = 0;
				sprop->CFD = 0;
			}
			sprop->CFlayer=(double) CFlayer;
			sprop->GWlayer=(double) GWlayer;

		
			/* 4. soil layers below the groundwater table are saturated - net water gain from soil system */
			for (layer = N_SOILLAYERS-1; layer >= 0; layer--)
			{
				GWboundL = sitec->soillayer_depth[layer];
			
				if (layer == 0)
					GWboundU  = -1*CRIT_PREC;
				else
					GWboundU  = sitec->soillayer_depth[layer-1];

				/* 4.1 BELOW the GW-table: saturation */
				if (layer > sprop->GWlayer)
				{
					soilw_sat = sprop->VWCsat[layer] * sitec->soillayer_thickness[layer] * water_density;

					/* GW fills the layer to saturation value */
					GWflux = soilw_sat - ws->soilw[layer];

					if (GWflux < 0)
					{
						if (fabs(GWflux > CRIT_PREC) && !errorCode)
						{
							printf("ERROR in groundwater_preproc.c: GWflux is invalid \n");
							errorCode = 1;
						}
						else
							GWflux = 0;
					}
				
					wf->GWmovchange[layer] = GWflux;

						
					ws->soilw[layer]        = soilw_sat;
					epv->VWC[layer]         = sprop->VWCsat[layer];
					sprop->VWCfc[layer]     = sprop->VWCsat[layer];	 

				}
				/* 4.2 IN and ABOVE the GW-table: saturation */
				else
				{
					/* 4.2.1 IN the GW-layer: saturation proportionally to the GW depth, higher FC-values, plus water from GW */
					if (sprop->GWlayer == (double)layer)
					{
						/* saturation proportionally to the GW depth */
						GWeff= (GWboundL - sprop->GWD) / sitec->soillayer_thickness[layer];
						if (sprop->CFlayer == (double)layer)
							CFeff = (sprop->GWD - sprop->CFD) / sitec->soillayer_thickness[layer];
						else
							CFeff = 1 - GWeff;

		

						/* higher FC */
						ratio = GWeff + CFeff;
						if ((ratio > 1 || ratio < 0) && !errorCode)
						{
							printf("ERROR in groundwater_preproc.c: invalid ratio in GWeff \n");
							errorCode = 1;
						}

						sprop->VWCfc[GWlayer] = sprop->VWCfc_base[GWlayer] * (1 - ratio) + sprop->VWCsat[GWlayer] * ratio;

						
						/* plus water or water loss because of  GW moving - on first simulation day - discharge, not movchange */
						VWCnew = GWeff * sprop->VWCsat[layer] + (1 - GWeff) * epv->VWC[layer];

						if (VWCnew > sprop->VWCfc[layer]) VWCnew = sprop->VWCfc[layer];

						if (VWCnew - epv->VWC[layer] > CRIT_PRECwater)
						{
							wf->GWmovchange[layer]  = (VWCnew - epv->VWC[layer]) * sitec->soillayer_thickness[layer] * water_density;
							ws->soilw[layer]       += (VWCnew - epv->VWC[layer]) * sitec->soillayer_thickness[layer] * water_density;
							epv->VWC[layer]         = VWCnew;
						} 

						/* if CF and GW are in the same layer: 3 layers in GWlayer - normal, capillary, saturated*/
						if (sprop->CFlayer == (double)layer)
						{
							sprop->dz_zoneCAPIL = sprop->GWD - sprop->CFD;	// capillary zone
							sprop->dz_zoneNORM = sprop->CFD - GWboundU;	// normal zone
						}
						/* if CF and GW are in the same layer: 2 layers in GWlayer - capillary, saturated*/
						else
						{
							sprop->dz_zoneCAPIL = sprop->GWD - GWboundU;
							sprop->dz_zoneNORM = 0;
						}

						sprop->dz_zoneSAT = GWboundL - sprop->GWD;        // saturated zone

						/* on the first day: initalization of VWC and concentration data of capillary and normalZone VWC-calculation and GWlayer */
						if (ctrl->firstsimday_flag == 1)
						{
							if (!errorCode && groundwater_firstday(md, sitec, GWS, sprop, soilInfo, ws, cs, ns))
							{
								printf("\n");
								printf("ERROR in groundwater_firstday.c for groundwater_preproc.c\n");
								errorCode = 1;
							}
						}
	
	
					}
					/* 4.2.2 ABOVE GW-table */
					else
					{
						/* 4.2.2.1 CF-layer */
						if (sprop->CFlayer <= (double) layer)
						{
							if (sprop->CFlayer == (double) layer)
								CFeff = (GWboundL - sprop->CFD)/sitec->soillayer_thickness[layer];
							else
								CFeff = 1;

							sprop->VWCfc[layer] = sprop->VWCfc_base[layer] * (1-CFeff) + sprop->VWCsat[layer] * CFeff;
							
						}
						/* 4.2.2.2 normal zone */
						else
						{
							sprop->VWCfc[layer] = sprop->VWCfc_base[layer];
						}
					}
				}
				wf->GWmovchange_total += wf->GWmovchange[layer];
			}
				
		}
		else
		{
			/* define GWlayer and actual CF from the previous day value */
			GWlayer=(int)sprop->GWlayer;
			CFlayer=(int)sprop->CFlayer;
		}

		/* ------------------------------------------------------------- */
		/* 5. determination of real capillary layers: potentially capillary zone is real CF if relSWC in layer below is higher than actual */
		/* special case: if capillary zone is in the bottom layer - always real capillary layer */

		if (CFlayer < N_SOILLAYERS-1)
		{ 
			for (layer = GWlayer - 1; layer >= CFlayer; layer--)
			{
				if (layer == sprop->CFlayer)
				{
					GWboundL = sitec->soillayer_depth[CFlayer];
					CFeff = (GWboundL - sprop->CFD) / sitec->soillayer_thickness[layer];
				}
				else
					CFeff = 1;


				if ((sprop->VWCsat[layer + 1] - epv->VWC[layer + 1]) > CRIT_PRECwater)
				{
					/* in only pontentially CF layers VWCfc is increasing gradually from VWCfc_base to saturation as the VWC of lower layer increases */
					if (epv->VWC[layer + 1] <= sprop->VWCfc_base[layer + 1])
						VWCfc_CF = sprop->VWCfc_base[layer];
					else
					{
						if (sprop->VWCsat[layer + 1] > sprop->VWCfc_base[layer + 1])
							ratio = (sprop->VWCsat[layer + 1] - epv->VWC[layer + 1]) / (sprop->VWCsat[layer + 1] - sprop->VWCfc_base[layer + 1]);
						else
							ratio = 1;

						VWCfc_CF = sprop->VWCfc_base[layer] + (1 - ratio) * (sprop->VWCsat[layer] - sprop->VWCfc_base[layer]);
					}
				}
				else
				{
					VWCfc_CF = sprop->VWCsat[layer];
				}


				sprop->VWCfc[layer] = sprop->VWCfc_base[layer] * (1 - CFeff) + VWCfc_CF * CFeff;
			} // end: for

		}
	
	}
	else
	{
		sprop->GWD = DATA_GAP;
		sprop->CFD = DATA_GAP;
	}
	

	sprop->preGWD     = sprop->GWD;
	sprop->preGWlayer = sprop->GWlayer;


	return (errorCode);
}


int groundwater_firstday(int md, const siteconst_struct* sitec, const groundwaterINIT_struct* GWS, soilprop_struct* sprop, soilInfo_struct* soilInfo,
	                     wstate_struct* ws, cstate_struct* cs, nstate_struct* ns)
{
	int errorCode = 0;
	double NH4ppmGW, NO3ppmGW, DOCppmGW;
	double soilw_zoneCAPIL, soilwGW, soilw_zoneNORM, soilC, GW_DOMconc, soilN;
	int GWlayer, layer, dm;
	double content_soilPRE[N_DISSOLVMATER][N_SOILLAYERS];
	double content_zoneNORM_PRE[N_DISSOLVMATER], content_zoneCAPIL_PRE[N_DISSOLVMATER];
	double conc_zoneNORM, conc_zoneCAPIL;

	GWlayer = (int)sprop->GWlayer;

	/*------------------------------------------------------------ */
	/* 1. calculation of water content of three layers of GWlayer (NORM, CF and SAT) and calculatuon in VWC of CF and NORM */

	soilwGW = sprop->VWCsat[GWlayer] * sprop->dz_zoneSAT * water_density;
	soilw_zoneCAPIL = (ws->soilw[GWlayer] - soilwGW) * (sprop->dz_zoneCAPIL / (sprop->dz_zoneCAPIL + sprop->dz_zoneNORM));
	soilw_zoneNORM = (ws->soilw[GWlayer] - soilwGW) * (sprop->dz_zoneNORM / (sprop->dz_zoneCAPIL + sprop->dz_zoneNORM));
	if (fabs((soilw_zoneCAPIL + soilwGW + soilw_zoneNORM) - ws->soilw[GWlayer]) > CRIT_PRECwater && !errorCode)
	{
		printf("\n");
		printf("ERROR in groundwater_preproc.c\n");
		errorCode = 1;
	}
	sprop->soilw_zoneNORM  = soilw_zoneNORM;
	sprop->soilw_zoneCAPIL = soilw_zoneCAPIL;
	sprop->soilw_zoneSAT   = soilwGW;

	sprop->VWC_zoneCAPIL = soilw_zoneCAPIL / (sprop->dz_zoneCAPIL * water_density);
	if (sprop->dz_zoneNORM)
		sprop->VWC_zoneNORM = soilw_zoneNORM / (sprop->dz_zoneNORM * water_density);
	else
		sprop->VWC_zoneNORM = 0;

	/*------------------------------------------------------------ */
   /* 2. Transfer of values  */

	/* call soil concentration calculation routine to calculate the concetration of soil */

	
	if (!errorCode && calc_soilconc(-1, 0, sprop, ws, cs, ns, soilInfo))
	{
		printf("ERROR in calc_soilconc.c for groundwater_preproc.c\n");
		errorCode = 1;
	}

	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		for (dm = 0; dm < N_DISSOLVMATER; dm++)
		{
			content_soilPRE[dm][layer] = soilInfo->content_soil[dm][layer];
		}
	}

	/*------------------------------------------------------------ */
	/* 3. calculation of concentration of groundwater  */

	if (GWS->GW_NH4ppm_array[md] == DATA_GAP)
	{
		NH4ppmGW = GW_NH4ppm;
		NO3ppmGW = GW_NO3ppm;
		DOCppmGW = GW_DOCppm;
	}
	else
	{
		NH4ppmGW = GWS->GW_NH4ppm_array[md];
		NO3ppmGW = GWS->GW_NO3ppm_array[md];
		DOCppmGW = GWS->GW_DOCppm_array[md];
	}



	soilInfo->conc_GW[0] = NH4ppmGW * 1e-6;
	soilInfo->conc_GW[1] = NO3ppmGW * 1e-6;
	GW_DOMconc = DOCppmGW * 1e-6;


	soilC = cs->soil1c_total + cs->soil2c_total + cs->soil3c_total + cs->soil4c_total;
	soilN = ns->soil1n_total + ns->soil2n_total + ns->soil3n_total + ns->soil4n_total;

	if (soilC && soilN)
	{
		soilInfo->conc_GW[2] = GW_DOMconc * (ns->soil1n_total / soilN);
		soilInfo->conc_GW[3] = GW_DOMconc * (ns->soil2n_total / soilN);
		soilInfo->conc_GW[4] = GW_DOMconc * (ns->soil3n_total / soilN);
		soilInfo->conc_GW[5] = GW_DOMconc * (ns->soil4n_total / soilN);
		soilInfo->conc_GW[6] = GW_DOMconc * (cs->soil1c_total / soilC);
		soilInfo->conc_GW[7] = GW_DOMconc * (cs->soil2c_total / soilC);
		soilInfo->conc_GW[8] = GW_DOMconc * (cs->soil3c_total / soilC);
		soilInfo->conc_GW[9] = GW_DOMconc * (cs->soil4c_total / soilC);
	}
	else
	{
		soilInfo->conc_GW[2] = GW_DOMconc * soil1ratio / sprop->soil1_CN;
		soilInfo->conc_GW[3] = GW_DOMconc * soil2ratio / sprop->soil2_CN;
		soilInfo->conc_GW[4] = GW_DOMconc * soil3ratio / sprop->soil3_CN;
		soilInfo->conc_GW[5] = GW_DOMconc * soil4ratio / sprop->soil4_CN;
		soilInfo->conc_GW[6] = GW_DOMconc * soil1ratio;
		soilInfo->conc_GW[7] = GW_DOMconc * soil2ratio;
		soilInfo->conc_GW[8] = GW_DOMconc * soil3ratio;
		soilInfo->conc_GW[9] = GW_DOMconc * soil4ratio;
	}
	
	
	/*------------------------------------------------------------ */
	/* 4. calculation of concentration and material content below GWlayer  */

	for (layer = GWlayer; layer < N_SOILLAYERS; layer++)
	{
		/*  GWlayer (3 zones)  */
		if (layer == GWlayer)
		{
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				content_zoneNORM_PRE[dm] = soilInfo->content_soil[dm][layer] * (sprop->dz_zoneNORM / sitec->soillayer_thickness[layer]);
				content_zoneCAPIL_PRE[dm] = soilInfo->content_soil[dm][layer] * (sprop->dz_zoneCAPIL / sitec->soillayer_thickness[layer]);
				content_soilPRE[dm][layer] = soilInfo->content_soil[dm][layer];

				conc_zoneNORM= soilInfo->conc_soil[dm][layer];
				conc_zoneCAPIL = soilInfo->conc_soil[dm][layer];

				soilInfo->content_zoneSAT[dm] = (soilInfo->conc_GW[dm] * sprop->soilw_zoneSAT) / soilInfo->dissolv_prop[dm];
				soilInfo->content_zoneCAPIL[dm] = (conc_zoneCAPIL * sprop->soilw_zoneCAPIL) / soilInfo->dissolv_prop[dm];
				soilInfo->content_zoneNORM[dm] = (conc_zoneNORM * sprop->soilw_zoneNORM) / soilInfo->dissolv_prop[dm];

				soilInfo->content_soil[dm][layer] = soilInfo->content_zoneNORM[dm] + soilInfo->content_zoneCAPIL[dm] + soilInfo->content_zoneSAT[dm];

				soilInfo->dismatGWmovchange_NORM[dm] = soilInfo->content_zoneNORM[dm] - content_zoneNORM_PRE[dm];
				soilInfo->dismatGWmovchange_CAPIL[dm] = soilInfo->content_zoneCAPIL[dm] - content_zoneCAPIL_PRE[dm];
				soilInfo->dismatGWmovchange[dm][layer] = soilInfo->content_soil[dm][layer] - content_soilPRE[dm][layer];
			}
		}
		else
		{ 
			/*  below GWlayer (all saturated)  */
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				content_soilPRE[dm][layer]           = soilInfo->content_soil[dm][layer];

				soilInfo->conc_soil[dm][layer]         = soilInfo->conc_GW[dm];
				soilInfo->content_soil[dm][layer]      = (soilInfo->conc_GW[dm] * ws->soilw[layer]) / soilInfo->dissolv_prop[dm];

				soilInfo->dismatGWmovchange[dm][layer] = soilInfo->content_soil[dm][layer] - content_soilPRE[dm][layer];
			}
		}

		/* src/snk variables*/
		for (dm = 0; dm < N_DISSOLVN; dm++)
		{
			if (soilInfo->dismatGWmovchange[dm][layer] > 0)
				ns->GWsrc_N += soilInfo->dismatGWmovchange[dm][layer];
			else
				ns->GWsnk_N += -1 * soilInfo->dismatGWmovchange[dm][layer];
		}
		for (dm = N_DISSOLVN; dm < N_DISSOLVMATER; dm++)
		{
			if (soilInfo->dismatGWmovchange[dm][layer] > 0)
				cs->GWsrc_C += soilInfo->dismatGWmovchange[dm][layer];
			else
				cs->GWsnk_C += -1 * soilInfo->dismatGWmovchange[dm][layer];
		}

		/* transfer value: content -> NH4, NO3, etc.*/
		if (!errorCode && calc_soilconc(layer, 1, sprop, ws, cs, ns, soilInfo))
		{
			printf("ERROR in calc_soilconc.c for groundwater_preproc.c\n");
			errorCode = 1;
		}
	}


	return (errorCode);
}

