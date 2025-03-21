/* 
groundwater_preproc.c
calculate GWD, CFD, GWeff, CFeff in function of GWD

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


int groundwater_preproc(control_struct* ctrl, const groundwaterINIT_struct* GWS, const siteconst_struct* sitec, soilprop_struct* sprop, soilInfo_struct* soilInfo,
	                    wstate_struct* ws, wflux_struct* wf, cstate_struct* cs, nstate_struct* ns)
{
	int layer,errorCode,md, year, GWlayer, CFlayer, dm;
	double CFact;
	double GWboundL, GWboundU, ratio;


	/*------------------------------------*/
	/* Initialization */
	
	errorCode=0;
	GWlayer=CFlayer = DATA_GAP;

	md = GWS->mgmdGW-1;
	year = ctrl->simstartyear + ctrl->simyr;


	/*------------------------------------*/
	/* Groundwater calculation only if GWD data is available from GW-file */

	if (GWS->GWnum != 0)	
	{	
		/*-----------------------------------------------------------------------------*/
		/* 1. Determination of groundwater days */
	
		if (md >= 0 && year == GWS->GWyear_array[md] && ctrl->month == GWS->GWmonth_array[md] && ctrl->day == GWS->GWday_array[md])
		{
			sprop->GWD = GWS->GWdepth_array[md];


			/* GW above the surface: pond water formation (limitation: maximal pond water height */
			if (sprop->GWD < 0)
			{	
				ws->GW_waterlogging = -1*sprop->GWD* m_to_mm;
				if (ws->pondw < ws->GW_waterlogging)
					wf->GW_to_pondw = -1 * sprop->GWD * m_to_mm - ws->pondw;
				else
					wf->GW_to_pondw = 0;

				sprop->GWD          = 0;

			}

			/* very deep GW-level */
			if (sprop->GWD >= sitec->soillayer_depth[N_SOILLAYERS - 1])
			{	
				
				sprop->CFD = sprop->GWD - sprop->CapillFringe[N_SOILLAYERS - 1];
				sprop->CFlayer = DATA_GAP;
				sprop->GWlayer = DATA_GAP;
				
			}


	

/*			else
			{


				sprop->soilw_NORMcf = ws->soilw[(int)sprop->CFlayer] * sprop->dz_NORMcf / sitec->soillayer_thickness[(int)sprop->CFlayer];
				sprop->soilw_CAPILcf = ws->soilw[(int)sprop->CFlayer] * sprop->dz_CAPILcf / sitec->soillayer_thickness[(int)sprop->CFlayer];
				sprop->soilw_NORMgw = ws->soilw[(int)sprop->GWlayer] * sprop->dz_NORMgw / sitec->soillayer_thickness[(int)sprop->GWlayer];
				sprop->soilw_CAPILgw = ws->soilw[(int)sprop->GWlayer] * sprop->dz_CAPILgw / sitec->soillayer_thickness[(int)sprop->GWlayer];
				sprop->soilw_SATgw = ws->soilw[(int)sprop->GWlayer] * sprop->dz_SATgw / sitec->soillayer_thickness[(int)sprop->GWlayer];

				sprop->VWC_NORMcf = sprop->soilw_NORMcf / (water_density * sprop->dz_NORMcf);
				sprop->VWC_CAPILcf = sprop->soilw_CAPILcf / (water_density * sprop->dz_CAPILcf);
				sprop->VWC_NORMgw = sprop->soilw_NORMgw / (water_density * sprop->dz_NORMgw);
				sprop->VWC_CAPILgw = sprop->soilw_CAPILgw / (water_density * sprop->dz_CAPILgw);
			}
			}*/
		}
		else
		{
			sprop->GWD = DATA_GAP;
			if (ctrl->firstsimday_flag != 1)
			{
				printf("\n");
				printf("ERROR in groundwater_preproc: groundwater data series must be continouos\n");
				errorCode = 1;
			}
			sprop->CFlayer = DATA_GAP;
			sprop->GWD_pre = DATA_GAP;
		}


		/*-----------------------------------------------------------------------------*/
		/* 2. Calculation of actual value of the depth of groundwater and capillary zone*/

		if (sprop->GWD != DATA_GAP && sprop->GWD < sitec->soillayer_depth[N_SOILLAYERS - 1] && sprop->GWD != sprop->GWD_pre)
		{
			/* initialization */
			for (layer = 0; layer < N_SOILLAYERS; layer++)
			{
				sprop->GWeff[layer] = -1;
				sprop->CFeff[layer] = -1;
			}

			/* pre values */	
			sprop->soilw_NORMcf_pre = sprop->soilw_NORMcf;
			sprop->soilw_CAPILcf_pre = sprop->soilw_CAPILcf;
			sprop->soilw_NORMgw_pre = sprop->soilw_NORMgw;
			sprop->soilw_CAPILgw_pre = sprop->soilw_CAPILgw;
			sprop->soilw_SATgw_pre = sprop->soilw_SATgw;

			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				soilInfo->content_NORMcf_pre[dm] = soilInfo->content_NORMcf[dm];
				soilInfo->content_CAPILcf_pre[dm] = soilInfo->content_CAPILcf[dm];
				soilInfo->content_NORMgw_pre[dm] = soilInfo->content_NORMgw[dm];
				soilInfo->content_CAPILgw_pre[dm] = soilInfo->content_CAPILgw[dm];
				soilInfo->content_SATgw_pre[dm] = soilInfo->content_SATgw[dm];
			}

		
			/* 2.1: GWlayer and CFlayer */
			if (sprop->GWD == 0)
			{
				GWboundU = 0;
				GWlayer = 0;
				CFlayer = 0;
				sprop->CFD = 0;
			}
			else
			{
				/*  lower and upper boundary of the laye */
				layer = N_SOILLAYERS - 1;
				while (GWlayer == DATA_GAP && layer >= 0)
				{
					GWboundL = sitec->soillayer_depth[layer];

					/* if groundwater table is in actual layer (above lower boundary): GWlayer   */
					if (layer > 0)
						GWboundU = sitec->soillayer_depth[layer - 1];
					else
						GWboundU = 0;

					if (sprop->GWD < GWboundL && sprop->GWD >= GWboundU)
					{
						/* if GW-table is under 10 meter */
						if (sprop->GWD > sitec->soillayer_depth[N_SOILLAYERS - 1])
						{
							GWlayer = N_SOILLAYERS;
							CFact = sprop->CapillFringe[N_SOILLAYERS - 1];
						}
						else
						{
							GWlayer = layer;
							CFact = sprop->CapillFringe[layer];
						}


						/* searching of upper boundary of capillary zone */
						if (sprop->GWD - CFact < GWboundU)
						{
							layer -= 1;
							CFact = sprop->CapillFringe[layer] - (sprop->GWD - GWboundU);
							while (CFlayer == DATA_GAP && layer >= 0)
							{
								GWboundL = sitec->soillayer_depth[layer];
								if (layer == 0)
									GWboundU = -1 * CRIT_PREC;
								else
									GWboundU = sitec->soillayer_depth[layer - 1];

								if (GWboundL - CFact <= GWboundU)
								{
									layer -= 1;
									CFact = sprop->CapillFringe[layer] - (sprop->GWD - GWboundU);
									if (CFact < 0)
									{
										CFlayer = layer + 1;
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
							if (sprop->CFD >= sitec->soillayer_depth[N_SOILLAYERS - 1]) CFlayer = N_SOILLAYERS;

						}

					}
					else
						layer -= 1;
				}

			}
		
			
			/* 2.2 define GWlayer and actual CF (maximum of capillary fringe upper layer is the top soil zone */
			if (CFlayer == DATA_GAP) 
			{
				CFlayer = 0;
				sprop->CFD = 0;
			}
			sprop->CFlayer=(double) CFlayer;
			sprop->GWlayer=(double) GWlayer;

		
			/* 2.3. soil layers below the groundwater table are saturated - net water gain from soil system */
			for (layer = N_SOILLAYERS-1; layer >= 0; layer--)
			{
				GWboundL = sitec->soillayer_depth[layer];
			
				if (layer == 0)
					GWboundU  = 0;
				else
					GWboundU  = sitec->soillayer_depth[layer-1];

				/* 4.1 BELOW the GW-table: saturation */
				if (layer > sprop->GWlayer)
				{
					sprop->CFeff[layer]     = 0;
					sprop->GWeff[layer]     = 1;

				}
				/* 2.3.1. IN and ABOVE the GW-table: saturation */
				else
				{
					/* IN the GW-layer: saturation proportionally to the GW depth, higher FC-values, plus water from GW */
					if (sprop->GWlayer == (double)layer)
					{
						sprop->dz_SATgw = GWboundL - sprop->GWD;        // saturated zone

						/* if CF and GW are in the same layer: 3 layers in GWlayer - normal, capillary, saturated*/
						if (sprop->CFlayer == (double)layer)
						{

							sprop->dz_CAPILgw = sprop->GWD - sprop->CFD; // capillary zone in GWlayer
							sprop->dz_NORMgw = sprop->CFD - GWboundU;	 // normal zone in GWlayer		
							sprop->dz_CAPILcf = 0;
							sprop->dz_NORMcf = 0;
						}
						/* if CF and GW are in the same layer: 2 layers in GWlayer - capillary, saturated*/
						else
						{
							sprop->dz_CAPILgw = sprop->GWD - GWboundU; // capillary zone in GWlayer
							sprop->dz_NORMgw = 0;
						}

					
						/* saturation proportionally to the GW depth */
						sprop->GWeff[layer] = sprop->dz_SATgw / sitec->soillayer_thickness[layer];
						sprop->CFeff[layer] = sprop->dz_CAPILgw / sitec->soillayer_thickness[layer];

						/* higher FC */
						ratio = sprop->GWeff[layer] + sprop->CFeff[layer];
						if ((ratio > 1 || ratio < 0) && !errorCode)
						{
							printf("ERROR in groundwater_preproc.c: invalid ratio in GWeff \n");
							errorCode = 1;
						}
		
			
					}
					/* 2.3.2 ABOVE GW-table */
					else
					{
						/* CF-layer */
						if (sprop->CFlayer <= (double) layer)
						{
							
							/* if CF and GW are in the same layer: 3 layers in GWlayer - normal, capillary, saturated*/
							if (sprop->CFlayer == (double)layer)
							{
								sprop->CFeff[layer] = (GWboundL - sprop->CFD) / sitec->soillayer_thickness[layer];
								sprop->dz_CAPILcf = GWboundL - sprop->CFD;	// capillary zone in unsat zone
								sprop->dz_NORMcf = sitec->soillayer_thickness[layer] - sprop->dz_CAPILcf;	// normal zone in unsat zone
							}
							/* if CF and GW are not in the same layer: 2 layers in GWlayer - capillary, saturated*/
							else
								sprop->CFeff[layer] = 1;

							
						}
						/* normal zone */
						else
							sprop->CFeff[layer] = -1;
					}
				}

			}

	
				
		}
		/* no change in GWD */
		else
		{
			/* define GWlayer and actual CF from the previous day value */
			GWlayer=(int)sprop->GWlayer;
			CFlayer=(int)sprop->CFlayer;


		}

		/* concentration of GW - from file or from contant value of MuSo */
		if (!errorCode && GWlayer != DATA_GAP && groundwater_concentration(md, GWS, sprop, soilInfo, ws, cs, ns))
		{
			printf("\n");
			printf("ERROR in groundwater_concentration.c for groundwater_preproc.c\n");
			errorCode = 1;
		}

		/* calculation of VWC of mixed layer (NORM and CAPIL) */
		if (sprop->dz_NORMgw + sprop->dz_CAPILgw != 0)
		{
			sprop->ratioNORMgw = sprop->dz_NORMgw / (sprop->dz_NORMgw + sprop->dz_CAPILgw);
			sprop->ratioCAPILgw = sprop->dz_CAPILgw / (sprop->dz_NORMgw + sprop->dz_CAPILgw);

			if (fabs(1 - sprop->ratioNORMgw - sprop->ratioCAPILgw) > CRIT_PREC)
			{
				printf("\n");
				printf("ERROR in ratio calculation in groundwater_preproc.c\n");
				errorCode = 1;
			}
		}
		else
		{
			sprop->ratioNORMgw = 0;
			sprop->ratioCAPILgw = 0;
		}

		if (sprop->dz_NORMcf + sprop->dz_CAPILcf)
		{ 
			sprop->ratioNORMcf = sprop->dz_NORMcf / (sprop->dz_NORMcf + sprop->dz_CAPILcf);
			sprop->ratioCAPILcf = sprop->dz_CAPILcf / (sprop->dz_NORMcf + sprop->dz_CAPILcf);
			if (fabs(1 - sprop->ratioNORMcf - sprop->ratioCAPILcf) > CRIT_PREC)
			{
				printf("\n");
				printf("ERROR in ratio calculation in groundwater_preproc.c\n");
				errorCode = 1;
			}
		}
		else
		{
			sprop->ratioNORMcf = 0;
			sprop->ratioCAPILcf = 0;
		}

	}
	else
	{
		if (ctrl->firstsimday_flag == 1)
		{
			sprop->GWD = DATA_GAP;
			sprop->CFD = DATA_GAP;
			sprop->GWlayer = DATA_GAP;
			sprop->CFlayer = DATA_GAP;
			sprop->dz_CAPILcf = 0;
			sprop->dz_NORMcf = 0;
			sprop->dz_CAPILgw = 0;
			sprop->dz_NORMgw = 0;
			sprop->dz_SATgw = 0;
			sprop->GWD_pre = DATA_GAP;
			sprop->realCAPILlayer = DATA_GAP;
			for (layer = 0; layer < N_SOILLAYERS; layer++)
			{
				sprop->GWeff[layer] = -1;
				sprop->CFeff[layer] = -1;
			}
			for (dm = 0; layer < N_DISSOLVMATER; dm++) soilInfo->GWconc[dm] = 0;

		}
	}
	

	

	return (errorCode);
}

int groundwater_concentration(int md, const groundwaterINIT_struct* GWS, soilprop_struct* sprop, soilInfo_struct* soilInfo,
	                          wstate_struct* ws, cstate_struct* cs, nstate_struct* ns)
{
	int dm, GWlayer;
	int errorCode = 0;
	double GWconc_fromFILE[N_DISSOLVMATER];
	double soilconc_array[N_DISSOLVMATER];      

	GWlayer = (int)sprop->GWlayer;

	/* calculation of soil concentration */
	
	soilconc_array[0] = (ns->NH4[GWlayer] * soilInfo->dissolv_prop[0]) / ws->soilw[GWlayer];
	soilconc_array[1] = (ns->NO3[GWlayer] * soilInfo->dissolv_prop[1]) / ws->soilw[GWlayer];
	soilconc_array[2] = (ns->soil1n[GWlayer] * soilInfo->dissolv_prop[2]) / ws->soilw[GWlayer];
	soilconc_array[3] = (ns->soil2n[GWlayer] * soilInfo->dissolv_prop[3]) / ws->soilw[GWlayer];
	soilconc_array[4] = (ns->soil3n[GWlayer] * soilInfo->dissolv_prop[4]) / ws->soilw[GWlayer];
	soilconc_array[5] = (ns->soil4n[GWlayer] * soilInfo->dissolv_prop[5]) / ws->soilw[GWlayer];
	soilconc_array[6] = (cs->soil1c[GWlayer] * soilInfo->dissolv_prop[6]) / ws->soilw[GWlayer];
	soilconc_array[7] = (cs->soil2c[GWlayer] * soilInfo->dissolv_prop[7]) / ws->soilw[GWlayer];
	soilconc_array[8] = (cs->soil3c[GWlayer] * soilInfo->dissolv_prop[8]) / ws->soilw[GWlayer];
	soilconc_array[9] = (cs->soil4c[GWlayer] * soilInfo->dissolv_prop[9]) / ws->soilw[GWlayer];

	GWconc_fromFILE[0] = GWS->GW_NH4ppm_array[md] * 1e-6;
	GWconc_fromFILE[1] = GWS->GW_NO3ppm_array[md] * 1e-6;
	GWconc_fromFILE[2] = GWS->GW_DON1ppm_array[md] * 1e-6;
	GWconc_fromFILE[3] = GWS->GW_DON2ppm_array[md] * 1e-6;
	GWconc_fromFILE[4] = GWS->GW_DON3ppm_array[md] * 1e-6;
	GWconc_fromFILE[5] = GWS->GW_DON4ppm_array[md] * 1e-6;
	GWconc_fromFILE[6] = GWS->GW_DOC1ppm_array[md] * 1e-6;
	GWconc_fromFILE[7] = GWS->GW_DOC2ppm_array[md] * 1e-6;
	GWconc_fromFILE[8] = GWS->GW_DOC3ppm_array[md] * 1e-6;
	GWconc_fromFILE[9] = GWS->GW_DOC4ppm_array[md] * 1e-6;

	/* if input data is avaialbe ins GWfile - concentration from file; else: concentration from soil (no effect of GW) */
	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		if (GWconc_fromFILE[dm] >= 0)
			soilInfo->GWconc[dm] = GWconc_fromFILE[dm];
		else
			soilInfo->GWconc[dm] = soilconc_array[dm];
	}
		
	

	return (errorCode);
}