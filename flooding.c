/* 
flooding.c
Calculate the effect of flooding (water, NH4 and NO3 from river)

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
#include "bgc_func.h"

int flooding(control_struct* ctrl, const siteconst_struct* sitec, const flooding_struct* FLS, soilprop_struct* sprop, epvar_struct* epv, 
	         wstate_struct* ws, wflux_struct* wf, cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf, soilInfo_struct* soilInfo, int* mondays)
{

	int errorCode=0;
	int md, year, layer, dm;
	double FL_to_soilw;
	double material_fromFL[N_DISSOLVMATER];
	int FLyday_start, FLyday_end;

	FL_to_soilw = 0;
	year = ctrl->simstartyear + ctrl->simyr;
	md = FLS->mgmdFL-1;

	/* 2. flooding affects soil water content (fills to saturation) */ 
	/* the flowing water also brings nutrients - in the first round, we assume that the concentration in the upper layer does not change */
	if (FLS->FLnum && md >= 0)
	{
		FLyday_start = date_to_doy(mondays, FLS->FLstart_month_array[md], FLS->FLstart_day_array[md]);
		FLyday_end   = date_to_doy(mondays, FLS->FLend_month_array[md], FLS->FLend_day_array[md]);

	
		if (year == FLS->FLstart_year_array[md] && ctrl->yday >= FLyday_start && ctrl->yday <= FLyday_end) 
		{	
			sprop->FLD = FLS->FLheight_array[md];

			/* concentration of GW - from file or from contant value of MuSo */
			if (!errorCode && flooding_concentration(md, FLS, soilInfo, ws, cs, ns))
			{
				printf("\n");
				printf("ERROR in flooding_concentration.c for flooding.c\n");
				errorCode = 1;
			}
	

			/* we assume, flooding wet the 0-10 cm soil layers */
			for (layer = 0; layer < 2; layer++)
			{
				if (epv->VWC[layer] < sprop->VWCsat[layer])
				{
					FL_to_soilw = (sprop->VWCsat[layer] - epv->VWC[layer]) * sitec->soillayer_thickness[layer] * water_density;


					/* calculation nutrients entering with flooding water */
					for (dm = 0; dm < N_DISSOLVMATER; dm++)
					{
						material_fromFL[dm]              = soilInfo->FLconc[dm] * FL_to_soilw;
						soilInfo->content_soil[dm][layer] += material_fromFL[dm];
					}
				

					ws->soilw[layer] += FL_to_soilw;
					epv->VWC[layer]   = ws->soilw[layer] / (sitec->soillayer_thickness[layer] * water_density);

					if (sprop->VWCsat[layer] - epv->VWC[layer] > CRIT_PREC_lenient)
					{
						printf("\n");
						printf("ERROR in water saturation calculation in flooding.c\n");
						errorCode = 1;
					}

					/* call soil concentration calculation routine to calculate the concetration of FL */
					if (!errorCode && check_soilcontent(layer, 1, sprop, cs, ns, soilInfo))
					{
						printf("ERROR in check_soilcontent.c for flooding.c\n");
						errorCode = 1;
					}
								
					
					wf->FL_to_soilw  += FL_to_soilw;
					for (dm = 0; dm < N_DISSOLVinorgN; dm++)                       nf->sminN_fromFL += material_fromFL[dm];
					for (dm = N_DISSOLVinorgN; dm < N_DISSOLVN; dm++)              nf->orgN_fromFL  += material_fromFL[dm];
					for (dm = N_DISSOLVN; dm < N_DISSOLVMATER; dm++)               cf->orgC_fromFL  += material_fromFL[dm];
					for (dm = 0; dm < N_DISSOLVN; dm++)                            ns->FLsrc_N      += material_fromFL[dm];
					for (dm = N_DISSOLVN; dm < N_DISSOLVMATER; dm++)               cs->FLsrc_C      += material_fromFL[dm];
				}

			}

			/* unit of FLheight is mm - 1 mm water on 1m2 surface is 1 kg */
			if (sprop->FLD > ws->pondw)
			{
				wf->FL_to_pondw = sprop->FLD - ws->pondw;
				ws->pondw       += wf->FL_to_pondw;
			}

			wf->FL_to_soilwTOTAL = wf->FL_to_pondw + wf->FL_to_soilw;

			if (ws->pondw > sprop->FLD)
			{
				wf->pondw_to_runoff += ws->pondw - sprop->FLD;
				ws->pondw             = sprop->FLD;
			}
		}
		else
		{
			sprop->FLD = DATA_GAP;
		}
		
		
		
	}
	else
	{
		for (dm = 0; dm < N_DISSOLVMATER; dm++) soilInfo->FLconc[dm] = 0;

	}
	return (errorCode);
}

int flooding_concentration(int md, const flooding_struct* FLS, soilInfo_struct* soilInfo, wstate_struct* ws, cstate_struct* cs, nstate_struct* ns)
{
	int dm, FLlayer;
	int errorCode = 0;
	double FLconc_fromFILE[N_DISSOLVMATER];
	double soilconc_array[N_DISSOLVMATER];

	FLlayer = 0;

	/* calculation of soil concentration */

	soilconc_array[0] = (ns->NH4[FLlayer] * soilInfo->dissolv_prop[0]) / ws->soilw[FLlayer];
	soilconc_array[1] = (ns->NO3[FLlayer] * soilInfo->dissolv_prop[1]) / ws->soilw[FLlayer];
	soilconc_array[2] = (ns->soil1n[FLlayer] * soilInfo->dissolv_prop[2]) / ws->soilw[FLlayer];
	soilconc_array[3] = (ns->soil2n[FLlayer] * soilInfo->dissolv_prop[3]) / ws->soilw[FLlayer];
	soilconc_array[4] = (ns->soil3n[FLlayer] * soilInfo->dissolv_prop[4]) / ws->soilw[FLlayer];
	soilconc_array[5] = (ns->soil4n[FLlayer] * soilInfo->dissolv_prop[5]) / ws->soilw[FLlayer];
	soilconc_array[6] = (cs->soil1c[FLlayer] * soilInfo->dissolv_prop[6]) / ws->soilw[FLlayer];
	soilconc_array[7] = (cs->soil2c[FLlayer] * soilInfo->dissolv_prop[7]) / ws->soilw[FLlayer];
	soilconc_array[8] = (cs->soil3c[FLlayer] * soilInfo->dissolv_prop[8]) / ws->soilw[FLlayer];
	soilconc_array[9] = (cs->soil4c[FLlayer] * soilInfo->dissolv_prop[9]) / ws->soilw[FLlayer];

	FLconc_fromFILE[0] = FLS->FL_NH4ppm_array[md] * 1e-6;
	FLconc_fromFILE[1] = FLS->FL_NO3ppm_array[md] * 1e-6;
	FLconc_fromFILE[2] = FLS->FL_DON1ppm_array[md] * 1e-6;
	FLconc_fromFILE[3] = FLS->FL_DON2ppm_array[md] * 1e-6;
	FLconc_fromFILE[4] = FLS->FL_DON3ppm_array[md] * 1e-6;
	FLconc_fromFILE[5] = FLS->FL_DON4ppm_array[md] * 1e-6;
	FLconc_fromFILE[6] = FLS->FL_DOC1ppm_array[md] * 1e-6;
	FLconc_fromFILE[7] = FLS->FL_DOC2ppm_array[md] * 1e-6;
	FLconc_fromFILE[8] = FLS->FL_DOC3ppm_array[md] * 1e-6;
	FLconc_fromFILE[9] = FLS->FL_DOC4ppm_array[md] * 1e-6;

	/* if input data is avaialbe ins FLfile - concentration from file; else: concentration from soil (no effect of FL) */
	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		if (FLconc_fromFILE[dm] >= 0)
			soilInfo->FLconc[dm] = FLconc_fromFILE[dm];
		else
			soilInfo->FLconc[dm] = soilconc_array[dm];
	}



	return (errorCode);
}