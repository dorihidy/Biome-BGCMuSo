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
	double conc1[N_DISSOLVMATER], material_fromFL[N_DISSOLVMATER];
	double soilN, soilC, FL_DOMconc;
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

			if (FLS->FL_NH4ppm_array[md] != DATA_GAP)
			{ 
				ctrl->FLconcFROMfile_flag = 1;
				

				/* flooding data */
				if (FLS->FL_NH4ppm_array[md] != -1)
				{
					/* NH4, NO3 and DOM */
					soilInfo->conc_FL[0] = FLS->FL_NH4ppm_array[md] * 1e-6;
					soilInfo->conc_FL[1] = FLS->FL_NO3ppm_array[md] * 1e-6;
					FL_DOMconc = FLS->FL_DOCppm_array[md] * 1e-6;


					soilC = cs->soil1c_total + cs->soil2c_total + cs->soil3c_total + cs->soil4c_total;
					soilN = ns->soil1n_total + ns->soil2n_total + ns->soil3n_total + ns->soil4n_total;

					if (soilC && soilN)
					{
						soilInfo->conc_FL[2] = FL_DOMconc * (ns->soil1n_total / soilN);
						soilInfo->conc_FL[3] = FL_DOMconc * (ns->soil2n_total / soilN);
						soilInfo->conc_FL[4] = FL_DOMconc * (ns->soil3n_total / soilN);
						soilInfo->conc_FL[5] = FL_DOMconc * (ns->soil4n_total / soilN);
						soilInfo->conc_FL[6] = FL_DOMconc * (cs->soil1c_total / soilC);
						soilInfo->conc_FL[7] = FL_DOMconc * (cs->soil2c_total / soilC);
						soilInfo->conc_FL[8] = FL_DOMconc * (cs->soil3c_total / soilC);
						soilInfo->conc_FL[9] = FL_DOMconc * (cs->soil4c_total / soilC);
					}
					else
					{
						for (dm = N_DISSOLVinorgN; dm < N_DISSOLVMATER; dm++) soilInfo->conc_FL[dm] = 0;
					}
				}
			}
			else
			{
				ctrl->FLconcFROMfile_flag = 0;
				
				/* call soil concentration calculation routine to calculate the concetration of soil (content2pool_flag=0: NH4, NO3 -> content). */
				if (!errorCode && calc_soilconc(-1, 0, sprop, ws, cs, ns, soilInfo))
				{
					printf("ERROR in calc_soilconc.c for flooding.c\n");
					errorCode = 1;
				}
			}
				

			/* we assume, flooding wet the 0-10 cm soil layers */
			for (layer = 0; layer < 2; layer++)
			{
				if (epv->VWC[layer] < sprop->VWCsat[layer])
				{
					FL_to_soilw = (sprop->VWCsat[layer] - epv->VWC[layer]) * sitec->soillayer_thickness[layer] * water_density;


					/* calculation ratio for save concentrations in the upper layer */
					if (ctrl->FLconcFROMfile_flag)
					{
						for (dm = 0; dm < N_DISSOLVMATER; dm++) conc1[dm] = soilInfo->conc_FL[dm];
					}
					else
					{ 
						for (dm = 0; dm < N_DISSOLVMATER; dm++) conc1[dm] = soilInfo->conc_soil[dm][layer];
					}

					/* calculation nutrients entering with flooding water */
					for (dm = 0; dm < N_DISSOLVMATER; dm++)
					{
						material_fromFL[dm]              = conc1[dm] * FL_to_soilw;
						soilInfo->content_soil[dm][layer] += material_fromFL[dm];
					}
				

					ws->soilw[layer] += FL_to_soilw;
					epv->VWC[layer]   = ws->soilw[layer] / (sitec->soillayer_thickness[layer] * water_density);

					if (sprop->VWCsat[layer] - epv->VWC[layer] > CRIT_PRECwater)
					{
						printf("\n");
						printf("ERROR in water saturation calculation in flooding.c\n");
						errorCode = 1;
					}

					/* call soil concentration calculation routine to calculate the concetration of FL */
					if (!errorCode && calc_soilconc(layer, 1, sprop, ws, cs, ns, soilInfo))
					{
						printf("ERROR in calc_soilconc.c for flooding.c\n");
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
		for (dm = 0; dm < N_DISSOLVMATER; dm++) soilInfo->conc_FL[dm] = 0;

	}
	return (errorCode);
}