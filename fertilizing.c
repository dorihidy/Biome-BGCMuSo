/* 
fertilizing.c
do fertilization  - increase the mineral soil nitrogen (sminn)

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

int fertilizing(const control_struct* ctrl, const siteconst_struct* sitec, const soilprop_struct* sprop, fertilizing_struct* FRZ,  
				cstate_struct* cs, nstate_struct* ns, wstate_struct* ws, cflux_struct* cf, nflux_struct* nf, wflux_struct* wf)
{


	/* fertilizing parameters .*/
	int errorCode=0;

	int layer;

	double FRZdepth;					    /* (m) actual  depth of fertilization */
	double fertilizer_DM;                   /* (kg/m2) dry matter content of fertilizer  */
	double fertilizer_WC;                   /* (kg/m2) water content of fertilizer  */
	double flab;                            /* (%) labile fraction of fertilizer */
	double fcel;                            /* (%) cellulose fraction of fertilizer */
	double fucel;                           /* (%) unshielded cellulose fraction of fertilizer */
	double fscel;                           /* (%) shielded cellulose fraction of fertilizer */
	double flig;                            /* (%) lignin fraction of fertilizer */
	double EFf_N2O;                         /* (kgN2O-N:kgN) fertilization emission factor for direct N2O emissions from synthetic fertililers */
	double FRZ_to_litrc, FRZ_to_litrn;		/* fertilizing carbon and nitrogen fluxes */
	int FRZlayer;                               /* (DIM) number of fertilization layer */
	
		
	double ratio, ratioSUM,ha_to_m2, diff;

	int md, year;

	year = ctrl->simstartyear + ctrl->simyr;
	md = FRZ->mgmdFRZ-1;

	FRZ_to_litrc=FRZ_to_litrn=ratio=ratioSUM=0;
	ha_to_m2 = 1./10000;

	
	/* On management days fertilizer is put on the field */
	if (FRZ->FRZ_num && md >= 0)
	{
		if (year == FRZ->FRZyear_array[md] && ctrl->month == FRZ->FRZmonth_array[md] && ctrl->day == FRZ->FRZday_array[md]) 
		{


			/* 2. input parameters in actual year from array */
			FRZdepth   = FRZ->FRZdepth_array[md]; 

			flab		= FRZ->litr_flab_array[md] / 100.;		/* from % to number */
			fcel		= FRZ->litr_fcel_array[md] / 100.;		/* from % to number */
			flig        = 1 - flab - fcel;

			/* calculate shielded and unshielded cellulose fraction */
			fscel = 0;
			fucel = 0;
			if (fcel)
			{
				ratio = flig/fcel;
				if (ratio <= 0.45)
				{
					fscel = 0.0;
					fucel = fcel;
				}
				else if (ratio > 0.45 && ratio < 0.7)
				{
					fscel = (ratio - 0.45)*3.2*fcel;
					fucel = fcel - fscel;
				}
				else
				{
					fscel = 0.8*fcel;
					fucel = 0.2*fcel;
				}
			}
		
	
			EFf_N2O = FRZ->EFfert_N2O[md];

			if ((flab + fucel + fscel + flig) > 0 && fabs(flab + fucel + fscel + flig - 1) > CRIT_PREC)
			{
				printf("ERROR in fertilizing parameters in management file: sum of labile/cellulose/ligning fraction parameters must equal to 1 (fertilizing.c)\n");
				errorCode=1;
			}

			

			/* DM and WC content of fertilizer: kg/m2 = kg fertilizer/ha * ha/m2 * (% to prop.) */
			fertilizer_DM = FRZ->fertilizer_array[md]  * ha_to_m2 * (FRZ->DM_array[md] / 100.);
			fertilizer_WC = FRZ->fertilizer_array[md]  * ha_to_m2 * (1 - FRZ->DM_array[md] / 100.);

			/* kgN/m2 = kg fertilizerDM/m2 * kgN/kg fertilizerDM */
			nf->FRZ_to_NH4 = fertilizer_DM * (FRZ->NH4content_array[md] / 100.);
			nf->FRZ_to_NO3 = fertilizer_DM * (FRZ->NO3content_array[md] / 100.);

			/* N2O emissions (kgN2O-N:kgN) */
			nf->N2OfluxFRZ_NH4 = nf->FRZ_to_NH4  * EFf_N2O;
			nf->N2OfluxFRZ_NO3 = nf->FRZ_to_NO3 * EFf_N2O;
			nf->N2OfluxFRZ = nf->N2OfluxFRZ_NH4 + nf->N2OfluxFRZ_NO3;



			/* due to N2O emission, FRZ_to_NH4 and FRZ_to_NO3 are decreasing */
			nf->FRZ_to_NH4 -= nf->FRZ_to_NH4 * EFf_N2O;
			nf->FRZ_to_NO3 -= nf->FRZ_to_NO3 * EFf_N2O;

			/* on fertilizing day a fixed amount of ammonium/nitrate/organic nitrogen/organic carbon/water enters into the soil */
			FRZ_to_litrn = fertilizer_DM * (FRZ->orgNcontent_array[md] / 100.);
			nf->FRZ_to_litr1n  = FRZ_to_litrn * flab;
			nf->FRZ_to_litr2n  = FRZ_to_litrn * fucel;
			nf->FRZ_to_litr3n  = FRZ_to_litrn * fscel;
			nf->FRZ_to_litr4n  = FRZ_to_litrn * flig;	
	
			FRZ_to_litrc = fertilizer_DM * (FRZ->orgCcontent_array[md] / 100.);
			cf->FRZ_to_litr1c  = FRZ_to_litrc * flab;
			cf->FRZ_to_litr2c  = FRZ_to_litrc * fucel;
			cf->FRZ_to_litr3c  = FRZ_to_litrc * fscel;
			cf->FRZ_to_litr4c  = FRZ_to_litrc * flig;	

			wf->FRZ_to_soilw   = fertilizer_WC;

		
			/* 4. fertilizing layer from depth */
			layer = 1;
			FRZlayer = 0;
			if (FRZdepth > sitec->soillayer_depth[0])
			{
				while (FRZlayer == 0 && layer < N_SOILLAYERS)
				{
					if ((FRZdepth > sitec->soillayer_depth[layer-1]) && (FRZdepth <= sitec->soillayer_depth[layer])) FRZlayer = layer;
					layer += 1;
				}
				if (FRZlayer == 0)
				{
					printf("ERROR in fertilizing depth calculation (fertilizing.c)\n");
					errorCode=1;
				}
			}


			/* 6. STATE UPDATE */
			ratio=0;
			for (layer = 0; layer <= FRZlayer; layer++)
			{
				if (FRZlayer == 0)
					ratio = 1;
				else	
				{
					if (sitec->soillayer_depth[layer] < FRZdepth)
						ratio = sitec->soillayer_thickness[layer] / FRZdepth;
					else
						ratio = (FRZdepth - sitec->soillayer_depth[layer-1]) / FRZdepth;
				}

				ratioSUM += ratio;

				cs->litr1c[layer] += cf->FRZ_to_litr1c * ratio;
				cs->litr2c[layer] += cf->FRZ_to_litr2c * ratio;
				cs->litr3c[layer] += cf->FRZ_to_litr3c * ratio;
				cs->litr4c[layer] += cf->FRZ_to_litr4c * ratio;

				ns->litr1n[layer] += nf->FRZ_to_litr1n * ratio;
				ns->litr2n[layer] += nf->FRZ_to_litr2n * ratio;
				ns->litr3n[layer] += nf->FRZ_to_litr3n * ratio;
				ns->litr4n[layer] += nf->FRZ_to_litr4n * ratio;

				ns->NH4[layer]  += nf->FRZ_to_NH4 * ratio;
				ns->NO3[layer]  += nf->FRZ_to_NO3 * ratio;

				/* water from fertilization -> soil layers, in case of oversaturation: pondw */	
				ws->soilw[layer] += wf->FRZ_to_soilw * ratio;

				diff = ws->soilw[layer] - sprop->VWCsat[layer] * sitec->soillayer_thickness[layer] * water_density;
				if (diff > CRIT_PRECwater) 
				{
					ws->pondw += diff;
					ws->soilw[layer] = sprop->VWCsat[layer] * sitec->soillayer_thickness[layer] * water_density;
				}

			}

			/* control */
			if (fabs(ratioSUM - 1) > CRIT_PREC)
			{
				printf("ERROR in fertilizing ratio calculation (fertilizing.c)\n");
				errorCode=1;
			}

			cs->FRZsrc_C += FRZ_to_litrc;
			ns->FRZsrc_N += nf->FRZ_to_NH4 + nf->FRZ_to_NO3 + FRZ_to_litrn + nf->N2OfluxFRZ;
			ws->FRZsrc_W += wf->FRZ_to_soilw;
			ns->Nvol_snk += nf->N2OfluxFRZ;
		

		} /* endif  */



	}
   return (errorCode);
}
	