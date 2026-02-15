/* 
multilayer_hydrolparams.c
calcultion of soil water potential, hydr. conductivity and hydr. diffusivity as a function of volumetric water content and constants related to texture
calculation of relative VWC data and critical VWC data for rootzone 

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

int multilayer_hydrolparams(siteconst_struct* sitec, soilprop_struct* sprop, wstate_struct* ws, epvar_struct* epv)
{
	/* given a list of site constants and the soil water mass (kg/m2),
	this function returns the soil water potential (MPa)
	inputs:
	ws.soilw				     (kg/m2) water mass per unit area
	sprop.rootzoneDepth_max     (m)     maximum depth of rooting zone               
	sprop.soilB				 (dimless)   slope of log(PSI) vs log(rwc)
	sprop.VWCsat				 (m3/m3)   volumetric water content at saturation
	sprop.PSIsat			   	(MPa)   soil matric potential at saturation

	uses the relation:
	PSI = PSIsat * (VWC/VWCsat)^(-b)
	hydrCONDUCT = hydrCONDUCTsat * (VWC/VWCsat)^(2b+3)
	hydrDIFFUS = b* hydrDIFFUS * (-1*PSI) * (VWC/VWCsat)^(b+2)

	For further discussion see:
	Cosby, B.J., G.M. Hornberger, R.B. Clapp, and T.R. Ginn, 1984.  A     
	   statistical exploration of the relationships of soil moisture      
	   characteristics to the physical properties of soils.  Water Res.
	   Res. 20:682-690.
	
	Saxton, K.E., W.J. Rawls, J.S. Romberger, and R.I. Papendick, 1986.
		Estimating generalized soil-water characteristics from texture.
		Soil Sci. Soc. Am. J. 50:1031-1036.

	Balsamo et al 2009 - 
	A Revised Hydrology for the ECMWF Model - Verification from Field Site to Water Storage IFS - JHydromet.pdf

	Chen and Dudhia 2001 - 
	Coupling an Advanced Land Surface-Hydrology Model with the PMM5 Modeling System Part I - MonWRev.pdf


	*/

	int errorCode=0;
	int layer;

	double VWCsat_RZ, VWCfc_RZ, VWCwp_RZ, VWChw_RZ, soilwAVAIL_RZ;
	double VWC_avg, VWC_RZmax, relVWCsat_fc_RZmax, relVWCfc_wp_RZmax, VWC_RZ, PSI_RZ, soilw_RZ, weight, weight_SUM, ratio;

	VWC_RZ = VWC_RZmax = VWCsat_RZ = VWCfc_RZ = VWCwp_RZ = VWChw_RZ = relVWCsat_fc_RZmax = relVWCfc_wp_RZmax = 0.0;
	VWC_avg = VWC_RZ = PSI_RZ = soilw_RZ = weight = weight_SUM = ratio = soilwAVAIL_RZ = 0;



	/* ***************************************************************************************************** */
	/* calculating VWC, PSI and hydr. cond. to every layer */


	for (layer=0; layer < N_SOILLAYERS; layer++)
	{
		
		/* convert kg/m2 --> m3/m2 --> m3/m3 */
		epv->VWC[layer]  = ws->soilw[layer] / (water_density * sitec->soillayer_thickness[layer]);

		if (sprop->VWCeq[layer] > sprop->VWCfc[layer])
			ws->soilwFCEQ[layer] = sprop->VWCeq[layer] * water_density * sitec->soillayer_thickness[layer];
		else
			ws->soilwFCEQ[layer] = sprop->VWCfc[layer] * water_density * sitec->soillayer_thickness[layer];

		epv->WFPS[layer] = epv->VWC[layer] / sprop->VWCsat[layer];	

   
		/* PSI, hydrCONDUCT and hydrDIFFUS ( Cosby et al.) from VWC ([1MPa=100m] [m/s] [m2/s] */
		epv->PSI[layer]  = sprop->PSIsat[layer] * pow( (epv->VWC[layer] /sprop->VWCsat[layer]), -1* sprop->soilB[layer]);
		
	
		/* pF from PSI: cm from MPa */
		epv->pF[layer] =log10(fabs(10000*epv->PSI[layer]));


		VWC_avg += epv->VWC[layer] * (sitec->soillayer_thickness[layer] / sitec->soillayer_depth[N_SOILLAYERS - 1]);
	
	
		epv->hydrCONDUCTact[layer]  = sprop->hydrCONDUCTsat[layer] * pow(epv->VWC[layer] / sprop->VWCsat[layer], 2 * sprop->soilB[layer] + 3);
		epv->hydrDIFFUSact[layer]   = (((sprop->soilB[layer] * sprop->hydrCONDUCTsat[layer] * (-100 * sprop->PSIsat[layer]))) / sprop->VWCsat[layer])
			                           * pow(epv->VWC[layer] / sprop->VWCsat[layer], sprop->soilB[layer] + 2);

		/*----------------------------------------------*/
		/* relative soil water content (FC-WP) */

		if (sprop->VWCfc[layer] - sprop->VWCwp[layer] > CRIT_PREC_lenient)
			epv->relVWCfc_wp[layer] = (epv->VWC[layer] - sprop->VWCwp[layer]) / (sprop->VWCfc[layer] - sprop->VWCwp[layer]);
		else
		{
			if (!errorCode)
			{	
				printf("ERROR with relative VWC calculation, multilayer_hydrolparams.c\n");
				errorCode = 1;
			}
		}
		if (fabs(epv->relVWCfc_wp[layer]) < CRIT_PREC_lenient) epv->relVWCfc_wp[layer] = 0;

		/* relative soil water content (SAT-FC) */
		if (sprop->VWCsat[layer] - sprop->VWCfc[layer] > CRIT_PREC_lenient)
			epv->relVWCsat_fc[layer] = (epv->VWC[layer] - sprop->VWCfc[layer]) / (sprop->VWCsat[layer]  - sprop->VWCfc[layer]);
		else
		{
			if (!errorCode)
			{
				printf("ERROR with relative VWC calculation, multilayer_hydrolparams.c\n");
				errorCode = 1;
			}
		}
		if (fabs(epv->relVWCsat_fc[layer]) < CRIT_PREC_lenient) epv->relVWCsat_fc[layer] = 0;

	
		/*----------------------------------------------*/
		/* calculation of rootzone variables - weight of the last layer depends on the depth of the root */
		if (epv->n_maxrootlayers && layer < epv->n_maxrootlayers)
		{
			VWC_RZmax += epv->VWC[layer] * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[epv->n_maxrootlayers - 1];
			relVWCsat_fc_RZmax += epv->relVWCsat_fc[layer] * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[epv->n_maxrootlayers - 1];
			relVWCfc_wp_RZmax += epv->relVWCfc_wp[layer] * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[epv->n_maxrootlayers - 1];
		}

		if (epv->n_rootlayers)
		{
			weight_SUM += epv->rootlengthProp[layer];

			VWC_RZ += epv->VWC[layer] * epv->rootlengthProp[layer];
			PSI_RZ += epv->PSI[layer] * epv->rootlengthProp[layer];
			VWCsat_RZ += sprop->VWCsat[layer] * epv->rootlengthProp[layer];
			VWCfc_RZ += sprop->VWCfc[layer] * epv->rootlengthProp[layer];
			VWCwp_RZ += sprop->VWCwp[layer] * epv->rootlengthProp[layer];
			VWChw_RZ += sprop->VWChw[layer] * epv->rootlengthProp[layer];
			soilw_RZ += ws->soilw[layer] * epv->rootlengthProp[layer];
			if (epv->VWC[layer] > sprop->VWCwp[layer])
				soilwAVAIL_RZ += ((epv->VWC[layer] - sprop->VWCwp[layer]) * sitec->soillayer_thickness[layer] * water_density) * epv->rootlengthProp[layer];
		}
		else
		{
			VWC_RZ = 0;
			PSI_RZ = 0;
			VWCsat_RZ = 0;
			VWCfc_RZ = 0;
			VWCwp_RZ = 0;
			VWChw_RZ = 0;
			soilw_RZ = 0;
			soilwAVAIL_RZ = 0;
		}



	}

	if (epv->rootDepth && fabs(1 - weight_SUM) > CRIT_PREC && !errorCode)
	{
		printf("ERROR in calculation of rootzone variables (multilayer_hydrolprocess.c) \n");
		errorCode = 1;
	}
	epv->VWCsat_RZ = VWCsat_RZ;
	epv->VWCfc_RZ = VWCfc_RZ;
	epv->VWCwp_RZ = VWCwp_RZ;
	epv->VWChw_RZ = VWChw_RZ;
	epv->VWC_RZmax = VWC_RZmax;
	epv->relVWCsat_fc_RZmax = relVWCsat_fc_RZmax;
	epv->relVWCfc_wp_RZmax = relVWCfc_wp_RZmax;

	epv->VWC_avg = VWC_avg;
	epv->VWC_RZ = VWC_RZ;
	epv->PSI_RZ = PSI_RZ;
	ws->soilw_RZ = soilw_RZ;
	ws->soilwAVAIL_RZ = soilwAVAIL_RZ;



	return(errorCode);
}

