 /* 
multilayer_Tsoil.c
calculation of soil temperature in the different soil layers based on the change of air temperature (direct connection)
to top soil layer and based on empirical function of temperature gradient in soil (Zheng et al.1993)

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

int multilayer_Tsoil(const control_struct* ctrl, const epconst_struct* epc, const siteconst_struct* sitec, const soilprop_struct* sprop, epvar_struct* epv, int yday, double snoww,
					 metvar_struct* metv)

{
	int errorCode=0;
	int layer;

	/* effect of air temperature change on top soil temperature */
	double effect_of_vegetation = 1;
	double heatcoeff_snow = 0.1;
	double heatcoeff_nosnow = 0.25;
	double heating_coefficient;


	/* daily averaged air tempreture on the given day (calculated from Tmax and Tmin), temp.gradient and local temperatures */
	double temp_diff_total, temp_diff, Tsoil, Tsoil_avg;

	double STv1, STv2, WC, FX, f1, ALX, TA, Td, ZD;

	/* 4M parameter (dimless) empirical parameter for Tsoil esimation */
	double c_param_Tsoil = 4;
	
	Tsoil_avg=0;

	/* *********************************************************** */
	/* 1. FIRST LAYER PROPERTIES */
	/* surface soil temperature change caused by air temp. change can be estimated from the air temperature using empirical equations */	

	
	if (snoww) 
		heating_coefficient = heatcoeff_snow;
	else
		heating_coefficient = heatcoeff_nosnow;
	
	/* shading effect of vegetation (if soil temperature is lower than air temperature the effect is zero) */
	if (metv->Tday > metv->Tsoil_surface && (epv->projLAI + epv->projLAI_STDB) > 0)
		effect_of_vegetation = exp(-1 * epc->ext_coef * (epv->projLAI+epv->projLAI_STDB));
	else 
		effect_of_vegetation = 1.0;

   if (effect_of_vegetation < 0.5) effect_of_vegetation = 0.5;

	/* empirical function for the effect of tair changing */
	metv->Tsoil_top_change = (metv->Tday - metv->Tsoil_surface) * heating_coefficient * effect_of_vegetation;


	epv->ES = heating_coefficient;
	epv->EV = effect_of_vegetation;
	epv->Tday = metv->Tday;

	/* ************************************************- */
	/* 2. TEMPERATURE OF DEEPER LAYER BASED ON TEMPERATURE GRADIENT BETWEEN SURFACE LAYER AND LOWERMOST LAYER (BELOW 3M) */

	/* on the first day the temperature of the soil layers are calculated based on the temperature of top and bottom layer */

	metv->Tsoil_surface += metv->Tsoil_top_change;
	
	temp_diff_total = metv->annTavgRA - metv->Tsoil_surface;
	
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		
	//	temp_diff = temp_diff_total * (0.1526 * log(depth) + 0.703);	
		temp_diff = temp_diff_total * (0.1448 * log(sitec->soillayer_midpoint[layer]) + 0.6667);
		metv->Tsoil[layer] = metv->Tsoil_surface + temp_diff;	

		/* unit change: rootzone_depth: m to cm, BD: kg/m3 to g/cm3 */
		STv1 = 1000 + 2500 * sprop->BD[layer]/(((sprop->BD[layer] / g_per_cm3_to_kg_per_m3) + 686 * exp(-5.63*sprop->BD[layer])));
		STv2 = log(500/STv1);
		WC = epv->VWC_avg / ((0.356-0.144*(sprop->BD[layer] / g_per_cm3_to_kg_per_m3)) * sitec->soillayer_depth[N_SOILLAYERS-2]*100); 
		FX = exp(STv2*pow((1-WC)/(1+WC),2));
		f1 = 1/(FX*STv1);

		ALX = 0.0174*(yday-200);
		TA = metv->annTavgRA + metv->annTrangeRA * cos(ALX)/2;
		Td = metv->tempradFra - TA;
	
		ZD = -1 * sitec->soillayer_midpoint[layer] * 1000 * f1 * c_param_Tsoil; // m to mm
	
		Tsoil = metv->annTavgRA + (metv->annTrangeRA/2 * cos(ALX + ZD) + Td)  * exp(ZD); // depth: m to cm 

		if (ctrl->STCM_flag) 
			metv->Tsoil[layer] = Tsoil;

		if (layer < N_SOILLAYERS-1) Tsoil_avg += metv->Tsoil[layer] * (sitec->soillayer_thickness[layer] / sitec->soillayer_depth[N_SOILLAYERS-2]);

	}

    metv->Tsoil_avg         = Tsoil_avg;

	return (errorCode);
}
	