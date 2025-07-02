 /*
calc_VWCeq.c
Calculation of equilibrium VWC from groudnwater depth (or from lower boundary condition of VWC)

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
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_func.h"
#include "bgc_constants.h"


int calc_VWCeq(control_struct* ctrl, const siteconst_struct* sitec, soilprop_struct* sprop)
{

	int layer;
	double Zbottom, Ztop, Dtop, Dbottom, Weq, hsat, h, Zcoeff, VWCeq;
	static double CFDact, GWDact;

	int errorCode = 0;

	/* --------------------------------------------------------------------------------------------------------------*/
	/* calculation of GWD (in case of no-GW simulation: GWD-estimation form PSI sat and lower boundary condition of VWC: ctrl->wstate_propFC (user can set in WSTATE_INIT) */
	if (sprop->GWD == DATA_GAP)
	{
		if (ctrl->firstsimday_flag)
		{
			/* Soil system is in equlibrium under bottom layer -> VWCbottom = VWCeq -> depth of GW can be calculated: h = hsat * (VWC / VWCsat) ^ -soilB.
				h : soil water potential (9distance from GWD)  in m
				hsat: distance of capillary zone from GWD (can be calculated from PSIsat: 1 MPa = 100 m) */
			hsat = -1 * sprop->PSIsat[N_SOILLAYERS - 1] * 100;
			h = hsat * pow(((sprop->VWCfc[N_SOILLAYERS - 1] * ctrl->wstate_propFC) / sprop->VWCsat[N_SOILLAYERS - 1]), -1 * sprop->soilB[N_SOILLAYERS - 1]);
			GWDact = sitec->soillayer_depth[N_SOILLAYERS - 1] + h;
			CFDact = GWDact - sprop->CapillFringe[N_SOILLAYERS - 1];


		}
	}
	else
	{
		GWDact = sprop->GWD;
		CFDact = sprop->CFD;
	}

	/* --------------------------------------------------------------------------------------------------------------*/
	/* calculation of GWD for VWCeq calculation (in case of GWD-simulation: input data, else: calculated data */

	if (ctrl->firstsimday_flag || sprop->GWD != sprop->GWD_pre)
	{
		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{

			/* distance from groundwater or boundary layer */
			if (layer == 0)
				Dtop = 0;
			else
				Dtop = sitec->soillayer_depth[layer - 1];

			Dbottom = sitec->soillayer_depth[layer];

			if (Dtop < CFDact)
			{
				if (layer > 0)
					Ztop = (GWDact - sitec->soillayer_depth[layer - 1]) * m_to_cm;
				else
					Ztop = GWDact * 100;

				if (Dbottom < CFDact)
				{
					Zbottom = (GWDact - sitec->soillayer_depth[layer]) * m_to_cm;

					Zcoeff = (pow(Ztop, 1 - 1. / sprop->soilB[layer]) - pow(Zbottom, 1 - 1. / sprop->soilB[layer])) / (1 - 1. / sprop->soilB[layer]);
					Weq = sprop->VWCsat[layer] * pow(-10000 * sprop->PSIsat[layer], 1. / sprop->soilB[layer]) * Zcoeff;
					sprop->VWCeq[layer] = Weq / (sitec->soillayer_thickness[layer] * m_to_cm);
				}
				else
				{
					Zbottom = (GWDact - CFDact) * m_to_cm;

					Zcoeff = (pow(Ztop, 1 - 1. / sprop->soilB[layer]) - pow(Zbottom, 1 - 1. / sprop->soilB[layer])) / (1 - 1. / sprop->soilB[layer]);
					Weq = sprop->VWCsat[layer] * pow(-10000 * sprop->PSIsat[layer], 1. / sprop->soilB[layer]) * Zcoeff;
					VWCeq = Weq / (Ztop - Zbottom);
					sprop->VWCeq[layer] = VWCeq;
				}
				 
				/* lower limitation of VWCeq: wilting point */
				if (sprop->VWCeq[layer] < 1.01 * sprop->VWCwp[layer]) sprop->VWCeq[layer] = 1.01 * sprop->VWCwp[layer];
		

			}
			else
			{
				sprop->VWCeq[layer] = sprop->VWCsat[layer];
			}

		}
	}

	return (errorCode);

}
