/*
calc_drainage.c
Calculation of percolation fluxes

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
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

int calc_drainage(int flagRAIN, double soilB, double INFILT, double VWC, double VWCsat, double VWCfc, double dz0, double DC, double conductSAT_cmday, double* DRN, double* EXCESS, double* VWCnew)
{

	int errorCode = 0;
	double HOLD, DRAIN, DRNact, DRNx, conduct_cmday;
	*EXCESS = 0;

	
	if (flagRAIN)
	{
		/* [cm] = m3/m3 * cm */
		HOLD = (VWCsat - VWC) * dz0;


		/* 1.1.2.  IF: INFILT > HOLD */
		if (INFILT > 0.0 && INFILT > HOLD)
		{
			/* drainage from soil profile [cm = m3/m3 * cm ] */
			DRAIN = DC * (VWCsat - VWCfc) * dz0;


			/* drainage rate throug soil layer (cm/day) */
			DRNact = INFILT - HOLD + DRAIN;


			/* drainage is limited: actual conductance with updated soil water content */
			conduct_cmday = conductSAT_cmday * pow((VWC + INFILT / dz0) / VWCsat, 3 + 2 * soilB);
			if (conduct_cmday > conductSAT_cmday) conduct_cmday = conductSAT_cmday;
		

			if ((DRNact - conduct_cmday) > 0.0)
			{
				DRNact = conduct_cmday;
				DRAIN = DRNact + HOLD - INFILT;
			}


			/* state update temporal varialbe */
			VWC = VWC + (INFILT - DRNact) / dz0;

			/* above saturation - */
			if (VWC >= VWCsat)
			{
				*EXCESS = (VWC - VWCsat) * dz0;
				VWC = VWCsat;
			}

			INFILT = DRNact;

		} /* END IF: INFILT > HOLD */
		else
		{ /* 1.1.3. BEGIN ELSE: INFILT < HOLD */

			VWC = VWC + INFILT / dz0;

			/* drainage is limited: actual conductance with updated soil water content */
			conduct_cmday = conductSAT_cmday * pow(VWC / VWCsat, 3 + 2 * soilB);
			if (conduct_cmday > conductSAT_cmday) conduct_cmday = conductSAT_cmday;

			/* BEGIN IF-ELSE: VWC > FC */
			if (VWC >= VWCfc)
			{

				DRAIN = (VWC - VWCfc) * DC * dz0;

				/* drainage rate throug soil layer (cm/day) */
				if (DRAIN > CRIT_PREC_lenient)
					DRNact = DRAIN;
				else
					DRNact = 0;




				if ((DRNact - conduct_cmday) > 0.0)
				{
					DRNact = conduct_cmday;
					DRAIN = DRNact;
				}

				VWC = VWC - DRNact / dz0;
			}
			else
			{
				INFILT = 0.0;
				DRNact = 0.0;

			} /* END IF-ELSE: VWC > FC */

		} /* END ELSE: INFILT < HOLD */
	}
	else
	{
		DRNx = MAX((VWC - VWCfc) * DC * dz0, 0);

		
		if (VWC < VWCfc)
			HOLD = (VWCfc- VWC) * dz0;
		else
			HOLD = 0.0;

		DRNact = MAX(INFILT + DRNx - HOLD, 0.0);
	
		/* drainage is limited: actual conductance with updated soil water content */
		conduct_cmday = conductSAT_cmday * pow((VWC + INFILT / dz0) / VWCsat, 3 + 2 * soilB);
		if (conduct_cmday > conductSAT_cmday) conduct_cmday = conductSAT_cmday;

		if ((DRNact - conduct_cmday) > 0.0) DRNact = conduct_cmday;
		if (DRNact < CRIT_PREC_lenient) DRNact = 0;

		VWC += INFILT / dz0 - DRNact / dz0;

		/* above saturation - */
		if (VWC - VWCsat > 0)
		{
			*EXCESS = (VWC - VWCsat) * dz0;
			VWC = VWCsat;
		}


	}

	*DRN = DRNact;
	*VWCnew = VWC;

	return (errorCode);

}
