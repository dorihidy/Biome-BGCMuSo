/* 
mulching.c
Calculation of mulching fluyes (cover surface with mulch layer)

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
#include "pointbgc_struct.h"
#include "bgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"

int mulching(control_struct* ctrl, const mulching_struct* MUL,  cstate_struct* cs, nstate_struct*ns, cflux_struct* cf, nflux_struct* nf)
{

	int errorCode=0;
	int md, year;
	double litrCabove_MUL, litrCNabove_MUL, cwdCabove_MUL, cwdCNabove_MUL;			            

	litrCabove_MUL = litrCNabove_MUL = cwdCabove_MUL = cwdCNabove_MUL = 0;
	year = ctrl->simstartyear + ctrl->simyr;
	md = MUL->mgmdMUL-1;


	if (MUL->MUL_num && md >= 0)
	{
		if (year == MUL->MULyear_array[md] && ctrl->month == MUL->MULmonth_array[md] && ctrl->day == MUL->MULday_array[md]) 
		{
			
			litrCabove_MUL = MUL->litrCabove_MUL[md];
			litrCNabove_MUL = MUL->litrCNabove_MUL[md];
			cwdCabove_MUL = MUL->cwdCabove_MUL[md];
			cwdCNabove_MUL = MUL->cwdCNabove_MUL[md];

			cf->litrc_from_MUL = litrCabove_MUL;
			cf->cwdc_from_MUL  = cwdCabove_MUL;
			if (cf->litrc_from_MUL > 0) nf->litrn_from_MUL = litrCabove_MUL / litrCNabove_MUL;
			if (cf->cwdc_from_MUL > 0)  nf->cwdn_from_MUL  = cwdCabove_MUL / cwdCNabove_MUL;

			cs->litr1c[0] += cf->litrc_from_MUL;
			ns->litr1n[0] += nf->litrn_from_MUL;
		
			cs->cwdc[0]  += cf->cwdc_from_MUL;	
			ns->cwdn[0]  += nf->cwdn_from_MUL;

			cs->litrCabove[0] += cf->litrc_from_MUL;
			cs->cwdCabove[0] += cf->cwdc_from_MUL;

			cs->litrCabove_total += cf->litrc_from_MUL;
			cs->cwdCabove_total += cf->cwdc_from_MUL;

	
			cs->MULsrc_C += cf->litrc_from_MUL + cf->cwdc_from_MUL;
			ns->MULsrc_N += nf->litrn_from_MUL + nf->cwdn_from_MUL;


		}
		else
		{
			cf->litrc_from_MUL = 0;
			cf->cwdc_from_MUL = 0;
		}
		
	
	}

	return (errorCode);
}