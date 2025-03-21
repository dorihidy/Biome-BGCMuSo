/* 
check_balance.c
daily test of mass balance (water, carbon, and nitrogen state variables)

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2022, D. Hidy [dori.hidy@gmail.com]
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

int check_water_balance(wstate_struct* ws, int first_balance)
{
	int errorCode=0;

	int layer;

	double balance, soilw_SUM, soilw_2m;
	balance=soilw_SUM=soilw_2m=0;
	static double old_balance;
	
	/* DAILY CHECK ON WATER BALANCE */

	/* control to avoid negative storage */
	for (layer=0; layer < N_SOILLAYERS; layer++)
	{
		if (ws->soilw[layer] < 0.0 && !errorCode)
		{
			printf("\n");
			printf("ERROR: negative soil water content\n");
			errorCode=1;
		}
		soilw_SUM += ws->soilw[layer];
		if (layer < 7) soilw_2m += ws->soilw[layer];

	}
	ws->soilw_SUM = soilw_SUM;
	ws->soilw_2m = soilw_2m;

	if (ws->snoww < 0.0 || ws->canopyw < 0  || ws->soilw_SUM < 0 || ws->pondw < 0)
	{
		printf("\n");
		printf("ERROR: negative water storage\n");
		errorCode=1;
	}
	
	/* sum of sources */
	ws->inW = ws->prcp_src + ws->GWsrc_W + ws->IRGsrc_W +  ws->FRZsrc_W + ws->FLsrc_W;
	

	
	/* sum of sinks */
	ws->outW = ws->soilEVP_snk + ws->snowSUBL_snk + ws->GWsnk_W +
		ws->EVPcanopyw_snk + ws->TRP_snk +
		ws->canopyw_THNsnk +			/* thinning */
		ws->canopyw_MOWsnk +			/* mowing */
		ws->canopyw_HRVsnk +			/* harvesting */
		ws->canopyw_PLGsnk +			/* ploughing */
		ws->canopyw_GRZsnk +			/* grazing */
		ws->runoff_snk     +			/* soil-water submodel */
		ws->deeppercolation_snk +		/* soil-water submodel .*/
		ws->pondEVP_snk;

	/* sum of current storage */
	ws->storeW = ws->soilw_SUM + ws->pondw + ws->snoww + ws->canopyw;
	
	/* calculate current balance */
	balance = ws->inW - ws->outW - ws->storeW;
	 
	/* calculate actual maximum balance error */
	if (!first_balance && (fabs(old_balance - balance) > ws->WbalanceERR))
	{
		ws->WbalanceERR = fabs(old_balance - balance);

	}
	old_balance = balance;



	return (errorCode);
}

int check_carbon_balance(cstate_struct* cs, int first_balance)
{
	int errorCode=0;
	int layer;

	static double old_balance;
	double balance;

	cs->soil1c_total = 0;
	cs->soil2c_total = 0;
	cs->soil3c_total = 0;
	cs->soil4c_total = 0;
	cs->litr1c_total = 0;
	cs->litr2c_total = 0;
	cs->litr3c_total = 0;
	cs->litr4c_total = 0;
	cs->cwdc_total = 0;

	/* control avoiding negative pools */
	if (cs->leafc < 0.0 ||  cs->leafc_storage < 0.0 || cs->leafc_transfer < 0.0 || 
		cs->frootc < 0.0 || cs->frootc_storage < 0.0 || cs->frootc_transfer < 0.0 || 
		cs->yieldc < 0.0 || cs->yieldc_storage < 0.0 || cs->yieldc_transfer < 0.0 || 
		cs->softstemc < 0.0 || cs->softstemc_storage < 0.0 || cs->softstemc_transfer < 0.0 || 
		cs->livestemc < 0.0 || cs->livestemc_storage < 0.0 || cs->livestemc_transfer < 0.0 || 
		cs->deadstemc < 0.0 || cs->deadstemc_storage < 0.0 || cs->deadstemc_transfer < 0.0 || 
		cs->livecrootc < 0.0 ||  cs->livecrootc_storage < 0.0 || cs->livecrootc_transfer < 0.0 || 
		cs->deadcrootc < 0.0 || cs->deadcrootc_storage < 0.0 || cs->deadcrootc_transfer < 0.0 || 
		cs->gresp_storage < 0.0 || cs->gresp_transfer < 0.0 ||
		cs->STDBc_leaf < 0 || cs->STDBc_froot < 0 || cs->STDBc_yield < 0|| cs->STDBc_softstem < 0 || 
		cs->CTDBc_leaf < 0 || cs->CTDBc_froot < 0 || cs->CTDBc_yield < 0|| cs->CTDBc_softstem < 0 || cs->CTDBc_cstem < 0  || cs->CTDBc_croot < 0)
		{
			printf("\n");
			printf("ERROR: negative carbon stock\n");
			errorCode=1;
		}

	
	/* summarizing cut-down and standing dead biomass */
	cs->CTDBc_above = cs->CTDBc_leaf + cs->CTDBc_yield + cs->CTDBc_softstem + cs->CTDBc_cstem;
	cs->CTDBc_below = cs->CTDBc_froot + cs->CTDBc_croot;

	cs->STDBc_above = cs->STDBc_leaf + cs->STDBc_yield + cs->STDBc_softstem;
	cs->STDBc_below = cs->STDBc_froot;



	/* summarizing soil and litter pools  */
	for (layer=0; layer < N_SOILLAYERS; layer++)
	{
	
		if ((cs->litr1c[layer] < 0.0 || cs->litr2c[layer] < 0.0 || cs->litr3c[layer] < 0.0 || cs-> litr4c[layer] < 0.0  || cs->cwdc[layer] < 0.0) && !errorCode)
		{	
			printf("ERROR: negative litter carbon stock\n");
			errorCode=1;
		}

		cs->soil1c_total += cs->soil1c[layer];
		cs->soil2c_total += cs->soil2c[layer];
		cs->soil3c_total += cs->soil3c[layer];
		cs->soil4c_total += cs->soil4c[layer];
		cs->soilC[layer] = cs->soil1c[layer] + cs->soil2c[layer] + cs->soil3c[layer] + cs->soil4c[layer];

		cs->litr1c_total += cs->litr1c[layer];
		cs->litr2c_total += cs->litr2c[layer];
		cs->litr3c_total += cs->litr3c[layer];
		cs->litr4c_total += cs->litr4c[layer];
		cs->litrC[layer] = cs->litr1c[layer] + cs->litr2c[layer] + cs->litr3c[layer] + cs->litr4c[layer];
		cs->cwdc_total += cs->cwdc[layer];


	}


	/* DAILY CHECK ON CARBON BALANCE */
	
	/* sum of sources: photosynthesis and managenet */
	cs->inC = cs->psnsun_src + cs->psnshade_src + cs->GWsrc_C + cs->FLsrc_C + cs->PLTsrc_C + cs->GRZsrc_C  + cs->FRZsrc_C + cs->MULsrc_C;

	
	
	/* sum of sinks: respiration, fire and management */
	cs->outC = cs->MRleaf_snk      + cs->GRleaf_snk      + cs->MRfroot_snk     + cs->GRfroot_snk + 
			   cs->MRyield_snk     + cs->GRyield_snk     + cs->MRsoftstem_snk  + cs->GRsoftstem_snk + 
			   cs->MRlivestem_snk  + cs->GRlivestem_snk  + cs->GRdeadstem_snk  + 
			   cs->MRlivecroot_snk + cs->GRlivecroot_snk + cs->GRdeadcroot_snk + 
			   cs->HRlitr1_snk + cs->HRlitr2_snk + cs->HRlitr4_snk + 
			   cs->HRsoil1_snk + cs->HRsoil2_snk + cs->HRsoil3_snk + cs->HRsoil4_snk + 
			   cs->FIREsnk_C +  cs->Cdeepleach_snk + 
		       cs->GWsnk_C +
			   cs->GRZsnk_C + cs->THN_transportC + cs->MOW_snkC + cs->HRV_snkC + cs->CWEsnk_C; 
		


		
	/* sum of current storage */
	cs->storeC = cs->leafc      + cs->leafc_storage      + cs->leafc_transfer +
				cs->frootc     + cs->frootc_storage     + cs->frootc_transfer + 
				cs->yieldc     + cs->yieldc_storage     + cs->yieldc_transfer +
				cs->softstemc  + cs->softstemc_storage  + cs->softstemc_transfer +
				cs->livestemc  + cs->livestemc_storage  + cs->livestemc_transfer + 
				cs->deadstemc  + cs->deadstemc_storage  + cs->deadstemc_transfer +
				cs->livecrootc + cs->livecrootc_storage + cs->livecrootc_transfer + 
				cs->deadcrootc + cs->deadcrootc_storage + cs->deadcrootc_transfer + 
				cs->gresp_storage + cs->gresp_transfer  + 
				cs->cwdc_total + 
				cs->litr1c_total + cs->litr2c_total + cs->litr3c_total + cs->litr4c_total + 
				cs->soil1c_total + cs->soil2c_total + cs->soil3c_total + cs->soil4c_total + 
				cs->cpool +
				cs->CTDBc_above + cs->STDBc_above + cs->CTDBc_below + cs->STDBc_below;   
	
	
	/* calculate current balance */
	balance = cs->inC - cs->outC - cs->storeC;
	 
	/* calculate actual maximum balance error */
 	if (!first_balance && (fabs(old_balance - balance) > cs->CbalanceERR))
	{
	 	cs->CbalanceERR = fabs(old_balance - balance);
	}
	old_balance = balance;


	return (errorCode);
}		

int check_nitrogen_balance(nstate_struct* ns, int first_balance)
{
	int errorCode=0;
	int layer=0;
	double balance;
	static double old_balance = 0.0;

	ns->NH4_total = 0;
	ns->NO3_total = 0;
	ns->litr1n_total = 0;
	ns->litr2n_total = 0;
	ns->litr3n_total = 0;
	ns->litr4n_total = 0;
	ns->cwdn_total = 0;
	ns->soil1n_total = 0;
	ns->soil2n_total = 0;
	ns->soil3n_total = 0;
	ns->soil4n_total = 0;
	
	/* CONTROL AVOIDING NITROGEN POOLS */
	if (ns->leafn < 0.0 ||  ns->leafn_storage < 0.0 || ns->leafn_transfer < 0.0 || 
		ns->frootn < 0.0 || ns->frootn_storage < 0.0 || ns->frootn_transfer < 0.0 || 
		ns->yieldn < 0.0 || ns->yieldn_storage < 0.0 || ns->yieldn_transfer < 0.0 || 
		ns->softstemn < 0.0 || ns->softstemn_storage < 0.0 || ns->softstemn_transfer < 0.0 || 
		ns->livestemn < 0.0 || ns->livestemn_storage < 0.0 || ns->livestemn_transfer < 0.0 || 
		ns->deadstemn < 0.0 || ns->deadstemn_storage < 0.0 || ns->deadstemn_transfer < 0.0 || 
		ns->livecrootn < 0.0 ||  ns->livecrootn_storage < 0.0 || ns->livecrootn_transfer < 0.0 || 
		ns->deadcrootn < 0.0 || ns->deadcrootn_storage < 0.0 || ns->deadcrootn_transfer < 0.0 ||
		ns->retransn < 0.0 ||
		ns->STDBn_leaf < 0 || ns->STDBn_froot < 0 || ns->STDBn_yield < 0|| ns->STDBn_softstem < 0 || 
		ns->CTDBn_leaf < 0 || ns->CTDBn_froot < 0 || ns->CTDBn_yield < 0|| ns->CTDBn_softstem < 0 || ns->CTDBn_cstem < 0  || ns->CTDBn_croot < 0)
	{	
		printf("\n");
		printf("ERROR: negative plant nitrogen pool\n");
		errorCode=1;
	}

	

	/* summarizing cut-down and standing dead biomass */
	ns->CTDBn_above = ns->CTDBn_leaf + ns->CTDBn_yield + ns->CTDBn_softstem + ns->CTDBn_cstem;
	ns->CTDBn_below = ns->CTDBn_froot + ns->CTDBn_croot;

	ns->STDBn_above = ns->STDBn_leaf + ns->STDBn_yield + ns->STDBn_softstem;
	ns->STDBn_below = ns->STDBn_froot;

	/* summarizing soil and litter pools  */
	for (layer=0; layer < N_SOILLAYERS; layer++)
	{
		if (ns->litr1n[layer] < 0.0 || ns->litr2n[layer] < 0.0 || ns->litr3n[layer] < 0.0 || ns-> litr4n[layer] < 0.0 || ns->cwdn[layer] < 0.0)
		{	
			printf("\n");
			printf("ERROR: negative litter nitrogen pool\n");
			errorCode=1;
		}


		ns->soil1n_total += ns->soil1n[layer];
		ns->soil2n_total += ns->soil2n[layer];
		ns->soil3n_total += ns->soil3n[layer];
		ns->soil4n_total += ns->soil4n[layer];
		ns->soilN[layer] = ns->soil1n[layer] + ns->soil2n[layer] + ns->soil3n[layer] + ns->soil4n[layer];

		ns->litr1n_total += ns->litr1n[layer];
		ns->litr2n_total += ns->litr2n[layer];
		ns->litr3n_total += ns->litr3n[layer];
		ns->litr4n_total += ns->litr4n[layer];
		ns->litrN[layer] = ns->litr1n[layer] + ns->litr2n[layer] + ns->litr3n[layer] + ns->litr4n[layer];
		ns->cwdn_total += ns->cwdn[layer];


		ns->NH4_total += ns->NH4[layer];
		ns->NO3_total += ns->NO3[layer];



	}
	

	/* DAILY CHECK ON NITROGEN BALANCE */
	
	/* sum of sources: fixation, deposition, spinup add, management */
	ns->inN = ns->Nfix_src + ns->Ndep_src + ns->SPINUPsrc + ns->GWsrc_N + ns->FLsrc_N +			            /* SUM of nitrogen plus from flooding */
			  ns->PLTsrc_N +  ns->GRZsrc_N +  ns->FRZsrc_N + ns->MULsrc_N;
	
	
	/* sum of sinks: volatilization, fire, deep leach, management */
	ns->outN = ns->Nvol_snk + ns->Nprec_snk + ns->FIREsnk_N + ns->Ndeepleach_snk + 
		       ns->GWsnk_N +
			   ns->GRZsnk_N + ns->THNsnk_N + ns->MOWsnk_N+  + ns->HRVsnk_N + ns->CWEsnk_N;


	/* sum of current storage */
	ns->storeN = ns->leafn      + ns->leafn_storage      + ns->leafn_transfer +
				 ns->frootn     + ns->frootn_storage     + ns->frootn_transfer + 			
				ns->yieldn     + ns->yieldn_storage     + ns->yieldn_transfer + 
				ns->softstemn  + ns->softstemn_storage  + ns->softstemn_transfer + 
				ns->livestemn  + ns->livestemn_storage  + ns->livestemn_transfer + 
				ns->deadstemn  + ns->deadstemn_storage  + ns->deadstemn_transfer + 
				ns->livecrootn + ns->livecrootn_storage + ns->livecrootn_transfer + 
				ns->deadcrootn + ns->deadcrootn_storage + ns->deadcrootn_transfer +
				ns->retransn + 
				ns->litr1n_total + ns->litr2n_total + ns->litr3n_total + ns->litr4n_total +
				ns->soil1n_total + ns->soil2n_total + ns->soil3n_total + ns->soil4n_total +
				ns->NH4_total  + ns->NO3_total  + ns->cwdn_total + 
				ns->npool        + 
				ns->CTDBn_above + ns->STDBn_above + ns->CTDBn_below + ns->STDBn_below;

	
	/* calculate current balance */
	balance = ns->inN - ns->outN - ns->storeN;
	 
	/* calculate actual maximum balance error */
	if (!first_balance && (fabs(old_balance - balance) > ns->NbalanceERR))
	{
		ns->NbalanceERR = fabs(old_balance - balance);

	}
	old_balance = balance;

	
	
	return (errorCode);
}

