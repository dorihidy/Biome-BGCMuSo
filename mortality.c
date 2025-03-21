/* 
mortality.c
daily mortality fluxes

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2022, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

Updated:
6 March 2000 (PET): Appears from tests with P.Pine data that combustion of deadstem should be higher: setting to 20% for testing.
9 February 1999, PET: Changed the treatment of fire mortality, so that instead of the entire deadstem and cwd pools being subject to the mortality rates,
only 5% of deadstem and 30% of cwd are subject to fire losses. Mike White found that there was a big difference between literature estimates of fire
emissions and model estimates for the boreal forest, and this is also a topic that Bob Keane had pointed out a while ago as a first-order improvement
to the fire treatment. Got the values for consumption efficiency from Agee, James K., 1993. Fire Ecology of Pacific Northwest Forests. Island Press, 
Washington, D.C. p 42.
The deadstem material that is not wilted (95%) is sent to CWD pools. CWD that is not burned (70%) stays in CWD pools.

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

int mortality(const control_struct* ctrl, const soilprop_struct* sprop, const epconst_struct* epc,
	          epvar_struct* epv, cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf, int simyr)
{
	int errorCode=0;
	int layer;
	double flux_from_carbon;
	

	/******************************************************************/
	/* 0. if no changing data constant EPC parameter are used - ATTENTION: WPM = WPMyr / nDAYS_OF_YEAR */

	if (ctrl->varWPM_flag)
		epv->WPM = epc->WPMyr_array[simyr]/nDAYS_OF_YEAR;
	else
		epv->WPM = epc->WPMyr/nDAYS_OF_YEAR;	


	if (ctrl->varFM_flag)
		epv->FM = epc->FMyr_array[simyr]/nDAYS_OF_YEAR;
	else
		epv->FM = epc->FMyr/nDAYS_OF_YEAR;	


	/******************************************************************/
	/* I. Non-fire mortality: these fluxes all enter litter or CWD pools */
	/******************************************************************/


	
	
	/* 2. daily fluxes due to mortality */

	if (epv->WPM)
	{
		/* in order to save the C:N ratio: N-fluxes are calculated from C-fluxes using C:N ratio parameters */
		if (epc->leaf_cn)
		{
			cf->m_leafc_to_litr1c  = epv->WPM * cs->leafc * epc->leaflitr_flab;  	 
			cf->m_leafc_to_litr2c  = epv->WPM * cs->leafc * epc->leaflitr_fucel;
			cf->m_leafc_to_litr3c  = epv->WPM * cs->leafc * epc->leaflitr_fscel;  	 
			cf->m_leafc_to_litr4c  = epv->WPM * cs->leafc * epc->leaflitr_flig;  	 

			nf->m_leafn_to_litr1n  = cf->m_leafc_to_litr1c / epc->leaf_cn;  	
			nf->m_leafn_to_litr2n  = cf->m_leafc_to_litr2c / epc->leaf_cn;  	
			nf->m_leafn_to_litr3n  = cf->m_leafc_to_litr3c / epc->leaf_cn;  	
			nf->m_leafn_to_litr4n  = cf->m_leafc_to_litr4c / epc->leaf_cn;  

			cf->m_leafc_storage_to_litr1c	= epv->WPM * cs->leafc_storage;
			cf->m_leafc_transfer_to_litr1c  = epv->WPM * cs->leafc_transfer;
		
			nf->m_leafn_storage_to_litr1n	= cf->m_leafc_storage_to_litr1c / epc->leaf_cn;
			nf->m_leafn_transfer_to_litr1n	= cf->m_leafc_transfer_to_litr1c / epc->leaf_cn;
	
		}
	
		if (epc->froot_cn)
		{
			cf->m_frootc_to_litr1c = epv->WPM * cs->frootc * epc->frootlitr_flab; 	 
			cf->m_frootc_to_litr2c = epv->WPM * cs->frootc * epc->frootlitr_fucel;
			cf->m_frootc_to_litr3c = epv->WPM * cs->frootc * epc->frootlitr_fscel; 	 
			cf->m_frootc_to_litr4c = epv->WPM * cs->frootc * epc->frootlitr_flig;

			nf->m_frootn_to_litr1n = cf->m_frootc_to_litr1c / epc->froot_cn; 	
			nf->m_frootn_to_litr2n = cf->m_frootc_to_litr2c / epc->froot_cn; 	
			nf->m_frootn_to_litr3n = cf->m_frootc_to_litr3c / epc->froot_cn; 	
			nf->m_frootn_to_litr4n = cf->m_frootc_to_litr4c / epc->froot_cn;

			cf->m_frootc_storage_to_litr1c	= epv->WPM * cs->frootc_storage;
			cf->m_frootc_transfer_to_litr1c	= epv->WPM * cs->frootc_transfer;
			nf->m_frootn_storage_to_litr1n	= cf->m_frootc_storage_to_litr1c / epc->froot_cn;
			nf->m_frootn_transfer_to_litr1n	= cf->m_frootc_transfer_to_litr1c / epc->froot_cn;

		}
	
		if (epc->yield_cn)
		{
			cf->m_yieldc_to_litr1c = epv->WPM * cs->yieldc * epc->yieldlitr_flab;  	 
			cf->m_yieldc_to_litr2c = epv->WPM * cs->yieldc * epc->yieldlitr_fucel;
			cf->m_yieldc_to_litr3c = epv->WPM * cs->yieldc * epc->yieldlitr_fscel;  	 
			cf->m_yieldc_to_litr4c = epv->WPM * cs->yieldc * epc->yieldlitr_flig;  	

			nf->m_yieldn_to_litr1n = cf->m_yieldc_to_litr1c / epc->yield_cn;  	
			nf->m_yieldn_to_litr2n = cf->m_yieldc_to_litr2c / epc->yield_cn;  	
			nf->m_yieldn_to_litr3n = cf->m_yieldc_to_litr3c / epc->yield_cn;  	
			nf->m_yieldn_to_litr4n = cf->m_yieldc_to_litr4c / epc->yield_cn;  	

			cf->m_yieldc_storage_to_litr1c	= epv->WPM * cs->yieldc_storage;
			cf->m_yieldc_transfer_to_litr1c	= epv->WPM * cs->yieldc_transfer;
		
			nf->m_yieldn_storage_to_litr1n	= cf->m_yieldc_storage_to_litr1c / epc->yield_cn;
			nf->m_yieldn_transfer_to_litr1n	= cf->m_yieldc_transfer_to_litr1c / epc->yield_cn;
	
		}

		if (epc->softstem_cn)
		{
			cf->m_softstemc_to_litr1c = epv->WPM * cs->softstemc * epc->softstemlitr_flab;  	 
			cf->m_softstemc_to_litr2c = epv->WPM * cs->softstemc * epc->softstemlitr_fucel;
			cf->m_softstemc_to_litr3c = epv->WPM * cs->softstemc * epc->softstemlitr_fscel;  	 
			cf->m_softstemc_to_litr4c = epv->WPM * cs->softstemc * epc->softstemlitr_flig;  	

			nf->m_softstemn_to_litr1n = cf->m_softstemc_to_litr1c / epc->softstem_cn;  	
			nf->m_softstemn_to_litr2n = cf->m_softstemc_to_litr2c / epc->softstem_cn;  	
			nf->m_softstemn_to_litr3n = cf->m_softstemc_to_litr3c / epc->softstem_cn;  	
			nf->m_softstemn_to_litr4n = cf->m_softstemc_to_litr4c / epc->softstem_cn;  	

			cf->m_softstemc_storage_to_litr1c	= epv->WPM * cs->softstemc_storage;
			cf->m_softstemc_transfer_to_litr1c	= epv->WPM * cs->softstemc_transfer;
		
			nf->m_softstemn_storage_to_litr1n	= cf->m_softstemc_storage_to_litr1c / epc->softstem_cn;
			nf->m_softstemn_transfer_to_litr1n	= cf->m_softstemc_transfer_to_litr1c / epc->softstem_cn;
	
		}
	
		cf->m_gresp_storage_to_litr1c		= epv->WPM * cs->gresp_storage;
		cf->m_gresp_transfer_to_litr1c		= epv->WPM * cs->gresp_transfer;	

		nf->m_retransn_to_litr1n			= epv->WPM * ns->retransn;


		/* woody biomass */
		if (epc->livewood_cn)
		{
			/* stem */
			cf->m_livestemc_to_cwdc    = epv->WPM * cs->livestemc;	
			nf->m_livestemn_to_cwdn    = cf->m_livestemc_to_cwdc  / epc->livewood_cn;	
			nf->m_livestemn_to_litr1n  = (epv->WPM * ns->livestemn) - nf->m_livestemn_to_cwdn;
			
			cf->m_livestemc_storage_to_litr1c	= epv->WPM * cs->livestemc_storage;
			cf->m_livestemc_transfer_to_litr1c	= epv->WPM * cs->livestemc_transfer;

			nf->m_livestemn_storage_to_litr1n	= cf->m_livestemc_storage_to_litr1c   / epc->livewood_cn;
			nf->m_livestemn_transfer_to_litr1n	= cf->m_livestemc_transfer_to_litr1c  / epc->livewood_cn;

			/* root */
			cf->m_livecrootc_to_cwdc   = epv->WPM * cs->livecrootc;   
			nf->m_livecrootn_to_cwdn   = cf->m_livecrootc_to_cwdc / epc->livewood_cn;  
			nf->m_livecrootn_to_litr1n = (epv->WPM * ns->livecrootn) - nf->m_livecrootn_to_cwdn;

			cf->m_livecrootc_storage_to_litr1c  = epv->WPM * cs->livecrootc_storage;
			cf->m_livecrootc_transfer_to_litr1c = epv->WPM * cs->livecrootc_transfer;
		
			nf->m_livecrootn_storage_to_litr1n  = cf->m_livecrootc_storage_to_litr1c  / epc->livewood_cn;
			nf->m_livecrootn_transfer_to_litr1n = cf->m_livecrootc_transfer_to_litr1c / epc->livewood_cn;
			
		}
		if (epc->deadwood_cn)
		{
			/* deadwood */
			cf->m_deadstemc_to_cwdc    = (epv->WPM + (1.0-epc->dscp)* epv->FM) * cs->deadstemc;	 
			cf->m_deadcrootc_to_cwdc   = (epv->WPM + (1.0-epc->dscp)* epv->FM) * cs->deadcrootc;   

			nf->m_deadstemn_to_cwdn    = cf->m_deadstemc_to_cwdc  / epc->deadwood_cn;	
			nf->m_deadcrootn_to_cwdn   = cf->m_deadcrootc_to_cwdc / epc->deadwood_cn;  

			cf->m_deadstemc_storage_to_litr1c	= epv->WPM * cs->deadstemc_storage;
			cf->m_deadstemc_transfer_to_litr1c	= epv->WPM * cs->deadstemc_transfer;
		
			nf->m_deadstemn_storage_to_litr1n	= cf->m_deadstemc_storage_to_litr1c   / epc->deadwood_cn;
			nf->m_deadstemn_transfer_to_litr1n	= cf->m_deadstemc_transfer_to_litr1c  / epc->deadwood_cn;

			cf->m_deadcrootc_storage_to_litr1c  = epv->WPM * cs->deadcrootc_storage;
			cf->m_deadcrootc_transfer_to_litr1c = epv->WPM * cs->deadcrootc_transfer;
		
			nf->m_deadcrootn_storage_to_litr1n  = cf->m_deadcrootc_storage_to_litr1c  / epc->deadwood_cn;
			nf->m_deadcrootn_transfer_to_litr1n = cf->m_deadcrootc_transfer_to_litr1c / epc->deadwood_cn;
			
		}
	
	
		/* 3. UPDATE STATE VARIALBES */

		/* Mortality is taken care of last and given special treatment for state update so that it doesn't interfere
		with the other fluxes that are based on proportions of state variables, especially the phenological fluxes */

		/* 3.1 ABOVEGROUND BIOMASS  - into the first layer in multilayer soil */
		cs->leafc       -= (cf->m_leafc_to_litr1c     + cf->m_leafc_to_litr2c     + cf->m_leafc_to_litr3c     + cf->m_leafc_to_litr4c);
		cs->yieldc      -= (cf->m_yieldc_to_litr1c    + cf->m_yieldc_to_litr2c    + cf->m_yieldc_to_litr3c    + cf->m_yieldc_to_litr4c);
		cs->softstemc   -= (cf->m_softstemc_to_litr1c + cf->m_softstemc_to_litr2c + cf->m_softstemc_to_litr3c + cf->m_softstemc_to_litr4c);
		cs->livestemc	-= cf->m_livestemc_to_cwdc;
		cs->deadstemc	-= cf->m_deadstemc_to_cwdc;


		/* litter turns into the first three soil layers  (non-woody biomass: proportion to soil layer thickness, woody-biomass: higher propotion in layer2 */

		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{
			cs->litr1c[layer] += (cf->m_leafc_to_litr1c + cf->m_yieldc_to_litr1c + cf->m_softstemc_to_litr1c) * sprop->PROPlayerDC[layer];
			cs->litr2c[layer] += (cf->m_leafc_to_litr2c + cf->m_yieldc_to_litr2c + cf->m_softstemc_to_litr2c) * sprop->PROPlayerDC[layer];
			cs->litr3c[layer] += (cf->m_leafc_to_litr3c + cf->m_yieldc_to_litr3c + cf->m_softstemc_to_litr3c) * sprop->PROPlayerDC[layer];
			cs->litr4c[layer] += (cf->m_leafc_to_litr4c + cf->m_yieldc_to_litr4c + cf->m_softstemc_to_litr4c) * sprop->PROPlayerDC[layer];
			cs->cwdc[layer] += (cf->m_livestemc_to_cwdc + cf->m_deadstemc_to_cwdc) * sprop->PROPlayerDC[layer];

			ns->litr1n[layer] += (nf->m_leafn_to_litr1n + nf->m_yieldn_to_litr1n + nf->m_softstemn_to_litr1n) * sprop->PROPlayerDC[layer];
			ns->litr2n[layer] += (nf->m_leafn_to_litr2n + nf->m_yieldn_to_litr2n + nf->m_softstemn_to_litr2n) * sprop->PROPlayerDC[layer];
			ns->litr3n[layer] += (nf->m_leafn_to_litr3n + nf->m_yieldn_to_litr3n + nf->m_softstemn_to_litr3n) * sprop->PROPlayerDC[layer];
			ns->litr4n[layer] += (nf->m_leafn_to_litr4n + nf->m_yieldn_to_litr4n + nf->m_softstemn_to_litr4n) * sprop->PROPlayerDC[layer];
			ns->cwdn[layer] += (nf->m_livestemn_to_cwdn + nf->m_deadstemn_to_cwdn) * sprop->PROPlayerDC[layer];

			ns->litr1n[layer] += nf->m_livestemn_to_litr1n * sprop->PROPlayerDC[layer];
			ns->litr1n[layer] += nf->m_livecrootn_to_litr1n * sprop->PROPlayerDC[layer];
		}


	
		ns->leafn       -= (nf->m_leafn_to_litr1n     + nf->m_leafn_to_litr2n     + nf->m_leafn_to_litr3n     + nf->m_leafn_to_litr4n);
		ns->yieldn      -= (nf->m_yieldn_to_litr1n    + nf->m_yieldn_to_litr2n    + nf->m_yieldn_to_litr3n    + nf->m_yieldn_to_litr4n);
		ns->softstemn   -= (nf->m_softstemn_to_litr1n + nf->m_softstemn_to_litr2n + nf->m_softstemn_to_litr3n + nf->m_softstemn_to_litr4n);
		ns->livestemn	-= nf->m_livestemn_to_cwdn;
		ns->deadstemn	-= nf->m_deadstemn_to_cwdn;


	
		/* special N-flux: live woody biomass to litter  */
		ns->livestemn  -= nf->m_livestemn_to_litr1n;
		ns->livecrootn -= nf->m_livecrootn_to_litr1n;
	
		/* 3.2 NON-STRUCTURED (transfer, storage, retrans  - into the first, labile layer in multilayer soil */
		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{
			cs->litr1c[layer] += (cf->m_leafc_storage_to_litr1c + cf->m_frootc_storage_to_litr1c + cf->m_yieldc_storage_to_litr1c + cf->m_softstemc_storage_to_litr1c +
				cf->m_livestemc_storage_to_litr1c + cf->m_deadstemc_storage_to_litr1c + cf->m_livecrootc_storage_to_litr1c + cf->m_deadcrootc_storage_to_litr1c +
				cf->m_leafc_transfer_to_litr1c + cf->m_frootc_transfer_to_litr1c + cf->m_yieldc_transfer_to_litr1c + cf->m_softstemc_transfer_to_litr1c +
				cf->m_livestemc_transfer_to_litr1c + cf->m_deadstemc_transfer_to_litr1c + cf->m_livecrootc_transfer_to_litr1c + cf->m_deadcrootc_transfer_to_litr1c +
				cf->m_gresp_storage_to_litr1c + cf->m_gresp_transfer_to_litr1c) * sprop->PROPlayerDC[layer];



			ns->litr1n[layer] += (nf->m_leafn_storage_to_litr1n + nf->m_frootn_storage_to_litr1n + nf->m_yieldn_storage_to_litr1n + nf->m_softstemn_storage_to_litr1n +
				nf->m_livestemn_storage_to_litr1n + nf->m_deadstemn_storage_to_litr1n + nf->m_livecrootn_storage_to_litr1n + nf->m_deadcrootn_storage_to_litr1n +
				nf->m_leafn_transfer_to_litr1n + nf->m_frootn_transfer_to_litr1n + nf->m_yieldn_transfer_to_litr1n + nf->m_softstemn_transfer_to_litr1n +
				nf->m_livestemn_transfer_to_litr1n + nf->m_deadstemn_transfer_to_litr1n + nf->m_livecrootn_transfer_to_litr1n + nf->m_deadcrootn_transfer_to_litr1n +
				nf->m_retransn_to_litr1n) * sprop->PROPlayerDC[layer];

		}


	
		cs->leafc_storage       -= cf->m_leafc_storage_to_litr1c;
		cs->frootc_storage      -= cf->m_frootc_storage_to_litr1c;
		cs->yieldc_storage      -= cf->m_yieldc_storage_to_litr1c;
		cs->softstemc_storage   -= cf->m_softstemc_storage_to_litr1c;
		cs->livestemc_storage   -= cf->m_livestemc_storage_to_litr1c;
		cs->deadstemc_storage   -= cf->m_deadstemc_storage_to_litr1c;
		cs->livecrootc_storage  -= cf->m_livecrootc_storage_to_litr1c;
		cs->deadcrootc_storage  -= cf->m_deadcrootc_storage_to_litr1c;
	
		cs->leafc_transfer      -= cf->m_leafc_transfer_to_litr1c;
		cs->frootc_transfer     -= cf->m_frootc_transfer_to_litr1c;
		cs->yieldc_transfer     -= cf->m_yieldc_transfer_to_litr1c;
		cs->softstemc_transfer  -= cf->m_softstemc_transfer_to_litr1c;
		cs->livestemc_transfer  -= cf->m_livestemc_transfer_to_litr1c;
		cs->deadstemc_transfer  -= cf->m_deadstemc_transfer_to_litr1c;
		cs->livecrootc_transfer -= cf->m_livecrootc_transfer_to_litr1c;
		cs->deadcrootc_transfer -= cf->m_deadcrootc_transfer_to_litr1c;
	
		cs->gresp_storage       -= cf->m_gresp_storage_to_litr1c;
		cs->gresp_transfer      -= cf->m_gresp_transfer_to_litr1c;

		
		ns->leafn_storage       -= nf->m_leafn_storage_to_litr1n;
		ns->frootn_storage      -= nf->m_frootn_storage_to_litr1n;
		ns->yieldn_storage      -= nf->m_yieldn_storage_to_litr1n;
		ns->softstemn_storage   -= nf->m_softstemn_storage_to_litr1n;
		ns->livestemn_storage   -= nf->m_livestemn_storage_to_litr1n;
		ns->deadstemn_storage   -= nf->m_deadstemn_storage_to_litr1n;
		ns->livecrootn_storage  -= nf->m_livecrootn_storage_to_litr1n;
		ns->deadcrootn_storage  -= nf->m_deadcrootn_storage_to_litr1n;
	
		ns->leafn_transfer      -= nf->m_leafn_transfer_to_litr1n;
		ns->frootn_transfer     -= nf->m_frootn_transfer_to_litr1n;
		ns->yieldn_transfer     -= nf->m_yieldn_transfer_to_litr1n;
		ns->softstemn_transfer  -= nf->m_softstemn_transfer_to_litr1n;
		ns->livestemn_transfer  -= nf->m_livestemn_transfer_to_litr1n;
		ns->deadstemn_transfer  -= nf->m_deadstemn_transfer_to_litr1n;
		ns->livecrootn_transfer -= nf->m_livecrootn_transfer_to_litr1n;
		ns->deadcrootn_transfer -= nf->m_deadcrootn_transfer_to_litr1n;
			
		ns->retransn       -= nf->m_retransn_to_litr1n;


	
		/* 3.3 BELOWGROUND BIOMASS  - DIVIDED BETWEEN THE DIFFERENT SOIL LAYERS in multilayer soil  */
	
		cs->frootc      -= (cf->m_frootc_to_litr1c    + cf->m_frootc_to_litr2c    + cf->m_frootc_to_litr3c    + cf->m_frootc_to_litr4c);
		cs->livecrootc	-=  cf->m_livecrootc_to_cwdc;
		cs->deadcrootc	-=  cf->m_deadcrootc_to_cwdc;

		ns->frootn     -= (nf->m_frootn_to_litr1n    + nf->m_frootn_to_litr2n    + nf->m_frootn_to_litr3n    + nf->m_frootn_to_litr4n);
		ns->livecrootn -=  nf->m_livecrootn_to_cwdn;
		ns->deadcrootn -=  nf->m_deadcrootn_to_cwdn;

		if (epv->rootDepth > CRIT_PREC)
		{
			for (layer = 0; layer < N_SOILLAYERS; layer++)
			{
				cs->litr1c[layer]  += (cf->m_frootc_to_litr1c) * epv->rootlengthProp[layer];	
				cs->litr2c[layer]  += (cf->m_frootc_to_litr2c) * epv->rootlengthProp[layer];	
				cs->litr3c[layer]  += (cf->m_frootc_to_litr3c) * epv->rootlengthProp[layer];	
				cs->litr4c[layer]  += (cf->m_frootc_to_litr4c) * epv->rootlengthProp[layer];
				cs->cwdc[layer]    += (cf->m_livecrootc_to_cwdc + cf->m_deadcrootc_to_cwdc) * epv->rootlengthProp[layer];

				ns->litr1n[layer]  += (nf->m_frootn_to_litr1n) * epv->rootlengthProp[layer];	
				ns->litr2n[layer]  += (nf->m_frootn_to_litr2n) * epv->rootlengthProp[layer];	
				ns->litr3n[layer]  += (nf->m_frootn_to_litr3n) * epv->rootlengthProp[layer];	
				ns->litr4n[layer]  += (nf->m_frootn_to_litr4n) * epv->rootlengthProp[layer];
				ns->cwdn[layer]    += (nf->m_livecrootn_to_cwdn + nf->m_deadcrootn_to_cwdn) * epv->rootlengthProp[layer];	
			}
		}
		else
		{
			for (layer = 0; layer < N_SOILLAYERS; layer++)
			{
				cs->litr1c[layer]  += (cf->m_frootc_to_litr1c) * epv->rootlengthLandD_prop[layer];	
				cs->litr2c[layer]  += (cf->m_frootc_to_litr2c) * epv->rootlengthLandD_prop[layer];	
				cs->litr3c[layer]  += (cf->m_frootc_to_litr3c) * epv->rootlengthLandD_prop[layer];	
				cs->litr4c[layer]  += (cf->m_frootc_to_litr4c) * epv->rootlengthLandD_prop[layer];
				cs->cwdc[layer]    += (cf->m_livecrootc_to_cwdc + cf->m_deadcrootc_to_cwdc) * epv->rootlengthProp[layer];

				ns->litr1n[layer]  += (nf->m_frootn_to_litr1n) * epv->rootlengthLandD_prop[layer];	
				ns->litr2n[layer]  += (nf->m_frootn_to_litr2n) * epv->rootlengthLandD_prop[layer];	
				ns->litr3n[layer]  += (nf->m_frootn_to_litr3n) * epv->rootlengthLandD_prop[layer];	
				ns->litr4n[layer]  += (nf->m_frootn_to_litr4n) * epv->rootlengthLandD_prop[layer];
				ns->cwdn[layer]    += (nf->m_livecrootn_to_cwdn + nf->m_deadcrootn_to_cwdn) * epv->rootlengthProp[layer];	
			}
		}
	

	}

	
	/************************************************************/
	/* II. Fire mortality: these fluxes all enter atmospheric sinks */
	/************************************************************/

	if (epv->FM)
	{
		
		/* 1. Daily fluxes due to mortality */

		/* 1.1 non-woody pools */
		if (epc->leaf_cn)
		{
			cf->m_leafc_to_fire          = epv->FM * cs->leafc; 
			cf->m_leafc_storage_to_fire  = epv->FM * cs->leafc_storage;
			cf->m_leafc_transfer_to_fire = epv->FM * cs->leafc_transfer;
		
			nf->m_leafn_to_fire	          = cf->m_leafc_to_fire          / epc->leaf_cn;  
			nf->m_leafn_storage_to_fire	  = cf->m_leafc_storage_to_fire  / epc->leaf_cn; 
			nf->m_leafn_transfer_to_fire  = cf->m_leafc_transfer_to_fire / epc->leaf_cn;

		}
	
		if (epc->froot_cn)
		{
			cf->m_frootc_to_fire          = epv->FM * cs->frootc; 
			cf->m_frootc_storage_to_fire  = epv->FM * cs->frootc_storage;
			cf->m_frootc_transfer_to_fire = epv->FM * cs->frootc_transfer;
		
			nf->m_frootn_to_fire	       = cf->m_frootc_to_fire          / epc->froot_cn;  
			nf->m_frootn_storage_to_fire   = cf->m_frootc_storage_to_fire  / epc->froot_cn; 
			nf->m_frootn_transfer_to_fire  = cf->m_frootc_transfer_to_fire / epc->froot_cn; 

		}

		if (epc->yield_cn)
		{
			cf->m_yieldc_to_fire          = epv->FM * cs->yieldc; 
			cf->m_yieldc_storage_to_fire  = epv->FM * cs->yieldc_storage;
			cf->m_yieldc_transfer_to_fire = epv->FM * cs->yieldc_transfer;
		
			nf->m_yieldn_to_fire	      = cf->m_yieldc_to_fire          / epc->yield_cn;  
			nf->m_yieldn_storage_to_fire  = cf->m_yieldc_storage_to_fire  / epc->yield_cn; 
			nf->m_yieldn_transfer_to_fire = cf->m_yieldc_transfer_to_fire / epc->yield_cn; 

		}

		if (epc->softstem_cn)
		{
			cf->m_softstemc_to_fire          = epv->FM * cs->softstemc; 
			cf->m_softstemc_storage_to_fire  = epv->FM * cs->softstemc_storage;
			cf->m_softstemc_transfer_to_fire = epv->FM * cs->softstemc_transfer;
		
			nf->m_softstemn_to_fire	          = cf->m_softstemc_to_fire          / epc->softstem_cn;  
			nf->m_softstemn_storage_to_fire	  = cf->m_softstemc_storage_to_fire  / epc->softstem_cn; 
			nf->m_softstemn_transfer_to_fire  = cf->m_softstemc_transfer_to_fire / epc->softstem_cn; 

		}

		/* 1.2 woody pools */
		if (epc->livewood_cn)
		{
			cf->m_livestemc_to_fire           = epv->FM * cs->livestemc;	
			cf->m_livestemc_storage_to_fire   = epv->FM * cs->livestemc_storage;	
			cf->m_livestemc_transfer_to_fire  = epv->FM * cs->livestemc_transfer;

			cf->m_livecrootc_to_fire          = epv->FM * cs->livecrootc;   
			cf->m_livecrootc_storage_to_fire  = epv->FM * cs->livecrootc_storage;	
			cf->m_livecrootc_transfer_to_fire = epv->FM * cs->livecrootc_transfer;	 

		    nf->m_livestemn_to_fire           = cf->m_livestemc_to_fire           / epc->livewood_cn;	
			nf->m_livestemn_storage_to_fire   = cf->m_livestemc_storage_to_fire   / epc->livewood_cn;		
			nf->m_livestemn_transfer_to_fire  = cf->m_livestemc_transfer_to_fire  / epc->livewood_cn;		
			
			nf->m_livecrootn_to_fire          = cf->m_livecrootc_to_fire          / epc->livewood_cn;	   
			nf->m_livecrootn_storage_to_fire  = cf->m_livecrootc_storage_to_fire  / epc->livewood_cn;		
			nf->m_livecrootn_transfer_to_fire = cf->m_livecrootc_transfer_to_fire / epc->livewood_cn;	 

			/*nf->m_livestemn_to_fire             = epv->FM * ns->livestemn;	
			nf->m_livestemn_storage_to_fire     = epv->FM * ns->livestemn_storage;
			nf->m_livestemn_transfer_to_fire    = epv->FM * ns->livestemn_transfer;
			
			nf->m_livecrootn_to_fire            = epv->FM * ns->livecrootn;
			nf->m_livecrootn_transfer_to_fire   = epv->FM * ns->livecrootn_transfer;
			nf->m_livecrootn_storage_to_fire    = epv->FM * ns->livecrootn_storage;*/	
			 
		}

		if (epc->deadwood_cn)
		{
			cf->m_deadstemc_to_fire           = epv->FM * cs->deadstemc * epc->dscp;	
			cf->m_deadstemc_storage_to_fire   = epv->FM * cs->deadstemc_storage;	
			cf->m_deadstemc_transfer_to_fire  = epv->FM * cs->deadstemc_transfer;	 
			cf->m_deadcrootc_to_fire          = epv->FM * cs->deadcrootc * epc->dscp;   
			cf->m_deadcrootc_storage_to_fire  = epv->FM * cs->deadcrootc_storage;	
			cf->m_deadcrootc_transfer_to_fire = epv->FM * cs->deadcrootc_transfer;	 

		 nf->m_deadstemn_to_fire           = cf->m_deadstemc_to_fire           / epc->deadwood_cn;	
			nf->m_deadstemn_storage_to_fire   = cf->m_deadstemc_storage_to_fire   / epc->deadwood_cn;		
			nf->m_deadstemn_transfer_to_fire  = cf->m_deadstemc_transfer_to_fire  / epc->deadwood_cn;		 
			nf->m_deadcrootn_to_fire          = cf->m_deadcrootc_to_fire          / epc->deadwood_cn;	   
			nf->m_deadcrootn_storage_to_fire  = cf->m_deadcrootc_storage_to_fire  / epc->deadwood_cn;		
			nf->m_deadcrootn_transfer_to_fire = cf->m_deadcrootc_transfer_to_fire / epc->deadwood_cn;	

			
			/* nf->m_deadstemn_to_fire           = epc->dscp * epv->FM * ns->deadstemn;
			nf->m_deadstemn_storage_to_fire   = epv->FM * ns->deadstemn_storage;
			nf->m_deadstemn_transfer_to_fire  = epv->FM * ns->deadstemn_transfer;
			
			nf->m_deadcrootn_to_fire          = epc->dscp * epv->FM * ns->deadcrootn;
			nf->m_deadcrootn_transfer_to_fire = epv->FM * ns->deadcrootn_transfer;
			nf->m_deadcrootn_storage_to_fire  = epv->FM * ns->deadcrootn_storage;*/

		}
	
	
		/* 1.3 transfer pools */
		cf->m_gresp_storage_to_fire     = epv->FM * cs->gresp_storage;
		cf->m_gresp_transfer_to_fire    = epv->FM * cs->gresp_transfer;
		
		nf->m_retransn_to_fire			= epv->FM * ns->retransn;


		/* 1.4 standing and cut-down dead biomass simulation */
		cf->m_STDBc_to_fire = epv->FM * (cs->STDBc_leaf + cs->STDBc_froot + cs->STDBc_yield + cs->STDBc_softstem);
		cf->m_CTDBc_to_fire = epv->FM * (cs->CTDBc_leaf + cs->CTDBc_froot + cs->CTDBc_yield + cs->CTDBc_softstem + cs->CTDBc_cstem + cs->CTDBc_croot);
	
		nf->m_STDBn_to_fire = epv->FM * (ns->STDBn_leaf + ns->STDBn_froot + ns->STDBn_yield + ns->STDBn_softstem);
		nf->m_CTDBn_to_fire = epv->FM * (ns->CTDBn_leaf + ns->CTDBn_froot + ns->CTDBn_yield + ns->CTDBn_softstem + ns->CTDBn_cstem + ns->CTDBn_croot);
	
	
	
		/* 1.5 litter and CWD fire fluxes: fire affects only the rootzone layers */

		cf->m_litr1c_to_fire_total=cf->m_litr2c_to_fire_total=cf->m_litr3c_to_fire_total=cf->m_litr4c_to_fire_total=0;
		nf->m_litr1n_to_fire_total=nf->m_litr2n_to_fire_total=nf->m_litr3n_to_fire_total=nf->m_litr4n_to_fire_total=nf->m_cwdn_to_fire_total=0;
		
		for (layer = 0; layer < epv->n_rootlayers; layer++)
		{
			cf->m_litr1c_to_fire[layer] = epv->FM   * cs->litr1c[layer];
			cf->m_litr2c_to_fire[layer] = epv->FM   * cs->litr2c[layer];
			cf->m_litr3c_to_fire[layer] = epv->FM   * cs->litr3c[layer];
			cf->m_litr4c_to_fire[layer] = epv->FM   * cs->litr4c[layer];
			cf->m_cwdc_to_fire[layer]   = epc->cwcp * epv->FM * cs->cwdc[layer];

			nf->m_litr1n_to_fire[layer] =  epv->FM   * ns->litr1n[layer];
			nf->m_litr2n_to_fire[layer] =  epv->FM   * ns->litr2n[layer];
			nf->m_litr3n_to_fire[layer] =  epv->FM   * ns->litr3n[layer];
			nf->m_litr4n_to_fire[layer] =  epv->FM   * ns->litr4n[layer];
			nf->m_cwdn_to_fire[layer]   =  epc->cwcp * epv->FM * ns->cwdn[layer];

			
			flux_from_carbon = cf->m_litr1c_to_fire[layer] * (cs->litr1c[layer] / ns->litr1n[layer]);
			if (nf->m_litr1n_to_fire[layer] > flux_from_carbon) nf->m_litr1n_to_fire[layer] = flux_from_carbon;
			flux_from_carbon = cf->m_litr2c_to_fire[layer] * (cs->litr2c[layer] / ns->litr2n[layer]);
			if (nf->m_litr2n_to_fire[layer] > flux_from_carbon) nf->m_litr2n_to_fire[layer] = flux_from_carbon;
			flux_from_carbon = cf->m_litr3c_to_fire[layer] * (cs->litr3c[layer] / ns->litr3n[layer]);
			if (nf->m_litr3n_to_fire[layer] > flux_from_carbon) nf->m_litr3n_to_fire[layer] = flux_from_carbon;
			flux_from_carbon = cf->m_litr4c_to_fire[layer] * (cs->litr4c[layer] / ns->litr4n[layer]);
			if (nf->m_litr4n_to_fire[layer] > flux_from_carbon) nf->m_litr4n_to_fire[layer] = flux_from_carbon;
			flux_from_carbon = cf->m_cwdc_to_fire[layer] * (cs->cwdc[layer] / ns->cwdn[layer]);
			if (nf->m_cwdn_to_fire[layer] > flux_from_carbon) nf->m_cwdn_to_fire[layer] = flux_from_carbon;

			cf->m_litr1c_to_fire_total  += cf->m_litr1c_to_fire[layer];
			cf->m_litr2c_to_fire_total  += cf->m_litr2c_to_fire[layer];
			cf->m_litr3c_to_fire_total  += cf->m_litr3c_to_fire[layer];
			cf->m_litr4c_to_fire_total  += cf->m_litr4c_to_fire[layer];
			cf->m_cwdc_to_fire_total    += cf->m_cwdc_to_fire[layer];
	
			nf->m_litr1n_to_fire_total  += nf->m_litr1n_to_fire[layer];
			nf->m_litr2n_to_fire_total  += nf->m_litr2n_to_fire[layer];
			nf->m_litr3n_to_fire_total  += nf->m_litr3n_to_fire[layer];
			nf->m_litr4n_to_fire_total  += nf->m_litr4n_to_fire[layer];
			nf->m_cwdn_to_fire_total    += nf->m_cwdn_to_fire[layer];


			
		}

	
		/* 2. update state variables for fire fluxes */
		/* this is the only place other than state_update.c routines wherestate variables get changed.  Mortality is taken care of last and  given special treatment for state update so that it doesn't interfere
		with the other fluxes that are based on proportions of state variables, especially the phenological fluxes */
	
		/* 2.1 ABOVEGROUND variables */	
	
		cs->leafc      -= cf->m_leafc_to_fire;
		cs->yieldc     -= cf->m_yieldc_to_fire;
		cs->softstemc  -= cf->m_softstemc_to_fire;
		cs->livestemc  -= cf->m_livestemc_to_fire;
		cs->deadstemc  -= cf->m_deadstemc_to_fire;

		cs->FIREsnk_C       += cf->m_leafc_to_fire + cf->m_yieldc_to_fire + cf->m_softstemc_to_fire + cf->m_livestemc_to_fire + cf->m_deadstemc_to_fire;

		ns->leafn      -= nf->m_leafn_to_fire;
		ns->yieldn     -= nf->m_yieldn_to_fire;
		ns->softstemn  -= nf->m_softstemn_to_fire;
		ns->livestemn  -= nf->m_livestemn_to_fire;
		ns->deadstemn  -= nf->m_deadstemn_to_fire;

		ns->FIREsnk_N       += nf->m_leafn_to_fire + nf->m_yieldn_to_fire + nf->m_softstemn_to_fire + nf->m_livestemn_to_fire + nf->m_deadstemn_to_fire;
	

		/* 2.2 BELOWGROUND variables */	
		
		cs->frootc     -= cf->m_frootc_to_fire;
		cs->livecrootc -= cf->m_livecrootc_to_fire;
		cs->deadcrootc -= cf->m_deadcrootc_to_fire;

		cs->leafc_storage       -= cf->m_leafc_storage_to_fire;
		cs->frootc_storage      -= cf->m_frootc_storage_to_fire;
		cs->yieldc_storage      -= cf->m_yieldc_storage_to_fire;
		cs->softstemc_storage   -= cf->m_softstemc_storage_to_fire;
		cs->livestemc_storage   -= cf->m_livestemc_storage_to_fire;
		cs->deadstemc_storage   -= cf->m_deadstemc_storage_to_fire;
		cs->livecrootc_storage  -= cf->m_livecrootc_storage_to_fire;
		cs->deadcrootc_storage  -= cf->m_deadcrootc_storage_to_fire;

		cs->leafc_transfer      -= cf->m_leafc_transfer_to_fire;
		cs->frootc_transfer     -= cf->m_frootc_transfer_to_fire;
		cs->yieldc_transfer     -= cf->m_yieldc_transfer_to_fire;
		cs->softstemc_transfer  -= cf->m_softstemc_transfer_to_fire;
		cs->livestemc_transfer  -= cf->m_livestemc_transfer_to_fire;
		cs->deadstemc_transfer  -= cf->m_deadstemc_transfer_to_fire;
		cs->livecrootc_transfer -= cf->m_livecrootc_transfer_to_fire;
		cs->deadcrootc_transfer -= cf->m_deadcrootc_transfer_to_fire;
		cs->gresp_storage       -= cf->m_gresp_storage_to_fire;
		cs->gresp_transfer      -= cf->m_gresp_transfer_to_fire;
	
		cs->FIREsnk_C   += cf->m_frootc_to_fire + cf->m_livecrootc_to_fire + cf->m_deadcrootc_to_fire + 
						  cf->m_leafc_storage_to_fire + cf->m_frootc_storage_to_fire + cf->m_yieldc_storage_to_fire + cf->m_softstemc_storage_to_fire + 
						  cf->m_livestemc_storage_to_fire + cf->m_deadstemc_storage_to_fire + cf->m_livecrootc_storage_to_fire + cf->m_deadcrootc_storage_to_fire + 
						  cf->m_leafc_transfer_to_fire + cf->m_frootc_transfer_to_fire + cf->m_yieldc_transfer_to_fire + cf->m_softstemc_transfer_to_fire + 
						  cf->m_livestemc_transfer_to_fire + cf->m_deadstemc_transfer_to_fire + cf->m_livecrootc_transfer_to_fire + cf->m_deadcrootc_transfer_to_fire +
						  cf->m_gresp_storage_to_fire + cf->m_gresp_transfer_to_fire;

		ns->frootn     -= nf->m_frootn_to_fire;
		ns->livecrootn -= nf->m_livecrootn_to_fire;
		ns->deadcrootn -= nf->m_deadcrootn_to_fire;

		ns->leafn_storage       -= nf->m_leafn_storage_to_fire;
		ns->frootn_storage      -= nf->m_frootn_storage_to_fire;
		ns->yieldn_storage      -= nf->m_yieldn_storage_to_fire;
		ns->softstemn_storage   -= nf->m_softstemn_storage_to_fire;
		ns->livestemn_storage   -= nf->m_livestemn_storage_to_fire;
		ns->deadstemn_storage   -= nf->m_deadstemn_storage_to_fire;
		ns->livecrootn_storage  -= nf->m_livecrootn_storage_to_fire;
		ns->deadcrootn_storage  -= nf->m_deadcrootn_storage_to_fire;

		ns->leafn_transfer      -= nf->m_leafn_transfer_to_fire;
		ns->frootn_transfer     -= nf->m_frootn_transfer_to_fire;
		ns->yieldn_transfer     -= nf->m_yieldn_transfer_to_fire;
		ns->softstemn_transfer  -= nf->m_softstemn_transfer_to_fire;
		ns->livestemn_transfer  -= nf->m_livestemn_transfer_to_fire;
		ns->deadstemn_transfer  -= nf->m_deadstemn_transfer_to_fire;
		ns->livecrootn_transfer -= nf->m_livecrootn_transfer_to_fire;
		ns->deadcrootn_transfer -= nf->m_deadcrootn_transfer_to_fire;
		ns->retransn             -= nf->m_retransn_to_fire;
	
		ns->FIREsnk_N   += nf->m_frootn_to_fire + nf->m_livecrootn_to_fire + nf->m_deadcrootn_to_fire + 
						  nf->m_leafn_storage_to_fire + nf->m_frootn_storage_to_fire + nf->m_yieldn_storage_to_fire + nf->m_softstemn_storage_to_fire + 
						  nf->m_livestemn_storage_to_fire + nf->m_deadstemn_storage_to_fire + nf->m_livecrootn_storage_to_fire + nf->m_deadcrootn_storage_to_fire + 
						  nf->m_leafn_transfer_to_fire + nf->m_frootn_transfer_to_fire + nf->m_yieldn_transfer_to_fire + nf->m_softstemn_transfer_to_fire + 
						  nf->m_livestemn_transfer_to_fire + nf->m_deadstemn_transfer_to_fire + nf->m_livecrootn_transfer_to_fire + nf->m_deadcrootn_transfer_to_fire +
						  nf->m_retransn_to_fire;

		/* 2.3. standing and cut-down dead biomass */
		cs->STDBc_leaf     -= epv->FM * cs->STDBc_leaf;
		cs->STDBc_froot    -= epv->FM * cs->STDBc_froot;
		cs->STDBc_yield    -= epv->FM * cs->STDBc_yield;
		cs->STDBc_softstem -= epv->FM * cs->STDBc_softstem;
	
		cs->FIREsnk_C     += cf->m_STDBc_to_fire;

		ns->STDBn_leaf     -= epv->FM * ns->STDBn_leaf;
		ns->STDBn_froot    -= epv->FM * ns->STDBn_froot;
		ns->STDBn_yield    -= epv->FM * ns->STDBn_yield;
		ns->STDBn_softstem -= epv->FM * ns->STDBn_softstem;

		ns->FIREsnk_N     += nf->m_STDBn_to_fire;

		cs->CTDBc_leaf     -= epv->FM * cs->CTDBc_leaf;
		cs->CTDBc_froot    -= epv->FM * cs->CTDBc_froot;
		cs->CTDBc_yield    -= epv->FM * cs->CTDBc_yield;
		cs->CTDBc_softstem -= epv->FM * cs->CTDBc_softstem;
		cs->CTDBc_cstem    -= epv->FM * cs->CTDBc_cstem;
		cs->CTDBc_croot    -= epv->FM * cs->CTDBc_croot;
	
		cs->FIREsnk_C     += cf->m_CTDBc_to_fire;

		ns->CTDBn_leaf     -= epv->FM * ns->CTDBn_leaf;
		ns->CTDBn_froot    -= epv->FM * ns->CTDBn_froot;
		ns->CTDBn_yield    -= epv->FM * ns->CTDBn_yield;
		ns->CTDBn_softstem -= epv->FM * ns->CTDBn_softstem;
		ns->CTDBn_cstem    -= epv->FM * ns->CTDBn_cstem;
		ns->CTDBn_croot    -= epv->FM * ns->CTDBn_croot;

		ns->FIREsnk_N     += nf->m_CTDBn_to_fire;

		/* 2.4 litter and CWD carbon state updates: fire affects only rootzone layer */

		for (layer = 0; layer < epv->n_rootlayers; layer++)
		{
			cs->litr1c[layer]   -= cf->m_litr1c_to_fire[layer];
			cs->litr2c[layer]   -= cf->m_litr2c_to_fire[layer];
			cs->litr3c[layer]   -= cf->m_litr3c_to_fire[layer];
			cs->litr4c[layer]   -= cf->m_litr4c_to_fire[layer];
			cs->cwdc[layer] 	-= cf->m_cwdc_to_fire[layer] ;

			cs->FIREsnk_C += cf->m_litr1c_to_fire[layer]  + cf->m_litr2c_to_fire[layer]  + cf->m_litr3c_to_fire[layer]  + cf->m_litr4c_to_fire[layer]  + cf->m_cwdc_to_fire[layer] ;

			ns->litr1n[layer]    -= nf->m_litr1n_to_fire[layer];
			ns->litr2n[layer]    -= nf->m_litr2n_to_fire[layer];
			ns->litr3n[layer]    -= nf->m_litr3n_to_fire[layer];
			ns->litr4n[layer]    -= nf->m_litr4n_to_fire[layer];
			ns->cwdn[layer]      -= nf->m_cwdn_to_fire[layer];

			ns->FIREsnk_N += nf->m_litr1n_to_fire[layer]  + nf->m_litr2n_to_fire[layer]  + nf->m_litr3n_to_fire[layer]  + nf->m_litr4n_to_fire[layer]  + nf->m_cwdn_to_fire[layer] ;
	
		}

	}
	
	return (errorCode);
}

