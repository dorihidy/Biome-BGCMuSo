/* 
zero_srcsnk.c
fill the source and sink variables with 0.0 at the start of the simulation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2025, D. Hidy [dori.hidy@gmail.com]
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
#include "bgc_struct.h"     /* structure definitions */
#include "bgc_func.h"       /* function prototypes */
#include "bgc_constants.h"

/* zero the source and sink state variables */
int zero_srcsnk(cstate_struct* cs, nstate_struct* ns, wstate_struct* ws, summary_struct* summary)
{
	int errorCode=0;

	/* zero the water sources and sinks  */
	ws->prcp_src = 0.0;
	ws->soilEVP_snk = 0.0;
	ws->snowSUBL_snk = 0.0;
	ws->EVPcanopyw_snk = 0.0;
	ws->TRP_snk = 0.0;
	ws->prcp_src = 0.0;
	ws->soilEVP_snk = 0.0;
	ws->canopyw_THNsnk = 0.0;
	ws->canopyw_MOWsnk = 0.0;
	ws->canopyw_HRVsnk = 0.0;	
	ws->canopyw_PLGsnk = 0.0;
	ws->canopyw_GRZsnk = 0.0;	
	ws->IRGsrc_W = 0.0;
	ws->condIRGsrc = 0;
	ws->FRZsrc_W = 0.0;
	ws->runoff_snk = 0.0;
	ws->pondEVP_snk = 0;
	ws->deeppercolation_snk = 0.0;  
	ws->GWsrc_W = 0.0;
	ws->GWsnk_W = 0.0;
	ws->FLsrc_W = 0;
	ws->WbalanceERR = 0.0;	
	ws->inW = 0.0;
	ws->outW = 0.0;
	ws->storeW = 0.0;
	ws->EVPsurface1cum = 0.0;
	ws->EVPsurface2cum = 0.0;
	ws->GW_waterlogging = 0.0;
	ws->Wprec_snk = 0.0;

	/* zero the carbon sources and sinks */
	cs->psnsun_src = 0.0;
	cs->psnshade_src = 0.0;
	cs->MRleaf_snk = 0.0;
	cs->GRleaf_snk = 0.0;
	cs->MRfroot_snk = 0.0;
	cs->GRfroot_snk = 0.0;
	cs->MRlivestem_snk = 0.0;
	cs->GRlivestem_snk = 0.0;
	cs->GRdeadstem_snk = 0.0;
	cs->MRlivecroot_snk = 0.0;
	cs->GRlivecroot_snk = 0.0;
	cs->GRdeadcroot_snk = 0.0;
	cs->HRlitr1_snk = 0.0;
	cs->HRlitr2_snk = 0.0;
	cs->HRlitr4_snk = 0.0;
	cs->HRsoil1_snk = 0.0;
	cs->HRsoil2_snk = 0.0;
	cs->HRsoil3_snk = 0.0;
	cs->HRsoil4_snk = 0.0;
	cs->Cdeepleach_snk = 0;
	cs->Cdeepleach_snk = 0;
	cs->flowHSsnk_C = 0;
	cs->FIREsnk_C = 0.0;

	cs->SNSCsnk_C = 0.0;
	cs->PLTsrc_C = 0.0; 
	cs->MULsrc_C = 0.0;
	cs->CWEsnk_C = 0.0;
	cs->THN_transportC = 0.0; 

	cs->MOW_snkC = 0;
	cs->GRZsnk_C = 0.0;  
	cs->GRZsrc_C = 0.0;
	cs->HRV_snkC = 0.0;
	cs->FRZsrc_C = 0.0;
	cs->yieldC_HRV = 0.0;
	cs->frootC_HRV = 0.0;
	cs->vegCabove_HRV = 0.0;
	cs->MRyield_snk = 0.0;
	cs->GRyield_snk = 0.0;
	cs->MRsoftstem_snk = 0.0;
	cs->GRsoftstem_snk = 0.0;
	cs->CbalanceERR = 0.0;
	cs->CNratioERR = 0.0;
	cs->inC = 0.0;
	cs->outC = 0.0;
	cs->storeC = 0.0;
	cs->NSCnw = 0.0;
	cs->NSCw = 0.0;
	cs->SCnw = 0.0;
	cs->SCw = 0.0;
	cs->GWsrc_C = 0;
	cs->GWsnk_C = 0;
	cs->FLsrc_C = 0;
	cs->litrCabove_total = 0;
	cs->litrCbelow_total = 0;
	cs->cwdCabove_total = 0;
	cs->cwdCbelow_total = 0;
	cs->Cprec_snk = 0;

	ns->Nfix_src = 0.0;
	ns->Ndep_src = 0.0;
	ns->Ndeepleach_snk = 0.0;
	ns->Nvol_snk = 0.0;
	ns->Nprec_snk = 0;
	ns->FIREsnk_N = 0.0;
	ns->SPINUPsrc = 0.0;
	ns->SNSCsnk_N = 0.0;
	ns->PLTsrc_N = 0.0; 
	ns->MULsrc_N = 0.0;
	ns->CWEsnk_N = 0.0;
	ns->THNsnk_N = 0.0; 
	ns->MOWsnk_N = 0;
	ns->GRZsnk_N = 0.0;  
	ns->GRZsrc_N = 0.0;
	ns->HRVsnk_N = 0.0; 
	ns->FRZsrc_N = 0.0; 
	ns->NbalanceERR = 0.0;
	ns->inN = 0.0;
	ns->outN = 0.0;
	ns->storeN = 0.0;
	ns->FLsrc_N = 0;
	ns->GWsrc_N = 0;
	ns->GWsnk_N = 0;

	/* zero the summary variables */
	summary->annprcp = 0.0;
	summary->anntavg = 0.0;
	summary->cumRunoff = 0;
	summary->cumRunoffH = 0;
	summary->cumRunoffD = 0;
	summary->cumWleach_RZmax = 0;
	summary->cumNleach_RZmax = 0;
	summary->cumNPP = 0;
	summary->cumNPPabove_w = 0;
	summary->cumNPPbelow_w = 0;
	summary->cumNPPabove_nw = 0;
	summary->cumNPPbelow_nw = 0;
	summary->cumNEP = 0;
	summary->cumNEE = 0;
	summary->cumGPP = 0;
	summary->cumNBP = 0;
	summary->cumNGB = 0;
	summary->cumMR = 0;
	summary->cumGR = 0;
	summary->cumHR = 0;
	summary->cumAR = 0;
	summary->cumTR = 0;
	summary->cumSR = 0;
	summary->cumN2Oflux = 0;
	summary->cumN2OfluxCeq = 0;
	summary->cumCH4flux = 0;
	summary->cumCH4fluxCeq = 0;
	summary->cumCloss_MGM = 0;
	summary->cumCplus_MGM = 0;
	summary->cumCloss_THN_w = 0;
	summary->cumCloss_THN_nw = 0;
	summary->cumCloss_MOW = 0;
	summary->cumCloss_HRV = 0;
	summary->cumYieldC_HRV = 0;
	summary->cumCloss_PLG = 0;
	summary->cumCloss_GRZ = 0;
	summary->cumCplus_GRZ = 0;
	summary->cumCplus_FRZ = 0;
	summary->cumCplus_PLT = 0;
	summary->cumCloss_PLT = 0;
	summary->cumNplus_PLT = 0;
	summary->cumNloss_PLT = 0;
	summary->cumNplus_GRZ = 0;
	summary->cumNplus_FRZ = 0;
	summary->cumNplus_FRZ_org = 0;
	summary->cumNplus_FRZ_NH4 = 0;
	summary->cumNplus_FRZ_NO3 = 0;
	summary->cumCloss_SNSC = 0;
	summary->cumCplus_STDB = 0;
	summary->cumCplus_CTDB = 0;
	summary->cumEVPsurface = 0;
	summary->cumETcanopy = 0;
	summary->cumET = 0;
	summary->cumPET = 0;
	summary->cumPETsurface = 0;
	summary->cumPETcanopy = 0;

	summary->cumTOTALchangeGW_orgN = 0;
	summary->cumTOTALchangeGW_NH4 = 0;
	summary->cumTOTALchangeGW_NO3 = 0;
	summary->cumTOTALplantUPto_orgN = 0;
	summary->cumTOTALecofunc_orgN = 0;
	summary->cumTOTALfertil_orgN = 0;
	summary->cumTOTALdischarge_orgN = 0;
	summary->cumTOTALrecharge_orgN = 0;

	summary->cumUNSATplantUPto_orgN = 0;
	summary->cumUNSATecofunc_orgN = 0;
	summary->cumUNSATfertil_orgN = 0;
	summary->cumUNSATdischarge_orgN = 0;
	summary->cumUNSATrecharge_orgN = 0;
	summary->cumUNSATecofunc_NH4 = 0;
	summary->cumUNSATecofunc_NO3 = 0;


	summary->cumUNSATfertil_NH4 = 0;
	summary->cumUNSATfertil_NO3 = 0;


	summary->cumUNSATdischarge_NH4 = 0;
	summary->cumUNSATdischarge_NO3 = 0;


	summary->cumUNSATrecharge_NH4 = 0;
	summary->cumUNSATrecharge_NO3 = 0;



	summary->cumMRdeficit_NSC = 0;
	summary->cumMRdeficit_SC = 0;
	summary->cumHRV_to_transpN = 0;
	summary->cumENVtoSMINN = 0;
	summary->cumFRZtoN = 0;
	summary->cumUNSATvolat = 0;
	summary->cumWflux_fromPRCP = 0;
	summary->cumWflux_fromFRZ = 0;
	summary->cumWflux_fromIRG = 0;
	summary->cumGWdischarge = 0;
	summary->cumGWrecharge = 0;
	summary->cumEVPfromGW = 0;
	summary->cumTRPfromGW = 0;
	summary->cumGWsrc = 0;
	summary->cumGWsnk = 0;
	summary->cumWinput = 0;
	summary->cumWoutput = 0;

	summary->cumCflux_lateral = 0;
	summary->harvestIndex = 0;
	summary->rootIndex = 0;
	summary->belowground_ratio = 0;
	summary->annmax_livingBIOMabove = 0.0;
	summary->annmax_livingBIOMbelow = 0.0;
	summary->annmax_BIOMaboveSUM = 0.0;
	summary->annmax_BIOMbelowSUM = 0.0;

	return (errorCode);
}
