/* 
annVARinit.c
initalization of annual, cumulative variables on first day of every simulation year

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2022, D. Hidy [dori.hidy@gmail.com]
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


#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

int annVARinit(summary_struct* summary, epvar_struct* epv, cstate_struct* cs, wstate_struct* ws)
{

	int errorCode=0;
	
	
	summary->annprcp        = 0;
	summary->anntavg        = 0;
	summary->cumRunoff      = 0;
	summary->cumRunoffH     = 0;
	summary->cumRunoffD     = 0;
	summary->cumWleach_RZmax    = 0;
	summary->cumNleach_RZmax    = 0;
	summary->Wleach_RZmax       = 0;
	summary->DOCleach_RZmax     = 0;
	summary->sminNleach_RZmax   = 0;
	summary->DONleach_RZmax     = 0;
	summary->cumNPP         = 0;
	summary->cumNPPabove_w  = 0;
	summary->cumNPPbelow_w  = 0;
	summary->cumNPPabove_nw = 0;
	summary->cumNPPbelow_nw = 0;
	summary->cumNEP  = 0;
	summary->cumNEE  = 0;
	summary->cumGPP  = 0;
	summary->cumNGB  = 0;
	summary->cumNBP  = 0;
	summary->cumMR  = 0;
	summary->cumGR  = 0;
	summary->cumHR  = 0;
	summary->cumAR = 0.0;
	summary->cumTR  = 0;
	summary->cumSR  = 0;
	summary->cumN2Oflux  = 0;
	summary->cumN2OfluxCeq = 0;
	summary->cumCH4fluxCeq = 0;
	summary->cumCflux_lateral = 0;
	summary->cumCH4flux  = 0;
	summary->cumCloss_MGM  = 0;
	summary->cumCplus_MGM  = 0;
	summary->cumCloss_THN_w = 0;
	summary->cumCloss_THN_nw = 0;
	summary->cumCloss_MOW  = 0;

	summary->cumCloss_GRZ  = 0;
	summary->cumCplus_GRZ  = 0;
	summary->cumCplus_FRZ  = 0;
		
	summary->cumCplus_PLT  = 0;
	summary->cumCloss_PLT  = 0;
	summary->cumCloss_HRV  = 0;
	summary->cumYieldC_HRV = 0;
	summary->cumCloss_PLG  = 0;

	summary->cumNplus_GRZ  = 0;
	summary->cumNplus_FRZ  = 0;
	summary->cumNplus_PLT  = 0;
	summary->cumNloss_PLT = 0;
	summary->cumNplus_FRZ_org = 0;
	summary->cumNplus_FRZ_NH4 = 0;
	summary->cumNplus_FRZ_NO3 = 0;
	summary->cumCloss_SNSC  = 0;
	summary->cumCplus_STDB  = 0;
	summary->cumCplus_CTDB  = 0;
	summary->cumEVPsurface  = 0;
	summary->cumETcanopy  = 0;
	summary->cumET  = 0;
	
	summary->cumMRdeficit_NSC  = 0;
	summary->cumMRdeficit_SC   = 0;

	cs->yieldC_HRV        = 0;
	cs->vegCabove_HRV     = 0;
	cs->frootC_HRV        = 0;
	cs->MRlimitSUM        = 0;

	epv->cumWS_anoxic = 0;
	epv->cumWS_drought = 0;
	epv->cumWS = 0;
	epv->cumNS = 0;
	epv->cumCalloc_plant = 0;
	epv->cumNalloc_plant = 0;

	ws->EVPsurface1cum       = 0;
	ws->EVPsurface2cum       = 0;



	epv->cumWS_anoxic     = 0;
	epv->cumWS_drought    = 0;
	epv->cumWS            = 0;
	epv->cumNS            = 0;
	epv->cumCalloc_plant  = 0;
	epv->cumNalloc_plant  = 0;


	epv->annmax_leafc = 0;
	epv->annmax_frootc = 0;
	epv->annmax_yieldc = 0;
	epv->annmax_softstemc = 0;
	epv->annmax_livestemc = 0;
	epv->annmax_livecrootc = 0;


	summary->cumUNSATchangeGW_orgN= 0;
	summary->cumUNSATplantUPto_orgN = 0;
	summary->cumUNSATecofunc_orgN = 0;
	summary->cumUNSATfertil_orgN = 0;
	summary->cumUNSATdischarge_orgN = 0;
	summary->cumUNSATrecharge_orgN = 0;

	summary->cumTOTALchangeGW_orgN = 0;
	summary->cumTOTALplantUPto_orgN = 0;
	summary->cumTOTALecofunc_orgN = 0;
	summary->cumTOTALfertil_orgN = 0;
	summary->cumTOTALdischarge_orgN = 0;
	summary->cumTOTALrecharge_orgN = 0;

	summary->cumUNSATchangeGW_NH4 = 0;
	summary->cumUNSATchangeGW_NO3 = 0;
	summary->cumUNSATecofunc_NH4= 0;
	summary->cumUNSATecofunc_NO3= 0;
	summary->cumUNSATfertil_NH4= 0;
	summary->cumUNSATfertil_NO3= 0;
	summary->cumUNSATdischarge_NH4= 0;
	summary->cumUNSATdischarge_NO3= 0;
	summary->cumUNSATrecharge_NH4= 0;
	summary->cumUNSATrecharge_NO3= 0;
	summary->cumUNSATdenitr = 0;
	summary->cumUNSATvolat = 0;

	summary->cumHRV_to_transpN = 0;
	summary->cumENVtoSMINN = 0;
	summary->cumFRZtoN = 0;


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
	summary->cumPET = 0;
	summary->cumPETcanopy = 0;
	summary->cumPETsurface = 0;
	summary->cumFLsoilw = 0;


	return (errorCode);

}

