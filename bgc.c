/*
bgc.c
Core BGC model logic

Includes in-line output handling routines that write to daily and annual
output files. This is the only library module that has external
I/O connections, and so it is the only module that includes bgc_io.h.

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2022, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Modified:
4/17/2000 (PET): Added annual average temperature and annual total precipitation
to the simple annual text output file.
Modified:
13/07/2000: Added input of Ndep from file. Changes are made by Galina Churkina. 
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>

#include "ini.h"    
#include "bgc_struct.h"     
#include "bgc_constants.h"
#include "bgc_func.h"
#include "pointbgc_struct.h"   /* data structures for point driver */
#include "bgc_io.h"
#include "pointbgc_func.h"


int bgc(bgcin_struct* bgcin, bgcout_struct* bgcout)
{
	
	/* variable declarations */
	int errorCode=0;
	int nyears, leap;

	/* iofiles and program control variables */
	control_struct     ctrl;

	/* meteorological variables */
	metarr_struct			metarr;
	metvar_struct			metv;
	co2control_struct		co2;
	NdepControl_struct		ndep;

	
	/* state and flux variables for water, carbon, and nitrogen */
	wstate_struct      ws;
	wflux_struct       wf;
	cinit_struct       cinit;
	cstate_struct      cs;
	cflux_struct       cf;
	nstate_struct      ns;
	nflux_struct       nf;


	/* primary ecophysiological variables */
	epvar_struct       epv;

	/* management variables */
	planting_struct     PLT;
	thinning_struct     THN;
	mowing_struct       MOW;
	grazing_struct      GRZ;
	harvesting_struct   HRV;
	ploughing_struct    PLG;
	fertilizing_struct  FRZ;
	irrigating_struct   IRG;
	mulching_struct     MUL;
	CWDextract_struct	CWE;
	flooding_struct     FLS;

	/* groundwater variables */
	groundwaterINIT_struct  GWS;

	/* site physical constants */
	siteconst_struct   sitec;

	/* soil proportion variables */
	soilprop_struct   sprop;

	/* ecophysiological constants */
	epconst_struct     epc;

	/* soil information struct */
	soilInfo_struct     soilInfo;


	/* phenological data */
	phenarray_struct   phenarr;
    phenology_struct   phen;


	/* photosynthesis structures */
	psn_struct         psn_sun, psn_shade;
	
	/* temporary nitrogen variables for decomposition and allocation */
	ntemp_struct       nt;
	
	/* summary variable structure */
	summary_struct     summary;

	/* miscelaneous variables for program control in main */
	int simyr = 0;
	int yday  = 0;
	int i     = 0;
	int first_balance;
	int annual_alloc;

	double dailyNdep; 

	double CbalanceERR = -100;
	double NbalanceERR = -100;
	double WbalanceERR = -100;
	double CNratioERR = -100;

	char* soiltype_name[] = { "sand", "loamySand", "sandyLoam", "loam", "siltLoam", "silt", "sandyClayLoam", "clayLoam", "siltyClayLoam", "sandyClay", "siltyClay", "clay" };

	/* local storage for daily and annual output variables and output mapping (array of pointers to double)  */
	double* dayarr=0;
	double* monavgarr = 0;
	double* annavgarr=0;
	double* annarr=0;
	double** output_map=0;

    /* variables used for monthly average output */
	int* mondays=0;
	int* enddays=0;
	

	/* copy the input structures into local structures */
	ws = bgcin->ws;
	cinit = bgcin->cinit;
	cs = bgcin->cs;
	ns = bgcin->ns;
	sitec = bgcin->sitec;
	sprop = bgcin->sprop;
	epc = bgcin->epc;
	GWS = bgcin->GWS;

	PLT = bgcin->PLT; 		/* planting variables */
	THN = bgcin->THN; 		/* thinning variables */
	MOW = bgcin->MOW; 		/* mowing variables */
	GRZ = bgcin->GRZ; 		/* grazing variables */
	HRV = bgcin->HRV;		/* harvesting variables */
	PLG = bgcin->PLG;		/* harvesting variables */
	FRZ = bgcin->FRZ;		/* fertilizing variables */
	IRG = bgcin->IRG; 		/* irrigating variables */
	MUL	= bgcin->MUL;		/* mulching variables */
	FLS = bgcin->FLS;		/* flooding variables */
	CWE = bgcin->CWE;		/* CWD-extract variables */

	/* note that the following three structures have dynamic memory elements,
	and so the notion of copying the input structure to a local structure
	value-by-value is not the same as above. In this case, the array pointers
	are being copied, so the local members use the same memory that was
	allocated in the calling function. Note also that bgc.c does not modify
	the contents of these structures. */
	ctrl = bgcin->ctrl;
	metarr = bgcin->metarr;
	co2 = bgcin->co2;
	ndep = bgcin->ndep;


	/********************************************************************************************************* */

	/* writing logfile */
	fprintf(bgcout->log_file.ptr, "NORMAL RUN\n");
	fprintf(bgcout->log_file.ptr, " \n");

	fprintf(bgcout->log_file.ptr, "VEGETATION TYPE\n");
	if (epc.woody == 0)
		fprintf(bgcout->log_file.ptr, "biome type           - NON-WOODY\n");
	else
		fprintf(bgcout->log_file.ptr, "biome type           - WOODY\n");
	
	if (epc.evergreen == 0)
		fprintf(bgcout->log_file.ptr, "woody type           - DECIDUOUS\n");
	else
		fprintf(bgcout->log_file.ptr, "woody type           - EVERGREEN\n");

	if (epc.c3_flag == 0)
		fprintf(bgcout->log_file.ptr, "photosyn.type        - C4 PSN\n");
	else
		fprintf(bgcout->log_file.ptr, "photosyn.type        - C3 PSN\n");
	fprintf(bgcout->log_file.ptr, " \n");

	fprintf(bgcout->log_file.ptr, "CALCULATION METHODS\n");
	
	if (ctrl.STCM_flag == 0)
		fprintf(bgcout->log_file.ptr, "temperature          - MuSo\n");
	else
		fprintf(bgcout->log_file.ptr, "temperature          - DSSAT\n");

	if (ctrl.photosynt_flag == 0)
		fprintf(bgcout->log_file.ptr, "photosynthesis       - Farquhar\n");
	else
		fprintf(bgcout->log_file.ptr, "photosynthesis       - DSSAT\n");

	if (ctrl.ET_flag == 0)
		fprintf(bgcout->log_file.ptr, "evapotranspiration   - Penman-Monteith\n");
	else
		fprintf(bgcout->log_file.ptr, "evapotranspiration   - Priestley-Taylor\n");


	if (ctrl.radiation_flag == 0)
		fprintf(bgcout->log_file.ptr, "radiation            - based on SWabs\n");
	else
	{
		if (ctrl.radiation_flag == 1)
			fprintf(bgcout->log_file.ptr, "radiation        - based on Rn based on Jiang et al (2015)\n");
		else
			fprintf(bgcout->log_file.ptr, "radiation        - based on Rn (difference between net short wawe and net long wave) \n");
	}

	if (ctrl.soilstress_flag == 0)
		fprintf(bgcout->log_file.ptr, "soilstress           - based on VWC\n");
	else
		fprintf(bgcout->log_file.ptr, "soilstress           - based on transpiration demand\n");

	if (ctrl.interception_flag == 0)
		fprintf(bgcout->log_file.ptr, "interception         - based on exponential (saturating) function of PRCP\n");
	else
		fprintf(bgcout->log_file.ptr, "interception         - based on linear function of PRCP\n");

	if (ctrl.MRdeficit_flag == 0)
		fprintf(bgcout->log_file.ptr, "MRdeficit            - based on from total NSC-SC pools\n");
	else
	{
		if (ctrl.MRdeficit_flag == 1) fprintf(bgcout->log_file.ptr, "MRdeficit            - based on from organ-separeted NSC-SC pools\n");
		if (ctrl.MRdeficit_flag == 2) fprintf(bgcout->log_file.ptr, "MRdeficit            - based on from woody/non-woody separeted NSC-SC pools\n");
		if (ctrl.MRdeficit_flag == 3) fprintf(bgcout->log_file.ptr, "MRdeficit            - cover from Cpool\n");
	}

	if (ctrl.Ksat_flag == 0)
		fprintf(bgcout->log_file.ptr, "Ksat                 - exponential function\n"); 
	else
		fprintf(bgcout->log_file.ptr, "Ksat                 - power function\n");


	if (ctrl.transferGDD_flag == 0)
		fprintf(bgcout->log_file.ptr, "transfer period      - EPC\n");
	else
		fprintf(bgcout->log_file.ptr, "transfer period      - GDD\n");

	if (ctrl.q10depend_flag == 0)
		fprintf(bgcout->log_file.ptr, "q10 value            - constant\n");
	else
		fprintf(bgcout->log_file.ptr, "q10 value            - temperature dependent\n");

	if (ctrl.phtsyn_acclim_flag == 0)
		fprintf(bgcout->log_file.ptr, "photosyn. acclim.    - no\n");
	else
	{
		if (ctrl.phtsyn_acclim_flag == 1)
			fprintf(bgcout->log_file.ptr, "photosyn. acclim.    - yes (based on Dyukarev, 2017)\n");
		else
			fprintf(bgcout->log_file.ptr, "photosyn. acclim.    - yes (based on Kattge and Knorr, 2007 \n");
	}
		

	if (ctrl.resp_acclim_flag == 0)
		fprintf(bgcout->log_file.ptr, "respiration acclim.  - no\n");
	else
		fprintf(bgcout->log_file.ptr, "respiration acclim.  - yes\n");

	if (ctrl.CO2conduct_flag == 0)
		fprintf(bgcout->log_file.ptr, "CO2 conduct. effect  - no effect\n");
	else
		fprintf(bgcout->log_file.ptr, "CO2 conduct. effect  - reduction\n");

	if (sprop.Tp1_decomp == DATA_GAP)
		fprintf(bgcout->log_file.ptr, "Decomposition Tresp. - Lloyd and Taylor\n");
	else
		fprintf(bgcout->log_file.ptr, "Decomposition Tresp. - Bell function\n");

	if (sprop.Tp1_nitrif == DATA_GAP)
		fprintf(bgcout->log_file.ptr, "Nitrification Tresp. - Lloyd and Taylor\n");
	else
		fprintf(bgcout->log_file.ptr, "Nitrification Tresp. - Bell function\n");

	fprintf(bgcout->log_file.ptr, " \n");

	fprintf(bgcout->log_file.ptr, "DATA SOURCES\n");

	if (ctrl.phenology_flag == 0) 
	{
		if (PLT.PLT_num != 0)
		{
			fprintf(bgcout->log_file.ptr, "SGS data - planting date from planting file\n");
		}
		else
		{
			if (ctrl.varSGS_flag == 0) 
				fprintf(bgcout->log_file.ptr, "SGS data - user-defined from EPC file\n");
			else
			{
				fprintf(bgcout->log_file.ptr, "SGS data - user-defined from annual varying SGS file\n");
				if (ctrl.onscreen) printf("INFORMATION: reading onday_normal.txt: annual varying SGS data\n");
			}
		}
		if (HRV.HRV_num)
		{
			fprintf(bgcout->log_file.ptr, "EGS data - harvesting date from harvesting file\n");
		}
		else
		{
			if (ctrl.varEGS_flag == 0) 
				fprintf(bgcout->log_file.ptr, "EGS data - user-defined from EPC file\n");
			else
			{
				fprintf(bgcout->log_file.ptr, "EGS data - user-defined from annual varying EGS file\n");
				if (ctrl.onscreen) printf("INFORMATION: reading onday_normal.txt: annual varying SGS data\n");
			}
		}
	}
	else
	{
		if (ctrl.GSI_flag == 0) 
		{
			if (PLT.PLT_num != 0) 
				fprintf(bgcout->log_file.ptr, "SGS data - planting date from planting file\n");
			else
				fprintf(bgcout->log_file.ptr, "SGS data - model estimation (with original method)\n");
			
			if (HRV.HRV_num)
				fprintf(bgcout->log_file.ptr, "EGS data - harvesting date from harvesting file\n");
			else
				fprintf(bgcout->log_file.ptr, "EGS data - model estimation (with original method)\n");
		}
		else
		{
			if (PLT.PLT_num != 0) 
				fprintf(bgcout->log_file.ptr, "SGS data - planting date from planting file\n");
			else
				fprintf(bgcout->log_file.ptr, "SGS data - model estimation (with GSI method)\n");
			
			if (HRV.HRV_num)
				fprintf(bgcout->log_file.ptr, "EGS data - harvesting date from harvesting file\n");
			else
				fprintf(bgcout->log_file.ptr, "EGS data - model estimation (with GSI method)\n");
		}
	}

	if (ctrl.varFM_flag == 0) 
		fprintf(bgcout->log_file.ptr, "FM data - constant \n");
	else
	{
		fprintf(bgcout->log_file.ptr, "FM data - annual varying\n");
		if (ctrl.onscreen) printf("INFORMATION: reading FM_normal.txt - annual varying WPM data\n");
	}

	if (ctrl.varWPM_flag == 0) 
		fprintf(bgcout->log_file.ptr, "WPM data - constant \n");
	else
	{
		fprintf(bgcout->log_file.ptr, "WPM data - annual varying\n");
		if (ctrl.onscreen) printf("INFORMATION: reading WPM_normal.txt - annual varying WPM data\n");
	}

	if (ctrl.varMSC_flag == 0) 
		fprintf(bgcout->log_file.ptr, "MSC data - constant\n");
	else
	{
		fprintf(bgcout->log_file.ptr, "MSC data - annual varying\n");
		if (ctrl.onscreen) printf("INFORMATION: reading conductance_normal.txt: annual varying MSC data\n");
	}

	if (PLT.PLT_num || THN.THN_num || MOW.MOW_num || GRZ.GRZ_num || HRV.HRV_num || PLG.PLG_num || FRZ.FRZ_num || IRG.IRG_num || 
		ctrl.condIRG_flag || ctrl.condMOW_flag)
		fprintf(bgcout->log_file.ptr, "management  - YES\n");
	else
		fprintf(bgcout->log_file.ptr, "management  - NO\n");

	if (GWS.GWnum == 0) 
		fprintf(bgcout->log_file.ptr, "groundwater - NO\n");
	else
	{
		fprintf(bgcout->log_file.ptr, "groundwater - YES\n");
		if (ctrl.onscreen) printf("INFORMATION: reading groundwater_normal.txt - daily GWD data\n");
	}

	if (FLS.FLnum == 0) 
		fprintf(bgcout->log_file.ptr, "flooding - NO\n");
	else
	{
		fprintf(bgcout->log_file.ptr, "flooding - YES\n");
		if (ctrl.onscreen) printf("INFORMATION: reading flooding_normal.txt - daily FLS data\n");
	}

	fprintf(bgcout->log_file.ptr, " \n");

	fprintf(bgcout->log_file.ptr, "SOIL PROPERTIES FOR 10 SOIL LAYERS (POTENTIALLY) ESTIMATED BY THE MODEL \n");
	fprintf(bgcout->log_file.ptr, "Soiltype (based on sand/silt content) :%14s%14s%14s%14s%14s%14s%14s%14s%14s%14s\n",soiltype_name[ctrl.soiltype_array[0]], soiltype_name[ctrl.soiltype_array[1]], soiltype_name[ctrl.soiltype_array[2]], soiltype_name[ctrl.soiltype_array[3]], soiltype_name[ctrl.soiltype_array[4]], soiltype_name[ctrl.soiltype_array[5]], soiltype_name[ctrl.soiltype_array[6]], soiltype_name[ctrl.soiltype_array[7]], soiltype_name[ctrl.soiltype_array[8]], soiltype_name[ctrl.soiltype_array[9]]);
	fprintf(bgcout->log_file.ptr, "Bulk density [kg/m3]                  :%14.1f%14.1f%14.1f%14.1f%14.1f%14.1f%14.1f%14.1f%14.1f%14.1f\n", sprop.BD[0], sprop.BD[1], sprop.BD[2], sprop.BD[3], sprop.BD[4], sprop.BD[5], sprop.BD[6], sprop.BD[7], sprop.BD[8], sprop.BD[9]);
	fprintf(bgcout->log_file.ptr, "Clapp-Hornberger b parameter [dimless]:%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f\n",sprop.soilB[0],sprop.soilB[1],sprop.soilB[2],sprop.soilB[3],sprop.soilB[4],sprop.soilB[5],sprop.soilB[6],sprop.soilB[7],sprop.soilB[8],sprop.soilB[9]);
	fprintf(bgcout->log_file.ptr, "VWC at saturation [m3/m3]:             %14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f\n",sprop.VWCsat[0],sprop.VWCsat[1],sprop.VWCsat[2],sprop.VWCsat[3],sprop.VWCsat[4],sprop.VWCsat[5],sprop.VWCsat[6],sprop.VWCsat[7],sprop.VWCsat[8],sprop.VWCsat[9]);
	fprintf(bgcout->log_file.ptr, "VWC at field capacity [m3/m3]:         %14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f\n",sprop.VWCfc[0],sprop.VWCfc[1],sprop.VWCfc[2],sprop.VWCfc[3],sprop.VWCfc[4],sprop.VWCfc[5],sprop.VWCfc[6],sprop.VWCfc[7],sprop.VWCfc[8],sprop.VWCfc[9]);
	fprintf(bgcout->log_file.ptr, "VWC at wilting point [m3/m3]:          %14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f\n",sprop.VWCwp[0],sprop.VWCwp[1],sprop.VWCwp[2],sprop.VWCwp[3],sprop.VWCwp[4],sprop.VWCwp[5],sprop.VWCwp[6],sprop.VWCwp[7],sprop.VWCwp[8],sprop.VWCwp[9]);
	fprintf(bgcout->log_file.ptr, "VWC at hygroscopic water [m3/m3]:      %14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f\n",sprop.VWChw[0],sprop.VWChw[1],sprop.VWChw[2],sprop.VWChw[3],sprop.VWChw[4],sprop.VWChw[5],sprop.VWChw[6],sprop.VWChw[7],sprop.VWChw[8],sprop.VWChw[9]);
	fprintf(bgcout->log_file.ptr, "PSI at saturation [MPa]:               %14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f\n",sprop.PSIsat[0],sprop.PSIsat[1],sprop.PSIsat[2],sprop.PSIsat[3],sprop.PSIsat[4],sprop.PSIsat[5],sprop.PSIsat[6],sprop.PSIsat[7],sprop.PSIsat[8],sprop.PSIsat[9]);
	fprintf(bgcout->log_file.ptr, "PSI at field capacity [MPa]:           %14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f\n",sprop.PSIfc[0],sprop.PSIfc[1],sprop.PSIfc[2],sprop.PSIfc[3],sprop.PSIfc[4],sprop.PSIfc[5],sprop.PSIfc[6],sprop.PSIfc[7],sprop.PSIfc[8],sprop.PSIfc[9]);
	fprintf(bgcout->log_file.ptr, "PSI at wilting point [MPa]:            %14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f\n",sprop.PSIwp[0],sprop.PSIwp[1],sprop.PSIwp[2],sprop.PSIwp[3],sprop.PSIwp[4],sprop.PSIwp[5],sprop.PSIwp[6],sprop.PSIwp[7],sprop.PSIwp[8],sprop.PSIwp[9]);
	fprintf(bgcout->log_file.ptr, "drainage coefficient [prop]:           %14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f\n",sprop.drainCoeff[0],sprop.drainCoeff[1],sprop.drainCoeff[2],sprop.drainCoeff[3],sprop.drainCoeff[4],sprop.drainCoeff[5],sprop.drainCoeff[6],sprop.drainCoeff[7],sprop.drainCoeff[8],sprop.drainCoeff[9]);
	fprintf(bgcout->log_file.ptr, "hydr. conduct. at saturation [m/day]:  %14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f\n",sprop.hydrCONDUCTsat[0]*nSEC_IN_DAY,sprop.hydrCONDUCTsat[1]*nSEC_IN_DAY,sprop.hydrCONDUCTsat[2]*nSEC_IN_DAY,sprop.hydrCONDUCTsat[3]*nSEC_IN_DAY,sprop.hydrCONDUCTsat[4]*nSEC_IN_DAY,sprop.hydrCONDUCTsat[5]*nSEC_IN_DAY,sprop.hydrCONDUCTsat[6]*nSEC_IN_DAY,sprop.hydrCONDUCTsat[7]*nSEC_IN_DAY,sprop.hydrCONDUCTsat[8]*nSEC_IN_DAY,sprop.hydrCONDUCTsat[9]*nSEC_IN_DAY);
	fprintf(bgcout->log_file.ptr, "capillary fringe [m]:                  %14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f%14.3f\n",sprop.CapillFringe[0],sprop.CapillFringe[1],sprop.CapillFringe[2],sprop.CapillFringe[3],sprop.CapillFringe[4],sprop.CapillFringe[5],sprop.CapillFringe[6],sprop.CapillFringe[7],sprop.CapillFringe[8],sprop.CapillFringe[9]);
	fprintf(bgcout->log_file.ptr, " \n");

	/********************************************************************************************************* */
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/* MEMORY ALLOCATION */
	
	/* in case of natural ecosystem, nyears = number of meteorological year, in case agricultural system: nyears = number of plantings */
	if (!PLT.PLT_num)
		nyears = ctrl.simyears;
	else
		nyears = PLT.PLT_num + ctrl.simyears;
	
	/* allocate memory for local output arrays */
	if (!errorCode && (ctrl.dodaily || ctrl.domonavg || ctrl.doannavg)) 
	{
		dayarr = (double*) malloc(ctrl.ndayout * sizeof(double));
		if (!dayarr)
		{
			printf("ERROR allocating for local daily output array in bgc.c\n");
			errorCode=301;
		}
		for (i=0 ; i<ctrl.ndayout ; i++) dayarr[i] = 0.0;
	}

	if (!errorCode && ctrl.domonavg) 
	{
		monavgarr = (double*) malloc(ctrl.ndayout * sizeof(double));
		if (!monavgarr)
		{
			printf("ERROR allocating for monthly average output array in bgc.c\n");
			errorCode=302;
		}

		for (i=0 ; i<ctrl.ndayout ; i++) monavgarr[i] = 0.0;
	}
	if (!errorCode && ctrl.doannavg) 
	{
		annavgarr = (double*) malloc(ctrl.ndayout * sizeof(double));
		if (!annavgarr)
		{
			printf("ERROR allocating for annual average output array in bgc.c\n");
			errorCode=303;
		}
		for (i=0 ; i<ctrl.ndayout ; i++) annavgarr[i] = 0.0;
	}
	if (!errorCode && ctrl.doannual)
	{
		annarr = (double*) malloc(ctrl.nannout * sizeof(double));
		if (!annarr)
		{
			printf("ERROR allocating for local annual output array in bgc.c\n");
			errorCode=304;
		}
		for (i=0 ; i<ctrl.nannout ; i++) annarr[i] = 0.0;

	}
	
	/* allocate space for the output map pointers */
	if (!errorCode) 
	{
		output_map = (double**) malloc(NMAP * sizeof(double*));
		if (!output_map)
		{
			printf("ERROR allocating for output map in bgc.c\n");
			errorCode=305;
		}
	}
	
	/* allocate space for the onday_arr and offday_arr: first column - year, second column: day*/
	if (!errorCode)
	{
		phenarr.onday_arr  = (int**) malloc(nyears*sizeof(int*));  
        phenarr.offday_arr = (int**) malloc(nyears*sizeof(int*));  
			
		for (i = 0; i<nyears; i++)
		{
			phenarr.onday_arr[i]  = (int*) malloc(2*sizeof(int));  
			phenarr.offday_arr[i] = (int*) malloc(2*sizeof(int));  
			phenarr.onday_arr[i]  = (int*) malloc(2*sizeof(int));  
			phenarr.offday_arr[i] = (int*) malloc(2*sizeof(int));  
		}


		if (!phenarr.onday_arr || !phenarr.offday_arr)
		{
			printf("ERROR allocating for onday_arr/offday_arr, bgc.c\n");
			errorCode=306;
		}
	}

	/* allocate space for enddays */
	if (!errorCode) 
	{
		enddays = (int*) malloc(nMONTHS_OF_YEAR * sizeof(int));
		if (!enddays)
		{
			printf("ERROR allocating for enddays in bgc.c\n");
			errorCode=307;
		}
	}


	/* allocate space for mondays */	
	if (!errorCode) 
	{
		mondays = (int*) malloc(nMONTHS_OF_YEAR * sizeof(int));
		if (!mondays)
		{
			printf("ERROR allocating for mondays in bgc.c\n");
			errorCode=308;
		}
	}
	
	
	/* initialize the output mapping array */
	if (!errorCode && output_map_init(output_map,&phen,&metv,&ws,&wf,&cs,&cf,&ns,&nf,&sprop,&epv,&soilInfo,&psn_sun,&psn_shade,&summary))
	{
		printf("ERROR in call to output_map_init.c from bgc.c\n");
		errorCode=401;
	}
	

	/* atmospheric pressure (Pa) as a function of elevation (m) */
	if (!errorCode && atm_pres(sitec.elev, &metv.pa))
	{
		printf("ERROR in atm_pres.c from bgc.c\n");
		errorCode=402;
	}


    /* calculate GSI to deterime onday and offday 	*/	
	if (ctrl.GSI_flag && !PLT.PLT_num)
	{
		if (!errorCode && GSI_calculation(&metarr, &sitec, &epc, &phenarr, &ctrl))
		{
			printf("ERROR in call to GSI_calculation.c, from bgc.c\n");
			errorCode=403;
		}
	}

	/* calculate conductance limitation factors */	
	if (!errorCode && conductLimit_calculations(&sprop, &epc, &epv))
	{
		printf("ERROR in call to conductLimit_calculations.c, from bgc.c\n");
		errorCode=404;
	}

	/* determine phenological signals */
 	if (!errorCode && prephenology(bgcout->log_file, &epc, &metarr, &PLT, &HRV, &ctrl, &phenarr))
	{
		printf("ERROR in call to prephenology.c, from bgc.c\n");
		errorCode=405;
	}

	/* if this simulation is using a restart file for its initialconditions, then copy restart info into structures */
	if (!errorCode && ctrl.read_restart)
	{
		if (!errorCode && restart_input(&ctrl, &epc, &sitec, &ws, &cs, &ns, &epv, &sprop, &(bgcin->restart_input)))
		{
			printf("ERROR in call to restart_input.c from bgc.c\n");
			errorCode=406;
		}
	}
	

	
	errorCode = firstday(&ctrl, &epc, &PLT, &sprop, &sitec, &cinit, &phen, &epv, &soilInfo, &ws, &cs, &ns, &summary, &psn_sun, &psn_shade);
	if (errorCode)
	{
		/* errorCode = 407 */ 
		printf("ERROR in call to firstday.c, from bgc.c\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}

	/* zero water, carbon, and nitrogen source and sink variables */
	if (!errorCode && zero_srcsnk(&cs,&ns,&ws,&summary))
	{
		printf("ERROR in call to zero_srcsnk.c, from bgc.c\n");
		errorCode=408;
	}


	/* initialize the indicator for first day of current simulation, so that the checks for mass balance can have two days for comparison */
	first_balance = 1;
		
	/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
	/* 1. BEGIN OF THE ANNUAL LOOP */

	for (simyr=0 ; !errorCode && simyr<ctrl.simyears ; simyr++)
	{


		/* set current month to 0 (january) at the beginning of each year */
		ctrl.curmonth = 0;

		 /* counters into control sturct */ 
		ctrl.simyr = simyr;
		ctrl.spinyears = 0;


		if (!errorCode && leapControl(ctrl.simstartyear+simyr, enddays, mondays, &leap))
		{
			printf("ERROR in call to leapControl.c from bgc.c\n");
			errorCode=500;
		}


		/* set the max lai, maturity and flowering variables, for annual diagnostic output */
		epv.annmax_lai = 0.0;
		epv.annmax_rootDepth = 0.0;
		epv.annmax_plantHeight = 0.0;
		summary.annmax_livingBIOMabove = 0.0;
		summary.annmax_livingBIOMbelow = 0.0;
		summary.annmax_BIOMaboveSUM = 0.0;
		summary.annmax_BIOMbelowSUM = 0.0;
	
		/* atmospheric CO2 handling */
		if (!(co2.varco2))
		{
			/* constant CO2 */
			metv.co2 = co2.co2ppm;
		}
		else 
		{
            /* CO2 from file */
			metv.co2 = co2.co2ppm_array[simyr];
		}
		
		/* atmospheric Ndep handling */
		if (!(ndep.varndep))
		{
			/*constant Ndep */
			dailyNdep = ndep.ndep / nDAYS_OF_YEAR;
		}
		else
		{	
			/* Ndep from file */
			dailyNdep = ndep.Ndep_array[simyr] / nDAYS_OF_YEAR;
		}
		
		if (ctrl.onscreen) printf("-------------------\n");
		if (ctrl.onscreen) printf("Year: %d\t\n",ctrl.simstartyear+simyr);

	

		/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
		/* 2. BEGIN OF THE DAILY LOOP */
		
		for (yday=0 ; !errorCode && yday<nDAYS_OF_YEAR ; yday++)
		{

			/* set the day index for meteorological and phenological arrays */
			ctrl.yday   = yday;
			ctrl.metday	= simyr*nDAYS_OF_YEAR + yday;		

			
			/* set fluxes to zero */
			if (!errorCode && make_zero_flux_struct(&ctrl, &wf, &cf, &nf, &soilInfo, &summary))
			{
				printf("ERROR in call to make_zero_flux_struct.c from bgc.c\n");
				errorCode=501;
			}

		
			/* initalizing annmax and cumulative variables */
			if (yday == 0)
			{
				if (!errorCode && annVARinit(&summary, &epv, &cs, &ws))
				{
					printf("ERROR in call to annVARinit.c from bgc.c\n");
					errorCode=502;
				}
			}
		

			/* set the day index for meteorological and phenological arrays */
			ctrl.metday = simyr*nDAYS_OF_YEAR + yday;
			

			/* nitrogen deposition and fixation */
			nf.ndep_to_sminn_total = dailyNdep;
			nf.nfix_to_sminn_total = epc.nfix / nDAYS_OF_YEAR;

	
			/* calculating actual onday and offday */
			if (!errorCode && dayphen(&ctrl, &epc, &phenarr, &PLT, &phen))
			{
				printf("ERROR in dayphen.c from bgc.c\n");
				errorCode=503;
			}
			

			/* setting MANAGEMENT and GW  DAYS based on input data */
			if (!errorCode && management(&ctrl, &FRZ, &GRZ, &HRV, &MOW, &PLT, &PLG, &THN, &IRG, &MUL, &CWE, &FLS, &GWS, mondays))
			{
				printf("ERROR in management days.c from bgc.c\n");
				errorCode=504;
			}
	

			/* GROUNDWATER calculations */
			if (!errorCode && groundwater_calculations(&ctrl, &sitec, &GWS, &sprop, &soilInfo, &epv, &ws, &wf, &cs, &ns))
			{
				printf("ERROR in groundwater_calculations.c from bgc.c\n");
				errorCode = 505;
			}


			/* daily meteorological variables from metarrays */
			if (!errorCode && daymet(&ctrl, &metarr, &epc, &metv, ws.snoww))
			{
				printf("ERROR in daymet.c from bgc.c\n");
				errorCode=506;
			}
			
			/* phenophases calculation */
			if (!errorCode && phenphase(bgcout->log_file, &ctrl, &epc, &sprop, &PLT, &phen, &metv, &epv, &cs))
			{
				printf("ERROR in phenphase.c from bgc.c\n");
				errorCode=507;
			}
	
			
			/* soil temperature calculations */
			if (!errorCode && multilayer_Tsoil(&ctrl,&epc, &sitec, &sprop, &epv, yday, ws.snoww, &metv))
			{
				printf("ERROR in multilayer_Tsoil.c from bgc.c\n");
				errorCode=508;
			}
			

			/* soilCover calculations */
			if (!errorCode && soilCover(&sitec, &sprop, &metv, &epv, &cs))
			{
				printf("ERROR in soilCover.c from bgc.c\n");
				errorCode=509;
			}
			

			/* phenology calculation */
			if (!errorCode && phenology(&ctrl, &epc, &cs, &ns, &phen, &metv, &epv, &cf, &nf))
			{
				printf("ERROR in phenology.c from bgc.c\n");
				errorCode=510;
			}
		
			/* calculate leaf area index, sun and shade fractions, and specific leaf area for sun and shade canopy fractions, then calculate canopy radiation interception and transmission */          
			if (!errorCode && radtrans(&ctrl, &phen, &cs, &epc, &sitec, &metv, &epv))
			{
				printf("ERROR in radtrans.c from bgc.c\n");
				errorCode=511;
			}
			
		
			/* update the annmax LAI/rootingDepth/plantHeight for annual diagnostic output */
			if (epv.projLAI > epv.annmax_lai)             epv.annmax_lai = epv.projLAI;
			if (epv.rootDepth > epv.annmax_rootDepth)      epv.annmax_rootDepth = epv.rootDepth;
			if (epv.plantHeight > epv.annmax_plantHeight)  epv.annmax_plantHeight = epv.plantHeight;

			
			/* IRRIGATING separately from other management routines*/
			if (!errorCode && irrigating(&ctrl, &IRG, &sitec, &sprop, &epv, &ws, &wf))
			{
				printf("ERROR in irrigating.c from bgc.c\n");
				errorCode=512;
			}
	
			/* precip routing (when there is precip) */
			if (!errorCode && metv.prcp && prcpANDrunoffH(&ctrl, &ws, &metv, &sprop, &epc, &epv, &wf))
			{
				printf("ERROR in prcpANDrunoffH.c from bgc.c\n");
				errorCode=513;
			}
			
			
			/* snowmelt (when there is a snowpack) */
			if (!errorCode && ws.snoww && snowmelt(&metv, &wf, &ws))
			{
				printf("ERROR in snowmelt.c from bgc.c\n");
				errorCode=514;
			}
			
			

			/* potential evaporation and transpiration */
			if (!errorCode && Elimit_and_PET(&ctrl, &epc, &sprop, &metv, &epv, &wf))
			{
				printf("ERROR in Elimit_and_PET.c from bgc.c\n");
				errorCode=515;
			}
			
	
			/* conductance calculation */
			if (!errorCode && conduct_calc(&ctrl, &metv, &epc, &epv, simyr))
			{
				printf("ERROR in conduct_calc.c from bgc.c\n");
				errorCode=516;
			}
			
	

			/* begin canopy bio-physical process simulation */
			/* do canopy ET calculations whenever there is leaf area displayed, since there may be intercepted water on the canopy that needs to be dealt with */
			if (!errorCode && epv.n_actphen >= epc.n_emerg_phenophase && metv.dayl)
			{
				/* evapotranspiration */
				if (!errorCode && cs.leafc && canopy_et(&ctrl, &metv, &epv, &wf))
				{
					printf("ERROR in canopy_et.c from bgc.c\n");
					errorCode=517;
				}
			}
		

			/* daily maintenance respiration */
			if (!errorCode && maint_resp(&ctrl, &PLT, &cs, &ns, &epc, &metv, &epv, &cf))
			{
				printf("ERROR in m_resp.c from bgc.c\n");
				errorCode=518;
			}
			

			/* photosynthesis calculation */
			if (!errorCode && cs.leafc && photosynthesis(&ctrl, &epc, &metv, &cs, &ws, &phen, &epv, &psn_sun, &psn_shade, &cf))
			{
				printf("ERROR in photosynthesis.c from bgc.c\n");
				errorCode=519;
			}
		
	
			/* daily litter and soil decomp and nitrogen fluxes */
			if (!errorCode && decomp(&metv,&epc, &sprop, &sitec,&cs,&ns,&epv,&cf,&nf,&nt))
			{
				printf("ERROR in decomp.c from bgc.c\n");
				errorCode=520;
			}
			

		
			/* Daily allocation gets called whether or not this is a current growth day, because the competition between decomp immobilization fluxes 
			and plant growth N demand is resolved here.  On days with no growth, no allocation occurs, but immobilization fluxes are updated normally */

			if (!errorCode && daily_allocation(&ctrl,&epc,&sprop,&metv,&ndep,&cs,&ns,&cf,&nf,&epv,&nt,0))
			{
				printf("ERROR in daily_allocation.c from bgc.c\n");
				errorCode=521;
			}
					
			
			/* heat stress during flowering can affect daily allocation of yield */
			if (epc.n_flowHS_phenophase > 0)
			{
				if (!errorCode && flowering_heatstress(&epc, &metv, &cs, &epv, &cf, &nf))
				{
					printf("ERROR in flowering_heatstress.c from bgc.c\n");
					errorCode=522;
				}
			}
			

			/* reassess the annual turnover rates for livewood --> deadwood, and for evergreen leaf and fine root litterfall. 
			This happens once each year, on the annual_alloc day (the last litterfall day - test for annual allocation day) */
			
			if (phen.remdays_litfall == 1) 
				annual_alloc = 1;
			else 
				annual_alloc = 0;
		

			/* litterfall rates once a year */
			if (!errorCode && annual_alloc)
			{
				if (!errorCode && annual_rates(&epc,&epv))
				{
					printf("ERROR in annual_rates.c from bgc.c\n");
					errorCode=523;
				}
				
			} 
			
			/* daily growth respiration */
			if (!errorCode && growth_resp(&epc, &cf))
			{
				printf("ERROR in growth_resp.c from bgc.c\n");
				errorCode=524;
			}
			
			
			/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
			/* 3. WATER CALCULATIONS WITH STATE UPDATE */

	
			/* EVAPORATION: calculation of actual evaporation from potential evaporation */
			if (!errorCode && potEVPsurface_to_actEVPsurface(&sprop, &epv, &ws, &wf))
			{
				printf("ERROR in potEVPsurface_to_actEVPsurface.c from bgc.c\n");
				errorCode=525;
			}

	
		    /* multilayer soil hydrology: percolation calculation based on PRCP, RUNOFF, EVP, TRP */
			if (!errorCode)
			{
				errorCode = multilayer_hydrolprocess(&ctrl, &sitec, &sprop, &soilInfo, &epc, &epv, &ws, &wf, &ns, &nf, &cs, &cf, &FLS, mondays);
				if (errorCode)
				{
					/* errorCode = 526 */
					printf("ERROR in call to multilayer_hydrolprocess.c, from bgc.c\n");
					writeErrorCode(errorCode);
					exit(errorCode);
				}
			}
	
			/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
			/* 4. STATE UPDATE */

			/* daily update of the water state variables */
			if (!errorCode && water_state_update(&wf, &ws))
			{
				printf("ERROR in water_state_update.c from bgc.c\n");
				errorCode=527;
			}

	

			/* daily update of carbon and nitrogen state variables */
   	    	if (!errorCode && CN_state_update(&sitec, &epc, &soilInfo, &sprop, &ctrl, &epv, &cf, &nf, &cs, &ns, annual_alloc, epc.evergreen))
			{
				printf("ERROR in CN_state_update.c from bgc.c\n");
				errorCode=528;
			}
			
  

			/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
			/* 5. MORTALITY AND NITROGEN FLUXES CALCULATION WITH OWN STATE UPDATE: 
			to insure that pools don't go negative due to mortality/leaching fluxes conflicting with other proportional fluxes */

			
			/* calculate daily senescence mortality fluxes and update state variables */
			if (!errorCode && senescence(&sprop, &epc, &GRZ, &metv, &ctrl, &cs, &cf, &ns, &nf, &epv))
			{
				printf("ERROR in senescence.c from bgc.c\n");
				errorCode=529;
			}
			
	
		    /* calculate daily mortality fluxes  and update state variables */
			if (!errorCode && mortality(&ctrl,  &sprop, &epc, &epv, &cs, &cf, &ns, &nf, simyr))
			{
				printf("ERROR in mortality.c from bgc.c\n");
				errorCode=530;
			}
			
			/* calculate the change of soil mineralized N in multilayer soil */
			if (!errorCode && multilayer_sminn(&ctrl, &metv, &sitec, &ndep, &cs, &cf, &ns, &nf, &sprop, &epv, &soilInfo))
			{
				printf("ERROR in multilayer_sminn.c from spinup_bgc.c\n");
				errorCode = 531;
			}

	
			/* calculate the leaching of N, DOC and DON from multilayer soil */
			if (!errorCode && multilayer_leaching( &sprop,  &soilInfo, &cs, &ns, &ws, &wf))
			{
				printf("ERROR in multilayer_leaching.c from bgc.c\n");
				errorCode=532;
			}		


			/* calculate summary variables */
			if (!errorCode && aboveANDbelow(&sprop, &epv, &cs, &cf))
			{
				printf("ERROR in aboveANDbelow.c from bgc.c\n");
				errorCode = 533;
			}



			/* ------------------------------------------------------------------*/
			/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
			/* 6. MANAGEMENT FLUXES */
			
			/* PLANTING */
			if (!errorCode && planting(&ctrl, &sitec, &sprop, &PLT, &epc, &epv, &phen, &cs, &ns, &cf, &nf))
			{
				printf("ERROR in planting.c from bgc.c\n");
				errorCode=534;
			}
			 
	
		   	/* THINNIG  */
			if (!errorCode && thinning(&ctrl, &epc, &THN, &cs, &ns, &ws, &cf, &nf, &wf, &epv))
			{
				printf("ERROR in thinning.c from bgc.c\n");
				errorCode=535;
			}
			

			/* MOWING  */
			if (!errorCode && mowing(&ctrl, &epc, &MOW, &epv, &cs, &ns, &ws, &cf, &nf, &wf))
			{
				printf("ERROR in mowing.c from bgc.c\n");
				errorCode=536;
			}
			

			/* grazing  */
			if (!errorCode && wf.FL_to_soilwTOTAL==0 && grazing(&ctrl, &epc, &sprop, &GRZ, &epv, &cs, &ns, &ws, &cf, &nf, &wf, mondays))
			{
				printf("ERROR in grazing.c from bgc.c\n");
				errorCode=537;
			}
			
	
		   	/* HARVESTING  */
			if (!errorCode && harvesting(bgcout->econout_file, &ctrl, &phen, &epc, &HRV, &IRG, &epv, &cs, &ns, &ws, &cf, &nf, &wf))
			{
				printf("ERROR in harvesting.c from bgc.c\n");
				errorCode=538;
			}
			
 
			/* PLOUGHING */
 			if (!errorCode && ploughing(&ctrl, &epc, &sitec, &sprop, &metv, &epv, &PLG, &cs, &ns, &ws, &cf, &nf, &wf))
			{
				printf("ERROR in ploughing.c from bgc.c\n");
				errorCode=539;
			}
			
	
	
			/* FERTILIZING  */
	        if (!errorCode && fertilizing(&ctrl, &sitec, &sprop, &soilInfo, &epv, &FRZ, &cs, &ns, &ws, &cf, &nf, &wf))
			{
				printf("ERROR in fertilizing.c from bgc.c\n");
				errorCode=540;
			}	
			

			/* MULCHING */
			if (!errorCode && mulching(&ctrl, &MUL, &cs, &ns, &cf, &nf))
			{
				printf("ERROR in mulching.c from bgc.c\n");
				errorCode=541;
			}
			

			/* CWD-extract */
			if (!errorCode && CWDextract(&ctrl, &CWE, &cs, &ns, &cf, &nf))
			{
				printf("ERROR in CWDextract.c from bgc.c\n");
				errorCode=542;
			}
			
				
			/* cut-down plant material (due to management) */
			if (!errorCode && cutdown2litter(&sprop, &epc, &epv, &cs, &cf, &ns, &nf))
			{
				printf("ERROR in cutdown2litter.c from bgc.c\n");
				errorCode=543;
			}
			

			/* calculating rooting depth, n_rootlayers, n_maxrootlayers, rootlengthProp */
 			 if (!errorCode && multilayer_rootDepth(&epc, &sprop, &cs, &sitec, &epv))
			 {
				printf("ERROR in multilayer_rootDepth.c from bgc.c\n");
				errorCode=544;
			 }
			 

			/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
			/* 7. ERROR CHECKING AND SUMMARY VARIABLES  */
	
			/* test for very low state variable values and force them to 0.0 to avoid rounding and floating point overflow errors */
			if (!errorCode && precision_control(&ws, &cs, &ns, &sprop, &soilInfo))
			{
				printf("ERROR in precision_control.c from bgc.c\n");
				errorCode=545;
			} 

	
			/* test for water balance*/
			if (!errorCode && sprop.GWlayer != DATA_GAP && check_virtualLayer_balance(&ctrl, &soilInfo, &sprop, &wf))
			{
				printf("ERROR in check_virtualLayer_balance.c from bgc.c\n");
				errorCode = 546;
			}


			/* test for water balance*/
 			if (!errorCode && check_water_balance(&ws, first_balance))
			{
				printf("ERROR in check_water_balance.c from bgc.c\n");
				errorCode=547;
			}
			
	
	        /* test for carbon balance */
			if (!errorCode && check_carbon_balance(&cs, first_balance))
			{
				printf("ERROR in check_carbon_balance.c from bgc.c\n");
				errorCode=548;
			}		

	
			/* test for nitrogen balance */
			if (!errorCode && check_nitrogen_balance(&ns, first_balance))
			{
				printf("ERROR in check_nitrogen_balance.c from bgc.c\n");
				errorCode=549;
			}	


			/* calculate summary variables */
 			if (!errorCode && cnw_summary(&epc, &sitec, &sprop, &metv, &cs, &cf, &ns, &nf, &wf, &soilInfo, &epv, &summary))
			{
				printf("ERROR in cnw_summary.c from bgc.c\n");
				errorCode=550;
			}
	
		
			/* output handling */
			if (!errorCode && output_handling(mondays, enddays, &ctrl, output_map, dayarr, monavgarr, annavgarr, annarr, 
				                            bgcout->dayout, bgcout->monavgout, bgcout->annavgout, bgcout->annout))
			{
				printf("ERROR in output_handling.c from bgc.c\n");
				errorCode=551;
			}
			

			/*  if no dormant period (e.g. evergreen): last day is the dormant day */
			if (phen.offday - phen.onday == 364 && phen.offday == phen.yday_total) 
			{
				epv.n_actphen = 0;
				phen.onday = -1;
				phen.offday = -1;
				phen.remdays_litfall =-1;
			}

		
			/* at the end of first day of simulation, turn off the first_balance switch */
			if (first_balance) first_balance = 0;


			/* if this is the last day of the current month: increment current month counter */
			if (yday == enddays[ctrl.curmonth]) ctrl.curmonth++;

			ctrl.firstsimday_flag = 0;

		}   /* end of daily model loop */
		

	}   /* end of annual model loop */


	/********************************************************************************************************* */
	/* wrinting log file */
	
	fprintf(bgcout->log_file.ptr, " \n");
	fprintf(bgcout->log_file.ptr, " \n");
	/* calculate conductance limitation factors 	*/
	if (!errorCode && conductLimit_information(bgcout->log_file, &ctrl, &epc, &epv))
	{
		printf("ERROR in call to conductLimit_calculations.c, from spinup_bgc.c\n");
		errorCode = 404;
	}

	/* most important informations - onscreen */
	if (ctrl.onscreen) 
	{
		printf("\n");
		printf("INFORMATION FROM THE LAST SIMULATION YEAR\n");
		printf("Maximum projected LAI [m2/m2]:                 %12.1f\n",epv.annmax_lai);
		printf("Aboveground litter carbon content [kgC/m2]:    %12.1f\n", cs.litrCabove_total);
		printf("Aboveground CWD carbon content [kgC/m2]:       %12.1f\n", cs.cwdCabove_total);
		printf("Total soil carbon content [kgC/m2]:            %12.1f\n",summary.soilC_total);
		printf("Total stable soil carbon content [kgC/m2/year]:%12.1f\n", cs.soil4c_total);
	}

	if (cs.CbalanceERR != 0) CbalanceERR = log10(cs.CbalanceERR);
	if (ns.NbalanceERR != 0) NbalanceERR = log10(ns.NbalanceERR);
	if (ws.WbalanceERR != 0) WbalanceERR = log10(ws.WbalanceERR);
	if (cs.CNratioERR  != 0) CNratioERR  = log10(cs.CNratioERR);

	fprintf(bgcout->log_file.ptr, " \n");
	fprintf(bgcout->log_file.ptr, "SOME IMPORTANT ANNUAL OUTPUTS FROM LAST SIMULATION YEAR\n");
    fprintf(bgcout->log_file.ptr, "Cumulative sum of GPP [gC/m2/year]:                      %12.1f\n",summary.cumGPP * kg_to_g);
	fprintf(bgcout->log_file.ptr, "Cumulative sum of TER [gC/m2/year]:                      %12.1f\n",summary.cumTR  * kg_to_g);
	fprintf(bgcout->log_file.ptr, "Cumulative sum of NEE [gC/m2/year]:                      %12.1f\n",summary.cumNEE * kg_to_g);
	fprintf(bgcout->log_file.ptr, "Cumulative sum of ET  [kgH2O/m2/year]:                   %12.1f\n",summary.cumET);
	fprintf(bgcout->log_file.ptr, "Cumulative sum of soil evaporation [kgH2O/m2/year]:      %12.1f\n",summary.cumEVPsurface);
	fprintf(bgcout->log_file.ptr, "Cumulative sum of transpiration [kgH2O/m2/year]:         %12.1f\n",summary.cumETcanopy);
	fprintf(bgcout->log_file.ptr, "Cumulative sum of N2O flux [gN/m2/year]:                 %12.2f\n",summary.cumN2Oflux* kg_to_g);
	fprintf(bgcout->log_file.ptr, "Maximum projected LAI [m2/m2]:                           %12.2f\n",epv.annmax_lai);
	fprintf(bgcout->log_file.ptr, "Maximum rooting depth [m2/m2]:                           %12.2f\n",epv.annmax_rootDepth);
	fprintf(bgcout->log_file.ptr, "Aboveground litter carbon content [kgC/m2/year]:         %12.2f\n",cs.litrCabove_total);
	fprintf(bgcout->log_file.ptr, "Aboveground CWD carbon content [kgC/m2/year]:            %12.2f\n",cs.cwdCabove_total);
	fprintf(bgcout->log_file.ptr, "Soil carbon content (in 0-30 cm soil layer) [%%]:         %12.2f\n",summary.SOCpercent_top30);
	fprintf(bgcout->log_file.ptr, "Total litter carbon content [kgC/m2/year]:               %12.2f\n",summary.litrC_total);
	fprintf(bgcout->log_file.ptr, "Total soil carbon content [kgC/m2/year]:                 %12.2f\n",summary.soilC_total);
	fprintf(bgcout->log_file.ptr, "Total stable soil carbon content [kgC/m2/year]:          %12.2f\n", cs.soil4c_total);
	fprintf(bgcout->log_file.ptr, "Averaged available soil ammonium content (0-30 cm) [ppm]:%12.2f\n",summary.NH4ppmAVAIL_top30);
	fprintf(bgcout->log_file.ptr, "Averaged available soil nitrate content (0-30 cm) [ppm]: %12.2f\n",summary.NO3ppmAVAIL_top30);
	fprintf(bgcout->log_file.ptr, "Averaged soil water content  [m3/m3]:                    %12.2f\n",epv.VWC_avg);
	fprintf(bgcout->log_file.ptr, " \n");
	fprintf(bgcout->log_file.ptr, "10-base logarithm of the maximum carbon balance diff.:   %12.1f\n",CbalanceERR);
	fprintf(bgcout->log_file.ptr, "10-base logarithm of the maximum nitrogen balance diff.: %12.1f\n",NbalanceERR);
	fprintf(bgcout->log_file.ptr, "10-base logarithm of the maximum water balance diff.:    %12.1f\n",WbalanceERR);
	fprintf(bgcout->log_file.ptr, "10-base logarithm of the C-N calc. numbering error:      %12.1f\n",CNratioERR);
	fprintf(bgcout->log_file.ptr, " \n");

	fprintf(bgcout->log_file.ptr, "WARNINGS \n");
	if (!ctrl.limitTRP_flag && !ctrl.limitEVP_flag && !ctrl.limitleach_flag && !ctrl.limitdiffus_flag &&
		!ctrl.limitSNSC_flag && !ctrl.limitMR_flag && !ctrl.limitDENIT_flag && !ctrl.noTRP_flag && !ctrl.pond_flag && !ctrl.grazingW_flag &&
		!ctrl.condMOWerr_flag && !ctrl.condIRGerr_flag && !ctrl.condIRGerr_flag && !ctrl.prephen1_flag && !ctrl.prephen2_flag && 
		!ctrl.bareground_flag && !ctrl.vegper_flag && !ctrl.allocControl_flag)
	{
		fprintf(bgcout->log_file.ptr, "no WARNINGS\n");
	}
	else
	{
		if (ctrl.limitTRP_flag)
		{
			fprintf(bgcout->log_file.ptr, "Limited transpiration due to dry soil\n");
			ctrl.limitTRP_flag = -1;
		}

		if (ctrl.limitEVP_flag)
		{
			fprintf(bgcout->log_file.ptr, "Limited evaporation due to dry soil\n");
			ctrl.limitEVP_flag = -1;
		}

		if (ctrl.limitleach_flag)
		{
			fprintf(bgcout->log_file.ptr, "Limited leaching\n");
			ctrl.limitleach_flag = -1;
		}

		if (ctrl.limitdiffus_flag)
		{
			fprintf(bgcout->log_file.ptr, "Limited diffusing\n");
			ctrl.limitdiffus_flag = -1;
		}

		if (ctrl.limitSNSC_flag)
		{
			fprintf(bgcout->log_file.ptr, "Limited genetically programmed leaf senescence\n");
			ctrl.limitSNSC_flag = -1;
		}

		if (ctrl.limitMR_flag)
		{
			fprintf(bgcout->log_file.ptr, "Limited maintanance respiration \n");
			ctrl.limitMR_flag = -1;
		}

		if (ctrl.limitDENIT_flag)
		{
			fprintf(bgcout->log_file.ptr, "Limited denitrification due high soil respiration\n");
			ctrl.limitDENIT_flag = -1;
		}

		if (ctrl.noTRP_flag)
		{
			fprintf(bgcout->log_file.ptr, "No transpiration (full limitation) due to dry soil\n");
			ctrl.noTRP_flag = -1;
		}

		if (ctrl.pond_flag)
		{
			fprintf(bgcout->log_file.ptr, "Pond water on soil surface\n");
			ctrl.noTRP_flag = -1;
		}

		if (ctrl.grazingW_flag)
		{
			fprintf(bgcout->log_file.ptr, "Not enough grass for grazing\n");
			ctrl.grazingW_flag = -1;
		}

		if (ctrl.condMOWerr_flag)
		{
			fprintf(bgcout->log_file.ptr, "If conditional MOWING flag is on, no MOWING is possible\n");
			ctrl.condMOWerr_flag = -1;
		}

	
		if (ctrl.condIRGerr_flag)
		{
			fprintf(bgcout->log_file.ptr, "If conditional IRRIGATING flag is on, no IRRIGATING is possible\n");
			ctrl.condIRGerr_flag = -1;
		}


		if (ctrl.prephen1_flag)
		{
			fprintf(bgcout->log_file.ptr, "In case of planting model-defined phenology is not possible (firstday:planting, lasTday:harvesting)\n");
			ctrl.prephen1_flag = -1;
		}

		if (ctrl.prephen2_flag)
		{
			fprintf(bgcout->log_file.ptr, "In case of user-defined phenology GSI calculation is not possible\n");
			ctrl.prephen2_flag = -1;
		}

		if (ctrl.bareground_flag)
		{
			fprintf(bgcout->log_file.ptr, "User-defined bareground run (onday and offday set to -9999 in EPC)\n");
			ctrl.bareground_flag = -1;
		}

		if (ctrl.vegper_flag)
		{
			fprintf(bgcout->log_file.ptr, "Vegetation period has not ended until the last day of year, the offday is equal to the last day of year\n");
			ctrl.vegper_flag = -1;
		}

		if (ctrl.allocControl_flag)
		{
			fprintf(bgcout->log_file.ptr, "Adjustment of allocation parameters due to small error (<10-4) in the setting of allocation parameters (see EPC file, lines 132-139)\n");
			ctrl.allocControl_flag = -1;
		}
	
	}


	/********************************************************************************************************* */

	
	/* 8. RESTART OUTPUT HANDLING */
	/* if write_restart flag is set, copy data to the output restart struct */
	if (!errorCode && ctrl.write_restart)
	{
		if (restart_output( &ws, &cs, &ns, &epv, &(bgcout->restart_output)))
		{
			printf("ERROR in call to restart_output.c from bgc.c\n");
			errorCode=600;
		}
	
	}

	
	/* free memory for local output arrays  */
	
    if ((!errorCode || errorCode > 301) && ctrl.dodaily) free(dayarr);
	if ((!errorCode || errorCode > 302) && ctrl.domonavg) free(monavgarr);
	if ((!errorCode || errorCode > 303) && ctrl.doannavg) free(annavgarr);
	if ((!errorCode || errorCode > 304) && ctrl.doannual) free(annarr); 
	if ((!errorCode || errorCode > 305)) free(output_map);
    if (!errorCode  || errorCode > 307) free(enddays);
	if (!errorCode  || errorCode > 308) free(mondays);

	if (((!errorCode || errorCode > 306) && !ctrl.GSI_flag) || ((!errorCode || errorCode > 405) && ctrl.GSI_flag))
	{
		free(phenarr.onday_arr);
		free(phenarr.offday_arr);
		if (ctrl.GSI_flag)
		{
			free(phenarr.Tmin_index);
			free(phenarr.vpd_index);
			free(phenarr.heatsum_index);
			free(phenarr.dayl_index);
			free(phenarr.gsi_indexAVG);
			free(phenarr.heatsum);

		}
	}

	

	/* print timing info if error */
	if (errorCode)
	{
		printf("ERROR at year %d\n",simyr-1);
		printf("ERROR at yday %d\n",yday-1);
	}
	


	/* return error status */	
	return (errorCode);
}
