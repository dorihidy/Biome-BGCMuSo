	/*
pointbgc.c
front-end to BIOME-BGC for single-point, single-biome simulations
Uses BBGC MuSo v6 library function

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
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <time.h>

#include "ini.h"              /* general file structure and I/O prototypes */
#include "bgc_struct.h"       /* data structures for bgc.c */
#include "pointbgc_struct.h"   /* data structures for point driver */
#include "pointbgc_func.h"     /* function prototypes for point driver */
#include "bgc_io.h"           /* bgc.c interface definition */
#include "bgc_constants.h"      /* array structure for epc-by-vegtype */

int main(int argc, char *argv[])
{
	int errorCode=0;
	int transient=0;

	/* bgc input and output structures */
	bgcin_struct bgcin;
	bgcout_struct bgcout;

	/* local control information */
	point_struct point;
	restart_ctrl_struct restart;
	climchange_struct scc;
	output_struct output;
	
	/* initialization file */
	file init;
	
	/* system time variables */
	struct tm *tm_ptr;
	time_t lt;
	
	/* get the system time at start of simulation  - WARNING: unsafe functions: localtime, asctime */
	lt = time(NULL);
	tm_ptr = localtime(&lt);
	strcpy(point.systime,asctime(tm_ptr));
	
	
	/* initialize the bgcin state variable structures before filling with
	values from ini file */
	if (presim_state_init(&bgcin.ws, &bgcin.cs, &bgcin.ns, &bgcin.cinit, &bgcin.ctrl, &bgcin.soilInfo))
	{
		printf("ERROR in call to presim_state_init.c from pointbgc.c ... Exiting.c\n");
		writeErrorCode(101);
		exit(101);
	}



	/******************************
	**                           **
	**  BEGIN READING INIT FILE  **
	**                           **
	******************************/
	
	/* wrinting on screen: model version */
	if(argc > 1)
	{
        if(!strcmp(argv[1],"-v"))
		{
           	printf("Model version: Biome-BGCMuSo7.0alfa\n");
			exit(0);
        }
    }
	

	/* read the name of the main init file from the command line and store as init.name */
	if (argc != 2)
	{
		printf("ERROR in reading the main init file from command line. Exiting\n");
		printf("Correct usage: <executable name>  <initialization file name>\n");
		writeErrorCode(102);
		exit(102);
	} 
	
	strcpy(init.name, argv[1]);

	/* open the main init file for ascii read and check for errors */
	if (file_open(&init,'i',1))
	{
		printf("ERROR opening init file, pointbg.c ... Exiting\n");
		writeErrorCode(103);
		exit(103);
	}

	/* read the header string from the init file */
	if (fgets(point.header, 100, init.ptr)==NULL)
	{
		printf("ERROR reading header string: pointbgc.c ... Exiting\n");
		writeErrorCode(201);
		exit(201);
	}

	/* open met file, discard header lines */
	errorCode = met_init(init, &point);
	if (errorCode)
	{
		printf("ERROR in call to met_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}

	/* read restart control parameters */
	errorCode = restart_init(init, &restart);
	if (errorCode)
	{
		printf("ERROR in call to restart_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}

	/* read simulation timing control parameters */
	errorCode = time_init(init, &point, &(bgcin.ctrl));
	if (errorCode)
	{
		printf("ERROR in call to time_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}
	
	
	/* read CO2 control parameters */
	errorCode = co2_init(init, &(bgcin.co2), &(bgcin.ctrl));
	if (errorCode)
	{
		printf("ERROR in call to co2_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}


	/* read varied nitrogen deposition block */
	errorCode = ndep_init(init, &bgcin.ndep, &(bgcin.ctrl));
	if (errorCode)
	{
		printf("ERROR in call to ndep_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}


	/* read site constants */
	errorCode = sitec_init(init, &bgcin.sitec, &(bgcin.ctrl));
	if (errorCode)
	{
		printf("ERROR in call to sitec_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}
	
	/* read soil properties */
	errorCode = sprop_init(init, &bgcin.sprop, &bgcin.ctrl);
	if (errorCode)
	{
		printf("ERROR in call to sprop_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}
	
	/* read ecophysiological constants */
	errorCode = epc_init(init, &bgcin.epc, &bgcin.ctrl, 1);
	if (errorCode)
	{
		printf("ERROR in call to epc_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}

	/* read management file with management information */
	errorCode = mgm_init(init, &bgcin.ctrl, &bgcin.FRZ, &bgcin.GRZ, &bgcin.HRV, &bgcin.MOW, &bgcin.PLT, &bgcin.PLG, &bgcin.THN, &bgcin.IRG, &bgcin.MUL, &bgcin.CWE);
	if (errorCode)
	{
		printf("ERROR in call to mgm_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}

	/* read simulation control flags */
	errorCode = simctrl_init(init, &bgcin.ctrl, &bgcin.epc, &bgcin.sprop, &bgcin.PLT);
	if (errorCode)
	{
		printf("ERROR in call to simctrl_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}

	/* initialize water state structure */
	errorCode = wstate_init(init, &bgcin.ctrl, &bgcin.sitec, &bgcin.sprop, &bgcin.ws);
	if (errorCode)
	{
		printf("ERROR in call to wstate_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}

	/* initialize carbon and nitrogen state structures */
	errorCode = cnstate_init(init, &bgcin.ctrl, &bgcin.epc, &bgcin.sprop, &bgcin.sitec, &bgcin.cs, &bgcin.cinit, &bgcin.ns);
	if (errorCode)
	{
		printf("ERROR in call to cstate_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(213);
		exit(213);
	}
  
	/* read scalar climate change parameters */
	errorCode = scc_init(init, &scc);
	if (errorCode)
	{
		printf("ERROR in call to scc_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}

	/* read conditional management strategies parameters */
	errorCode = conditionalMGM_init(init, &bgcin.ctrl, &bgcin.IRG, &bgcin.MOW);
	if (errorCode)
	{
		printf("ERROR in call to scc_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}

	/* read the output control information */
	if ((bgcin.co2.varco2 == 1 || bgcin.ndep.varndep == 1) && bgcin.ctrl.spinup == 1) transient = 1;
	errorCode = output_init(init, transient, &bgcin.HRV, &output);
	if (errorCode)
	{
		printf("ERROR in call to output_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}
	

	
	/* -------------------------------------------------------------------------*/

	/* read final line out of init file to test for proper file structure */
	errorCode = end_init(init);
	if (errorCode)
	{
		printf("ERROR in call to end_init.c from pointbgc.c... exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}
	fclose(init.ptr);

	/* read meteorology file, build metarr arrays, compute running avgs */
	errorCode = metarr_init(&point, &bgcin.metarr, &scc, &bgcin.sitec, &bgcin.ctrl);
	if (errorCode)
	{
		printf("ERROR in call to metarr_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}
	fclose(point.metf.ptr);

	/* read groundwater depth if it is available */
	errorCode = groundwater_init(&bgcin.GWS, &bgcin.ctrl);
	if (errorCode)
	{
		printf("ERROR in call to groundwater_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}

	/* read flooding height if it is available */
	errorCode = flooding_init(&bgcin.FLS,&bgcin.ctrl);
	if (errorCode)
	{
		printf("ERROR in call to flooding_init.c from pointbgc.c... Exiting\n");
		writeErrorCode(errorCode);
		exit(errorCode);
	}
	


	/* copy some of the info from input structure to bgc simulation control structure */
 	bgcin.ctrl.onscreen = output.onscreen;
	bgcin.ctrl.dodaily = output.dodaily;
	bgcin.ctrl.domonavg = output.domonavg;
	bgcin.ctrl.doannavg = output.doannavg;
	bgcin.ctrl.doannual = output.doannual;
	bgcin.ctrl.ndayout = output.ndayout;
	bgcin.ctrl.nannout = output.nannout;
	bgcin.ctrl.daycodes = output.daycodes;
	bgcin.ctrl.daynames = output.daynames;
	bgcin.ctrl.anncodes = output.anncodes;
	bgcin.ctrl.annnames = output.annnames;
	bgcin.ctrl.read_restart = restart.read_restart;
	bgcin.ctrl.write_restart = restart.write_restart;
	bgcin.ctrl.condIRG_flag = bgcin.IRG.condIRG_flag;  
	bgcin.ctrl.condMOW_flag = bgcin.MOW.condMOW_flag;   

	/* copy the output file structures into bgcout */
	bgcout.dayout  = output.dayout;
	bgcout.dayoutT = output.dayoutT;
	
	bgcout.monavgout  = output.monavgout;
	bgcout.monavgoutT = output.monavgoutT;
	
	bgcout.annavgout  = output.annavgout;
	bgcout.annavgoutT = output.annavgoutT;
	
	bgcout.annout  = output.annout;
	bgcout.annoutT = output.annoutT;

	bgcout.log_file = output.log_file;
	bgcout.econout_file = output.econout_file;
	
	
	
	/* if using an input restart file, read a record */
	if (restart.read_restart)
	{
		fread(&(bgcin.restart_input),sizeof(restart_data_struct),1,restart.in_restart.ptr);
	}

	/*********************
	**                  **
	**  CALL BIOME-BGC  **
	**                  **
	*********************/

	
	/* all initialization complete, call model */
	/* either call the spinup code or the normal simulation code */
	if (bgcin.ctrl.spinup)
	{
		errorCode = spinup_bgc(&bgcin, &bgcout);
	 	if (errorCode)
		{
			fprintf(output.log_file.ptr, "\n");
			fprintf(output.log_file.ptr, "ERROR in spinup run\n");
			fprintf(output.log_file.ptr, "error code: %i\n", errorCode);
			fprintf(output.log_file.ptr, "SIMULATION STATUS [0 - failure; 1 - success]\n");
			fprintf(output.log_file.ptr, "0\n");
			writeErrorCode(errorCode);
			exit(errorCode);
		}
		else
		{
			fprintf(output.log_file.ptr, "\n");
			fprintf(output.log_file.ptr, "SIMULATION STATUS [0 - failure; 1 - success]\n");
			fprintf(output.log_file.ptr, "1\n");
		}
	
	}
	else
	{   
		if (bgcin.ctrl.read_restart == 0)
		{
			fprintf(output.log_file.ptr, "\n");
			fprintf(output.log_file.ptr, "WARNING: in normal run: \n");
			fprintf(output.log_file.ptr, " read_restart=0 means only using W_STATE block of INI file instead of soilw walues of restart file \n");
			bgcin.ctrl.read_restart = 2;
		}
		errorCode = bgc(&bgcin, &bgcout);
		if (errorCode)
		{
			fprintf(output.log_file.ptr, "\n");
			fprintf(output.log_file.ptr, "ERROR in normal run\n");
			fprintf(output.log_file.ptr, "error code: %i\n", errorCode);
			fprintf(output.log_file.ptr, "\n");
			fprintf(output.log_file.ptr, "SIMULATION STATUS [0 - failure; 1 - success]\n");
			fprintf(output.log_file.ptr, "0\n");
			writeErrorCode(errorCode);
			exit(errorCode);
		}
		else
		{
			fprintf(output.log_file.ptr, "\n");
			fprintf(output.log_file.ptr, "SIMULATION STATUS [0 - failure; 1 - success]\n");
			fprintf(output.log_file.ptr, "1\n");
		}
	}
		

	/* if using an output restart file, write a record */
	if (restart.write_restart)
	{
		fwrite(&(bgcout.restart_output),sizeof(restart_data_struct),1,
			restart.out_restart.ptr);
	}

	/* post-processing output handling, if any, goes here */
	
	/* free memory */
	free(bgcin.metarr.Tmax_array);
	free(bgcin.metarr.Tmin_array);
	free(bgcin.metarr.prcp_array);
	free(bgcin.metarr.vpd_array);
	free(bgcin.metarr.Tday_array);
	free(bgcin.metarr.Tavg_array);
	free(bgcin.metarr.TavgRA11_array);
	free(bgcin.metarr.TavgRA30_array);
	free(bgcin.metarr.TavgRA10_array);
	free(bgcin.metarr.tempradF_array);
	free(bgcin.metarr.tempradFra_array);
	free(bgcin.metarr.swavgfd_array);
	free(bgcin.metarr.par_array);
	free(bgcin.metarr.dayl_array);
	free(bgcin.metarr.annTavg_array);
	free(bgcin.metarr.monTavg_array);
    free(bgcin.metarr.annTrange_array);
	free(bgcin.metarr.annTavgRA_array);
    free(bgcin.metarr.annTrangeRA_array);

	if(bgcin.PLT.PLT_num)
	{
		free(bgcin.PLT.PLTyear_array);  
		free(bgcin.PLT.PLTmonth_array); 
		free(bgcin.PLT.PLTday_array); 
		free(bgcin.PLT.germDepth_array); 
		free(bgcin.PLT.n_seedlings_array); 
		free(bgcin.PLT.weight_1000seed_array); 
		free(bgcin.PLT.seed_carbon_array); 
		free(bgcin.PLT.filename_array);
	}

	if (bgcin.THN.THN_num)
	{
		free(bgcin.THN.THNyear_array);  
		free(bgcin.THN.THNmonth_array); 
		free(bgcin.THN.THNday_array); 
		free(bgcin.THN.thinningRate_w_array); 
		free(bgcin.THN.thinningRate_nw_array); 
		free(bgcin.THN.transpCoeff_w_array); 
		free(bgcin.THN.transpCoeff_nw_array); 
	}

	if (bgcin.MOW.MOW_num)
	{
		free(bgcin.MOW.MOWyear_array);  
		free(bgcin.MOW.MOWmonth_array); 
		free(bgcin.MOW.MOWday_array); 
		free(bgcin.MOW.LAI_limit_array); 
		free(bgcin.MOW.transportMOW_array); 
	}

	if (bgcin.GRZ.GRZ_num)
	{
		free(bgcin.GRZ.GRZstart_year_array);  
		free(bgcin.GRZ.GRZstart_month_array); 
		free(bgcin.GRZ.GRZstart_day_array); 
		free(bgcin.GRZ.GRZend_year_array);  
		free(bgcin.GRZ.GRZend_month_array); 
		free(bgcin.GRZ.GRZend_day_array); 
		free(bgcin.GRZ.trampling_effect); 
		free(bgcin.GRZ.weight_LSU); 
		free(bgcin.GRZ.stocking_rate_array); 
		free(bgcin.GRZ.DMintake_array);
		free(bgcin.GRZ.DMintake2excr_array); 
		free(bgcin.GRZ.excr2litter_array); 
		free(bgcin.GRZ.DM_Ccontent_array); 
		free(bgcin.GRZ.EXCR_Ncontent_array);
		free(bgcin.GRZ.EXCR_Ccontent_array); 
		free(bgcin.GRZ.Nexrate); 
		free(bgcin.GRZ.EFman_N2O); 
		free(bgcin.GRZ.EFman_CH4);
		free(bgcin.GRZ.EFfer_CH4);
	}

	if (bgcin.HRV.HRV_num)
	{
		free(bgcin.HRV.HRVyear_array);  
		free(bgcin.HRV.HRVmonth_array); 
		free(bgcin.HRV.HRVday_array); 
		free(bgcin.HRV.propStemResidue_array); 
		free(bgcin.HRV.propRootResidue_array);
		free(bgcin.HRV.transportHRVleaf_array); 
		free(bgcin.HRV.transportHRVstem_array);
		free(bgcin.HRV.transportHRVyield_array);
	}

	if (bgcin.PLG.PLG_num)
	{
		free(bgcin.PLG.PLGyear_array);  
		free(bgcin.PLG.PLGmonth_array); 
		free(bgcin.PLG.PLGday_array); 
		free(bgcin.PLG.PLGdepths_array); 
	}

	if (bgcin.FRZ.FRZ_num)
	{
		free(bgcin.FRZ.FRZyear_array);  
		free(bgcin.FRZ.FRZmonth_array); 
		free(bgcin.FRZ.FRZday_array); 
		free(bgcin.FRZ.FRZdepth_array); 
		free(bgcin.FRZ.fertilizer_array); 
		free(bgcin.FRZ.DM_array); 
		free(bgcin.FRZ.NO3content_array); 
		free(bgcin.FRZ.NH4content_array); 
		free(bgcin.FRZ.UREAcontent_array); 
		free(bgcin.FRZ.orgNcontent_array);
		free(bgcin.FRZ.orgCcontent_array); 
		free(bgcin.FRZ.litr_flab_array); 
		free(bgcin.FRZ.litr_fcel_array); 
		free(bgcin.FRZ.EFfert_N2O);
	}

	if (bgcin.IRG.IRG_num)
	{
		free(bgcin.IRG.IRGyear_array);  
		free(bgcin.IRG.IRGmonth_array); 
		free(bgcin.IRG.IRGday_array); 
		free(bgcin.IRG.IRGquantity_array);
		free(bgcin.IRG.IRGheight_array); 
	}
	
	if (bgcin.MUL.MUL_num)
	{
		free(bgcin.MUL.MULyear_array);  
		free(bgcin.MUL.MULmonth_array); 
		free(bgcin.MUL.MULday_array); 
		free(bgcin.MUL.litrCabove_MUL);
		free(bgcin.MUL.litrCNabove_MUL); 
		free(bgcin.MUL.cwdCabove_MUL);
		free(bgcin.MUL.cwdCNabove_MUL); 
	}

	if (bgcin.CWE.CWE_num)
	{
		free(bgcin.CWE.CWEyear_array);  
		free(bgcin.CWE.CWEmonth_array); 
		free(bgcin.CWE.CWEday_array); 
		free(bgcin.CWE.removePROP_CWE);
	}

	if (bgcin.GWS.GWnum)
	{
		free(bgcin.GWS.GWyear_array);	
		free(bgcin.GWS.GWmonth_array);	
		free(bgcin.GWS.GWday_array);	
        free(bgcin.GWS.GWdepth_array);	
		free(bgcin.GWS.GW_NH4ppm_array);
		free(bgcin.GWS.GW_NO3ppm_array);
		free(bgcin.GWS.GW_DOC1ppm_array);
		free(bgcin.GWS.GW_DOC2ppm_array);
		free(bgcin.GWS.GW_DOC3ppm_array);
		free(bgcin.GWS.GW_DOC4ppm_array);
		free(bgcin.GWS.GW_DON1ppm_array);
		free(bgcin.GWS.GW_DON2ppm_array);
		free(bgcin.GWS.GW_DON3ppm_array);
		free(bgcin.GWS.GW_DON4ppm_array);
	}

	if (bgcin.FLS.FLnum)
	{
		free(bgcin.FLS.FLstart_year_array);	
		free(bgcin.FLS.FLstart_month_array);	
		free(bgcin.FLS.FLstart_day_array);	
		free(bgcin.FLS.FLend_year_array);	
		free(bgcin.FLS.FLend_month_array);	
		free(bgcin.FLS.FLend_day_array);
        free(bgcin.FLS.FLheight_array);
		free(bgcin.FLS.FL_NH4ppm_array);
		free(bgcin.FLS.FL_NO3ppm_array);
		free(bgcin.FLS.FL_DOC1ppm_array);
		free(bgcin.FLS.FL_DOC2ppm_array);
		free(bgcin.FLS.FL_DOC3ppm_array);
		free(bgcin.FLS.FL_DOC4ppm_array);
		free(bgcin.FLS.FL_DON1ppm_array);
		free(bgcin.FLS.FL_DON2ppm_array);
		free(bgcin.FLS.FL_DON3ppm_array);
		free(bgcin.FLS.FL_DON4ppm_array);
	}

	if (bgcin.co2.varco2) free(bgcin.co2.co2ppm_array);
	if (bgcin.co2.varco2) free(bgcin.co2.co2yrs_array);
	if (bgcin.ndep.varndep) free(bgcin.ndep.Ndep_array);
	if (bgcin.ndep.varndep) free(bgcin.ndep.Nyrs_array);

	if (bgcin.ctrl.varSGS_flag) free(bgcin.epc.SGS_array);
	if (bgcin.ctrl.varEGS_flag) free(bgcin.epc.EGS_array);
	if (bgcin.ctrl.varWPM_flag) free(bgcin.epc.WPMyr_array);
	if (bgcin.ctrl.varFM_flag) free(bgcin.epc.FMyr_array);
	if (bgcin.ctrl.varMSC_flag) free(bgcin.epc.MSC_array);

	if (output.ndayout != 0) free(output.daycodes);
	if (output.ndayout != 0) free(output.daynames);
	if (output.nannout != 0) free(output.anncodes);
	if (output.nannout != 0) free(output.annnames);
	
	/* close files */
	if (restart.read_restart) fclose(restart.in_restart.ptr);
	if (restart.write_restart) fclose(restart.out_restart.ptr);
	if (output.dodaily) 
	{
		fclose(output.dayout.ptr);
		if (transient) fclose(output.dayoutT.ptr);
	}
	if (output.domonavg) 
	{
		fclose(output.monavgout.ptr);
		if (transient) fclose(output.monavgoutT.ptr);
	}
	if (output.doannavg) 
	{
		fclose(output.annavgout.ptr);
		if (transient) fclose(output.annavgoutT.ptr);
	}
	if (output.doannual) 
	{
		fclose(output.annout.ptr);
		if (transient) fclose(output.annoutT.ptr);
	}

	fclose(output.log_file.ptr);
	if (bgcin.HRV.HRV_num) fclose(output.econout_file.ptr);

	return (errorCode); 

    /* end of main */	
 } 
	
