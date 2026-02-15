/* 
simctrl_init.c
read simulation control flags for pointbgc simulation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2025, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
20.03.01 Galina Churkina added variable "sum" substituting  t1+t2+t3 in IF statement,
which gave an error.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ini.h"
#include "bgc_struct.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"



int simctrl_init(file init, control_struct* ctrl, epconst_struct* epc, soilprop_struct* sprop, planting_struct* PLT)
{
	int errorCode=0;
	int layer;
	char key[] = "SIMULATION_CONTROL";
	char keyword[STRINGSIZE];
	double hydrCONDUCTsat_cmday;

	
	/********************************************************************
	**                                                                 **
	** Begin reading initialization file block starting with keyword:  **
	** SIMULATION_CONTROL                                                      ** 
	**                                                                 **
	********************************************************************/
	
	
	/* scan for the SOI file keyword, exit if not next */
	if (!errorCode && scan_value(init, keyword, 's'))
	{
		printf("ERROR reading keyword for control data\n");
		errorCode=211;
	}
	if (!errorCode && strcmp(keyword, key))
	{
		printf("Expecting keyword --> %s in file %s\n",key,init.name);
		errorCode=211;
	}


	if (!errorCode && scan_value(init, &ctrl->phenology_flag, 'i'))
	{
		printf("ERROR reading phenology flag, simctrl_init.c\n");
		errorCode=21101;
	}

	/* get flag of GSI flag */
	if (!errorCode && scan_value(init, &ctrl->GSI_flag, 'i'))
	{
		printf("ERROR reading flag indicating usage of GSI file: simctrl_init.c\n");
		errorCode=21102;
	}
	
	if (!errorCode && scan_value(init, &ctrl->transferGDD_flag, 'i'))
	{
		printf("ERROR reading transferGDD_flag, simctrl_init.c\n");
		errorCode=21103;
	}

	/* control of phenophase number */
	if (!errorCode)
	{
		if (ctrl->transferGDD_flag && epc->n_emerg_phenophase < 1)
		{
			printf("ERROR in phenophase parametrization: if transferGDD_flag = 1 -> n_emerg_phenophase must be specified in EPC file.c\n");
			errorCode=2110301;
		}
	}

	/* temperature dependent q10 value */
	if (!errorCode && scan_value(init, &ctrl->q10depend_flag, 'i'))
	{
		printf("ERROR reading q10depend_flag, simctrl_init.c\n");
		errorCode=21104;
	}
	
	/* acclimation */
	if (!errorCode && scan_value(init, &ctrl->phtsyn_acclim_flag, 'i'))
	{
		printf("ERROR reading acclimation flag of photosynthesis, simctrl_init.c\n");
		errorCode=21105;
	}
	if (!errorCode && scan_value(init, &ctrl->resp_acclim_flag, 'i'))
	{
		printf("ERROR reading acclimation flag of respiration, simctrl_init.c\n");
		errorCode=21106;
	}
	
	/* get flag of CO2 conductance reduction */
	if (!errorCode && scan_value(init, &ctrl->CO2conduct_flag, 'i'))
	{
		printf("ERROR reading CO2conduct_flag, simctrl_init.c\n");
		errorCode=21107;
	}

	/* soil temperature calculation flag */
	if (!errorCode && scan_value(init, &ctrl->STCM_flag, 'i'))
	{
		printf("ERROR reading soil temperature calculation flag: simctrl_init.c\n");
		errorCode=21108;
	}
	
	
	/*  photosynthesis calculation method flag */
	if (!errorCode && scan_value(init, &ctrl->photosynt_flag, 'i'))
	{
		printf("ERROR reading photosynthesis calculation method flag: simctrl_init.c\n");
		errorCode=21109;
	}


	/*  evapotranspiration calculation method flag */
	if (!errorCode && scan_value(init, &ctrl->ET_flag, 'i'))
	{
		printf("ERROR reading evapotranspiration calculation method flag: simctrl_init.c\n");
		errorCode=21110;
	}

	if (ctrl->ET_flag && epc->PT_ETcritT == DATA_GAP)
	{
		printf("ERROR in evaporatranspiration calculation: if ET_flag = 1 -> PT_ETcritT must be specified in EPC file.c\n");
		errorCode=2111001;
	}

	/*  radiation calculation method flag */
	if (!errorCode && scan_value(init, &ctrl->radiation_flag, 'i'))
	{
		printf("ERROR reading radiation calculation method flag: simctrl_init.c\n");
		errorCode=21111;
	}

	/*  soilstress calculation method flag */
	if (!errorCode && scan_value(init, &ctrl->soilstress_flag, 'i'))
	{
		printf("ERROR reading soilstress calculation method flag: simctrl_init.c\n");
		errorCode=21112;
	}

	/*  interception calculation method flag */
	if (!errorCode && scan_value(init, &ctrl->interception_flag, 'i'))
	{
		printf("ERROR reading interception calculation method flag: simctrl_init.c\n");
		errorCode=21113;
	}

	/*  MR-deficit calculation method flag */
	if (!errorCode && scan_value(init, &ctrl->MRdeficit_flag, 'i'))
	{
		printf("ERROR reading MR-deficit calculation method flag: simctrl_init.c\n");
		errorCode=21114;
	}

	/* soil water calculation flag */
	if (!errorCode && scan_value(init, &ctrl->Ksat_flag, 'i'))
	{
		printf("ERROR reading Ksat estimation method flag: simctrl_init.c\n");
		errorCode = 21115;
	}


	/* control: in case of planting/harvesting, model-defined phenology is not possible: first day - planting day, last day - harvesting day */
	if (ctrl->phenology_flag == 1 && PLT->PLT_num && ctrl->spinup != 1) 
	{
		ctrl->prephen1_flag = 1;
		ctrl->phenology_flag = 0;
	}

	/* control: in case of user-defined phenology, GSI-method is not possible */
	if (ctrl->phenology_flag == 0 && ctrl->GSI_flag) 
	{
		ctrl->prephen2_flag = 1;
		ctrl->GSI_flag = 0;
	}
	
	/* if measured data is availabe (in cm/day) */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (sprop->hydrCONTUCTsatMES_cmPERday[layer] == (double)DATA_GAP)
		{
			if (ctrl->Ksat_flag == 1)
				sprop->hydrCONDUCTsat[layer] = 7.05556 * 1e-6 * pow(10, (-0.6 + 0.0126 * sprop->sand[layer] - 0.0064 * sprop->clay[layer]));
			else
				sprop->hydrCONDUCTsat[layer] = (50 * exp(-0.075 * sprop->clay[layer]) + 200 * exp(-0.075 * (100 - sprop->sand[layer]))) / 100 / nSEC_IN_DAY;
		}
		else
			sprop->hydrCONDUCTsat[layer] = sprop->hydrCONTUCTsatMES_cmPERday[layer] / (m_to_cm * nSEC_IN_DAY);

		sprop->hydrDIFFUSsat[layer] = (sprop->soilB[layer] * sprop->hydrCONDUCTsat[layer] * (-100 * sprop->PSIsat[layer])) / sprop->VWCsat[layer];
		sprop->hydrCONDUCTfc[layer] = sprop->hydrCONDUCTsat[layer] * pow(sprop->VWCfc[layer] / sprop->VWCsat[layer], 2 * sprop->soilB[layer] + 3);
		sprop->hydrDIFFUSfc[layer] = (((sprop->soilB[layer] * sprop->hydrCONDUCTsat[layer] * (-100 * sprop->PSIsat[layer]))) / sprop->VWCsat[layer])
			* pow(sprop->VWCfc[layer] / sprop->VWCsat[layer], sprop->soilB[layer] + 2);

		hydrCONDUCTsat_cmday = sprop->hydrCONDUCTsat[layer] * m_to_cm * nSEC_IN_DAY; // saturated hydraulic conductivity (cm/day = m/s * 100 * sec/day)
		if (sprop->drainCoeff_mes[layer] != DATA_GAP)
			sprop->drainCoeff[layer] = sprop->drainCoeff_mes[layer];
		else
			sprop->drainCoeff[layer] = 0.1122 * pow(hydrCONDUCTsat_cmday, 0.339);

	}



  return (errorCode);
  
}