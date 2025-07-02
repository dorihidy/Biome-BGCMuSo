/* 
multilayer_sminn.c
Calculating the content of soil mineral nitrogen in multilayer soil (plant N upate, soil processes, nitrification and denitrification, depostion and fixing). 

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2025, D. Hidy [dori.hidy@gmail.com]
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

int multilayer_sminn(control_struct* ctrl, const metvar_struct* metv,const siteconst_struct* sitec, const NdepControl_struct* ndep, 
	                 cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf, soilprop_struct* sprop, epvar_struct* epv, soilInfo_struct* soilInfo)
{
	int errorCode=0;
	int layer=0;
	int GWlayer, CFlayer, dm;
	double NH4_prop,SR_layer, NO3avail, NO3avail_ppm;
	double pH, Tsoil, WFPS, net_miner,NH4dissolv, N2OfluxNITRIF,NH4_to_nitrif;
	double weight;
	double sminn_layer[N_SOILLAYERS];
	double sminn_to_soilCTRL, sminn_to_npoolCTRL, ndep_to_sminnCTRL, nfix_to_sminnCTRL;
	double SR_total,NO3_to_denitr,ratioN2_N2O;
	double NH4_to_soilSUM_total, NO3_to_soilSUM_total;
	double ratioNORM, ratioCAPIL, ratioSAT, diffNORM, diffCAPIL, diffSAT, diff, VWCcritDENIT;
	
	NH4_prop=net_miner=SR_layer=sminn_to_soilCTRL=sminn_to_npoolCTRL=ndep_to_sminnCTRL=nfix_to_sminnCTRL= NH4_to_soilSUM_total=NO3_to_soilSUM_total= VWCcritDENIT=0;



	GWlayer = (int)sprop->GWlayer;
	CFlayer = (int)sprop->CFlayer;
	/* *****************************************************************************************************

	1.Deposition and fixation - INPUT
	2.Plant N uptake from SMINN: sminn_to_npool 
	3.Immobilization-mineralization: sminn_to_soilSUM - due microbial soil processes SMINN is changing in the soil (determined in allocation).
	4. Denitrification and Nitrification */


	/* SR_total calculation - unit: kgC/ha */
	SR_total = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++) 
	{  
		SR_total += (cf->soil1_hr[layer] + cf->soil2_hr[layer] + cf->soil3_hr[layer] + cf->soil4_hr[layer]) * kg_to_g;
		for (dm = 0; layer < N_DISSOLVN; layer++) soilInfo->dismatTOTALecofunc[dm][layer] = 0;

	}
	/*-----------------------------------------------------------------------------*/
	/* Calculations layer by layer */

	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		/*********************************************/
		/* 1. Deposition and fixation */
		/*************************************/

		/* calculation of sminn and NH4 prop in soil layers */

		NH4dissolv = ns->NH4[layer] * sprop->NH4_mobilen_prop;
		sminn_layer[layer] = (NH4dissolv + ns->NO3[layer]);

		if (sminn_layer[layer])
			NH4_prop = NH4dissolv / sminn_layer[layer];
		else
			NH4_prop = ndep->NdepNH4_coeff;


		/* Deposition  only in 0-30 cm */
		if (layer < 3)
		{
			weight = sitec->soillayer_thickness[layer] / sitec->soillayer_depth[2];
			nf->ndep_to_NH4[layer] = (nf->ndep_to_sminn_total * weight) * ndep->NdepNH4_coeff;
			nf->ndep_to_NO3[layer] = (nf->ndep_to_sminn_total * weight) * (1 - ndep->NdepNH4_coeff);
		}
		else
		{
			nf->ndep_to_NH4[layer] = 0;
			nf->ndep_to_NO3[layer] = 0;
		}

		/* N-fixation in rootlayer (based on rootlength proportion) or in absence of root: top 30 cm */
		if (epv->n_rootlayers == 0)
		{
			if (layer < 3)
			{
				weight = sitec->soillayer_thickness[layer] / sitec->soillayer_depth[2];
				nf->nfix_to_NH4[layer] = nf->nfix_to_sminn_total * weight;
			}
			else
				nf->nfix_to_NH4[layer] = 0;
		}
		else
		{
			nf->nfix_to_NH4[layer] = nf->nfix_to_sminn_total * epv->rootlengthProp[layer];
			/* no root is assumed in the top soil layer, but N fixation should occur also in top soil layer  */
			if (epv->rootlengthProp[0] == 0 && layer < 2)
			{
				nf->nfix_to_NH4[layer] = nf->nfix_to_sminn_total * epv->rootlengthProp[1] * (sitec->soillayer_thickness[layer] / sitec->soillayer_depth[1]);
			}
			else
				nf->nfix_to_NH4[layer] = nf->nfix_to_sminn_total * epv->rootlengthProp[layer];
		}

		nf->environment_to_sminn[layer] = nf->ndep_to_NH4[layer] + nf->ndep_to_NO3[layer] + nf->nfix_to_NH4[layer];

		ndep_to_sminnCTRL += nf->ndep_to_NH4[layer] + nf->ndep_to_NO3[layer];
		nfix_to_sminnCTRL += nf->nfix_to_NH4[layer];

		/*********************************************/
		/* 2. Immobilization-mineralization fluxes */

		nf->NH4_to_soilSUM[layer] = nf->sminn_to_soilSUM[layer] * NH4_prop;
		nf->NO3_to_soilSUM[layer] = nf->sminn_to_soilSUM[layer] * (1 - NH4_prop);

		nf->NH4_to_soilSUM_total += nf->NH4_to_soilSUM[layer];
		nf->NO3_to_soilSUM_total += nf->NO3_to_soilSUM[layer];

		sminn_to_soilCTRL += nf->NH4_to_soilSUM[layer] + nf->NO3_to_soilSUM[layer];

		/*********************************************/
		/* 3. Plant N uptake  in rootzone */

		nf->NH4_to_npool[layer] = nf->sminn_to_npool[layer] * NH4_prop;
		nf->NO3_to_npool[layer] = nf->sminn_to_npool[layer] * (1 - NH4_prop);
		sminn_to_npoolCTRL += nf->NH4_to_npool[layer] + nf->NO3_to_npool[layer];

		/*********************************************/
		/* 4. Nitrification */

		if (nf->NH4_to_soilSUM[layer] > 0)
			net_miner = nf->sminn_to_soiln_s4[layer];
		else
			net_miner = nf->sminn_to_soiln_s4[layer] - nf->NH4_to_soilSUM[layer];

		WFPS = epv->WFPS[layer];
		NH4dissolv = ns->NH4[layer] * sprop->NH4_mobilen_prop;
		Tsoil = metv->Tsoil[layer];
		pH = sprop->pH[layer];

		if (!errorCode && nitrification(layer, sprop, net_miner, Tsoil, pH, WFPS, NH4dissolv, epv, &N2OfluxNITRIF, &NH4_to_nitrif))
		{
			printf("\n");
			printf("ERROR in nitrification.c for multilayer_sminn.c \n");
			errorCode = 1;
		}


		nf->N2OfluxNITRIF[layer] = N2OfluxNITRIF;
		nf->NH4_to_nitrif[layer] = NH4_to_nitrif;

		/*********************************************/

		NO3avail = ns->NO3[layer] - nf->NO3_to_npool[layer];
		if (nf->NO3_to_soilSUM[layer] > 0) NO3avail -= nf->NO3_to_soilSUM[layer];

		NO3avail_ppm = ns->NO3[layer] / (sprop->BD[layer]  * sitec->soillayer_thickness[layer]) * multi_ppm;


		SR_layer = (cf->soil1_hr[layer] + cf->soil2_hr[layer] + cf->soil3_hr[layer] + cf->soil4_hr[layer]) * kg_to_g;
		
		if (sprop->critWFPS_denitr == DATA_GAP)
			VWCcritDENIT = sprop->VWCfc[layer];
		else
			VWCcritDENIT = sprop->VWCsat[layer] * sprop->critWFPS_denitr;

		if (!errorCode && denitrification(ctrl->soiltype_array[layer], NO3avail_ppm, pH, WFPS, SR_layer, &NO3_to_denitr, &ratioN2_N2O))
		{
			printf("\n");
			printf("ERROR in denitrification.c for multilayer_sminn.c \n");
			errorCode = 1;
		}


		if (epv->VWC[layer] > VWCcritDENIT)
			nf->NO3_to_denitr[layer] = sprop->denitr_coeff * SR_layer * NO3avail * (epv->VWC[layer] - VWCcritDENIT) / (sprop->VWCsat[layer] - VWCcritDENIT);
		else
			nf->NO3_to_denitr[layer] = 0;

		if (nf->NO3_to_denitr[layer] > NO3avail)
		{
			nf->NO3_to_denitr[layer] = NO3avail;
			ctrl->limitDENIT_flag = 1;
		}


		nf->N2OfluxDENITR[layer] = nf->NO3_to_denitr[layer] / (1 + ratioN2_N2O * sprop->N2Oratio_denitr);
		nf->N2fluxDENITR[layer] = nf->NO3_to_denitr[layer] - nf->N2OfluxDENITR[layer];



		/*********************************************/
		/* 7. STATE UPDATE */

		/* NH4 */
		soilInfo->dismatTOTALecofunc[0][layer] = (nf->ndep_to_NH4[layer] + nf->nfix_to_NH4[layer] +
			               (nf->litr1n_to_release[layer] + nf->litr2n_to_release[layer] + nf->litr4n_to_release[layer]) * 0.5 -
			                nf->NH4_to_soilSUM[layer] - nf->NH4_to_npool[layer] - nf->NH4_to_nitrif[layer]);
		/* NO3 */
		soilInfo->dismatTOTALecofunc[1][layer] = (nf->ndep_to_NO3[layer] + (nf->NH4_to_nitrif[layer] - nf->N2OfluxNITRIF[layer]) +
			               (nf->litr1n_to_release[layer] + nf->litr2n_to_release[layer] + nf->litr4n_to_release[layer]) * 0.5 -
			                nf->NO3_to_soilSUM[layer] - nf->NO3_to_npool[layer] - nf->NO3_to_denitr[layer]);


		/* soilN */
		soilInfo->dismatTOTALecofunc[2][layer] = nf->sminn_to_soil1n_l1[layer];
		soilInfo->dismatTOTALecofunc[3][layer] = nf->sminn_to_soil2n_l2[layer] + nf->sminn_to_soil2n_s1[layer];
		soilInfo->dismatTOTALecofunc[4][layer] = nf->sminn_to_soil3n_l4[layer] + nf->sminn_to_soil3n_s2[layer];
		soilInfo->dismatTOTALecofunc[5][layer] = nf->sminn_to_soil4n_s3[layer] + nf->sminn_to_soiln_s4[layer];

		if (layer < GWlayer)
		{
			soilInfo->dismatUNSATecofunc[0][layer] = soilInfo->dismatTOTALecofunc[0][layer];
			soilInfo->dismatUNSATecofunc[1][layer] = soilInfo->dismatTOTALecofunc[1][layer];
			soilInfo->dismatUNSATecofunc[2][layer] = soilInfo->dismatTOTALecofunc[2][layer];
			soilInfo->dismatUNSATecofunc[3][layer] = soilInfo->dismatTOTALecofunc[3][layer];
			soilInfo->dismatUNSATecofunc[4][layer] = soilInfo->dismatTOTALecofunc[4][layer];
			soilInfo->dismatUNSATecofunc[5][layer] = soilInfo->dismatTOTALecofunc[5][layer];
		}

		/* state update */
		ns->NH4[layer] += soilInfo->dismatTOTALecofunc[0][layer];
		ns->NO3[layer] += soilInfo->dismatTOTALecofunc[1][layer];

		ns->soil1n[layer] += soilInfo->dismatTOTALecofunc[2][layer];
		ns->soil2n[layer] += soilInfo->dismatTOTALecofunc[3][layer];
		ns->soil3n[layer] += soilInfo->dismatTOTALecofunc[4][layer];
		ns->soil4n[layer] += soilInfo->dismatTOTALecofunc[5][layer];

		ns->npool += (nf->NH4_to_npool[layer] + nf->NO3_to_npool[layer]);
		ns->Nfix_src += nf->nfix_to_NH4[layer];
		ns->Ndep_src += nf->ndep_to_NH4[layer] + nf->ndep_to_NO3[layer];

		ns->Nvol_snk += nf->N2OfluxNITRIF[layer];
		ns->Nvol_snk += nf->NO3_to_denitr[layer];


		/*********************************************/
		/* 9. Total flux calculation */

		nf->N2fluxDENITR_total += nf->N2fluxDENITR[layer];
		nf->N2OfluxDENITR_total += nf->N2OfluxDENITR[layer];
		nf->N2OfluxNITRIF_total += nf->N2OfluxNITRIF[layer];

		nf->NO3_to_denitr_total += nf->NO3_to_denitr[layer];
		nf->NH4_to_nitrif_total += nf->NH4_to_nitrif[layer];

		nf->NH4_to_npool_total += nf->NH4_to_npool[layer];
		nf->NO3_to_npool_total += nf->NO3_to_npool[layer];

		nf->environment_to_NH4_total += nf->ndep_to_NH4[layer] + nf->nfix_to_NH4[layer];
		nf->environment_to_NO3_total += nf->ndep_to_NO3[layer];



	}

	if (!errorCode && (fabs(sminn_to_soilCTRL - nf->NH4_to_soilSUM_total - nf->NO3_to_soilSUM_total) > CRIT_PREC || fabs(sminn_to_npoolCTRL - nf->sminn_to_npool_total) > CRIT_PREC ||
		fabs(ndep_to_sminnCTRL - nf->ndep_to_sminn_total) > CRIT_PREC || fabs(nfix_to_sminnCTRL - nf->nfix_to_sminn_total) > CRIT_PREC))
	{
		printf("\n");
		printf("ERROR: in calculation of nitrogen state update (multilayer_sminn.c)\n");
		errorCode=1;
	}


	
	/*********************************************/
	/* II. CONTROL */

	/* transfer value: NH4, NO3, DOC, DON - > content_array  etc.*/
	if (!errorCode && check_soilcontent(-1, 0, sprop, cs, ns, soilInfo))
	{
		printf("ERROR in check_soilcontent.c for multilayer_sminn.c\n");
		errorCode = 1;
	}


	/******************************************************************************************/
	/* III. Groundwater transport */


	if (sprop->GWlayer != DATA_GAP)
	{

		/* if capillary zone exists in unsaturated zone (not in GWlayer) */
		if (sprop->dz_CAPILcf)
		{

			for (dm = 0; dm < N_DISSOLVN; dm++)
			{
				/* calculation of NORM and CAPIL ratio */
				if (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm])
				{
					ratioNORM = soilInfo->content_NORMcf[dm] / (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm]);
					ratioCAPIL = soilInfo->content_CAPILcf[dm] / (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm]);
				}
				else
				{
					ratioNORM = sprop->dz_NORMcf / (sprop->dz_NORMcf + sprop->dz_CAPILcf);
					ratioCAPIL = sprop->dz_CAPILcf / (sprop->dz_NORMcf + sprop->dz_CAPILcf);
				}
				if (fabs(1 - ratioNORM - ratioCAPIL) > CRIT_PREC)
				{
					printf("\n");
					printf("ERROR in ratio calculation in multilayer_sminn.c\n");
					errorCode = 1;
				}
		
				diffNORM = soilInfo->dismatTOTALecofunc[dm][CFlayer] * ratioNORM;
				diffCAPIL = soilInfo->dismatTOTALecofunc[dm][CFlayer] * ratioCAPIL;

				/* NORM zone: negative storage value is not possible - covered by GW */
				if (diffNORM + soilInfo->content_NORMcf[dm] < 0) diffNORM = -1 * soilInfo->content_NORMcf[dm];

				soilInfo->content_NORMcf[dm] += diffNORM;

				/* CAPIL zone: negative storage value is not possible - covered by GW */
				diffCAPIL = soilInfo->dismatTOTALecofunc[dm][CFlayer] - diffNORM;


				if (soilInfo->content_CAPILcf[dm] + diffCAPIL < 0)
				{
					if (fabs(soilInfo->content_CAPILcf[dm] + diffCAPIL) > CRIT_PREC)
					{
						printf("ERROR in content_CAPILgw calculation in multilayer_sminn.c\n");
						errorCode = 1;
					}
					else
						diffCAPIL = soilInfo->content_CAPILcf[dm];
				}
				soilInfo->content_CAPILcf[dm] += diffCAPIL;

				if (fabs(soilInfo->content_soil[dm][CFlayer] - (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm])) > CRIT_PREC_lenient*10)
				{
					printf("ERROR in content_CAPILcf calculation in multilayer_sminn.c\n");
			        errorCode = 1;
				}

			}
		}
		for (layer = GWlayer; layer < N_SOILLAYERS; layer++)
		{
			/* in GWlayer: soilInfo->dismatTOTALecofunc is distributed proportionally - soilInfo->dismatTOTALecofunc in saturated zone is accounted for  GWecofunc: -: sink, +: source*/

			if (layer == GWlayer)
			{


				for (dm = 0; dm < N_DISSOLVN; dm++)
				{
					if (soilInfo->content_NORMgw[dm] + soilInfo->content_CAPILgw[dm] + soilInfo->content_SATgw[dm])
					{
						ratioNORM = soilInfo->content_NORMgw[dm] / (soilInfo->content_NORMgw[dm] + soilInfo->content_CAPILgw[dm] + soilInfo->content_SATgw[dm]);
						ratioCAPIL = soilInfo->content_CAPILgw[dm] / (soilInfo->content_NORMgw[dm] + soilInfo->content_CAPILgw[dm] + soilInfo->content_SATgw[dm]);
						ratioSAT = soilInfo->content_SATgw[dm] / (soilInfo->content_NORMgw[dm] + soilInfo->content_CAPILgw[dm] + soilInfo->content_SATgw[dm]);
					}
					else
					{ 

						ratioNORM = sprop->dz_NORMgw / (sprop->dz_CAPILgw + sprop->dz_NORMgw + sprop->dz_SATgw);
						ratioCAPIL = sprop->dz_CAPILgw / (sprop->dz_CAPILgw + sprop->dz_NORMgw + sprop->dz_SATgw);
						ratioSAT = sprop->dz_SATgw / (sprop->dz_CAPILgw + sprop->dz_NORMgw + sprop->dz_SATgw);
					}
					if (fabs(1 - ratioNORM - ratioCAPIL - ratioSAT) > CRIT_PREC)
					{
						printf("\n");
						printf("ERROR in ratio calculation in multilayer_sminn.c\n");
						errorCode = 1;
					}

					diffNORM = soilInfo->dismatTOTALecofunc[dm][layer] * ratioNORM;
					diffCAPIL = soilInfo->dismatTOTALecofunc[dm][layer] * ratioCAPIL;
					diffSAT = soilInfo->dismatTOTALecofunc[dm][layer] * ratioSAT;


					/* NORM zone: negative storage value is not possible - covered by GW */
					if (soilInfo->content_NORMgw[dm] + diffNORM < 0)
						diff = -1 * diffNORM - soilInfo->content_NORMgw[dm];
					else
						diff = 0;

					soilInfo->dismatGWecofunc_NORM[dm] = diff;
					soilInfo->content_NORMgw[dm] += diffNORM + soilInfo->dismatGWecofunc_NORM[dm];

					/* CAPIL zone: negative storage value is not possible - covered by GW */
					if (soilInfo->content_CAPILgw[dm] + diffCAPIL < 0)
						diff = -1 * diffCAPIL - soilInfo->content_CAPILgw[dm];
					else
						diff = 0;

					soilInfo->dismatGWecofunc_CAPIL[dm] = diff;
					soilInfo->content_CAPILgw[dm] += diffCAPIL + soilInfo->dismatGWecofunc_CAPIL[dm];

					/* SAT zone: no soilInfo->dismatTOTALecofunc is possible (constant conc) - covered by GW  */
					soilInfo->dismatGWecofunc[dm][layer] = -1 * diffSAT;
					soilInfo->content_SATgw[dm] += diffSAT + soilInfo->dismatGWecofunc[dm][layer];

	
			
					/* UNSAT soilInfo->dismatTOTALecofunc */
					soilInfo->dismatUNSATecofunc[dm][layer] = diffNORM + diffCAPIL;

	                soilInfo->content_soil[dm][layer] = soilInfo->content_NORMgw[dm] + soilInfo->content_CAPILgw[dm] + soilInfo->content_SATgw[dm];
				}

			}
			else
			{
				for (dm = 0; dm < N_DISSOLVN; dm++)
				{
					soilInfo->dismatGWecofunc[dm][layer] = -1 * soilInfo->dismatTOTALecofunc[dm][layer];
					soilInfo->content_soil[dm][layer]   += soilInfo->dismatGWecofunc[dm][layer];
				}
			}

			/* transfer value: content_array ->NH4, NO3, DOC, DON  etc.*/
			if (!errorCode && check_soilcontent(layer, 1, sprop, cs, ns, soilInfo))
			{
				printf("ERROR in check_soilcontent.c for multilayer_sminn.c\n");
				errorCode = 1;
			}

		}

		/* src/snk variables*/
		for (dm = 0; dm < N_DISSOLVN; dm++)
		{
			if (soilInfo->dismatGWecofunc_NORM[dm] > 0)
				ns->GWsrc_N += soilInfo->dismatGWecofunc_NORM[dm];
			else
				ns->GWsnk_N += -1 * soilInfo->dismatGWecofunc_NORM[dm];

			if (soilInfo->dismatGWecofunc_CAPIL[dm] > 0)
				ns->GWsrc_N += soilInfo->dismatGWecofunc_CAPIL[dm];
			else
				ns->GWsnk_N += -1 * soilInfo->dismatGWecofunc_CAPIL[dm];

			for (layer = 0; layer < N_SOILLAYERS; layer++)
			{
				if (soilInfo->dismatGWecofunc[dm][layer] > 0)
					ns->GWsrc_N += soilInfo->dismatGWecofunc[dm][layer];
				else
					ns->GWsnk_N += -1 * soilInfo->dismatGWecofunc[dm][layer];
			}
		}


		
		
	}



	return (errorCode);
}

int nitrification(int layer, const soilprop_struct* sprop, double net_miner, double Tsoil, double pH, double WFPS, double NH4dissolv, 
	                  epvar_struct* epv, double* N2OfluxNITRIF, double* NH4_to_nitrif)
{
	int errorCode = 0;
	double NH4_to_nit,N2O_flux_NIT;
	
	/* calculation of scalar functions: Tsoil response function, ps_nitrif and WFPS_scalar */

	epv->ps_nitrif[layer] = sprop->pHp2_nitrif + (sprop->pHp1_nitrif-sprop->pHp2_nitrif)/(1 + exp((pH-sprop->pHp3_nitrif)/sprop->pHp4_nitrif));	

				
	if (sprop->Tp1_nitrif == DATA_GAP)
	{
		/* no decomp processes for Tsoil < -10.0 C */
		if (Tsoil < sprop->Tmin_decomp)	
				epv->ts_nitrif[layer] = 0.0;
		else
			epv->ts_nitrif[layer] = exp(sprop->Tp2_nitrif*((1.0/sprop->Tp3_nitrif)-(1.0/((Tsoil+Celsius2Kelvin)-sprop->Tp4_nitrif))));
	}
	else
	{
		/* no decomp processes for Tsoil < -10.0 C */
		if (Tsoil < sprop->Tp1_nitrif)	
				epv->ts_nitrif[layer] = 0.0;
		else
			epv->ts_nitrif[layer] = sprop->Tp1_nitrif/(1+pow(fabs((Tsoil-sprop->Tp4_nitrif)/sprop->Tp2_nitrif),sprop->Tp3_nitrif));
			
	}

	if (WFPS < sprop->minWFPS_nitrif)
	{
		epv->ws_nitrif[layer] = 0;
	}
	else
	{
		if (WFPS < sprop->opt1WFPS_nitrif)
		{
			epv->ws_nitrif[layer] = (WFPS - sprop->minWFPS_nitrif) / (sprop->opt1WFPS_nitrif - sprop->minWFPS_nitrif);	
		}
		else
		{
			if (WFPS < sprop->opt2WFPS_nitrif)
			{
				epv->ws_nitrif[layer] = 1;
			}
			else
			{
				epv->ws_nitrif[layer] = (1 - WFPS) / (1 - sprop->opt2WFPS_nitrif) + (WFPS-sprop->opt2WFPS_nitrif)/(1-sprop->opt2WFPS_nitrif)* sprop->scalarWFPSmin_nitrif;
			}
		}
	}
		
	if (net_miner > 0)
	{
		NH4_to_nit= sprop->netMiner_to_nitrif * net_miner + 
										sprop->maxNitrif_rate * NH4dissolv * epv->ts_nitrif[layer] * epv->ws_nitrif[layer] * epv->ps_nitrif[layer];
	}
	else
	{
		NH4_to_nit = sprop->maxNitrif_rate * NH4dissolv *  epv->ts_nitrif[layer] * epv->ws_nitrif[layer] * epv->ps_nitrif[layer];
	}
		
	N2O_flux_NIT   = NH4_to_nit * sprop->N2Ocoeff_nitrif;


	if (NH4_to_nit < 0)
	{
		printf("\n");
		printf("ERROR: in nitrification calulation (multilayer_sminn.c)\n");
		errorCode=1;
	}
	
	*NH4_to_nitrif = NH4_to_nit;
	*N2OfluxNITRIF   =  N2O_flux_NIT;

	return (errorCode);
}
int denitrification(int soiltype, double NO3avail_ppm, double pH, double WFPS, double SR_total, double* NO3_to_denitr, double* ratioN2_N2O)
{
	int errorCode = 0;
	double pDEN1, FrNO3, FrCO2, FrWFPS, FrPH;
	double FdNO3, FdWFPS, FdCO2, FdPH;
	double a, b, c, d, e;
	double denitr_flux, denitr_ratio;

	FrNO3 = FrCO2 = FrWFPS = FrPH = 0;
	FdNO3 = FdWFPS = FdCO2 = FdPH = 0;

	/* coarse: sand, loamy sand, sandy loam */
	if (soiltype <= 2)
	{
		a = 1.56;
		b = 12;
		c = 16;
		d = 2.01;
	}
	else
	{
		/* medium: loam types */
		if (soiltype <= 8)
		{
			a = 4.82;
			b = 14;
			c = 16;
			d = 1.39;
		}
		/* fine: clay types */
		else
		{
			a = 60;
			b = 18;
			c = 22;
			d = 1.06;
		}
	}

	FdNO3 = 11000 + (40000 * atan(PI * 0.002 * (NO3avail_ppm - 180))) / PI;

	e = c / pow(b, d * WFPS);

	FdWFPS = a / pow(b, e);

	FdCO2 = 24000 / (1 + (200 / exp(0.35 * SR_total)));


	if (pH < 6.5)
	{
		if (pH <= 3.5)
			FdPH = 0.001;
		else
			FdPH = (pH - 3.5) / 3;
	}
	else
		FdPH = 1;

	denitr_flux = MIN(FdNO3, FdCO2) * FdWFPS * FdPH;

	/* unit of denitr_flux: gN/ha/d -> kgN/m2/day */
	*NO3_to_denitr = denitr_flux * 0.001 * 0.0001;


	FrNO3 = (1 - (0.5 + ((atan(PI * 0.01 * (NO3avail_ppm - 190))) / PI))) * 25;

	FrCO2 = 13 + ((30.78 * atan(PI * 0.07 * (SR_total - 13))) / PI);

	pDEN1 = pow(13, 2.2 * WFPS);
	FrWFPS = 1.4 / pow(13, 17 / pDEN1);

	FrPH = 1. / (1470 * exp(-1.1 * pH));

	denitr_ratio = MIN(FrNO3, FrCO2) * FrWFPS * FrPH;

	*ratioN2_N2O = denitr_ratio;

	return (errorCode);
}
