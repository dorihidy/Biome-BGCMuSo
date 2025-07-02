/* 
soilstress_calculation.c
calculation of drought and anoxic soil water content stress layer by layers 

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
#include "bgc_constants.h"
#include "bgc_func.h"    

int soilstress_calculation(const control_struct* ctrl, const epconst_struct* epc, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf)

{
	int layer;
	double m_vwcR_layer, m_WS_avg, m_NS_avg, m_WSanoxic_avg, m_WSdrought_avg;
	double m_WSanoxic[N_SOILLAYERS], m_WSdrought[N_SOILLAYERS];
		

	int errorCode = 0;
	m_WS_avg = m_NS_avg = m_WSanoxic_avg = m_WSdrought_avg = 0;
	epv->nlayer_fullWS = 0;
	
	/* 1. Calculation SWC-stress (can occur only after emergence  - if root is presence) */
	if (epv->rootDepth)
	{		
		/*--------------------------------------------*/
		/* 1.1 SWC-stress calculation based on VWC   */
		if (ctrl->soilstress_flag == 0)
		{
			for (layer = 0; layer < N_SOILLAYERS; layer++)
			{
				/* DROUGHT STRESS */
				if (epv->VWC[layer] < epv->VWC_WScrit1[layer])
				{
					/* if VWCratio_crit is set to -9999: no soilstress calculation -> m_WS = 1 */
					if (epc->VWCratio_WScrit1 != DATA_GAP)
					{
						if (epv->VWC[layer] <= sprop->VWCwp[layer])
							m_vwcR_layer = 0;	
						else
							m_vwcR_layer = pow((epv->VWC[layer] - sprop->VWCwp[layer]) / (epv->VWC_WScrit1[layer] - sprop->VWCwp[layer]), sprop->curvatureWS[layer]);

							
					}	
					else
						m_vwcR_layer = 1;

					m_WSdrought[layer] = m_vwcR_layer;
					m_WSanoxic[layer] = 1;
										
				}
				else
				{
					/* ANOXIC CONDITION STRESS */
					if (epv->VWC[layer] >= epv->VWC_WScrit2[layer])
					{
						/* if VWCratio_crit is set to -9999: no soilstress calculation -> m_WS = 1 */
						if (epc->VWCratio_WScrit2 != DATA_GAP)
						{
							if (fabs(sprop->VWCsat[layer] - epv->VWC_WScrit2[layer]) > CRIT_PREC)
								m_vwcR_layer  = (sprop->VWCsat[layer] - epv->VWC[layer])/(sprop->VWCsat[layer] - epv->VWC_WScrit2[layer]);	
							else
								m_vwcR_layer  = 0;
						}
						else
							m_vwcR_layer  = 1;

						m_WSdrought[layer] = 1;
						m_WSanoxic[layer] = m_vwcR_layer;

					}
					/* NO STRESS */
					else
					{ 
						m_vwcR_layer       = 1;
						m_WSdrought[layer] = 1;
						m_WSanoxic[layer]  = 1;
					}
						
				}


				/* if all layers are saturated or below WP -> full stress */
				if ((epv->VWC[layer] >= epv->VWC_WScrit2[layer] || fabs(epv->VWC[layer] - epv->VWC_WScrit2[layer]) < CRIT_PREC_lenient || epv->VWC[layer] < sprop->VWCwp[layer]) && epv->rootlengthProp[layer] > 0)
				{
					epv->nlayer_fullWS += 1;
				}
				
				epv->m_WS_layer[layer] =  m_vwcR_layer;
				m_WS_avg	    += epv->m_WS_layer[layer] * epv->rootlengthProp[layer];
				m_WSdrought_avg += m_WSdrought[layer]     * epv->rootlengthProp[layer];
				m_WSanoxic_avg  += m_WSanoxic[layer]      * epv->rootlengthProp[layer];
			}
	

		}
		/* 1.2 SWC-stress calculation based on transpiration demand-possibitiy  */
		else
		{
			if (epv->n_rootlayers)
			{
				for (layer = 0; layer < epv->n_rootlayers; layer++)
				{	
					if (wf->TRPsoilw_demand[layer])
						epv->m_WS_layer[layer] = wf->TRPsoilw[layer] / wf->TRPsoilw_demand[layer];
					else
						epv->m_WS_layer[layer] = 1;

					/* control */
					if ((epv->m_WS_layer[layer] < 0 || epv->m_WS_layer[layer] > 1) && !errorCode)
					{
						printf("\n");
						printf("ERROR in soilstress_calculation.c\n");
						errorCode=1;
					}

					/* if all layers are saturated or below WP -> full stress */
					if ((fabs(epv->VWC[layer] - epv->VWC_WScrit2[layer]) < CRIT_PREC_lenient || epv->VWC[layer] < sprop->VWCwp[layer]) && epv->rootlengthProp[layer] > 0)
					{
						epv->nlayer_fullWS += 1;
					}

		
					m_WS_avg	 += epv->m_WS_layer[layer] * epv->rootlengthProp[layer];
				}
			}
			else
			{
				for (layer = 0; layer < N_SOILLAYERS; layer++) epv->m_WS_layer[layer] = 1;
				m_WS_avg = 1;
			}
			m_WSanoxic_avg  = 1;
			m_WSdrought_avg = m_WS_avg;
		}
	}
	else
	{
		for (layer = 0; layer < N_SOILLAYERS; layer++) 
		{
				epv->m_WS_layer[layer] = 1;
				m_WSdrought[layer]     = 1;
				m_WSanoxic[layer]      = 1;
				ws->soilwAVAIL[layer] = 0;
		}

		m_WS_avg = 1;
		m_WSdrought_avg = 1;
		m_WSanoxic_avg = 1;
	}
	/* avoid rounding and floating point overflow errors */
	if (fabs(m_WS_avg - 1.) < CRIT_PREC) m_WS_avg = 1;
	if (fabs(m_WSanoxic_avg - 1.) < CRIT_PREC) m_WSanoxic_avg = 1;
	if (fabs(m_WSdrought_avg - 1.) < CRIT_PREC) m_WSdrought_avg = 1;

	/* if all layers are saturated: full limitation */
	if (epv->nlayer_fullWS > 0 && epv->nlayer_fullWS >= epv->n_rootlayers - epv->germ_layer)
	{
		m_WS_avg = 0;
		for (layer = 0; layer < N_SOILLAYERS; layer++) epv->m_WS_layer[layer] = 0;
		
		if (epv->VWC[epv->n_rootlayers] <= sprop->VWCwp[epv->n_rootlayers])
		{
			m_WSanoxic_avg = 1;
			m_WSdrought_avg = 0;
		}
		else
		{
			m_WSanoxic_avg = 0;
			m_WSdrought_avg = 1;
		}
	}

		
	epv->m_WS = m_WS_avg;
	epv->m_WSanoxic = m_WSanoxic_avg;
	epv->m_WSdrought = m_WSdrought_avg;
	/****************************************************************************************/
	/* 2. calculating cumulative SWC stress and extreme temperature effect */

		
	/* 2.1 calculating WATER STRESS DAYS regarding to the average soil moisture conditions */
	
	if (epv->m_WS == 1) epv->WSlenght = 0;
	if (epv->n_actphen > 0 && epv->n_actphen >= epc->n_emerg_phenophase)
	{
		epv->WSlenght      += (1 - epv->m_WS);
		epv->cumWS         += (1 - epv->m_WS);
		epv->cumWS_anoxic  += (1 - epv->m_WSanoxic);
		epv->cumWS_drought += (1 - epv->m_WSdrought);
	}

	
	if (epv->WSlenght < epc->WSlenght_crit)
		epv->m_WSlenght = 1;
	else
		epv->m_WSlenght = epc->WSlenght_crit/epv->WSlenght;

	
	

	/****************************************************************************************/
	/* 3. N-stress based on immobilization ratio */

	if (epv->rootDepth)
		for (layer = 0; layer < N_SOILLAYERS; layer++) m_NS_avg += epv->IMMOBratio[layer] * epv->rootlengthProp[layer];
	else
		m_NS_avg = 1;

	if (1-m_NS_avg < CRIT_PREC) m_NS_avg = 1;
	
	if (epv->n_actphen >= epc->n_emerg_phenophase) epv->cumNS += (1-m_NS_avg);

	epv->m_NS = m_NS_avg;

	
	

	return (errorCode);
}