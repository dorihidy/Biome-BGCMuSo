/*
bgc_func.h
header file for function prototypes

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2025, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/


int leapControl(int year, int* enddays, int* mondays, int* leap);

int soilb_estimation(double sand, double silt, double* soilB, double* VWCsat,double* VWCfc, double* VWCwp,  
	                 double* BD, double* RCN, double* p1diffus, double* p2diffus, double* p3diffus, double* curvatureWS, int* soiltype);

int multilayer_soilcalc(control_struct* ctrl,  soilprop_struct* sprop);


int output_map_init(double** output_map, phenology_struct* phen, metvar_struct* metv, wstate_struct* ws,
	wflux_struct* wf, cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf, 
	soilprop_struct* sprop, epvar_struct* epv, soilInfo_struct* soilInfo, psn_struct* psn_sun, psn_struct* psn_shade, summary_struct* summary);

int make_zero_flux_struct(const control_struct* ctrl, wflux_struct* wf, cflux_struct* cf, nflux_struct* nf, soilInfo_struct* soilInfo, summary_struct* summary);

int annVARinit(summary_struct* summary, epvar_struct* epv, cstate_struct* cs, wstate_struct* ws);

int atm_pres(double elev, double* pa);

int restart_input(const control_struct* ctrl, const epconst_struct* epc, const siteconst_struct* sitec, 
	              wstate_struct* ws, cstate_struct* cs, nstate_struct* ns, epvar_struct* epv, soilprop_struct* sprop, restart_data_struct* restart);

int firstday(const control_struct* ctrl,  const epconst_struct* epc, const planting_struct* PLT, 
	         soilprop_struct* sprop, siteconst_struct* sitec, cinit_struct* cinit, phenology_struct* phen, epvar_struct* epv, soilInfo_struct* soilInfo,
	         wstate_struct* ws, cstate_struct* cs, nstate_struct* ns, summary_struct* summary, psn_struct* psn_sun, psn_struct* psn_shade);

int multilayer_hydrolparams(siteconst_struct* sitec, soilprop_struct* sprop, wstate_struct* ws, epvar_struct* epv);

int zero_srcsnk(cstate_struct* cs, nstate_struct* ns, wstate_struct* ws, summary_struct* summary);

int dayphen(control_struct* ctrl, const epconst_struct* epc, const phenarray_struct* phenarr, const planting_struct* PLT, phenology_struct* phen);

int management(control_struct* ctrl, fertilizing_struct* FRZ, grazing_struct* GRZ, harvesting_struct* HRV, mowing_struct* MOW, 
			   planting_struct* PLT, ploughing_struct* PLG, thinning_struct* THN, irrigating_struct* IRG, 
			   mulching_struct* MUL, CWDextract_struct* CWE, flooding_struct* FLS, groundwaterINIT_struct* GWS, int* mondays);

int groundwater_calculations(control_struct* ctrl, const siteconst_struct* sitec, const groundwaterINIT_struct* GWS, soilprop_struct* sprop, soilInfo_struct* soilInfo, epvar_struct* epv,
	                         wstate_struct* ws, wflux_struct* wf, cstate_struct* cs, nstate_struct* ns);
	int groundwater_preproc(control_struct* ctrl, const groundwaterINIT_struct* GWS, const siteconst_struct* sitec, soilprop_struct* sprop, soilInfo_struct* soilInfo, wstate_struct* ws, wflux_struct* wf);
	int groundwater_concentration(int md, const groundwaterINIT_struct* GWS, soilprop_struct* sprop, soilInfo_struct* soilInfo, wstate_struct* ws);
	int calc_VWCeq(control_struct* ctrl, const siteconst_struct* sitec, soilprop_struct* sprop);
	int groundwater_movement(const siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, soilInfo_struct* soilInfo,
		                     wstate_struct* ws, wflux_struct* wf, cstate_struct* cs, nstate_struct* ns);

int daymet(const control_struct* ctrl,const metarr_struct* metarr, const epconst_struct* epc, metvar_struct* metv, double snoww);

int phenphase(file logfile, const control_struct* ctrl, const epconst_struct* epc, const soilprop_struct* sprop, const planting_struct* PLT, 
	          phenology_struct* phen, metvar_struct* metv, epvar_struct* epv, cstate_struct* cs);
	int vernalization(const epconst_struct* epc, const metvar_struct* metv, phenology_struct* phen);
	int photoslow(const epconst_struct* epc, const metvar_struct* metv, phenology_struct* phen);

int multilayer_Tsoil(const control_struct* ctrl, const epconst_struct* epc, const siteconst_struct* sitec, const soilprop_struct* sprop, epvar_struct* epv, int yday, double snoww,
					 metvar_struct* metv);
int soilCover(siteconst_struct* sitec, soilprop_struct* sprop, metvar_struct* metv,  epvar_struct* epv, cstate_struct* cs);

int phenology(const control_struct* ctrl, const epconst_struct* epc, const cstate_struct* cs, const nstate_struct* ns,
	          phenology_struct* phen, metvar_struct* metv,epvar_struct* epv, cflux_struct* cf, nflux_struct* nf);
	int leaf_litfall(const epconst_struct* epc, double litfallc, cflux_struct* cf, nflux_struct* nf);
	int froot_litfall(const epconst_struct* epc, double litfallc, cflux_struct* cf, nflux_struct* nf);
	int yield_litfall(const epconst_struct* epc, double litfallc, cflux_struct* cf, nflux_struct* nf);
	int softstem_litfall(const epconst_struct* epc, double litfallc, cflux_struct* cf, nflux_struct* nf);
	int transfer_fromGDD(const epconst_struct* epc, const cstate_struct* cs, const nstate_struct* ns, 
	                     phenology_struct* phen, metvar_struct* metv, epvar_struct* epv, cflux_struct* cf, nflux_struct* nf);

int multilayer_rootDepth(const epconst_struct* epc, const soilprop_struct* sprop, const cstate_struct* cs, siteconst_struct* sitec,  epvar_struct* epv);
	int calc_nrootlayers(int flag, double maxRD, double frootc, const siteconst_struct* sitec, epvar_struct* epv);

int radtrans(const control_struct* ctrl, const phenology_struct* phen, const cstate_struct* cs, const epconst_struct* epc, const siteconst_struct* sitec,
	         metvar_struct* metv, epvar_struct* epv);

int prcpANDrunoffH(const control_struct* ctrl, const wstate_struct* ws, const metvar_struct* metv, const soilprop_struct* sprop,  const epconst_struct* epc, 
	               epvar_struct* epv, wflux_struct* wf);

int snowmelt(const metvar_struct* metv, wflux_struct* wf, wstate_struct* ws);

int Elimit_and_PET(const control_struct* ctrl, const epconst_struct* epc, const soilprop_struct* sprop, const metvar_struct* metv, epvar_struct* epv, wflux_struct* wf);

int soilEVP_calc(control_struct* ctrl, const siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf);

int conduct_calc(const control_struct* ctrl, const metvar_struct* metv, const epconst_struct* epc, epvar_struct* epv, int simyr);

int maint_resp(const control_struct* ctrl, const planting_struct* PLT, const cstate_struct* cs, const nstate_struct* ns, const epconst_struct* epc, const metvar_struct* metv,
	           epvar_struct* epv, cflux_struct* cf);

int canopy_et(const control_struct* ctrl, const metvar_struct* metv, epvar_struct* epv, wflux_struct* wf);

int penmon(const pmet_struct* in, int out_flag,	double* et);

int photosynthesis(const control_struct* ctrl, const epconst_struct* epc, const metvar_struct* metv, const cstate_struct* cs, const wstate_struct* ws, const phenology_struct* phen,
	               epvar_struct* epv, psn_struct* psn_sun, psn_struct* psn_shade, cflux_struct* cf); 
	int farquhar(const control_struct* ctrl, const metvar_struct* metv, psn_struct* psn);

int decomp(const metvar_struct* metv, const epconst_struct* epc, soilprop_struct* sprop, const siteconst_struct* sitec, const cstate_struct* cs, const nstate_struct* ns, 
	       epvar_struct* epv, cflux_struct* cf, nflux_struct* nf, ntemp_struct* nt);
	int CH4flux_estimation(const soilprop_struct* sprop, int layer, double VWC, double T, double* CH4flux);

int daily_allocation(const control_struct* ctrl, const epconst_struct* epc, const soilprop_struct* sprop, const metvar_struct* metv, const NdepControl_struct* ndep,
	                 cstate_struct*cs,  nstate_struct* ns, cflux_struct* cf, nflux_struct* nf, epvar_struct* epv, ntemp_struct* nt, double naddfrac);


int flowering_heatstress(const epconst_struct* epc, const metvar_struct* metv, cstate_struct* cs, epvar_struct* epv, cflux_struct* cf, nflux_struct* nf);

int annual_rates(const epconst_struct* epc, epvar_struct* epv);

int growth_resp(const epconst_struct* epc, cflux_struct* cf);

int multilayer_transpiration(control_struct* ctrl, const siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, 
	                         wstate_struct* ws, wflux_struct* wf);

int irrigating(const control_struct* ctrl, const irrigating_struct* IRG, const siteconst_struct* sitec, const soilprop_struct* sprop,
	           epvar_struct* epv, wstate_struct* ws, wflux_struct* wf);

int multilayer_hydrolprocess(control_struct* ctrl, siteconst_struct* sitec, soilprop_struct* sprop, soilInfo_struct* soilInfo, const epconst_struct* epc, epvar_struct* epv,
	                         wstate_struct* ws, wflux_struct* wf, nstate_struct* ns, nflux_struct* nf, cstate_struct* cs, cflux_struct* cf, flooding_struct* FLS, int* mondays);
	int infiltANDpond(wstate_struct* ws, wflux_struct* wf);
	int pondANDrunoffD(control_struct* ctrl, siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf);
	int tipping(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf);
		int calc_drainage(int flagRAIN, double soilB, double INFILT, double VWC, double VWCsat, double VWCfc, double dz0, double DC, double conductSAT_cmday, double* DRN, double* EXCESS, double* VWCnew);
	int capillary_tipping(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf);
	int groundwater_tipping(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf);
	int diffusion(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf);
	int capillary_diffusion(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf);
	int groundwater_diffusion(siteconst_struct* sitec, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf);
	int calc_diffus(int layer, const soilprop_struct* sprop, double dz0, double VWC0, double VWC0_sat, double VWC0_EqFC, double VWC0_wp, double VWC0_hw,
		            double dz1, double VWC1, double VWC1_sat, double VWC1_EqFC, double VWC1_wp, double VWC1_hw, double fluxLimit, double dLk, double* DBAR, double* soilwDiffus);
	int soilstress_calculation(const control_struct* ctrl, const epconst_struct* epc, soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf);
	int flooding(control_struct* ctrl, const siteconst_struct* sitec, const flooding_struct* FLS, soilprop_struct* sprop, epvar_struct* epv,
		         wstate_struct* ws, wflux_struct* wf, cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf, soilInfo_struct* soilInfo, int* mondays);
		int flooding_concentration(int md, const flooding_struct* FLS, soilInfo_struct* soilInfo, wstate_struct* ws);
	int potEVPsurface_to_actEVPsurface(soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf);
		int EVPphase1TOphase2(const soilprop_struct* sprop, epvar_struct* epv, wstate_struct* ws, wflux_struct* wf);
	int hydrol_control(siteconst_struct* sitec, soilprop_struct* sprop, wstate_struct* ws, wflux_struct* wf, epvar_struct* epv);

int water_state_update(wflux_struct* wf, wstate_struct* ws);

int CN_state_update(const siteconst_struct* sitec, const epconst_struct* epc, soilInfo_struct* soilInfo, soilprop_struct* sprop, control_struct* ctrl, epvar_struct* epv,
	                      cflux_struct* cf, nflux_struct* nf, cstate_struct* cs, nstate_struct* ns, int alloc, int evergreen);
	int MRdeficit_calculation(const epconst_struct* epc, control_struct* ctrl, cflux_struct* cf, nflux_struct* nf, cstate_struct* cs, nstate_struct* ns);
	int CNratio_control(cstate_struct* cs, double CNratio, double cpool, double npool, double cflux, double nflux, double CNratio_flux);

int senescence( const soilprop_struct* sprop, const epconst_struct* epc, const grazing_struct* GRZ, const metvar_struct* metv,
			   control_struct* ctrl, cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf, epvar_struct* epv);
	int genprog_senescence(const epconst_struct* epc, const metvar_struct* metv, epvar_struct* epv, cflux_struct* cf, nflux_struct* nf);

int mortality(const control_struct* ctrl, const soilprop_struct* sprop, const epconst_struct* epc,
	          epvar_struct* epv, cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf, int simyr);	

int multilayer_sminn(control_struct* ctrl, const metvar_struct* metv, const siteconst_struct* sitec, const NdepControl_struct* ndep,
	                 cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf, soilprop_struct* sprop, epvar_struct* epv, soilInfo_struct* soilInfo);
	int nitrification(int layer, const soilprop_struct* sprop, double net_miner, double Tsoil, double pH, double WFPS, double NH4dissolv, epvar_struct* epv, double* N2OfluxNITRIF, double* NH4_to_nitrif);
	int denitrification(int soiltype, double NO3avail_ppm, double pH, double WFPS, double SR_total,  double* NO3_to_denitr, double* ratioN2_N2O);

int multilayer_leaching(const siteconst_struct* sitec, soilprop_struct* sprop, soilInfo_struct* soilInfo, cstate_struct* cs, nstate_struct* ns,  wstate_struct* ws, wflux_struct* wf);
	int capillary_leaching(int dm, const siteconst_struct* sitec, soilprop_struct* sprop, soilInfo_struct* soilInfo, wstate_struct* ws, wflux_struct* wf,
		                   double* dismatLeachNORM, double* dismatLeachCAPIL, double* dischargeNORM, double* dischargeCAPIL, double* rechargeNORM, double* rechargeCAPIL);
	int groundwater_leaching(int dm, const siteconst_struct* sitec, soilprop_struct* sprop, soilInfo_struct* soilInfo, wstate_struct* ws, wflux_struct* wf,
		                     double* dismatLeachNORM, double* dismatLeachCAPIL, double* dischargeNORM, double* dischargeCAPIL, double* rechargeNORM, double* rechargeCAPIL);
	int check_soilcontent(int layerFlag, int pool2content_flag, const soilprop_struct* sprop, cstate_struct* cs, nstate_struct* ns, soilInfo_struct* soilInfo);
	int calc_DISSOLVandBOUND(int firsttime_flag, int partlyORtotal_flag, const soilprop_struct* sprop, soilInfo_struct* soilInfo);
	int check_virtualLayer_balance(control_struct* ctrl, soilInfo_struct* soilInfo, soilprop_struct* sprop, wflux_struct* wf);

int planting(control_struct* ctrl, const siteconst_struct* sitec, const soilprop_struct* sprop, const planting_struct* PLT, epconst_struct* epc,
	         epvar_struct* epv, phenology_struct* phen, cstate_struct* cs, nstate_struct*ns, cflux_struct* cf, nflux_struct* nf);
	int planttype_determination(control_struct* ctrl, phenology_struct* phen);
	int conductLimit_calculations(const soilprop_struct* sprop, const epconst_struct* epc, epvar_struct* epv);

int thinning(const control_struct* ctrl, const epconst_struct* epc, const thinning_struct* THN, 
	         cstate_struct* cs, nstate_struct* ns, wstate_struct* ws, cflux_struct* cf, nflux_struct* nf, wflux_struct* wf, epvar_struct* epv);

int mowing(const control_struct* ctrl, const epconst_struct* epc, const mowing_struct* MOW, epvar_struct* epv, 
		   cstate_struct* cs, nstate_struct* ns, wstate_struct* ws, cflux_struct* cf, nflux_struct* nf, wflux_struct* wf);

int grazing(control_struct* ctrl, const epconst_struct* epc, const soilprop_struct* sprop, grazing_struct* GRZ, epvar_struct* epv,
	        cstate_struct* cs, nstate_struct* ns, wstate_struct* ws, cflux_struct* cf, nflux_struct* nf, wflux_struct* wf, int* mondays);

int harvesting(file econout, control_struct* ctrl, phenology_struct* phen, const epconst_struct* epc, const harvesting_struct* HRV, const irrigating_struct* IRG,
	           epvar_struct* epv, cstate_struct* cs, nstate_struct* ns, wstate_struct* ws, cflux_struct* cf, nflux_struct* nf, wflux_struct* wf);

int ploughing(const control_struct* ctrl, const epconst_struct* epc, siteconst_struct* sitec, soilprop_struct* sprop, soilInfo_struct* soilInfo, metvar_struct* metv, epvar_struct* epv,
	          ploughing_struct* PLG, cstate_struct* cs, nstate_struct* ns, wstate_struct* ws, cflux_struct* cf, nflux_struct* nf, wflux_struct* wf);

int fertilizing(const control_struct* ctrl, const siteconst_struct* sitec, soilprop_struct* sprop, soilInfo_struct* soilInfo, epvar_struct* epv, fertilizing_struct* FRZ,
				cstate_struct* cs, nstate_struct* ns, wstate_struct* ws, cflux_struct* cf, nflux_struct* nf, wflux_struct* wf);

int mulching(control_struct* ctrl, const mulching_struct* MUL,  cstate_struct* cs, nstate_struct*ns, cflux_struct* cf, nflux_struct* nf);

int CWDextract(control_struct* ctrl, const CWDextract_struct* CWE, cstate_struct* cs, nstate_struct*ns, cflux_struct* cf, nflux_struct* nf);

int cutdown2litter(const soilprop_struct* sprop, const epconst_struct* epc, const epvar_struct* epv, cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf);

int precision_control(control_struct* ctrl, soilprop_struct* sprop, wstate_struct* ws, cstate_struct* cs, nstate_struct* ns, soilInfo_struct* soilInfo);
	int precision_reset(int dm, int layer, control_struct* ctrl, soilprop_struct* sprop, cstate_struct* cs, nstate_struct* ns, soilInfo_struct* soilInfo, double* cpool, double* npool, int* errorCode);

int check_water_balance (wstate_struct* ws, int first_balance);
int check_carbon_balance(cstate_struct* cs, int first_balance);
int check_nitrogen_balance(nstate_struct* ns, int first_balance);

int aboveANDbelow(soilprop_struct* sprop, epvar_struct* epv, cstate_struct* cs, cflux_struct* cf);

int cnw_summary(const epconst_struct* epc, const siteconst_struct* sitec, const soilprop_struct* sprop, const metvar_struct* metv,
	            const cstate_struct* cs, const cflux_struct* cf, const nstate_struct* ns, const nflux_struct* nf, const wflux_struct* wf, const soilInfo_struct* soilInfo,
				epvar_struct* epv, summary_struct* summary);

int restart_output(const wstate_struct* ws, const cstate_struct* cs, const nstate_struct* ns, const epvar_struct* epv, restart_data_struct* restart);