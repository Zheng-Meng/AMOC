
% Compute_fluxes_osnap_master_2023.m
% =====================================================
%
% CALL change_dt.m
% CALL fun_osnap_sst_adt_winds_ave.m
% CALL fun_osnap_mooring_ave_instru.m
% CALL fun_osnap_mooring_vinterp_random.m or fun_osnap_mooring_vinterp_random2.m
% CALL fun_HT_FWT_decomposition
% CALL vinterp_single_prof.m
%
% =====================================================
%
% THIS ROUTINE IS THE MASTER SCRIPT TO COMPUTE THE OVERTURNING FLUXES
% ACROSS THE OSNAP SECTIONS.
%
%
% -------------------------------------------------------------------------
%
% INPUT
%       Monthly T/S/U/V climatology from FLAME;
%       Monthly T/S climatology from WOA13;
%       monthly gridded analysis of T/S/rho/ptmp/pden;
%       The statistical distribution of t/s/u/v for each monthly period;
%       The statisttical distribution of WIND, and SSH for each month;
%
% OUTPUT
%   data fields:
%       temperature (x,z,t)
%       potential temperature (x,z,t)
%       practical salinity (x,z,t)
%       density (x,z,t)
%       potential density (x,z,t)
%       velocity normal to the section (x,z,t)
%       
%   fluxes:
%       volume flux (x,z,t) and (x,sigma,t)
%       heat flux (x,z,t) and (x,sigma,t)
%       salt flux (x,z,t) and (x,sigma,t)
%       MOC(t)
%       MHT(t)
%       MFT(t)
%
% -------------------------------------------------------------------------
% 2017-03-19  Modified from Compute_moc_osnap_wMC_1122.m;
% 2017-03-30  Set up mask_land_v and _i using ETOPO2 on the V/I-Grid;
%             The velo. from RTWB1 are copied into the wedge west of the
%               mooring;
% 2017-03-31  Topograph-following vinterp is double checked;
%             Topography-following velocity interp for the bottom 300m;
%             The velo. from RTEB1 are copied into the wedge east of the
%               mooring;
% 2017-04-06  For the bottom triangles in the CM array, copy values from
%               the horizontally interpolated  one (above the common maximum
%               depth of the two bounding moorings) into the triangles;
%             Apply topography-following interp. only at depths > 1000m; 
% 2017-04-08  Add an option of new model climatology (u/v/t/s) above the
%               Labrador shelf (from Brad) (see Section A2);
% 2017-04-09  The Ekman layer depth= 7.6*U_wind/sqrt(sind(lat_i));
% 2017-04-10  Add an option to apply the surface reference velocity at each
%               grid cell #flag_full_grid_baro_corr#;
% 2017-06-27  Update mooring indices for 53 moorings;
%             Update the part to randomly draw observed T/S/U profiles for
%               Monte Carlo iterations;
% 2017-07-07  Add an option to test of using no OUC glider data;
% 2017-07-12  Update the setup of mask_obs_uv to allow a direct velocity 
%               interpolation between NOCM3-M4 when there is no OOI
%               data;
%             Update the velocity interpolation to fill the WB&EB wedges
%               in the Rockall Trough;
%             Update the reference salinity;
% 2017-07-16  Copy velocities from OM3 to OM2 since OM2 is west of OM3 and
%               has only T/S records;
% 2017-07-20  Add an option to call vinterp_osnap_mooring_random.m instead
%               of using regridded profiles directly;
%             Add an option to use pre-regridded daily profiles 
%               #flag_test_use_regridded_profiles#;
% 2017-07-26  Add an option for the flux estimates across the full array
%               with a net flow of -1.6 Sv at West and 1.6 Sv at East;
% 2017-07-31  Add an option to run MC with STE or STD;
% 2017-08-21  Add a test using only geostrophic velocties from CF6 to M1,
%               and also apply the topography-following interpolation;
% 2017-08-23  Extend the test first introduced on 8/21 to all other current
%               meter mooring arrays;
% 2017-09-18  Add two tests: (a) apply time-varying surface altimetry
%               reference only in areas without deep moorings and keep 
%               the deep-mooring-referenced geostrophic profiles otherwise; 
%               (b) use a temporal mean barotropic component to provide the
%               surface reference in areas without deep moorings; This
%               mean component is obtained by averging the time-varying
%               barotropic velocities obtained from test (a);
% 2017-10-13  Test of using alternative ADT products provided by WHOI;
% 2017-10-16  Test the use the AGV product from AVISO to provide 
%               the surface reference;
% 2017-10-25  Add an option on how to distribute Vcomp
%               #flag_vcomp_everywhere#;
% 2017-11-15  Add an option for calculating fluxes 
%               at different time interval #flag_data_interval#;
%             Add an option to make sure that full-grid velo. shears
%               only calcualted where they have to, i.e., above the western 
%               flank of the Hatton Bank #flag_full_grid_shears== 2#;
%             Modify #flag_test_noOUCglider_for_OA# to load an OA product
%               that incoporates no OUC glider data; 
%             Add a test to use an OA product that incoporates no glider 
%               data at all #flag_test_OA_noglider_for_OA#;
% 2017-11-17  Use pp_osnap_sst_adt_winds_ave.m to provide temporally 
%               averaged surface observations;
% 2017-11-19  Add an option to test using no OOI data in the calculations 
%               #flag_test_no_OOI_data#;
% 2017-11-20  The calculation period is determined by the latest deployment
%               of the LS moorings (8/23/2014) and the earliest end date 
%               of the RTEB1 CM (3/28/2016);
% 2017-11-29  Uncertainty estimation incorporates errors in the mean
%               barotropic component;
% 2018-02-23  Apply topography-following property interpolation;
% 2018-03-06  Modify #flag_keep_vrefshortmoorings_umm12# to a new test
%               option that can keep bottom-referenced geostrophic profiles
%               in any selected areas 
%               #flag_keep_vrefshortmoorings_selected#;
%             Modify #flag_vcomp_selected# to provide two options on how to
%               distribute Vcomp: (a) in areas without deep-mooring
%               reference; (b) in basins kown to be barotropic;
%             Revert the change on 2/23 so not to apply topo-following
%               property interpolation;
% 2018-03-16  By default, not to use gridded daily profiles for calculating
%               fluxes at time interval > 1 day;
%               #flag_regridded_profiles = 0#;
% 2018-03-18  Add to test of filling the time-mean velo. in all the areas
%               without direct measurements;
% 2018-03-27  Rewrite the part for drawing a profile from 30 of
%               daily regridded profiles - YET TO TEST FOR WEST/ALL;
%             Provide 4 MOC definitions;
% 2018-04-04  Add to test of filling the time-mean velo. only in the SAMS
%               glider survey domain;
% 2018-04-22  For MC simulation with loaded daily profiles, test of drawing
%               random rho profile for each time interval as well -- 
%               TESTING ONLY, NOTE USED;
% 2018-04-24  Change to flag_regridded_profiles:
%               tested of using averaged T/S/U/V profiles and then derive
%               dens/pden/ptmp profiles -- TESTING ONLY, NOT USED;
% 2018-05-15  Update flag_vcomp_constraint==1: +0.6 Sv at OSNAP-East 
%               to match the long-term mean transport at the BS (1 Sv) 
%               reported by Woodgate 2018 PiO;
% 2018-05-18  Update SE in 20180315_test1_v_baro_mean.mat;
% 2018-05-19  Call updated fun_osnap_sst_adt_winds_ave.m and
%               fun_osnap_mooring_ave_instru.m based on change_dt;
% 2018-05-21  Output the Ekman transport from each MC run;
% 2018-05-22  Output the total compensation transports from each MC
%               run;
%             Update long-term observations with uncertainty for MC:
%               vnet_BS= -1 +/- 0.05 Sv (Woodgate 2018 PiO),
%               vnet_DS= 1.6 +/- 0.2 Sv (Curry et al. JPO);
%             Update flag_moc_def to provide a combination of any
%               definitions;
%             Add flag_test_gamma_n to test of using neutral density;
% 2018-05-24  Add flag_test_v_baro_mean and flag_MC_MOC;
%             Update epsilons using 1/1000 for MOC and its 1/50 for MHT and
%               MFT;
%             The script can be ran simultaneously in diffferent working
%               direcotries.
% 2018-05-25  Update flag_MC_MOC to specify a fixed number iterations;
% 2018-05-26  Update epsilons for OSNAP-West, using 1/10000 for MOC and 
%               its 1/100 for MHT and MFT;
%             Maximum #iterations changed to 5000;
% 2018-05-27  Update flag_test_direct_velo_only with a new time-mean
%               velocity;
% 2018-05-30  Small epsilons for OSNAP-West are only used for the 
%               testing purpose;
%             For OSNAP-All, add flag_HT_FWT_decomposition to perform 
%               the MHT and MFT decompositions by calling 
%               fun_HT_FWT_decomposition();
%             New S0 based on the 0524T results;
% 2018-05-31  Output Text for West and East respectively; Not to output
%               vflux_ext and vflux_ek anymore;
% 2018-06-11  Update flag_long_term_mean_vbaro and rename to 
%               flag_vbaro_mean; Update flag_test_v_baro_mean and rename it
%               to flag_test_vbaro_mean;
% 2018-06-26  For the inshore Labrador Current, only use model velocity
%               climatology (not property);
% 2018-08-17  Spatial smooth AVISO products (Gourcuff et al. 2010) when
%               using full-grid altimetry reference (tentative); 
% 2018-09-25  Add a test of withholding the RTADCP data;
% 2018-09-26  Update flag_test_gamma_n for sigma2 or gamma_n;
% 2018-10-05  Updated fun_osnap_mooring_vinterp_random_v2 to deal with data
%               loss from RTADCP;
% 2018-10-05  Add an option to deal with RTADCP missing data: either using
%               EB1 data (CALL fun_osnap_mooring_vinterp_random.m) or using
%               time-mean RTADCP data (CALL fun_osnap_mooring_vinterp_random_v2.m);
% 2018-10-24  Implement the same RTADCP data processing when using gridded
%               daily profiles (flag_regridded_profiles= 1);
% 2018-12-03  Add flag_test_noRR to test of removing the RR arrays (both
%               NIOZ and UM);
% 2018-12-04  Add flagt_test_fewermr to test of reducing the density of the
%               current meter array;
% 2018-12-11  Modify flag_test_direct_velo_only to add option 3 for filling
%               the areas having direct velocity measurements with the 
%               time-mean velocities; 
%             Remove redundant codes regarding flag_vcomp_everywhere for
%               simplicity;
% 2019-03-20  Fill *_v land points with NaNs for all variables except v_v;
%             Always use mask_land_v together with v_v; 
% 2019-03-29  Update with all 2014-2018 data;
% 2019-03-31  Fill the 2016-2018 data gap at OM5 (C2M1954) and OM6 (C3M1955);
%             Fill the 2018 data gap at OM42 (UMM1 ADCP) with the time-mean
%             data from that instrument - need to double check with Bill;
% 2019-04-01  Updated the reference salinity (s0);
%             When no data, don't use FLAME uv climatology. VOID this step
%               instead;
% 2019-07-22  Updated all data path/files;
%             Updated most of mooring indices except for a few test runs-
%             note: mooring indices can be all substituted with longitudes
%             for simplicity, which however are now both used;
% 2019-07-23  Updated flags - more details to add!!!!!
% 2019-07-25  Updated (one more) input for fun_HT_FWT_decomposition;
% 2019-08-01  Test option: Using no C3 data at all;
% 2019-08-05  Enabled the use of flag_test_fewermr;
% 2019-08-06  Updated regarding OM3 and OM5 velocity loss;
% 2019-08-19  Add flag_C3_recon to reconstruct C3 profile from K7 when none
%               of the C3 measurements were returned;
% 2019-09-24  By default, NOT to reference the shears to D4 and D5;
% 2019-11-10  Add the calculation of v_baro_mean for flag_data_interval=1;
% 2019-11-11  Add flag_test_glider_shear - remove abnormal shears right
%               next to UMM4 and set the shears to zero below 2000m;
% 2019-11-12  Add flags to remove RREX and D5 moorings;
% 2019-11-17  Add GloSea5 V output on the Labrador Shelf;
% 2019-11-21  Corrected FLAME_clim_mm data;
% 2020-01-01  Calculate MFT using boundary-mean salinity (updated s0)
%               and add to output MST (not used as part of MC);
% 2020-01-15  Add tests about SAMS glider output (flag_test_OA_noSAMSgdr,
%               flag_test_SAMSgdr_geov);
% 2020-01-22  Adjust flag_test_D5_nouse to also remove added D4 instrument;
%             Add flag_test_OA_depth to loaad depOA product;
% 2020-02-06  Test of new model UV climatology the Labrador Shelf
%               (flag_test_LC_UV_clim);
% 2020-02-12  Test only using geostrophy for M2-M3 (flag_test_M2M3_geo); 
% 2020-02-14  Test for the time priod when the data are available;
% 2020-02-27  Disable the change of i_t_begin and i_t_end related to
%               flag_test_M2M3_geo, so to run over the full 47 months;
% 2020-03-05  BS-OSNAP boundary salinity is now the same for All, West and
%               East when calculating MFT;
%             Updated OA datasets;
% 2020-03-06  Updated s0 using the latest OA products;
% 2020-03-11  Add flag_test_pd_MOCwest to specifiy a density range when determining 
%               the MOC at OSNAP-W: [27.0 27.8];
% 2020-04-05  Output uncertainty as a function of transport in each density
%               bin;
%             Calculating fluxes for West and East, All at the same time;
% 2020-08-22  Updated fun_HT_FWT_decomposition.m;
% 2020-08-24  Updated to additionally ouput the net throughflow component
%               in the HF and FT in density space;
% 2021-06-01  changed from 30-day mean to monthly mean
% 2021-07-08  Updated to include IB5
% 2022-01-25  updated LSA and LSB, fixed bugs, added switch


clear;
% clc;
warning off;

% startup


disp('+' );
disp('+' );
disp(['+ Now in Compute_fluxes_osnap_2023:',datestr(now)]);
disp('+' );



%-----------------+
%% SET FLAGS      |
%-----------------+

% #flag_section#
% flag defines which section for the calculations
%   'east' -- OSNAP-East section
%   'west' -- OSNAP-West section
%   'all'  -- OSNAP-All
flag_section = 'all';


% #flag_data_interval#
% Time interval for the flux estimates;
% Call fun_osnap_sst_adt_winds_ave.m 
%   to provide temporally averaged surface
%   products at the section over each 
%   specified time interval;
flag_data_interval = 'm'; 
if flag_data_interval=='m'
    flag_data_interval_4calculation=30;
elseif flag_data_interval==1
    flag_data_interval_4calculation=1;
end


% #flag_topo_interp#
% 1-- [Default] apply topography-following
%       velocity interpolation for the bottom 200 meters;
% 0-- do nothing;
flag_topo_interp = 1;
                    

% #flag_OA#
% 1-- [Default] Use OA T/S other than climatology field for
%       computing the MHT and MFT;
% 0-- Use the WOA13 monthly climatology fields;
flag_OA = 1;


% #flag_regridded_profiles#
% 1-- Create mean or random regridded T/S/rho/U/V profile 
%       from each mooring based on daily regridded profiles;
% 0-- [Default] Call #fun_osnap_mooring_vinterp_random*.m# to provide 
%       regridded profiles for each time step at each mooring;
flag_regridded_profiles = 0;


% #flag_MC#
% 0-- Use the mean values;
% 1-- Use a Monte Carlo simulation to estimate the uncertainty in the mean;
flag_MC = 0;


% #flag_MC_MOC#
% 0-- [Default] Use epsilon_MOC epsilon_MHT epsilon_MFT at the same time;
% 1-- Use epsilon_MOC and epsilon_MHT to stop the interation;
% 2-- Use only epsilon_MOC to stop MC iteration;
% [#iterations] -- Run a fixed number of iterations, must >= 10;
flag_MC_MOC = 0;


% #flag_MC_SE#
% 0-- Draw random values from mean+/-SD;
% 1-- [Default] Draw random values from mean+/-SE;
%
% NOTE
%   This flag only works when #flag_MC= 1# & #flag_regridded_profiles= 0#
flag_MC_SE = 1;


% #flag_moc_def# Choose from the following options:
%   0-- [Default] Maximum of the overturning streamfucntion [sigma0];
%   1-- Sum of all northward transport [sigma];
%   2-- [Default] Maximum of the overturning streamfucntion [depth];
%   3-- Sum of all northward transport [depth];
flag_moc_def = [0 2];


% #flag_full_grid_shears#
% 0-- Dens profiles between UMM3-UMM4 and between UMM4-SAMS_glider_endpoint 
%       are used for the thermal wind calculations;
% 1-- All pre-defined dens profiles are used for the thermal wind 
%       calculations;
% 2-- [Default] Only dens profiles between UMM4 and 14.7W (SAMS glider 
%       survey eastern endpoint) are used for 
%       the thermal wind calculatios;
flag_full_grid_shears = 2;


% #flag_vcomp_constraint#
% 0-- Zero-net-flow constraint at the whole section;
% 1-- -1.6 Sv net transport at West and 0.6 Sv (BS-DS) at East;
% 2-- [Default] -1.6 Sv net transport at West and +1.6 at East;
flag_vcomp_constraint = 2;


% #flag_vcomp_everywhere#
% 0-- [Default] Apply the corrections to the geostrophic profiles only, 
%       specified by mask_flux_i;
% 1-- Apply the correction velocity to the whole section;
%
% NOTE
%   Flag 0 accords to mask_flux_i for distributing Vcomp;
flag_vcomp_everywhere = 0; 


% #flag_vcomp_selected#
% 0-- OFF;
% 1-- [Default] Apply Vcomp to where the surface reference is being used;
% 2-- [Test] Apply Vcomp to the Labrador and Irminger Sea interior;
%
% NOTE
%   This flag works only when #flag_vcomp_everywhere== 0#;
%   Flag option 1 and 2 make chagnes to mask_flux_i;
flag_vcomp_selected = 1;


% #flag_vref_shortmoorings#
% 0-- Don't reference to deep moorings by assuming a
%       level-of-no-motion at the bottom of each 
%       geostrophic segment;
% 1-- [Default] Reference to the velocities from deep
%       moorings wherever applicable;
flag_vref_shortmoorings= 1;



% #flag_keep_vrefshortmoorings_all#
% 1-- Keep *all* the deep-mooring referenced velocities;
%     So not to add the barotropic velocity if the 
%       shears have been referenced to the bottom velocities 
%       from deep moorings; 
% 0-- [Default] Apply the barotropic correction for all the 
%       geostrophic profiles; 
%
% NOTE
%   This flag makes changes to mask_v_baro;
%   This flag works when #flag_vref_shortmoorings= 1#;
flag_keep_vrefshortmoorings_all = 0;



% #flag_keep_vrefshortmoorings_selected#
% 1-- [Default] 
%       [updated on 9/24/2019] Keep the deep-mooring reference 
%           only in the designated areas, i.e., in all the geostrophic 
%           segments except DSOW5, UM-D4 and UM-D5; 
%       [updated on 7/23/2019] Keep the deep-mooring reference 
%           only in the designated areas, i.e., in all the geostrophic 
%           segments except for DSOW5; 
%       [old] Keep the deep-mooring reference only in the
%           designated areas, i.e., in all the geostrophic segments 
%           except for DSOW5 and UM-D4; 
% 0-- Apply the barotropic correction as indicated by
%       #flag_keep_vrefshortmoorings_all#;
%
% NOTE
%   This flag makes changes to mask_v_baro;
%   This flag works when #flag_vref_shortmoorings#= 1;
%   This flag is valid only when #flag_keep_vrefshortmoorings_all#== 0;
flag_keep_vrefshortmoorings_selected = 1;



% #flag_baro_corr#
% 1-- [Default] Add a barotropic component to the referenced geostrophic 
%       profiles; 
%     The barotropic velocity is obtained from the difference 
%       between the surface velocity from altimetry and the one 
%       referenced to the bottom;
% 0-- Do nothing to upset the geostrophic profiles that have been
%       referenced to the bottom;
%
% NOTE
%   This flag calculates v_ssh and v_baro and sets up mask_v_baro;
flag_baro_corr = 1; 



% #flag_ssh_mean#
% 1-- [default] Use the time-mean ssh to provide a surface reference for
%       the thermal winds;
% 0-- Use time-varying ssh;
flag_ssh_mean = 0; % using the time-mean vssh


% #flag_vbaro_mean# 
% 1-- [Default] Use long-term mean v_baro 
%       to provide a mean barotropic component;
% 0-- Apply time-varying surface reference velocities;
%
% NOTE
%   This flag works when #flag_baro_corr# = 1;
%   This flag is to load and use a pre-defined v_baro_mean;
flag_vbaro_mean = 1;



% #flag_full_grid_baro_corr#  
% 1-- Use the altimetry surface velocity across each grid cell;
% 0-- [Default] Use the altimetry surface velocity derived from the SSH at
%       density profiles that are used in #flag_full_grid_shears#;
%
% NOTE
%   This flag only works when #flag_baro_corr# = 1;
%   This flag is to update v_ssh, which will be used later 
%       for calculating v_baro;
%	This flag is void when #flag_vbaro_mean# = 1;
flag_full_grid_baro_corr = 0;



% #flag_aviso_smoothing#
% 1-- Spatially average altimetry velocity over length scales ~100km;
% 0-- OFF;
%
% NOTE
%   This flag works to update v_ssh;
%   This flag works on every grid point when flag_full_grid_baro_corr= 1; 
% [!!!!!TODO] This flag works on where the full-grid shears are calculated 
%       when flag_baro_corr= 1 (e.g., in the glider domain);
flag_aviso_smoothing = 0;
    % !!! need to be updated !!!

    

% #flag_vary_dek#
% 1-- [Default] Time&location-varying Ekman layer depth;
% 0-- A constant Ekman layer depth everywhere (50m);
flag_vary_dek = 1;



% #flag_OOI_for_tw#
% 1-- [Default] Use the OOI data for the thermal wind calculation whenever
%       available;
% 0-- Don't use the OOI data for the thermal wind calculation;
flag_OOI_for_tw = 1;



% #flag_opa_clim#
% 0-- [Default] No changes to make;
% 1-- Use OPA monthly climatology for LC;
flag_opa_clim = 0; 



% #flag_LC_UV#
% 0-- [Default] Use FLAME monthly climatology for LC;
% 1-- Use NEMO 5-day output for LC;
% 2-- Use GloSea5 monthly output for LC;
flag_LC_UV = 0;



% #flag_LC_TS#
% 0-- [Default] Use WOA18 monthly climatology for LC;
% 1-- Use GLOSEA5 monthly output for LC;
flag_LC_TS = 0;



% #flag_mean_RTADCP#
% 1-- [Default] using time-mean RTADCP when no ADCP data returned
% 0-- Do nothing to leave blank at RTADCP when no data
%       returned, and EB1 data will be used to fill the wedge;
flag_mean_RTADCP = 1;


% #flat_C3_recon#
% 1 -- [Default] Use K7 to reconstruct C3 only if no C3 data are returned;
% 0 -- Use the time-mean C3 when no data are returned;
% *** NOTE
%       This flag will be eliminated later, after the data file has
%       incoporated this reconstruction.  As of now, the data file has
%       the mean C3 data when no measurements are recorded.
flag_C3_recon = 1;



% #flag_HT_FWT_decomposition#
% 1-- Decomposing MHT and MFT in both density and depth spaces;
% 0-- OFF;
flag_HT_FWT_decomposition = 1;




% -------------------------
% Testing modules
% -------------------------


% #flag_test_LC_UV_clim#
% 0-- Use FLAME climatology;
% 1-- Use NEMOnew climatology;
% 2-- Use GLORYS climatology;
% 3-- Use VITALS climatology;
% 4-- [Default] Use multiple-model-mean climatology;
flag_test_LC_UV_clim = 4;



% #flag_test_mr_nouse#
% 1-- Exclude data from D4/D5/RREX;
% 0-- OFF;
flag_test_mr_nouse = 0;
    flag_test_d5_nouse = 0;
    flag_test_rrex_nouse = 0;

    flag_test_d5_mean = 0; % using mean data instead of time-varying measurements

    
    
% #flag_test_glider_shear#
% 1-- [Default] Update the vertical velo shears in the SAMSgdr domain by
%       eliminating abnormal values right next to UMM4 and below 2000m;
% 0-- OFF;
flag_test_glider_shear = 1;


% % #flag_test_no_OOI_data# [Disabled on 7/23/2019]
% % 1-- Remove all of the OOI data from the calculations;
% % 0-- [Default] Load the OOI data whenever available;
% flag_test_no_OOI_data = 0;


% #flag_test_model_clim_case1#
% 1- Test with filling uv (or ts) climatology at the whole section;
% 0-- OFF;
flag_test_model_clim_case1 = 0;


% #flag_test_model_clim_case2#
% 1- Test with filling property climatology where mode_d== 1;
% 0-- OFF;
flag_test_model_clim_case2 = 0;


% #flag_test_OA_argoonly#
% 1-- Use OA with Argo input only;
% 0-- OFF;
flag_test_OA_argoonly = 0;

% % #flag_test_OA_noWHOIOUCglider# [Disabled on 7/23/2019]
% % 1-- Load an OA product that incoporates no data from OUC gliders;
% % 0-- OA product incorporating Argo, mooring and glider data;
% flag_test_OA_noWHOIOUCglider = 0;

% #flag_test_OA_noglider_for_OA#
% 1-- Load an OA product that incoporates no data from any gliders; 
% 0-- OFF;
flag_test_OA_noglider = 0;


% #flag_test_SAMSgdr_TS#
% 0-- [default] OA products incorporating SAMS glider TS;
% 1-- [test] TS information from SAMS glider survey not included in OA;
flag_test_OA_noSAMSgdr = 0;


% #flag_test_OA_depth#
% 0-- [default] OFF;
% 1-- Product created using depth_level methods;
flag_test_OA_depth= 0;


% #flag_test_Vgeo_only_in_CMarray#
% 1-- Calculate geostrophic velocities ONLY in areas with
%       current meter moorings (i.e., CF6-NOCM5),
%       with a follow-on step to apply the
%       topogrpahy-following velocity interpolation if necessary;
% 0-- OFF;
flag_test_Vgeo_only_in_CMarray = 0;


% #flag_test_AGV#
% 0-- [OFF] Use ADT data for calculating the surface geo. velocity;
% 1-- Use AGV product to provide the surface reference;
flag_test_AGV = 0;


% #flag_test_direct_velo_only# - need to update!
% 0-- OFF;
% 1-- Fill the areas WITHOUT direct velocity measurments with 
%       the time-mean velocities;
% 2-- Fill the SAMS' glider survey domain with the time-mean velocities;
% 3-- Fill the areas having direct velocity measurements with 
%       the time-mean velocities;
flag_test_direct_velo_only = 0;


% #flag_test_pden#
% 1-- neutral density;
% 2-- sigma2;
% 0-- OFF; sigma_theta;
flag_test_pden = 0;


% #flag_test_vbaro_mean#   
% 1-- Test of using a different v_baro_mean;
% 0-- OFF;
%
% NOTE
%   flag_test_vbaro_mean= 1 only works if flag_vbaro_mean= 1;
flag_test_vbaro_mean = 0;


% #flag_test_MC_West_epsilon#
% 1-- Test of using smaller epsilons to stop MC iteration at OSNAP-West;
% 0-- OFF;
flag_test_MC_West_epsilon = 0;


% % #flag_test_RTADCP# [Disabled on 7/23/2019]
% % 1-- Not to use any data from RTADCP;
% % 2-- Use the time-mean from RTADCP; 
% %     This requires modifications to   
% %       fun_osnap_mooring_vinterp_random_v2.m;
% % 0-- OFF;
% %
% % NOTE
% %   This flag is not working when using regridded profiles
% %       (flag_regridded_profiles= 1);
% %   For 2014-2018 runs, not meaningful; gaps have been filled with 
% %       the time-mean; KEEP it to 0!
% %
% flag_test_RTADCP = 0;


% % #flag_test_noRR# [Disabled on 7/23/2019]
% % 1-- No data will be used from the moorings at two sides of the RR;
% % 2-- In addition to (1), ignore the RR when calculating the shears;
% % 0-- OFF;
% flag_test_noRR = 0;


% #flag_test_fewermr# 
% 1-- Withhold the following moorings from the calculation:
%       C3, K8, LS6, LS4, LS2, CF1, CF3, CF5, CF7, IC2, IC4;
%     So to keep all the deep moorings and dynamic height moorings (no need
%     to modify the time-mean barotropic velocities);
% 0-- OFF;
flag_test_fewermr= 0;


% #flag_test_C3#
% 1-- Use time-mean observations;
% 2-- Copied from K7;
% 3-- Simply drop C3 data;
% 4-- Reconstructed based on K7;
% 0-- OFF;
flag_test_C3= 0;


% #flag_test_SAMSgdr_geov#
% 0-- [default] Not to use velocity information (DAC) from SAMS glider
% survey;
% 1-- [test] In the glider domain, replacing the surface altimetry velocity
% with the mean surface geostrophic velocity derived from multiple glider
% transects during 2014-2016 (Loic et al. 2018 JGR);
flag_test_SAMSgdr_geov = 0;



% #flag_test_M2M3_geo#
% 0- OFF;
% 1- [Default] Removing D4,D5; Add the deepest MCTD on D4 to M2; 
%       Calculating the geostrophic velocities between M2-M3;
flag_test_M2M3_geo = 1;


% #flag_test_pd_MOCwest#
% 1- Add a density range when searching for the maximum of the
%       streamfunction at OSNAP West;
flag_test_pd_MOCwest = 0;



% #flag_test_AllWE#
% 1- Calculate fluxes for West and East as well for 'all';
flag_test_AllWE = 1;

% #flag_test_IB5#
% 1- [Default] include IB5 in thermal wind calculation by shifting the
%           western extent of all grid shear position for glider to IB5
% 0- use UMM4/IB4 as the western extent of all grid shear calculation
flag_test_IB5=1;

% #flag_test_UMD123_geo
% 0- [Default] change nothing. D123 are tall moorings for 2018-2020 period,
%    and directly measured velocity is used in the calculation
% 1- remove the velocity observation above 1200 m for 2018-2020, and fill
%    the gap between M1 and M2 with geostrophy
flag_test_UMD123_geo=0;

flag_test_remove_NOCM5=0;

flag_test_LC_TS_clim=0;


% #flag_test_LSAB
% 0- [Default] include LSAB in both velocity field and T/S for MHT/MFT in
%       201809-202006
% 1- remove LSAB in 201809-202006
flag_test_remove_LSAB=0;


% LOG
disp('+');
disp('+');
disp('+-------------+');
disp('| Set flags   |');
disp('+-------------+');
disp('+');
disp(['+    flag_section= ',flag_section]);
disp(['+    flag_data_interval= ',num2str(flag_data_interval)]);
disp(['+    flag_vcomp_constraint= ',num2str(flag_vcomp_constraint)]);
disp(['+    flag_vcomp_everywhere= ',num2str(flag_vcomp_everywhere)]);
disp(['+    flag_vcomp_selected= ',num2str(flag_vcomp_selected)]);
disp('+');
disp(['+    flag_regridded_profiles= ',num2str(flag_regridded_profiles)]);
disp(['+    flag_MC= ',num2str(flag_MC)]);
disp(['+    flag_MC_MOC= ',num2str(flag_MC_MOC)]);
disp(['+    flag_MC_SE= ',num2str(flag_MC_SE)]);
disp('+');
disp(['+    flag_topo_interp= ',num2str(flag_topo_interp)]);
disp(['+    flag_OA= ',num2str(flag_OA)]);
disp(['+    flag_moc_def= ',num2str(flag_moc_def)]);
disp('+');
disp(['+    flag_full_grid_shears= ',num2str(flag_full_grid_shears)]);
disp(['+    flag_OOI_for_tw= ',num2str(flag_OOI_for_tw)]);
% disp(['+    flag_test_no_OOI_data= ',num2str(flag_test_no_OOI_data)]);
disp('+');
disp(['+    flag_baro_corr= ',num2str(flag_baro_corr)]);
disp(['+    flag_ssh_mean= ',num2str(flag_ssh_mean)]);
disp(['+    flag_vbaro_mean= ',num2str(flag_vbaro_mean)]);
disp(['+    flag_full_grid_baro_corr= ',num2str(flag_full_grid_baro_corr)]);
disp(['+    flag_aviso_smoothing= ',num2str(flag_aviso_smoothing)]);
disp('+');
disp(['+    flag_vref_shortmoorings= ',num2str(flag_vref_shortmoorings)]);
disp(['+    flag_keep_vrefshortmoorings_all= ',num2str(flag_keep_vrefshortmoorings_all)]);
disp(['+    flag_keep_vrefshortmoorings_selected= ',num2str(flag_keep_vrefshortmoorings_selected)]);
disp('+');
disp(['+    flag_vary_dek= ',num2str(flag_vary_dek)]);
disp(['+    flag_opa_clim= ',num2str(flag_opa_clim)]);
disp(['+    flag_LC_UV= ',num2str(flag_LC_UV)]);
disp(['+    flag_LC_TS= ',num2str(flag_LC_TS)]);
disp(['+    flag_C3_recon= ',num2str(flag_C3_recon)]);
disp(['+    flag_mean_RTADCP= ',num2str(flag_mean_RTADCP)]);
disp('+');
disp(['+    flag_HT_FWT_decomposition= ',num2str(flag_HT_FWT_decomposition)]);
disp('+');
disp(['+    flag_test_glider_shear= ',num2str(flag_test_glider_shear)]);
disp(['+    flag_test_mr_nouse= ',num2str(flag_test_mr_nouse)]);
disp(['+    flag_test_d5_nouse= ',num2str(flag_test_d5_nouse)]);
disp(['+    flag_test_d5_mean= ',num2str(flag_test_d5_mean)]);
disp(['+    flag_test_rrex_nouse= ',num2str(flag_test_rrex_nouse)]);
disp(['+    flag_test_model_clim_case1= ',num2str(flag_test_model_clim_case1)]);
disp(['+    flag_test_model_clim_case2= ',num2str(flag_test_model_clim_case2)]);
disp(['+    flag_test_OA_argoonly= ',num2str(flag_test_OA_argoonly)]);
disp(['+    flag_test_OA_noglider= ',num2str(flag_test_OA_noglider)]);
disp(['+    flag_test_Vgeo_only_in_CMarray= ',num2str(flag_test_Vgeo_only_in_CMarray)]);
disp(['+    flag_test_AGV= ',num2str(flag_test_AGV)]);
disp(['+    flag_test_direct_velo_only= ',num2str(flag_test_direct_velo_only)]);
disp(['+    flag_test_pden= ',num2str(flag_test_pden)]);
disp(['+    flag_test_vbaro_mean= ',num2str(flag_test_vbaro_mean)]);
disp(['+    flag_test_MC_West_epsilon= ',num2str(flag_test_MC_West_epsilon)]);
% disp(['+    flag_test_RTADCP= ',num2str(flag_test_RTADCP)]);
% disp(['+    flag_test_C3= ',num2str(flag_test_C3)]);
disp(['+    flag_test_OA_noSAMSgdr= ',num2str(flag_test_OA_noSAMSgdr)]);
disp(['+    flag_test_SAMSgdr_geov= ',num2str(flag_test_SAMSgdr_geov)]);
disp(['+    flag_test_OA_depth= ',num2str(flag_test_OA_depth)]);
disp(['+    flag_test_LC_UV_clim= ',num2str(flag_test_LC_UV_clim)]);
disp(['+    flag_test_M2M3_geo= ',num2str(flag_test_M2M3_geo)]);
disp(['+    flag_test_pd_MOCwest= ',num2str(flag_test_pd_MOCwest)]);
disp(['+    flag_test_AllWE= ',num2str(flag_test_AllWE)]);
disp(['+    flag_test_IB5= ',num2str(flag_test_IB5)]);

disp('+');
% disp(['+    [Array optimization] flag_test_noRR= ',num2str(flag_test_noRR)]);
% disp(['+    [Array optimization] flag_test_fewermr= ',num2str(flag_test_fewermr)]);
disp('+');



% -------------------------
% Checking flags
% -------------------------

if flag_test_M2M3_geo== 1 && flag_test_d5_nouse== 1
    error('check for flag_test_d5_nouse= 0!')
end

if flag_aviso_smoothing== 1 && flag_full_grid_baro_corr== 0
    error('check flag_aviso_smoothing!');
end



if flag_test_direct_velo_only== 1 && flag_vcomp_selected== 1
    error('flag_vcomp_selected= 0????');
end


% flag_vcomp_everywhere must
% be 0 for flag_vcomp_selected to be effective;
if flag_vcomp_everywhere == 1 && flag_vcomp_selected== 1
    error('flag_vomp_everywhere= 0!');
end


% One of three sections should be used;
if strcmp(flag_section,'east') || ...
        strcmp(flag_section,'west') || ...
        strcmp(flag_section,'all')
else
    error('SEC 180!');
end


% Make sure thermal wind shears are referenced to 
% the same reference as they were when calculating 
% the mean barotropic component;
if flag_vbaro_mean == 1 && flag_vref_shortmoorings == 0
    error('flag_vref_shortmoorings= 1 ???');
end



% flag_test*glider* are used to modify the OA input that excludes glider
% data and the correponding flag_full_grid_shears;
% if flag_test_OA_noWHOIOUCglider== 1 && flag_full_grid_shears~= 0
%     
%     error('flag_full_grid_shears= 0 ???');
% end
% if flag_test_OA_noglider== 1 && flag_full_grid_shears~= 2
%     
%     error('flag_full_grid_shears= 2 ???');
%     
% end



% Must use the bottom mooring reference before choosing to keep the
% mooring-velocity referenced profiles;
if flag_vref_shortmoorings == 0 && (flag_keep_vrefshortmoorings_all == 1 ...
        || flag_keep_vrefshortmoorings_selected == 1)
        
    disp('Reference to NO deep moorings?');
    
end



% Make sure only selected areas are usd for keeping the deep mooring
% reference ...
if flag_keep_vrefshortmoorings_selected== 1 && flag_keep_vrefshortmoorings_all== 1
    error('flag_keep_vrefshortmoorings_all must be 0 if flag_keep_vrefshortmoorings_selected is 1!');
end



% If no bartropic correction applied ...
if flag_baro_corr == 0 && (flag_keep_vrefshortmoorings_all== 1 || ...
        flag_keep_vrefshortmoorings_selected== 1)
      
    disp('NOT USED: flag_keep_vrefshortmoorings_all/flag_keep_vrefshortmoorings_selected ...');
    
end


% MC simulations cannot operate on daily interval
if flag_MC == 1 && flag_data_interval == 1
    error('flag_data_interval must > 1 for MC simulations!');
    
end


% MC simulations draw random observations from each instrument
if flag_MC == 1 && flag_regridded_profiles == 1
    disp('NOTE #flag_regridded_profiles= 1# randomly drawn daily profiles used for MC simulations!');
end


% For mean MOC calculations (w/o MC), use averaged profiles...
if flag_MC == 0 && flag_regridded_profiles == 0
    
    disp(' ');
    disp('NOTE: Use pre-calculated daily regriddd profiles may save some time!');
end



% When using v_baro_mean, all geostrophic profiles should be
% referened to velocities from the deep moorings whenever available;
if flag_vbaro_mean== 1 && flag_vref_shortmoorings== 0
    error('Incompatability found!  flag_vbaro_mean vs. flag_vref_shortmoorings');
end




%------------------------+
%% PREPARE DATA          |
%------------------------+

% Constants
er = 6370e3; % Earth radius [m]
gr = 9.8; % gravity acceleration [m^2/s]
rho0 = 1027;  % reference density [kg/m^3]
% rhoCp = 4.1e6; % [J/m^3/C^1]

disp('+');
disp('+');
disp('+---------------+');
disp('| Prepare input |');
disp('+---------------+');


% Main working directory
dirio = './';
dirout= './outputs/results/';




% -------------------------------------------------------------------------
%% (A) Set up coordinates
% -------------------------------------------------------------------------

% =========================================================================
% flag_test_pden ==========================================================
% =========================================================================
switch flag_test_pden
    case 0
        % sigma-axis
        disp('+ Setting up sigma-axis of [23.3:0.01:28.1]');
        dens_bins = (1023.3:0.01:1028.1)';
    
    case 1
        % neutral density
        disp('+ Setting up sigma-axis of [23.3:0.01:28.1]');
        dens_bins = (1023.3:0.01:1028.1)';
    
	case 2
        % sigma2
        disp('+ Setting up sigma-axis of [30:0.01:38]');
        dens_bins = (1030:0.01:1038)';
        
    otherwise
        error('DEN 013!');
    
end
% =========================================================================
% =========================================================================
sz_d = length(dens_bins);

% x-axis (osnap grid)
[lon,lat,~,depth] = create_osnap_grid_2019(flag_section,0.25);
sz_n = length(lon);
sz_z = length(depth);
disp('+ Setting up x- and z-axis using create_osnap_grid.m');


% time-axies @ 00:00:00
disp('+ Setting up time-axis...');
% time = datenum(2014,7,1,0,0,0):flag_data_interval:datenum(2020,6,28,0,0,0); 
time = datenum(2014,1:84,1,12,0,0); % 20140101-20201201
sz_t_all = length(time);



% -------------------------------------------------------------------------
%% (B) Load ADT, WIND, SST data
% -------------------------------------------------------------------------
% Load data...
% ** SST ** 
    data_sect_SST = [dirio,'data/SST/20210511_SST_OSNAP_All_20140101-20201231.mat'];  
%     data_sect_SST = [dirio,'data/SST/0716_SST_OSNAP_All_20140701-20180630.mat'];
    disp(['+ data_sect_SST= ',data_sect_SST]);
    load(data_sect_SST);

% ** ADT ** 
    data_sect_ADT = [dirio,'data/ADT/20210512_ADT_OSNAP_All_20140101-20201231.mat']; 
%     data_sect_ADT = [dirio,'data/ADT/0628_ADT_OSNAP_All_20140701-20180630.mat']; 
    disp(['+ data_sect_ADT= ',data_sect_ADT]);
    load(data_sect_ADT,'osnap_uvh');   
    
% ** AGV ** 
    data_sect_AGV = [dirio,'data/ADT/20210512_ADT_OSNAP_All_20140101-20201231.mat']; 
%     data_sect_AGV = [dirio,'data/ADT/0628_ADT_OSNAP_All_20140701-20180630.mat']; 
        
    if flag_test_AGV== 1
        % LOG
        disp(['+ data_sect_AGV= ',data_sect_AGV]);
    end
    
    load(data_sect_AGV);    
    osnap_uvh.agv = osnap_uvh.v; % dummy;
    
        
% ** WIND ** 
%     data_sect_WIND = [dirio,'data/wind/20210511_WIND_OSNAP_All_20140101-20201231.mat'];
%     data_sect_WIND = [dirio,'data/wind/20211111_WIND_OSNAP_All_20140101-20201231.mat'];
    data_sect_WIND = [dirio,'data/wind/20211220_WINDstress_OSNAP_All_20140101-20201231.mat'];
%     data_sect_WIND = [dirio,'data/wind/0628_WIND_OSNAP_All_20140701-20180630.mat'];
    disp(['+ data_sect_WIND= ',data_sect_WIND]);
    load(data_sect_WIND);
    
   
%  Temporally averaging when needed ...
if flag_data_interval == 1
    data_adt = osnap_uvh;
    data_agv = osnap_uvh; 
    data_wind = osnap_winds;  data_wind.time = data_wind.time - 0.5; % change to 00:00:00 each day
    data_sst = osnap_sst;
    
else
    [data_adt, data_agv, data_wind, data_sst] = fun_osnap_sst_adt_winds_ave(time(1),time(end)+flag_data_interval_4calculation,flag_data_interval,osnap_uvh,osnap_uvh,osnap_winds,osnap_sst);
%     save ./temp/temp_adt_agv_wind_sst data_adt data_agv data_wind data_sst
% load temp/temp_adt_agv_wind_sst
    % //check point// all output time-axis are at 00:00:00 every day
end
% return

% =========================================================================
% FLAG_SSH_MEAN ===========================================================
% =======================================================================
% Use time-mean SSH
if flag_ssh_mean== 1 && flag_baro_corr== 1

    % SSH over the OSNAP time period
    ipt1 = osnap_uvh.time>= time(1) & osnap_uvh.time<= time(end)+ flag_data_interval_4calculation; 
    time_in = osnap_uvh.time(ipt1);
    ssh_in = osnap_uvh.h(:,ipt1);
    
    % initializatiion output
    ssh_out = nan(sz_n,sz_t_all,4);

    
% LOOP SECTION    
    for i_n = 1:sz_n
        
        % time series
        tmp = ssh_in(i_n,:);
        
        if all(isnan(tmp))
            continue;
        end
            
        
        % INTEGRAL TIME SCALE
        % determining the integral time scale
        % by integrating the autocorrelation function
        % to the first zero corssing [Thomson and Emery, 2014];
        ii = isnan(tmp);
        if ~isempty(ii)
            tmp(ii) = interp1(time_in(~ii),tmp(~ii),time_in(ii),'linear');
        end
        
        % size of one time step [days]
        delta_t = nanmean(diff(time_in)); 

        % autocorrelations
        acorr = xcov(tmp,tmp,'coeff');

        % number of valid observations
        N = length(tmp);

        % integrating the autocorrelation 
        % function to the
        % first zero crossing
        acorr1 = acorr(N:end); % one-sided autocorrelation function

        ip1 = find(acorr1< 0, 1, 'first') - 1;

        if ip1== 0
            tau_i = delta_t/2*acorr1(1);
        elseif ip1> 0
            tau_i = delta_t/2*sum((acorr1(1:ip1-1) + acorr1(2:ip1)));
        end
            
        ssh_out(i_n,:,4) = tau_i;
        
        
        % STD and STE
        ssh_out(i_n,:,1) = nanmean(tmp); % averages
        tmp_sd = nanstd(tmp); 
        ssh_out(i_n,:,2) = tmp_sd;

        % SE = SD/sqrt(DOFeff), at compatible time scale
        ssh_out(i_n,:,3) = tmp_sd./sqrt(length(time_in)/tau_i/2);
  
% END LOOP SECTION
    end
   
    
    % upate
    data_adt.time = time;
    data_adt.lon = osnap_uvh.lon;
    data_adt.lat = osnap_uvh.lat;
    data_adt.h  = ssh_out;
    
    clear ipn ssh_in ssh_out time_in dofeff tau_i tmp_sd delta_t
    % //// check-point ////
end
% =========================================================================
% =========================================================================



% Prepare for subsections...

% *** SST ***
% =========================================================================
% flag_regridded_profiles (1/2) ===========================================
% =========================================================================
% SST data used for regridding temperature
% profiles;
if flag_regridded_profiles== 0
    lon_sst = data_sst.lon;
    switch flag_section
        case 'east'
            ipx = find(lon_sst>=-44);
        case 'west'
            ipx = find(lon_sst<=-45);
        case 'all'
            ipx = 1:length(lon_sst);
    end
    data_sst.lon = lon_sst(ipx);
    data_sst.sst = data_sst.sst(ipx,:,:);
    save('TemporaryFile_sst_ave.mat','data_sst'); % to be used in fun_osnap_mooring_vinterp_random_v*.m
end
% =========================================================================
% =========================================================================


% *** ADT ***
time_adt = data_adt.time;  
lon_ssh = data_adt.lon;
switch flag_section
    case 'east'
        ipx = find(lon_ssh>=-44);
    case 'west'
        ipx = find(lon_ssh<=-45);
    case 'all'
        ipx = 1:length(lon_ssh);
end
lon_ssh = lon_ssh(ipx);
ssh = data_adt.h(ipx,:,:);
   

% *** AGV ***
% =========================================================================
% FLAG_TEST_AGV ===========================================================
% =========================================================================
% Load AVISO surface geostrophic velocities
if flag_test_AGV== 1
    
    % data_agv
    time_agv = data_agv.time;  
    lon_agv = data_agv.lon;
    switch flag_section
        case 'east'
            ipx = find(lon_agv>=-44);
        case 'west'
            ipx = find(lon_agv<=-45);
        case 'all'
            ipx = 1:length(lon_agv);
    end
    lon_agv = lon_agv(ipx);
    agv = data_agv.agv(ipx,:,:);

end
% =========================================================================


% *** WIND ***
time_wind = data_wind.time;
lon_wind = data_wind.lon;
switch flag_section
    case 'east'
        ipx = find(lon_wind>=-44);
    case 'west'
        ipx = find(lon_wind<=-45);
    case 'all'
        ipx = 1:length(lon_wind);
end
lon_wind = lon_wind(ipx);
taux = data_wind.taux(ipx,:,:);
tauy = data_wind.tauy(ipx,:,:);
u10 = data_wind.u10(ipx,:,:);
v10 = data_wind.v10(ipx,:,:);



% -------------------------------------------------------------------------
%% (C) Load model data [UVTS]
% -------------------------------------------------------------------------
data_sect_model_clim = [dirio,'data/climatology/1122_FLAME_clim_mm_OSNAP_All.mat'];
load(data_sect_model_clim);
disp(['+ data_sect_MODEL_CLIM= ',data_sect_model_clim]);

% For individual section
lon_clim_model = flame_clim.lon;
switch flag_section
    case 'east'
        ipx = find(lon_clim_model>=-44);
    case 'west'
        ipx = find(lon_clim_model<=-45);
	case 'all'
        ipx = 1:length(lon_clim_model); % keep all grid points
end
        
% Copy data for later use
lon_clim_model = lon_clim_model(ipx);
t_clim_model = flame_clim.t(ipx,:,:); % temperature (C)
s_clim_model = flame_clim.s(ipx,:,:);
u_clim_model = flame_clim.u(ipx,:,:);
v_clim_model = flame_clim.v(ipx,:,:);



% -------------------------------------------------------------------------
%% (E) Load monthly WOA18 climatology [TSRHO]
% -------------------------------------------------------------------------
% This is to use for the MHT/MFT calculations 
%   only if no OA product is provided;
% data_sect_TS_clim= [dirio,'outputs/OA/1202_WOA18andCTD_TS_OSNAP_All_200levels.mat'];
data_sect_TS_clim= [dirio,'data/WOA/20220125_WOA18andCTD_TS_OSNAP_All_200levels.mat'];
load(data_sect_TS_clim);
disp(['+ data_sect_TS_CLIM= ',data_sect_TS_clim]);

% longitude
lon_clim = woa18_sect_updated.lon;

% FLAG_SECTION
switch flag_section
    case 'west'
        ipx = lon_clim <= -45;
    case 'east'
        ipx = lon_clim >= -44;
    case 'all'
        ipx = 1:length(lon_clim); % keep all grid points
% END FLAG_SECTION            
end

lon_clim = woa18_sect_updated.lon(ipx); % [1,grid]
lat_clim = woa18_sect_updated.lat(ipx); % [1,grid]
t_clim = woa18_sect_updated.temp(ipx,:,:); % [grid,depth,month]
s_clim = woa18_sect_updated.salt(ipx,:,:); % [grid,depth,month]
depth_clim = woa18_sect_updated.depth;




% =========================================================================
% FLAG_OA =================================================================
% =========================================================================
if flag_OA == 1

    % *********************************************************************
    % (D) Load 30-day OA T/S/dens/pden/ptmp fields
    % *********************************************************************
    % Latest OA product using Argo, OSNAP moorings and gliders data
%     data_sect_OA = [dirio,'/OA/1202_OSNAP_All_OA_201407-201806_ArgoMooringGlider.mat']; 

    data_sect_OA = [dirio,'outputs/OA/20211205_OSNAP_All_OA_201401-202012_ArgoMooringGlider_vmm.mat'];
%     data_sect_OA = [dirio,'outputs/OA/20210616_OSNAP_All_OA_201401-202012_ArgoMooringGlider_v30day_para.mat'];  % updated Argo floats after 3/5/2021
%     data_sect_OA = [dirio,'data/intermediate_output/20200305_OSNAP_All_OA_201407-201806_ArgoMooringGlider.mat'];  % updated Argo floats after 3/5/2020

    
    
% =======================================================================
% flag_test_OA_depth =================================================
% =======================================================================
    % Testing of using OA product created on depth levels
    if flag_test_OA_depth== 1
        
%         data_sect_OA = [dirio,'/OA/20191121_OSNAP_All_OA_201407-201604_ArgoMooringGlider.mat'];  
        data_sect_OA = [dirio,'outputs/OA/20200317_OSNAP_All_d epOA_201407-201806_ArgoMooringGlider.mat'];
        
    end
% ========================================================================
    
    
% =======================================================================
% flag_test_OA_argoonly =================================================
% =======================================================================
    % Testing of using OA product with Argo input only
    if flag_test_OA_argoonly== 1
        
        % For testing purpose only (ArgoOnly OA)
        data_sect_OA = [dirio,'outputs/OA/20210601_OSNAP_All_OA_201401-202012_ArgoOnly_vmm.mat'];  
%  data_sect_OA = ['/Volumes/TOSHIBA/Work/osnap_2019_test/outputs/OA/1106_OSNAP_All_OA_201407-201806_ArgoOnly_v2.mat'];  
    end
% ========================================================================


% ========================================================================
% flag_test_OA_noSAMSgdr =================================================
% ========================================================================
    if flag_test_OA_noSAMSgdr == 1

        % An OA product using Argo + Moorings + ONLY WHOI/OUC glider data
        data_sect_OA = [dirio,'/OA/0115_OSNAP_All_OA_201407-201806_ArgoMooringOUCGlider.mat'];    

    end
% ========================================================================


% ========================================================================
% FLAG_TEST_OA_NOGLIDER ==================================================
% ========================================================================
    if flag_test_OA_noglider == 1

        % An OA product using Argo + Mooring data, and no glider data
        % incorporated;
        data_sect_OA = [dirio,'outputs/OA/20211105_OSNAP_All_OA_201401-202012_ArgoMooring_vmm.mat'];      

    end
% ========================================================================



    disp(['+ data_sect_OA= ',data_sect_OA]);
    load(data_sect_OA);
    lon_oa = ts_analysis_mm.lon;

    switch flag_section
        case 'east'
            ipx = find(lon_oa>=-44);
        case 'west'
            ipx = find(lon_oa<=-45);
        case 'all'
            ipx = 1:length(lon_oa);
    end


    % save
    time_analysis = ts_analysis_mm.time - 0.5; % changed to 00:00:00 each day to match all other time axis
    t_analysis = ts_analysis_mm.t(ipx,:,:);
    s_analysis = ts_analysis_mm.s(ipx,:,:);
    d_analysis = ts_analysis_mm.d(ipx,:,:);
    pd_analysis = ts_analysis_mm.pd(ipx,:,:);
    pt_analysis = ts_analysis_mm.pt(ipx,:,:);
    
    if nanmedian(pd_analysis(:))< 100
        pd_analysis = pd_analysis+1000;
    end


% END IF FLAG_OA
end
% =========================================================================



% -------------------------------------------------------------------------
%% (F) Load data for the unmeasured Labrador Current componnet
% -------------------------------------------------------------------------

% =========================================================================
% FLAG_OPA_CLIM (1/2) =====================================================
% =========================================================================
% Fill the inshore unmeasured LC componenet
% with OPA model climatology;

% IF FLAG_OPA_CLIM
if flag_opa_clim== 1 
    
    % (A2) Load model climatology of u/v/t/s
    % above the Labrador shelf and west of OM2;
    if strcmp(flag_section,'west') || strcmp(flag_section,'all')

        data_sect_model_clim = [dirio,'data/climatology/0408_OPA_model_clim_Labrador_Current.mat'];
        load(data_sect_model_clim);
        disp(['+ data_sect_UV for Labrador Current= ',data_sect_model_clim]);

        % Copy data for later use
        lon_LC_UV_clim= opa_clim.lon;
        lat_LC_UV_clim= opa_clim.lat;

        u_LC_model_clim = opa_clim.u;
        v_LC_model_clim = opa_clim.v;

%         t_clim_model_OPA = opa_clim.t; 
%         s_clim_model_OPA = opa_clim.s;
%         sz1= length(lon_LC_UV_clim);
% 
%         % Derive dens, pden and ptmp
%         % from climatological T/S t
%         pt_clim_model_OPA = nan(sz1,200,12);
%         d_clim_model_OPA = nan(sz1,200,12);
%         pd_clim_model_OPA = nan(sz1,200,12);
% 
%         for it = 1:12
%             pres = gsw_p_from_z(-repmat(opa_clim.depth',sz1,1),lat_LC_UV_clim); % pres
%             [SA,~,~] = gsw_SA_Sstar_from_SP(s_clim_model_OPA(:,:,it),pres,lon_LC_UV_clim,lat_LC_UV_clim); % absolute salinity 
%             pt_clim_model_OPA(:,:,it) = gsw_pt_from_t(SA,t_clim_model_OPA(:,:,it),pres,0); % potential temperature  
%             d_clim_model_OPA(:,:,it) = gsw_rho_t_exact(SA,t_clim_model_OPA(:,:,it),pres); % density 
%             pd_clim_model_OPA(:,:,it) = gsw_pot_rho_t_exact(SA,t_clim_model_OPA(:,:,it),pres,0); % potential density
%         end

        clear sz1

    end
    
% END IF FLAG_OPA_CLIM    
end
% =========================================================================



% =========================================================================
% FLAG_LC_UV (1/2) =====================================================
% =========================================================================
% Fill the inshore unmeasured LC componenet UV with NEMO output
if flag_LC_UV == 1 
    
    % (A2) Load model climatology of u/v/t/s
    % above the Labrador shelf and west of OM2;
    if strcmp(flag_section,'west') || strcmp(flag_section,'all')

        % Load data
        data_sect_model = [dirio,'data/climatology/0723_NEMO_5d_OSNAP_West_2014010-20171226_200levels.mat'];
        load(data_sect_model,'nemo_labsea');
        disp(['+ data_sect_UV for Labrador Current= ',data_sect_model]);

        % Copy data for later use
%         lon_LC_uv= nemo_labsea.lon;
%         lat_LC_uv= nemo_labsea.lat;
        time_LC_uv = nemo_labsea.time;

        u_LC_model = nemo_labsea.u;
        v_LC_model = nemo_labsea.v;
    end
    
    
    
elseif flag_LC_UV == 2
    
    if strcmp(flag_section,'west') || strcmp(flag_section,'all')

        % Load data
        data_sect_model = [dirio,'/climatology/1117_GloSea5_mm_TSV_LabSeaShelf_201401-201812_200levels.mat'];
        load(data_sect_model,'glosea5_labshelf');
        disp(['+ data_sect_UV for Labrador Current= ',data_sect_model]);

        % Copy data for later use
%         lon_LC_uv= glosea5_labshelf.lon;
%         lat_LC_uv= glosea5_labshelf.lat;
        time_LC_uv = glosea5_labshelf.time;

        v_LC_model = glosea5_labshelf.v;
        u_LC_model = nan(size(v_LC_model)); % note, v_LC_model itself to use
    end
    
end
% =========================================================================
% =========================================================================

    
% =========================================================================
% FLAG_LC_TS(1/2) =========================================================
% =========================================================================
% Fill the inshore unmeasured LC componenet TS from GloSea5
if flag_LC_TS == 1 
    
    % Load reanalysis TS
    % above the Labrador shelf 
    % and west of OM2;
    if strcmp(flag_section,'west') || strcmp(flag_section,'all')

        % Load data
        data_sect_model = [dirio,'/climatology/1117_GloSea5_mm_TSV_LabSeaShelf_201401-201812.mat'];
        load(data_sect_model,'glosea5_labshelf');
        disp(['+ data_sect_TS for Labrador Current = ',data_sect_model]);

        % Copy data for later use
%         lon_LC_rea= glosea5_labshelf.lon;
%         lat_LC_rea= glosea5_labshelf.lat;
        time_LC_rea = glosea5_labshelf.time;

        t_LC_rea = glosea5_labshelf.t;
        s_LC_rea = glosea5_labshelf.s;
    end
end
% =========================================================================
% =========================================================================



% =========================================================================
% flag_vbaro_mean (1/2) ===================================================
% =========================================================================
% Load long-term mean barotropic velocities on I-Grid
%
% NOTE
%   v_baro_mean should be obtained by referencing
%   geostrophic velocities to the bottom (a bottom level of no 
%   motion or velocity at the top of deep moorings).  
%   So, for compatibility, when using this version, 
%   flag_vref_shotmoorings sould be 1.
if flag_vbaro_mean == 1


% IF FLAG_VCOMP_CONSTRAINT
        if flag_vcomp_constraint == 1
            % d1, W= -1.6 & E= 0.6;
%             data_sect_mean_barotropic = [dirio,'/20180520_test1_v_baro_mean.mat']; 

        elseif flag_vcomp_constraint == 2
            % d1, W= -1.6 & E= 1.6;
            
%             data_sect_mean_barotropic = [dirio,'/20191112_test1_v_baro_mean_first47mon.mat'];  % 2014-2018
%             data_sect_mean_barotropic = [dirio,'/20200117_test1_v_baro_mean_first47mon.mat']; % including meanD5, d1
%             data_sect_mean_barotropic = [dirio,'/20200117_test2_v_baro_mean_first47mon.mat']; % excluding D5, d1
%             data_sect_mean_barotropic = [dirio,'/20200210_test1i__v_baro_mean_first47mon.mat']; % including meanD5 & multLC, d1
%               data_sect_mean_barotropic = [dirio,'/20200212_test1i__v_baro_mean_first47mon.mat']; % Post-OS2020 configuration
%             data_sect_mean_barotropic = [dirio,'/20200305_test1i__v_baro_mean.mat']; % Post-OS2020 configuration, updated OA_20200305

%             data_sect_mean_barotropic = [dirio,'data/mean_barotropic_velocity/20200324_test1i__v_baro_mean.mat']; % Post-OS2020 configuration, updated OA_20200305
%             data_sect_mean_barotropic = [dirio,'outputs/results/20211111_test_v_baro_mean.mat']; % updated OA_20211110 from monthly mean
              data_sect_mean_barotropic = [dirio,'outputs/results/20220215_test_v_baro_mean_M5corrected.mat']; % updated routine to adopt NOCM5 removal after 201806
            
            
% =======================================================================
% flag_full_grid_baro_corr ==============================================
% =======================================================================
            if flag_full_grid_baro_corr== 1
                % d1, W= -1.6 & E= 1.6;
                % & full-grid baro corr applied;
                
                error('EJK 190!');
                
%                 if flag_aviso_smoothing== 1
%                     % altimetry data are averaged over 100km
% %                     data_sect_mean_barotropic = [dirio,'/20180817_test1_v_baro_mean.mat']; 
%                     
%                 elseif flag_aviso_smoothing== 0
%                     % altimetry data are kept at original resolution
% %                     data_sect_mean_barotropic = [dirio,'/20180823_test1_v_baro_mean.mat']; 
%                 end
                
            end
% =======================================================================      

% % =======================================================================
% % flag_test_noRR ========================================================
% % =======================================================================
%             if flag_test_noRR== 1
%                 % d1 & RR removed (large bottom triangles)
% %                 data_sect_mean_barotropic = [dirio,'/20181203_test1_v_baro_mean.mat']; 
%                 
%             elseif flag_test_noRR== 2
%                 % d1 && RR removed and ignored when calculating the shears
% %                 data_sect_mean_barotropic = [dirio,'/20181203_test1a_v_baro_mean.mat']; 
%                 
%             end
% % =======================================================================


% END IF FLAG_VCOMP_CONSTRAINT
        end
            

% =========================================================================
% flag_test_vbaro_mean ====================================================
% =========================================================================
% Alternative long-term mean barotropic velocities
%   used for test runs;
%
    if flag_test_vbaro_mean == 1   % 0

        % For the first 21 months
%         data_sect_mean_barotropic = [dirio,'/20200122_test1i_v_baro_mean_first21mon.mat'];  % 2014-2016; CONTROL
%         data_sect_mean_barotropic = [dirio,'/20200122_test2i_v_baro_mean_first21mon.mat'];  % 2014-2016; TEST OA
%         data_sect_mean_barotropic = [dirio,'/20200122_test3i_v_baro_mean_first21mon.mat'];  % 2014-2016; TEST D4.D5
        


        if flag_test_SAMSgdr_geov== 1
            data_sect_mean_barotropic = [dirio,'/20200116_test2_F_v_baro_mean_first47mon.mat'];
        end

        if flag_test_OA_depth== 1
            data_sect_mean_barotropic = [dirio,'/20200331_test1i__v_baro_mean.mat'];
        end


%             if flag_test_LC_UV_clim~= 0
%                 data_sect_mean_barotropic = [dirio,'/20200206_test1ir_v_baro_mean_first21mon.mat'];
%             end


        if flag_test_d5_nouse== 1
            data_sect_mean_barotropic = [dirio,'/20200331_test2i__v_baro_mean.mat'];
        end


        if flag_test_d5_mean== 1
%                 data_sect_mean_barotropic = [dirio,'/20200210_test1i__v_baro_mean_first47mon.mat'];
            data_sect_mean_barotropic = [dirio,'/20200213_test1i__v_baro_mean.mat'];
        end


% END IF flag_test_vbaro_mean        
    end

    
    
    % LOG
    disp(['+ data_sect_mean_barotropic= ',data_sect_mean_barotropic]);
    
    % data
    load(data_sect_mean_barotropic);
    data_tmp = v_baro_mean;
    clear v_baro_mean
    
    

    switch flag_section
    case 'east'
        ipx = find(data_tmp.lon >= -44);
    case 'west'
        ipx = find(data_tmp.lon <= -45);
    case 'all'
        ipx = 1:length(data_tmp.lon);
    end
    
    v_baro_mean(:,1) = data_tmp.mean(ipx);
    v_baro_mean(:,2) = data_tmp.sd(ipx);
    v_baro_mean(:,3) = data_tmp.se(ipx);

end
% =========================================================================
% =========================================================================  



% -------------------------------------------------------------------------
%% Set up the flux grid (I-Grid)
% -------------------------------------------------------------------------
% display(['+ Creating the flux grid ...']);


% longitude and latitude
lon_i = 0.5*(lon(1:end-1)+lon(2:end))'; % [nx1]
lat_i = 0.5*(lat(1:end-1)+lat(2:end))'; % [nx1]


% angle of each grid
d_x = (lon(2:end)-lon(1:end-1))'.*(pi/180).*er.*cosd(lat_i); % [nx1]
d_y = (lat(2:end)-lat(1:end-1))'.*(pi/180).*er; % [nx1]
dist_i = sqrt(d_x.^2 + d_y.^2); % [nx1]
ang_i= atan(d_y./d_x); % [nx1]


% depth
dz_i = depth(2:end) - depth(1:end-1); % [1xz]
depth_i = depth(1:end-1) + 0.5*dz_i; % [1xz]


% area
area_i = repmat(dist_i,1,sz_z-1).*repmat(dz_i,sz_n-1,1);


% land mask on T-Grid,
% modified from OA temp field (updated on 7/24);
mask_land = t_analysis(:,:,1);
mask_land(~isnan(mask_land)) = 1;
mask_land(isnan(mask_land)) = 0;


% land masks for V/I-Grids
mask_land_v = zeros(sz_n-1,sz_z);
mask_land_i = zeros(sz_n-1,sz_z-1);
        

% bathymetry on V/I-grids
load([dirio,'data/Topography/0626_ETOPO2_OSNAP_All_I-Grid.mat']);
lon_i_topo = topo_I.lon;

% FLAG_SECTION
switch flag_section
    case 'west'
        ipx = lon_i_topo <= -45;
    case 'east'
        ipx = lon_i_topo >= -44;
    case 'all'
        ipx = 1:length(lon_i_topo); % keep all grid points
    otherwise
        error('SEC 198!');
% END FLAG_SECTION            
end
lon_i_topo = lon_i_topo(ipx);
bdep_sect_i = -topo_I.z(ipx);
bdep_sect_i(bdep_sect_i<0) = 0;


% land masks updated with the ETOPO2 bathymetry
for i_n = 1:sz_n-1
    ipz = depth<= bdep_sect_i(i_n);
    mask_land_v(i_n,ipz) = 1;
    mask_land_v(i_n,~ipz) = 0;
    
    ipz = depth_i<= bdep_sect_i(i_n);
    mask_land_i(i_n,ipz) = 1;
    mask_land_i(i_n,~ipz) = 0;
end


% area updated accordingly
area_i(mask_land_i==0) = 0;


% coriolis paramter at each grid point
f_i = gsw_f(lat_i); 



%----------------------+
%% COMPUTE FLUXES      |
%----------------------+
disp('+');
disp('+');
disp('+------------------+');
disp('| Compute fluxes   |');
disp('+------------------+');



% Reference salinity for computing 
% the freshwater fluxes;
% Area-weighted salinity;
% e.g.,
% s0 = nansum(nansum((area.*nanmean(salt,3))))/nansum(nansum(area));

% Last updated on 3.14.2018 based on 20180314T1 results:
%   s0 = 34.9189 for OSNAP All
%   s0 = 34.8507 for OSNAP West
%   s0 = 34.9682 for OSNAP East

% Last updated on 5.30.2018 based on 20180524T1,2,3 results:
%   s0 = 34.9187 for OSNAP All
%   s0 = 34.8501 for OSNAP West
%   s0 = 34.9681 for OSNAP East

% Last updated on 4/1/2019 based on 20190331T1 results:
%   s0_All  = 34.93;
%   s0_West = 34.85;
%   s0_East = 34.99;

% Last updated on 7/25/2019 based on 20190525T1 results:
%   s0 = 34.94for OSNAP All
%   s0 = 34.86 for OSNAP West
%   s0 = 34.99 for OSNAP East

% Last updated on 11/11/2019 based on 20191111T1 results:
%   s0 = 34.94 for OSNAP All
%   s0 = 34.86 for OSNAP West
%   s0 = 35.01 for OSNAP East


% Last updated on 3/6/2020 based on 20200305T1i results (updated densOA):
%   s0 = 34.9342 for OSNAP All
%   s0 = 34.8499 for OSNAP West
%   s0 = 34.9950 for OSNAP East


% OSNAP section mean salinity
s0 = 34.93;


% boundary-mean salinity for the wider Arctic bounded by the Bering Strait
% and OSNAP
Area_osnap = 0.7e10;
S_bs = 32.47;
Area_bs = 4.2e6;
Sbar = (s0*Area_osnap + S_bs*Area_bs)/(Area_osnap+Area_bs);


% update s0, from OSNAP section-mean salinity to OSNAP-BS boundary mean
s0 = Sbar;


% LOG
disp(['+ OSNAP-BS boundary-mean salinity s0= ',num2str(s0)]);
disp('+');


% Set the time range to calculate the fluxes
% i_t_begin = find(time<= datenum(2014,8,23,0,0,0),1,'last');  % LS mooring deployed
% i_t_end = find(time<= datenum(2020,6,30,0,0,0),1,'last');  % GEOMAR mooring recovered

% /////////////////// TEST BEGIN //////////////////////////////////////////
i_t_begin= find(time<= datenum(2014,8,23,0,0,0),1,'last'); % 2015-03-01
i_t_end = find(time<= datenum(2016,6,30,0,0,0),1,'last'); % End of the first 21month
% i_t_end = find(time<= datenum(2014,8,31,0,0,0),1,'last'); % 1 mon
% /////////////////// TEST END ////////////////////////////////////////////



% =========================================================================
% FLAG_TEST_RTADCP || FLAG_MEAN_RTADCP ====================================
% =========================================================================
% Test of not using the ADCP data ...
if flag_mean_RTADCP== 0 % || flag_test_RTADCP~= 0
    disp('+ TESTING RTADCP, running time period changed to 10/30/2014-6/19/2015...');
    i_t_begin = find(time<= datenum(2014,10,30,1,0,0),1,'last');  % RTADCP data begin
    i_t_end = find(time<= datenum(2015,6,19,23,0,0),1,'last');  % RTADCP data end
end
% =========================================================================
% 
% 
% 
% % =========================================================================
% % FLAG_TEST_D5 FLAG_TEST_D5_MEAN FLAG_TEST_M2M3_GEO =======================
% % =========================================================================
% if flag_test_d5_nouse== 1 || flag_test_d5_mean== 1 || flag_test_M2M3_geo== 1
%     disp('+ TESTING D5, running time period changed to 7/10/2016-7/13/2018...');
%     i_t_begin = find(time<= datenum(2016,7,10,0,0,0),1,'last');  
%     i_t_end = find(time<= datenum(2018,7,13,0,0,0),1,'last');
% end
% % =========================================================================
% 
% 

% =========================================================================
% FLAG_TEST_RREX ==========================================================
% =========================================================================
if flag_test_rrex_nouse== 1
    disp('+ TESTING RREX, running time period changed to 6/29/2015-7/26/2017...');
    i_t_begin = find(time<= datenum(2015,6,29,0,0,0),1,'last');  
    i_t_end = find(time<= datenum(2017,7,26,0,0,0),1,'last');
end
% =========================================================================



% LOG
disp(['+ Running from time#',num2str(i_t_begin),' (',datestr(time(i_t_begin),'yyyy-mm-dd HH:MM:SS'),') to time#',num2str(i_t_end),' (',datestr(time(i_t_end),'yyyy-mm-dd HH:MM:SS'),')']);
disp('+');


% Determine the number/range of input files
sz_t_max = i_t_end-i_t_begin+1; 


% LOG
disp(['+ Number of output sz_t_max= ',num2str(sz_t_max)]);
disp('+');


% INITILIZATION
per         = nan(sz_t_max,1); % time
velo_per    = nan(sz_n-1,sz_z-1,sz_t_max);
temp_per    = nan(sz_n-1,sz_z-1,sz_t_max);
ptmp_per    = nan(sz_n-1,sz_z-1,sz_t_max);
salt_per    = nan(sz_n-1,sz_z-1,sz_t_max);
dens_per    = nan(sz_n-1,sz_z-1,sz_t_max);
pden_per    = nan(sz_n-1,sz_z-1,sz_t_max);
vflux_per   = nan(sz_n-1,sz_z-1,sz_t_max);
hflux_per   = nan(sz_n-1,sz_z-1,sz_t_max);
sflux_per   = nan(sz_n-1,sz_z-1,sz_t_max);
vflux_d_per = nan(sz_n-1,sz_d,sz_t_max);
hflux_d_per = nan(sz_n-1,sz_d,sz_t_max);
sflux_d_per = nan(sz_n-1,sz_d,sz_t_max);
% vflux_ek_per = nan(sz_t_max,1);

if flag_baro_corr==1
    v_ssh_per = nan(sz_n-1,sz_t_max);
    v_baro_per = nan(sz_n-1,sz_t_max);
end
    

% Output flux time series
% from each MC run
MOC_mc = nan(sz_t_max,2);
MOC1_mc = nan(sz_t_max,2);
MOC2_mc = nan(sz_t_max,2);
MOC3_mc = nan(sz_t_max,2);

MHT_mc = nan(sz_t_max,2);
MFT_mc = nan(sz_t_max,2);
MST_mc = nan(sz_t_max,2);

Tek_mc = nan(sz_t_max,2);
Text_mc = nan(sz_t_max,2);
Text_mc_w = nan(sz_t_max,2);
Text_mc_e = nan(sz_t_max,2);

Tprof_mc = nan(sz_t_max,sz_z-1,2);
Tprof_d_mc = nan(sz_t_max,sz_d,2);


if strcmp(flag_section,'all') && flag_test_AllWE== 1
    % *** West ***
    MOC_mc_w = nan(sz_t_max,2);
    MOC1_mc_w = nan(sz_t_max,2);
    MOC2_mc_w = nan(sz_t_max,2);
    MOC3_mc_w = nan(sz_t_max,2);

    MHT_mc_w = nan(sz_t_max,2);
    MFT_mc_w = nan(sz_t_max,2);
    MST_mc_w = nan(sz_t_max,2);

    Tek_mc_w = nan(sz_t_max,2);
    
    Tprof_mc_w = nan(sz_t_max,sz_z-1,2);
    Tprof_d_mc_w = nan(sz_t_max,sz_d,2);


    % *** East ***
    MOC_mc_e = nan(sz_t_max,2);
    MOC1_mc_e = nan(sz_t_max,2);
    MOC2_mc_e = nan(sz_t_max,2);
    MOC3_mc_e = nan(sz_t_max,2);

    MHT_mc_e = nan(sz_t_max,2);
    MFT_mc_e = nan(sz_t_max,2);
    MST_mc_e = nan(sz_t_max,2);

    Tek_mc_e = nan(sz_t_max,2);
    
    Tprof_mc_e = nan(sz_t_max,sz_z-1,2);
    Tprof_d_mc_e = nan(sz_t_max,sz_d,2);

end


% =========================================================================
% flag_HT_FWT_decomposition ===============================================
% =======================================================================
% MHT and MFT decompsition in both density and depth space
if flag_HT_FWT_decomposition == 1

    % output from each MC run
    % *** Heat ***   
    MHT_mc_dia = nan(sz_t_max,2);
    MHT_mc_iso = nan(sz_t_max,2);
    MHT_mc_over = nan(sz_t_max,2);
    MHT_mc_gyre = nan(sz_t_max,2);
    MHT_mc_net = nan(sz_t_max,2);
    MHT_mc_netz = nan(sz_t_max,2);
    
    % *** Freshwater ***
    MFT_mc_dia = nan(sz_t_max,2);
    MFT_mc_iso = nan(sz_t_max,2);
    MFT_mc_over = nan(sz_t_max,2);
    MFT_mc_gyre = nan(sz_t_max,2);
    MFT_mc_net = nan(sz_t_max,2);
    MFT_mc_netz = nan(sz_t_max,2);
    

    if strcmp(flag_section,'all') && flag_test_AllWE== 1
        % *** West ***
            % *** Heat ***   
        MHT_mc_dia_w = nan(sz_t_max,2);
        MHT_mc_iso_w = nan(sz_t_max,2);
        MHT_mc_over_w = nan(sz_t_max,2);
        MHT_mc_gyre_w = nan(sz_t_max,2);
        MHT_mc_net_w = nan(sz_t_max,2);
        MHT_mc_netz_w = nan(sz_t_max,2);

        % *** Freshwater ***
        MFT_mc_dia_w = nan(sz_t_max,2);
        MFT_mc_iso_w = nan(sz_t_max,2);
        MFT_mc_over_w = nan(sz_t_max,2);
        MFT_mc_gyre_w = nan(sz_t_max,2);
        MFT_mc_net_w = nan(sz_t_max,2);
        MFT_mc_netz_w = nan(sz_t_max,2);
    
        % *** East ***
            % *** Heat ***   
        MHT_mc_dia_e = nan(sz_t_max,2);
        MHT_mc_iso_e = nan(sz_t_max,2);
        MHT_mc_over_e = nan(sz_t_max,2);
        MHT_mc_gyre_e = nan(sz_t_max,2);
        MHT_mc_net_e = nan(sz_t_max,2);
        MHT_mc_netz_e = nan(sz_t_max,2);

        % *** Freshwater ***
        MFT_mc_dia_e = nan(sz_t_max,2);
        MFT_mc_iso_e = nan(sz_t_max,2);
        MFT_mc_over_e = nan(sz_t_max,2);
        MFT_mc_gyre_e = nan(sz_t_max,2);
        MFT_mc_net_e = nan(sz_t_max,2);
        MFT_mc_netz_e = nan(sz_t_max,2);
    end
    
end
% =========================================================================
% =========================================================================





% Counter
i_per = 1; % the number of output




%% LOOP TIMESTEP
for i_t = i_t_begin:i_t_end

    
    % LOG
    disp('+');
    display(['+ Process i_t= ',num2str(i_t),' time= ',datestr(time(i_t),'yyyy-mm-dd HH:MM:SS')]);
    disp('+');
    
    
    % variable for MC-
    % running mean MOC, MHT and MFT;
    delta_MOC = 100;
    delta_MOC1 = 100;
    delta_MOC2 = 100;
    delta_MOC3 = 100;
    delta_MHT = 100;
    delta_MFT = 100;
    
    
    % variable for MC-
    % output from one MC run;
    MOC_n = nan(1,5000);
    MOC1_n = nan(1,5000);
    MOC2_n = nan(1,5000);
    MOC3_n = nan(1,5000);
    MHT_n = nan(1,5000);
    MFT_n = nan(1,5000);
    MST_n = nan(1,5000);
    Tek_n = nan(1,5000);
    Text_n = nan(1,5000);
    Text_n_w = nan(1,5000);
    Text_n_e = nan(1,5000);
    
    Tprof_d_n = nan(sz_d,5000);
    Tprof_n = nan(sz_z-1,5000);
    
    
    if strcmp(flag_section,'all') && flag_test_AllWE== 1
        % variable for MC-
        % output from one MC run;
        MOC_n_w = nan(1,5000);
        MOC1_n_w = nan(1,5000);
        MOC2_n_w = nan(1,5000);
        MOC3_n_w = nan(1,5000);
        MHT_n_w = nan(1,5000);
        MFT_n_w = nan(1,5000);
        MST_n_w = nan(1,5000);
        Tek_n_w = nan(1,5000);
        Tprof_d_n_w = nan(sz_d,5000);
        Tprof_n_w = nan(sz_z-1,5000);
        
        MOC_n_e = nan(1,5000);
        MOC1_n_e = nan(1,5000);
        MOC2_n_e = nan(1,5000);
        MOC3_n_e = nan(1,5000);
        MHT_n_e = nan(1,5000);
        MFT_n_e = nan(1,5000);
        MST_n_e = nan(1,5000);
        Tek_n_e = nan(1,5000);
        Tprof_d_n_e = nan(sz_d,5000);
        Tprof_n_e = nan(sz_z-1,5000);
    end
    
% =========================================================================
% flag_HT_FWT_decomposition ===============================================
% =======================================================================
% MHT and MFT decompsition in both density and depth space
    if flag_HT_FWT_decomposition == 1

        % initialization
        
        % *** Heat ***
        MHT_n_dia = nan(1,5000);
        MHT_n_iso = nan(1,5000);
        MHT_n_over = nan(1,5000);
        MHT_n_gyre = nan(1,5000);
        MHT_n_net = nan(1,5000);

  
        % *** Freshwater ***
        MFT_n_dia = nan(1,5000);
        MFT_n_iso = nan(1,5000);
        MFT_n_net = nan(1,5000);
        MFT_n_over = nan(1,5000);
        MFT_n_gyre = nan(1,5000);
        MFT_n_netz = nan(1,5000);

%         MFT_n_wBS = nan(1,5000);
%         MFT_n_div = nan(1,5000);
%         MFT_n_net_div = nan(1,5000);

        if strcmp(flag_section,'all') && flag_test_AllWE== 1
            % *** West ***
                % *** Heat ***
            MHT_n_dia_w = nan(1,5000);
            MHT_n_iso_w = nan(1,5000);
            MHT_n_net_w = nan(1,5000);
            MHT_n_over_w = nan(1,5000);
            MHT_n_gyre_w = nan(1,5000);
            MHT_n_netz_w = nan(1,5000);


            % *** Freshwater ***
            MFT_n_dia_w = nan(1,5000);
            MFT_n_iso_w = nan(1,5000);
            MFT_n_net_w = nan(1,5000);
            MFT_n_over_w = nan(1,5000);
            MFT_n_gyre_w = nan(1,5000);
            MFT_n_netz_w = nan(1,5000);

            
            % *** East ***
                % *** Heat ***
            MHT_n_dia_e = nan(1,5000);
            MHT_n_iso_e = nan(1,5000);
            MHT_n_net_e = nan(1,5000);
            MHT_n_over_e = nan(1,5000);
            MHT_n_gyre_e = nan(1,5000);
            MHT_n_netz_e = nan(1,5000);


            % *** Freshwater ***
            MFT_n_dia_e = nan(1,5000);
            MFT_n_iso_e = nan(1,5000);
            MFT_n_net_e = nan(1,5000);
            MFT_n_over_e = nan(1,5000);
            MFT_n_gyre_e = nan(1,5000);
            MFT_n_netz_e = nan(1,5000);
            
        end
        
    end
% =========================================================================
% =========================================================================
    
    
    % variable for MC-
    % number of Monte Carlo iterations;
    cnt_mc = 1;
    
    
    % variable for MC-
    % stopping criteria based on:
    %   * 1/1000 of O(~stdev) for MOC;
    %   * 1 MOC change is about 1/50 MHT 
    %       and MFT change;
    epsilon_MOC = 1e-3;
    epsilon_MHT = epsilon_MOC/50;
    epsilon_MFT = epsilon_MOC/50;
            
    
    
% =========================================================================
% FLAG_TEST_MC_WEST_EPSILON ===============================================
% =======================================================================
% Test of using epsilons adjusted for OSNAP W
% The results show negligible diff. but with
% increased hours of computing;
%
% Epsilons are based on:
%   * 1/1000 of O(~stdev) for MOC_West;
%   * 1 MOC change is about 1/50 MHT 
%       and MFT change;

    if flag_test_MC_West_epsilon== 1 && strcmp(flag_section,'west')
        epsilon_MOC = 1e-4;
        epsilon_MHT = epsilon_MOC/100;
        epsilon_MFT = epsilon_MOC/100;
    end
% =========================================================================
% =========================================================================  
    
    

% =========================================================================
% FLAG_MC_MOC =============================================================
% =======================================================================  
    if flag_MC_MOC >= 10
        % run a number of iterations specified by flag_MC_MOC
        epsilon_MOC = 1e-10;
        epsilon_MHT = epsilon_MOC;
        epsilon_MFT = epsilon_MOC;
    end
% =========================================================================
% =========================================================================  

    
% LOOP MONTE CARLO SIMULATION    
    while (abs(delta_MOC) >= epsilon_MOC ||...
           abs(delta_MHT) >= epsilon_MHT ||...
           abs(delta_MFT) >= epsilon_MFT) && ...
           cnt_mc <= 5000    % maximum 5000 iterations

       
% =========================================================================
% FLAG_MC & FLAG_MC_MOC ===================================================
% =======================================================================       
        if flag_MC == 1

            switch flag_MC_MOC
                case 0 
                    % do nothing;

                    if i_t == i_t_begin && cnt_mc == 1
                        disp(['>>> MC iteration stops using epsilon_MOC= ',num2str(epsilon_MOC), ' and epsilon_MHT= ',num2str(epsilon_MHT),' and epsilon_MFT= ',num2str(epsilon_MFT)]);
                    end


                case 1

                    if i_t == i_t_begin && cnt_mc == 1
                        disp(['>>> MC iteration stops using epsilon_MOC= ',num2str(epsilon_MOC), ' and epsilon_MHT= ',num2str(epsilon_MHT)]);
                    end


                    % Consider MOC and MHT ...
                    if abs(delta_MOC) < epsilon_MOC && ...
                            abs(delta_MHT) < epsilon_MHT
                         break;
                    end

                case 2

                    if i_t == i_t_begin && cnt_mc == 1
                        disp(['>>> MC iteration stops using epsilon_MOC= ',num2str(epsilon_MOC)]);
                    end

                    % Consider MOC only ...
                    if abs(delta_MOC) < epsilon_MOC
                        disp('delta_MOC < epsilon_MOC, will quick MC...');
                         break;
                    end

                otherwise

                    if i_t == i_t_begin && cnt_mc == 1
                        disp(['>>> MC iteration stops after #',num2str(flag_MC_MOC), ' iterations <<<']);
                    end

                    % Consider a fixed number of iterations...
                    if cnt_mc > flag_MC_MOC
                        disp(['cnt_mc_mc> ',num2str(flag_MC_MOC),', will quick MC...']);
                         break;
                    end

            end


        end
% =========================================================================
% =========================================================================
 

        % ----------------------------------------------------------
        %% (1) Prepare background TSUV fields
        % ----------------------------------------------------------
 
        % ------------------------------
        % (1.1) CLIMATOLOGICAL T/S/U/V
        % ------------------------------
        % Get climatology t/s/u/v for the current month
        im = str2double(datestr(time(i_t)+ flag_data_interval_4calculation/2,'mm')); % which month
        
        % *** U/V ***
        u_t_c = u_clim_model(:,:,im); % from FLAME
        v_t_c = v_clim_model(:,:,im); % from FLAME

        % Set velocity on the dry nodes to 0
        u_t_c(isnan(u_t_c)) = 0;
        v_t_c(isnan(v_t_c)) = 0;

        
        % *** T/S ***
        if flag_test_LC_TS_clim==1
            t_t_c = t_clim(:,:,im); % from updated WOA18 
            s_t_c = s_clim(:,:,im); % from updated WOA18
            t_t_c(1:21,:)=nanmean(t_clim(1:21,:,:),3);
            s_t_c(1:21,:)=nanmean(s_clim(1:21,:,:),3);
        else
            t_t_c = t_clim(:,:,im); % from updated WOA18
            s_t_c = s_clim(:,:,im); % from updated WOA18
        end
        
% %         t_t_c = t_clim_model(:,:,im); % from FLAME
% %         s_t_c = s_clim_model(:,:,im); % from FLAME
        

        % Derive dens, pden and ptmp
        % from climatological T/S 
        pres = gsw_p_from_z(-repmat(depth,sz_n,1),lat'); % pres
        [SA,~,~] = gsw_SA_Sstar_from_SP(s_t_c,pres,lon',lat'); % absolute salinity 
        pt_t_c = gsw_pt_from_t(SA,t_t_c,pres,0); % potential temperature  
        d_t_c = gsw_rho_t_exact(SA,t_t_c,pres); % density 
        pd_t_c = gsw_pot_rho_t_exact(SA,t_t_c,pres,0); % potential density

        

        % ---------------------------------------------------
        %% (2) Fill up t-grid and prepare masks
        % ---------------------------------------------------

        % -----------------------
        %% (2.1) OBSERVATIONS
        % -----------------------
        
        % ********************************
        % T/S/U/V regridded profiles
        % ********************************
        % Get OSNAP U/V/T/S data from the distributions
        % averaged over the calculation time interval;  

        
        % Two ways of providing random profiles from each mooring ...
        % NOTE 
        %   These two methods provide slighly different temporal 
        %   mean profiles;
        %



% =========================================================================
% flag_regridded_profiles (2/2) ===========================================
% ====================================================================  
        if flag_regridded_profiles== 1
% IF FLAG_REGRIDDED_PORFILES
            % *** Method1 ***
            % Work with daily regridded profiles 
            %   at each mooring, which are prepared 
            %   beforehand using pp_osnap_mooring_vinterp_2019.m;
            
            if ~exist('osnap','var')
                % load regridded daily profiles
                data_mr_osnap= [dirio,'outputs/regridded/20210616_mooring_data_OSNAP_All_20140701-20200630_REGRIDDED_all_4noMC.mat']; 
%                 data_mr_osnap= [dirio,'data/intermediate_output/0101_mooring_data_OSNAP_All_20140701-20180630_REGRIDDED.mat'];
%                 data_mr_osnap= [dirio,'/mooring/regridded/1117_mooring_data_OSNAP_All_daily_regridded_combined.mat'];  %0313* for 2014-2016


                load(data_mr_osnap,'osnap');
                
                % LOG
                disp(['+ data_OSNAP_mooring_daily_regridded_profiles= ',data_mr_osnap]);
            end
            
            
            
            
            
% % =========================================================================
% % FLAG_TEST_C3 ============================================================
% % =========================================================================
% % for 'West' or 'All' only...
% if strcmp(flag_section,'west') || strcmp(flag_section,'all') 
%     
%     % setting flag
%     if flag_test_C3== 0
%         % do nothing;
% 
%     elseif flag_test_C3== 1
%         % LOG
%         if i_t == i_t_begin
%             disp('*******C3 testing underway...');
%         end
%     
%         % found time after the 2016 deployment
%         ipt = osnap.time> datenum(2016,5,18,12,0,0);
% 
%         % time-mean
%         osnap.t(5,:,:) = repmat(nanmean(osnap.t(5,:,ipt),3),[1 1 length(osnap.time)]);
%         osnap.s(5,:,:) = repmat(nanmean(osnap.s(5,:,ipt),3),[1 1 length(osnap.time)]);
%         osnap.u(5,:,:) = repmat(nanmean(osnap.u(5,:,ipt),3),[1 1 length(osnap.time)]);
%         osnap.v(5,:,:) = repmat(nanmean(osnap.v(5,:,ipt),3),[1 1 length(osnap.time)]);
% 
%     elseif flag_test_C3== 2 
%         % LOG
%         if i_t == i_t_begin
%             disp('*********C3 testing underway...');
%         end
%     
%         % copy from K7;
%         %   TODO: consider some sort of 
%         %   coefficients determiend 
%         %   by the measurements at C3 and K7?
%         osnap.t(5,:,:) = osnap.t(6,:,:);
%         osnap.s(5,:,:) = osnap.s(6,:,:);
%         osnap.u(5,:,:) = osnap.u(6,:,:);
%         osnap.v(5,:,:) = osnap.v(6,:,:);
%         
%         
% 	elseif flag_test_C3== 3 
%         % LOG
%         if i_t == i_t_begin
%             disp('*********C3 testing underway...');
%         end
%         
%         % delete C3 data completely...
%         osnap.t(5,:,:) = nan;
%         osnap.s(5,:,:) = nan;
%         osnap.u(5,:,:) = nan;
%         osnap.v(5,:,:) = nan;
%         
%         
%         
%     elseif flag_test_C3== 4 
%         % LOG
%         if i_t == i_t_begin
%             disp('*********C3 testing underway...');
%         end
%         
%         
%         %% *** TS ***
% %         % method1, copying from K7
% %         osnap.t(5,:,:) = osnap.t(6,:,:); % K7
% %         osnap.s(5,:,:) = osnap.s(6,:,:); % K7
% 
% %         % method2, time-mean C3
% %         % found time after the 2016 deployment
% %         ipt0 = osnap.time> datenum(2016,5,18,12,0,0);
% % 
% %         % time-mean
% %         osnap.t(5,:,:) = repmat(nanmean(osnap.t(5,:,ipt0),3),[1 1 length(osnap.time)]);
% %         osnap.s(5,:,:) = repmat(nanmean(osnap.s(5,:,ipt0),3),[1 1 length(osnap.time)]);
% 
%         % method3, recons. from K7
%         t1 = squeeze(osnap.t(5,:,:)); 
%         s1 = squeeze(osnap.s(5,:,:));
%         t2 = squeeze(osnap.t(6,:,:)); % K7
%         s2 = squeeze(osnap.s(6,:,:)); % K7
%         
%         
%         % found time for the 2014-2016 deployment
%         ipt = osnap.time< datenum(2016,5,18,12,0,0);
%         
%         
%         % method 1, based on the mean profiles to reconstruct C3 TS
%         % mean profiles
%         tp1 = nanmean(t1(:,ipt),2);
%         tp2 = nanmean(t2(:,ipt),2);
%         
% %         t1_adj1 = nan(size(t1));
% %         s1_adj1 = nan(size(t1));
% 
%         ipx = ~isnan(tp1);
%         
%         stats_t1 = regstats(tp1(ipx),tp2(ipx),'linear',{'beta','rsquare'});
%         t1(ipx,:) = t2(ipx,:).*stats_t1.beta(2) + stats_t1.beta(1);
%         
%         sp1 = nanmean(s1(:,ipt),2);
%         sp2 = nanmean(s2(:,ipt),2);
%         
%         ipx = ~isnan(sp1);
%         
%         stats_s1 = regstats(sp1(ipx),sp2(ipx),'linear',{'beta','rsquare'});
%         s1(ipx,:) = s2(ipx,:).*stats_s1.beta(2) + stats_s1.beta(1);
%         
%         %%
%         clear stats_t1 stats_s1 tp1 tp2 sp1 sp2
%         
%         
%         
%         % *** UV ***
%         u1 = squeeze(osnap.u(5,:,:)); 
%         v1 = squeeze(osnap.v(5,:,:));
%         u2 = squeeze(osnap.u(6,:,:)); % K7
%         v2 = squeeze(osnap.v(6,:,:)); % K7
%         
%         
%         % found time for the 2014-2016 deployment
%         ipt = osnap.time< datenum(2016,5,18,12,0,0);
%         
%         
%         % method 1, based on the mean profiles to reconstruct C3 UV
%         % mean profiles
%         vp1 = nanmean(v1(:,ipt),2);
%         vp2 = nanmean(v2(:,ipt),2);
%         
%         % variables for reconstruction
%         v1_adj1 = nan(size(v1));
%         u1_adj1 = nan(size(v1));
% 
%         % valid water depth
%         ipx = ~isnan(vp1);
%         
%         % regression
%         stats_v1 = regstats(vp1(ipx),vp2(ipx),'linear',{'beta','rsquare'});
%         v1_adj1(ipx,:) = v2(ipx,:).*stats_v1.beta(2) + stats_v1.beta(1);
%         
%         % repeat above for ucur
%         up1 = nanmean(u1(:,ipt),2);
%         up2 = nanmean(u2(:,ipt),2);
%         
%         stats_u1 = regstats(up1(ipx),up2(ipx),'linear',{'beta','rsquare'});
%         u1_adj1(ipx,:) = u2(ipx,:).*stats_u1.beta(2) + stats_u1.beta(1);
%         
%         
%         clear stats_u1 stats_v1 up1 up2 vp1 vp2
%         
%         
% %         % method 2, regressing time-series at each level, tested -
% %         yielding a good mean profile, but rsquare is way too small;
% %         v1_adj = nan(size(v1));
% %         u1_adj = nan(size(v1));
% %         
% %         coef_v = nan(42,3);
% %         coef_u = nan(42,3);
% %         
% %         for xz = 1:42
% %             vs1 = v1(xz,ipt);
% %             vs2 = v2(xz,ipt);
% %             ipx = ~isnan(vs1) & ~isnan(vs2);
% %             stats_v = regstats(vs1(ipx),vs2(ipx),'linear',{'beta','rsquare'});
% %             coef_v(xz,1:2) = stats_v.beta;
% %             coef_v(xz,3) = stats_v.rsquare;
% %             v1_adj(xz,:) = v2(xz,:).*stats_v.beta(2) + stats_v.beta(1);
% %             
% %             us1 = u1(xz,ipt);
% %             us2 = u2(xz,ipt);
% %             ipx = ~isnan(us1) & ~isnan(us2);
% %             stats_u = regstats(us1(ipx),us2(ipx),'linear',{'beta','rsquare'});
% %             coef_u(xz,1:2) = stats_u.beta;
% %             coef_u(xz,3) = stats_u.rsquare;
% %             u1_adj(xz,:) = u2(xz,:).*stats_u.beta(2) + stats_u.beta(1);
% %         end
%         
% 
%         % rewrite
%         osnap.u(5,:,:) = u1_adj1; 
%         osnap.v(5,:,:) = v1_adj1;
%         
%         clear u1 v1 u2 v2 u1_adj1 v1_adj1
%         
%     else
%         error('LUI 018!');
% 
%     end
%     
% else
%     error('CCC 333!');
%     
% end
% % =========================================================================
% % =========================================================================

            
            
            % find indices for profiles within
            % the current time interval
            time_mr = osnap.time;
            sz_mr = length(osnap.lon);
            if flag_data_interval=='m'
                ipt = find(time_mr>= time(i_t) & time_mr< time(i_t+1));
            else
                ipt = find(time_mr>= time(i_t) & time_mr< time(i_t)+ flag_data_interval);
            end
            
            
%             % LOG
%             for i0 = 1:length(ipt)
%                 disp([' Found time_mr= ',datestr(time_mr(ipt(i0)),'yyyy-mm-dd HH:MM:SS')]);
%             end


            % initialize variables used for flux estimates
            t_mr = nan(sz_mr,sz_z);
            s_mr = nan(sz_mr,sz_z);
            d_mr = nan(sz_mr,sz_z);
            u_mr = nan(sz_mr,sz_z);
            v_mr = nan(sz_mr,sz_z);


% LOOP MOORING     
            for i_mr = 1:sz_mr

                % select profiles within this time interval
                % and convert to [depth,time]
                t1 = shiftdim(osnap.t(i_mr,:,ipt));
                s1 = shiftdim(osnap.s(i_mr,:,ipt));
                d1 = shiftdim(osnap.d(i_mr,:,ipt));
                u1 = shiftdim(osnap.u(i_mr,:,ipt));
                v1 = shiftdim(osnap.v(i_mr,:,ipt));

% =========================================================================
% FLAG_MC =================================================================
% ====================================================================            
                if flag_MC== 1 && flag_data_interval> 1
                    % Randomly select T/S/U/V profiles;
                    % Calculate density from T/S;
% IF FLAG_MC                       

                    % *** Random T profile ***
                    
                    % check for valid profiles
                    ip_m = nan(1,length(ipt));
                    for ii = 1:length(ipt)
                        ip_m(ii) = any(~isnan(t1(:,ii)));
                    end

%                     % LOG
%                     disp('xxx');
%                     disp(num2str(sum(ip_m)));

                    % select from valid profiles
                    if all(ip_m== 0)
                        % no valid profiles
                        % do nothing
                        
                        % LOG
%                         disp(['mr= ',num2str(i_mr),', no available T data!']);
                        
                    else 
                        % empty profiles with no data (all nans)
                        t1(:,ip_m== 0) = [];
                        ip_m(ip_m== 0) = [];
                        
                        % randomly pick out one T profile
                        rndt = randi([1 length(ip_m)]);
                        t_mr(i_mr,:) = t1(:,rndt);
                        
                        
%                         % LOG
%                         disp(['   i_mr= ',num2str(i_mr),'   T-total= ',num2str(length(ip_m)),' T-selected= ',num2str(rndt)]);
                    
                    end

                    
                    % *** Random S profile ***
                    
                    % check for valid profiles
                    ip_m = nan(1,length(ipt));
                    for ii = 1:length(ipt)
                        ip_m(ii) = any(~isnan(s1(:,ii)));
                    end


%                     % LOG
%                     disp(num2str(sum(ip_m)));
                    
                    
                    % select from valid profiles
                    if all(ip_m== 0)
                        % no valid profiles
                        % do nothing
                        
                        % LOG
%                         disp(['mr= ',num2str(i_mr),', no available T data!']);
                        
                    else 
                        % empty profiles with no data (all nans)
                        s1(:,ip_m== 0) = [];
                        ip_m(ip_m== 0) = [];
                        
                        % randomly draw a profile from each mooring;
                        rndt = randi([1 length(ip_m)]); 
                        s_mr(i_mr,:) = s1(:,rndt);
                        
                        
%                         % LOG
%                         disp(['   i_mr= ',num2str(i_mr),'   S-total= ',num2str(length(ip_m)),' S-selected= ',num2str(rndt)]);
                    
                    end
                    
                    
% %                     % -------- FOR TESTING PURPOSE ONLY ---------------
% %                     % *** Random Dens profile ***
% %                     
% %                     % check for valid profiles
% %                     ip_m = nan(1,length(ipt));
% %                     for ii = 1:length(ipt)
% %                         ip_m(ii) = any(~isnan(d1(:,ii)));
% %                     end
% % 
% % 
% % %                     % LOG
% % %                     disp(num2str(sum(ip_m)));
% %                     
% %                     
% %                     % select from valid profiles
% %                     if all(ip_m== 0)
% %                         % no valid profiles
% %                         % do nothing
% %                         
% %                         % LOG
% % %                         disp(['mr= ',num2str(i_mr),', no available density data!']);
% %                         
% %                     else 
% %                         % empty profiles with no data (all nans)
% %                         d1(:,ip_m== 0) = [];
% %                         ip_m(ip_m== 0) = [];
% %                         
% %                         % randomly draw a profile from each mooring;
% %                         rndt = randi([1 length(ip_m)]); 
% %                         d_mr(i_mr,:) = d1(:,rndt);
% %                         
% %                         
% % %                         % LOG
% % %                         disp(['   i_mr= ',num2str(i_mr),'   S-total= ',num2str(length(ip_m)),' S-selected= ',num2str(rndt)]);
% %                     
% %                     end
                    
                    
                    
                    % *** Derived density profile ***
                    pres = gsw_p_from_z(-depth,osnap.lon(i_mr))';
                    [SA,~,~] = gsw_SA_Sstar_from_SP(s_mr(i_mr,:),pres,osnap.lon(i_mr),osnap.lat(i_mr)); % absolute salinity 
                    d_mr(i_mr,:) = gsw_rho_t_exact(SA,t_mr(i_mr,:),pres); % in-situ density 

                    
                    
                    % *** Random U and V profiles ***
                    
                    % check for valid profiles
                    ip_m = nan(1,length(ipt));
                    for ii = 1:length(ipt)
                        ip_m(ii) = any(~isnan(u1(:,ii))) && any(~isnan(v1(:,ii)));
                    end


%                     % LOG
%                     disp(num2str(sum(ip_m)));
                    
                    
                    % select from valid profiles
                    if all(ip_m== 0)
                        % no valid profiles
                        % do nothing
                        
                        % LOG
%                         disp(['mr= ',num2str(i_mr),', no available U data!']);
                        
                    else 
                        % empty profiles with no data (all nans)
                        u1(:,ip_m== 0) = [];
                        v1(:,ip_m== 0) = [];
                        ip_m(ip_m== 0) = [];
                        
                        % randomly draw a profile from each mooring
                        % ucur
                        rndt = randi([1 length(ip_m)]); 
                        u_mr(i_mr,:) = u1(:,rndt);
                        
%                         % LOG
%                         disp(['   i_mr= ',num2str(i_mr),'   U-total= ',num2str(length(ip_m)),' U-selected= ',num2str(rndt)]);
                        
                        % vcur
                        rndt = randi([1 length(ip_m)]);  % not from the same day
                        v_mr(i_mr,:) = v1(:,rndt);
                        
                        
%                         % LOG
%                         disp(['   i_mr= ',num2str(i_mr),'   V-total= ',num2str(length(ip_m)),' V-selected= ',num2str(rndt)]);
                        
                    end


                else % for non-MC runs

                    % averages from daily profiles 
                    t_mr(i_mr,:) = nanmean(t1,2);
                    s_mr(i_mr,:) = nanmean(s1,2);
                    d_mr(i_mr,:) = nanmean(d1,2);
                    u_mr(i_mr,:) = nanmean(u1,2);
                    v_mr(i_mr,:) = nanmean(v1,2);
                              
                    
                    
%                     % FOR TESTING PURPOSE ONLY ....
%                     % alternatively, derive density profile
%                     pres = gsw_p_from_z(-depth,osnap.lon(i_mr))';
%                     [SA,~,~] = gsw_SA_Sstar_from_SP(s_mr(i_mr,:),pres,osnap.lon(i_mr),osnap.lat(i_mr)); % absolute salinity 
%                     d_mr(i_mr,:) = gsw_rho_t_exact(SA,t_mr(i_mr,:),pres); % in-situ density 
                    
                    
% END IF FLAG_MC
                end
% ====================================================================
% =========================================================================

                

% END LOOP MOORING                
            end

            
            % save to a structure, with 
            % an increasing longitude...
            osnap_mr.lon = osnap.lon;
            osnap_mr.lat = osnap.lat;
            osnap_mr.t = t_mr;
            osnap_mr.s = s_mr;
            osnap_mr.d = d_mr;
            osnap_mr.u = u_mr;
            osnap_mr.v = v_mr;
            
            % ///check-point/// loaded mooring data

            
% /////////////////////////////////////////////////////////////////////////
% Specific processing for RTADCP (OM58)
% /////////////////////////////////////////////////////////////////////////
% % Use the time-mean data when no data were returned from RTADCP;  
% % May not need anymore;
%             if flag_mean_RTADCP== 1 && all(isnan(v_mr(53,:)))
% 
%                 osnap_mr.u(53,:) = shiftdim(nanmean(osnap.u(53,:,:),3));
%                 osnap_mr.v(53,:) = shiftdim(nanmean(osnap.v(53,:,:),3));
%             end
% /////////////////////////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
           


% FLAG_SECTION
            switch flag_section

                case 'west'
                % 22 moorings
                    osnap_mr.lon = osnap_mr.lon(1:22);
                    osnap_mr.lat = osnap_mr.lat(1:22);
                    osnap_mr.t = osnap_mr.t(1:22,:);
                    osnap_mr.s = osnap_mr.s(1:22,:);
                    osnap_mr.d = osnap_mr.d(1:22,:);
                    osnap_mr.u = osnap_mr.u(1:22,:);
                    osnap_mr.v = osnap_mr.v(1:22,:);

                case 'east'
                % the remaining moorings
                    osnap_mr.lon = osnap_mr.lon(23:end);
                    osnap_mr.lat = osnap_mr.lat(23:end);
                    osnap_mr.t = osnap_mr.t(23:end,:);
                    osnap_mr.s = osnap_mr.s(23:end,:);
                    osnap_mr.d = osnap_mr.d(23:end,:);
                    osnap_mr.u = osnap_mr.u(23:end,:);
                    osnap_mr.v = osnap_mr.v(23:end,:);

                case 'all'
                % do nothing 
                % to keep all moorings
% END FLAG_SECTION
            end
        
        
            
        elseif flag_regridded_profiles == 0

            % *** Method2 ***
            % Work with daily moored data from each instrument, using
            %   fun_osnap_mooring_ave_instru to create 
            %   TemporaryFile_mooring_data_ave.mat;
            % Regrid vertically using
            %   fun_osnap_mooring_vinterp_random_v*.m together with
            % 	TemporaryFile_sst_ave.mat;
        
            % calculate averages (w/ or w/o uncertainty) for once
            if i_t== i_t_begin && cnt_mc== 1

                % directory to daily merged and ungridded moored data
                data_osnap_mr_daily = [dirio,'outputs/merged'];
% data_osnap_mr_daily = ['/Volumes/TOSHIBA/Work/osnap_2019_test/outputs/merged'];
                % LOG
                disp(['+ data_OSNAP_moorings_raw= ',data_osnap_mr_daily]);

                
                % prepare averaged data, saved as
                %   TemporaryFile_mooring_data_ave.mat;
                fun_osnap_mooring_ave_instru_mm(data_osnap_mr_daily,flag_data_interval);

            end
            
        
 %%           
% ========================================================================= 
% FLAG_TEST_D5_MEAN =======================================================
% ========================================================================= 
% Use the time-mean data instead of actual data...
            if flag_test_d5_mean== 1 && ~strcmp(flag_section,'west') && i_t== i_t_begin
                
                % load mooring datat at instrument levels
                dir_osnap_obs = 'TemporaryFile_mooring_data_ave.mat'; % data_mr_ave.om[number]
                load(dir_osnap_obs);


                clear tmp tmp1 tmp2

                % found M2, D4, D5
%                 tmp1 = [-28.0200 -26.9673 -25.6778]; % M2 (OM51), D4 (OM52), D5 (OM62) 
                
  

                
                % *** reconstructing the added D4 (OM52)
                ipt = find(data_mr.om62.time> datenum(2016,7,10,0,0,0)); % added D4 cm/microcat available afterwards
                
                
                % For T and S
                % Fill in the first D4 mcat (1800m) data with a weighted
                %   average of the data from level 1600m and 1800m from M2
                data_mr.om52.t(1,ipt) = 0.5*data_mr.om51.t(9,ipt)+0.5*data_mr.om51.t(10,ipt);
                data_mr.om52.s(1,ipt) = 0.7*data_mr.om51.s(9,ipt)+0.3*data_mr.om51.s(10,ipt);


                % For TS pres
                regcoeff = regstats(data_mr.om52.t_pres(1,ipt), data_mr.om52.t_pres(2,ipt),'linear',{'beta','rsquare'});
                data_mr.om52.t_pres(1,ipt) = data_mr.om52.t_pres(2,ipt).*regcoeff.beta(2) + regcoeff.beta(1);
                data_mr.om52.s_pres(1,ipt) = data_mr.om52.t_pres(1,ipt);

                
                % For U and V
                % Use a linear regression of the currents observed on 1800m and 
                %   2300m during the 2016-2018 records to fill in the currents 
                %   on 1800m for the first two years

                regcoeff = regstats(data_mr.om52.u(1,ipt,1), data_mr.om52.u(2,ipt,1),'linear',{'beta','rsquare'});
                data_mr.om52.u(1,ipt,1) = data_mr.om52.u(2,ipt,1).*regcoeff.beta(2) + regcoeff.beta(1);

                regcoeff = regstats(data_mr.om52.v(1,ipt,1), data_mr.om52.v(2,ipt,1),'linear',{'beta','rsquare'});
                data_mr.om52.v(1,ipt,1) = data_mr.om52.v(2,ipt,1).*regcoeff.beta(2) + regcoeff.beta(1);

                % copy pressure record from mcat
                data_mr.om52.v_pres(1,ipt) = data_mr.om52.t_pres(1,ipt);
                
                
                
                
                % *** using the mean at D5 (OM62)
                tmp = length(data_mr.om62.time);
                
                data_mr.om62.v_pres = repmat(nanmean(data_mr.om62.v_pres,2),1,tmp);
                data_mr.om62.t_pres = repmat(nanmean(data_mr.om62.t_pres,2),1,tmp);
                data_mr.om62.s_pres = repmat(nanmean(data_mr.om62.s_pres,2),1,tmp);
                data_mr.om62.u(:,:,1) = repmat(nanmean(data_mr.om62.u(:,:,1),2),1,tmp);
                data_mr.om62.v(:,:,1) = repmat(nanmean(data_mr.om62.v(:,:,1),2),1,tmp);
                data_mr.om62.t(:,:,1) = repmat(nanmean(data_mr.om62.t(:,:,1),2),1,tmp);
                data_mr.om62.s(:,:,1) = repmat(nanmean(data_mr.om62.s(:,:,1),2),1,tmp);
                
                
                
                % update the mat file
                save('TemporaryFile_mooring_data_ave.mat','data_mr');
                
                
            end
% ========================================================================= 

            
% ========================================================================= 
% FLAG_TEST_M2M3_GEO (1/4) ================================================
% ========================================================================= 
% Extending UMM2 with the deepest instrument from D4, 
% and excluding D4 and D5 entirely;
%
            if flag_test_M2M3_geo== 1 && ~strcmp(flag_section,'west') && i_t== i_t_begin && cnt_mc== 1
                
                
                % load mooring datat at instrument levels
                dir_osnap_obs = 'TemporaryFile_mooring_data_ave.mat'; % data_mr_ave.om[number]
                load(dir_osnap_obs);


                clear tmp tmp1 tmp2

                % found M2, D4, D5
%                 tmp1 = [-28.0200 -26.9673 -25.6778]; % M2 (OM51), D4 (OM52), D5 (OM62) 
                

                % add deepest MCTD from D4 to M2
                data_mr.om51.t_pres(end+1,:) = data_mr.om52.t_pres(end,:);
                data_mr.om51.t(end+1,:,:) = data_mr.om52.t(end,:,:);
                data_mr.om51.s_pres(end+1,:) = data_mr.om52.s_pres(end,:);
                data_mr.om51.s(end+1,:,:) = data_mr.om52.s(end,:,:);
                
                
                % remove D4
                data_mr.om52.v_pres = nan;
                data_mr.om52.t_pres = nan;
                data_mr.om52.s_pres = nan;
                data_mr.om52.u = nan;
                data_mr.om52.v = nan;
                data_mr.om52.t = nan;
                data_mr.om52.s = nan;
                
                % remove D5
                data_mr.om62.v_pres = nan;
                data_mr.om62.t_pres = nan;
                data_mr.om62.s_pres = nan;
                data_mr.om62.u = nan;
                data_mr.om62.v = nan;
                data_mr.om62.t = nan;
                data_mr.om62.s = nan;
                
                
                % update the mat file
                disp('+ Mooring data: Modified data_mr by constructing UMM2_extended, and removing D4, D5...');
                save('TemporaryFile_mooring_data_ave.mat','data_mr');
                
            end
%% =========================================================================  
            if flag_test_IB5==0 && ~strcmp(flag_section,'west') && i_t== i_t_begin && cnt_mc== 1
                
                 % load mooring datat at instrument levels
                dir_osnap_obs = 'TemporaryFile_mooring_data_ave.mat'; % data_mr_ave.om[number]
                load(dir_osnap_obs);
                
                % remove IB5
                data_mr.om63.v_pres = nan;
                data_mr.om63.t_pres = nan;
                data_mr.om63.s_pres = nan;
                data_mr.om63.u = nan;
                data_mr.om63.v = nan;
                data_mr.om63.t = nan;
                data_mr.om63.s = nan;
                
                % update the mat file
                disp('+ Mooring data: IB5 removed with flag_test_IB5=0');
                save('TemporaryFile_mooring_data_ave.mat','data_mr');
                
            end
%             return
            if flag_test_UMD123_geo==1
                
                % load mooring datat at instrument levels
                dir_osnap_obs = 'TemporaryFile_mooring_data_ave.mat'; % data_mr_ave.om[number]
                load(dir_osnap_obs);
                
                % remove D1
                data_mr.om48.u(data_mr.om48.v_pres(:,70)<1200,:,:) = nan;
                data_mr.om48.v(data_mr.om48.v_pres(:,70)<1200,:,:) = nan;
                data_mr.om48.t(data_mr.om48.t_pres(:,70)<1200,:,:) = nan;
                data_mr.om48.s(data_mr.om48.s_pres(:,70)<1200,:,:) = nan;
                data_mr.om48.v_pres(data_mr.om48.v_pres(:,70)<1200,:,:) = nan;
                data_mr.om48.t_pres(data_mr.om48.t_pres(:,70)<1200,:,:) = nan;
                data_mr.om48.s_pres(data_mr.om48.s_pres(:,70)<1200,:,:) = nan;
                
                % remove D2
                data_mr.om49.u(data_mr.om49.v_pres(:,70)<1200,:,:) = nan;
                data_mr.om49.v(data_mr.om49.v_pres(:,70)<1200,:,:) = nan;
                data_mr.om49.t(data_mr.om49.t_pres(:,70)<1200,:,:) = nan;
                data_mr.om49.s(data_mr.om49.s_pres(:,70)<1200,:,:) = nan;
                data_mr.om49.v_pres(data_mr.om49.v_pres(:,70)<1200,:,:) = nan;
                data_mr.om49.t_pres(data_mr.om49.t_pres(:,70)<1200,:,:) = nan;
                data_mr.om49.s_pres(data_mr.om49.s_pres(:,70)<1200,:,:) = nan;
                
                % remove D3
                data_mr.om50.u(data_mr.om50.v_pres(:,70)<1200,:,:) = nan;
                data_mr.om50.v(data_mr.om50.v_pres(:,70)<1200,:,:) = nan;
                data_mr.om50.t(data_mr.om50.t_pres(:,70)<1200,:,:) = nan;
                data_mr.om50.s(data_mr.om50.s_pres(:,70)<1200,:,:) = nan;
                data_mr.om50.v_pres(data_mr.om50.v_pres(:,70)<1200,:,:) = nan;
                data_mr.om50.t_pres(data_mr.om50.t_pres(:,70)<1200,:,:) = nan;
                data_mr.om50.s_pres(data_mr.om50.s_pres(:,70)<1200,:,:) = nan;
                
                % update the mat file
                disp('+ Mooring data: UMD1-3 data above 1200 m removed with flag_test_UMD123_geo=1');
                save('TemporaryFile_mooring_data_ave.mat','data_mr');
                
            end
%             return
            if flag_test_remove_NOCM5==1
                
                % load mooring datat at instrument levels
                dir_osnap_obs = 'TemporaryFile_mooring_data_ave.mat'; % data_mr_ave.om[number]
                load(dir_osnap_obs);
                % remove NOCM5
                data_mr.om41.v_pres = nan;
                data_mr.om41.t_pres = nan;
                data_mr.om41.s_pres = nan;
                data_mr.om41.u = nan;
                data_mr.om41.v = nan;
                data_mr.om41.t = nan;
                data_mr.om41.s = nan;
                
                % update the mat file
                disp('+ Mooring data: NOCM5 removed with flag_test_remove_NOCM5=1');
                save('TemporaryFile_mooring_data_ave.mat','data_mr');
                
            end
            
            if flag_test_remove_LSAB==1
                 % load mooring datat at instrument levels
                dir_osnap_obs = 'TemporaryFile_mooring_data_ave.mat'; % data_mr_ave.om[number]
                load(dir_osnap_obs);
                % remove NOCM5
                data_mr.om64.v_pres = nan;
                data_mr.om64.t_pres = nan;
                data_mr.om64.s_pres = nan;
                data_mr.om64.u = nan;
                data_mr.om64.v = nan;
                data_mr.om64.t = nan;
                data_mr.om64.s = nan;
                
                data_mr.om65.v_pres = nan;
                data_mr.om65.t_pres = nan;
                data_mr.om65.s_pres = nan;
                data_mr.om65.u = nan;
                data_mr.om65.v = nan;
                data_mr.om65.t = nan;
                data_mr.om65.s = nan;
                
                % update the mat file
                disp('+ Mooring data: NOCM5 removed with flag_test_remove_NOCM5=1');
                save('TemporaryFile_mooring_data_ave.mat','data_mr');
            end

% IF FLAG_MC        
            if flag_MC == 1 && flag_data_interval > 1
                  
                % random
                [osnap_mr] = fun_osnap_mooring_vinterp_random_test_IB5(flag_section,1,flag_MC_SE,time(i_t),dirio);
            
            else
                % regridding daily obs.
                [osnap_mr] = fun_osnap_mooring_vinterp_random_test_IB5(flag_section,0,1,time(i_t),dirio);
                
% END IF FLAG_MC            
            end


% END IF FLAG_REGRIDDED_PROFILES   
        end
% =======================================================================
% ========================================================================= 



% ========================================================================= 
% FLAG_C3_RECON ===========================================================
% =========================================================================  
        if flag_C3_recon== 1 && (strcmp(flag_section, 'west') || strcmp(flag_section, 'all'))

            % Determine the regression coefficient
            if ~exist('stats_v1','var')

                % use daily-regridded data, for only once,
                % to calculate regression coefficients
                if ~exist('osnap','var')
                    data_mr_osnap= [dirio,'outputs/regridded/20210616_mooring_data_OSNAP_All_20140701-20200630_REGRIDDED_all_4noMC.mat'];
%                     data_mr_osnap= [dirio,'data/intermediate_output/0101_mooring_data_OSNAP_All_20140701-20180630_REGRIDDED.mat']; 
                    load(data_mr_osnap,'osnap');
                end

                % *** TS ***
                % Regress C3 TS on K7 based on the mean profiles
                t1 = squeeze(osnap.t(5,:,:)); 
                s1 = squeeze(osnap.s(5,:,:));
                t2 = squeeze(osnap.t(6,:,:)); % K7
                s2 = squeeze(osnap.s(6,:,:)); % K7


                % found time with good C3 data (the 2014-2016 deployment)
                ipt = osnap.time< datenum(2016,5,18,12,0,0);


                % mean profiles
                tp1 = nanmean(t1(:,ipt),2);
                tp2 = nanmean(t2(:,ipt),2);
                ipzs = ~isnan(tp1);
                stats_t1 = regstats(tp1(ipzs),tp2(ipzs),'linear',{'beta','rsquare'});

                sp1 = nanmean(s1(:,ipt),2);
                sp2 = nanmean(s2(:,ipt),2);
                ipzs = ~isnan(sp1);
                stats_s1 = regstats(sp1(ipzs),sp2(ipzs),'linear',{'beta','rsquare'});


                clear tp1 tp2 sp1 sp2


                % *** UV ***
                % Regress C3 TS on K7 based on the mean profiles
                u1 = squeeze(osnap.u(5,:,:)); 
                v1 = squeeze(osnap.v(5,:,:));
                u2 = squeeze(osnap.u(6,:,:)); % K7
                v2 = squeeze(osnap.v(6,:,:)); % K7


                % found time with good C3 data (the 2014-2016 deployment)
                ipt = osnap.time< datenum(2016,5,18,12,0,0);


                % mean profiles
                vp1 = nanmean(v1(:,ipt),2);
                vp2 = nanmean(v2(:,ipt),2);
                ipzv = ~isnan(vp1); % valid water depth
                stats_v1 = regstats(vp1(ipzv),vp2(ipzv),'linear',{'beta','rsquare'});

                up1 = nanmean(u1(:,ipt),2);
                up2 = nanmean(u2(:,ipt),2);
                stats_u1 = regstats(up1(ipzv),up2(ipzv),'linear',{'beta','rsquare'});

                clear vp1 vp2 up1 up2

            end


            % Reconstruct UVTS at C3
            %  -When no good C3 data were returned (the 2016-2018 deployment);
            %  -Overwrite the time-mean data in the original C3 data file for
            %       that time period;
            if time(i_t) >= datenum(2016,5,18,12,0,0)
                osnap_mr.t(5,ipzs) = osnap_mr.t(6,ipzs).*stats_t1.beta(2) + stats_t1.beta(1);
                osnap_mr.s(5,ipzs) = osnap_mr.s(6,ipzs).*stats_s1.beta(2) + stats_s1.beta(1);

                osnap_mr.v(5,ipzv) = osnap_mr.v(6,ipzv).*stats_v1.beta(2) + stats_v1.beta(1);
                osnap_mr.u(5,ipzv) = osnap_mr.u(6,ipzv).*stats_u1.beta(2) + stats_u1.beta(1);
            end

            % clear regridded mooring data
            clear osnap
            
        end
% ========================================================================= 
% =========================================================================  


            
% % ========================================================================= 
% % FLAG_TEST_NORR (1/2) ====================================================
% % =========================================================================  
% % Remove the data from the RR arrays ...
% if flag_test_noRR== 1 && strcmp(flag_section,'all')
%     
%     % moorings to remove
%     %   e.g., between IC0 and UM-D4
%     tmp = [-035.1258 +059.2147 2938;... % IC0 (OM42)
%            -026.9678 +058.0097 2670]; % UM-D4 (OM52)
%        
%     % found grid indices
%     tmp1 = osnap_mr.lon>= tmp(1,1) & osnap_mr.lon<= tmp(2,1);
%     
%     % nan the data
%     osnap_mr.t(tmp1,:) = nan;
%     osnap_mr.s(tmp1,:) = nan;
%     osnap_mr.d(tmp1,:) = nan;
%     osnap_mr.u(tmp1,:) = nan;
%     osnap_mr.v(tmp1,:) = nan;
%         
% end
% % =========================================================================   




% ========================================================================= 
% FLAG_TEST_FEWERMR =======================================================
% ========================================================================= 
        if flag_test_fewermr== 1 && strcmp(flag_section,'all')

            % moorings to remove
            %   C3,K8,LS6,4,2,CF1,3,5,7,IC2,4;
            tmp = [5 7 17 19 21 23 25 27 29 39 41];

            % nan the data
            osnap_mr.t(tmp,:) = nan;
            osnap_mr.s(tmp,:) = nan;
            osnap_mr.d(tmp,:) = nan;
            osnap_mr.u(tmp,:) = nan;
            osnap_mr.v(tmp,:) = nan;

            clear tmp

        end
% ========================================================================= 



% =========================================================================
% FLAG_MEAN_RTADCP ========================================================
% =========================================================================
% The mooring data have been pre-filled with the time-mean velocities when
% no data are returned.  
        if flag_mean_RTADCP== 0 

            if i_t== i_t_begin && cnt_mc== 1
                disp('- ');
                disp('- Mooring data: Copy EB1 data to RTADCP when no data returned from RTADCP...');
            end

            % found mooring
            tmp = -9.3386; % RTADCP
            tmp1= find(osnap_mr.lon== tmp);

            % overwrite with EB1 data
            if time(i_t)<= datenum(2014,10,30) || time(i_t)>= datenum(2015,6,19)
                osnap_mr.u(tmp1,1:38) = osnap_mr.u(tmp1-1,1:38); % for RTADCP range only
                osnap_mr.v(tmp1,1:38) = osnap_mr.v(tmp1-1,1:38); % for RTADCP range only
            end
            
            clear tmp tmp1 tmp2
            
        end
% =========================================================================



% ========================================================================= 
% FLAG_TEST_MR_NOUSE ======================================================
% =========================================================================
% Test of excluding certain moorings from the calculations...
        if flag_test_mr_nouse== 1 && ~strcmp(flag_section,'west')
            
            
            % ************** EDIT BELOW THIS LINE *************
            tmp = nan; % mooring
            tmpz = nan; % depth range
            
            if flag_test_rrex_nouse== 1
                
                % *** RREX ***
                if i_t== i_t_begin && cnt_mc== 1
                    disp('- Mooring data: RREX mooring data excluded...');
                end            

                % list of moorings to remove
                tmp(end+1:end+3) = [-33.2592 -32.1593 -31.5585]; % RREX moorings (OM59,60,61)
                tmpz(end+1:end+3) = sz_z;
            end
            
            if flag_test_d5_nouse== 1
                
                % *** D5 ***
                if i_t== i_t_begin && cnt_mc== 1
                    disp('- Mooring data: D5 [all] and D4 [top] data excluded...');
                end            

                % list of moorings to remove
                tmp(end+1:end+2) = [-26.9673 -25.6778]; % D4 (OM52) D5 (OM62)
                tmpz(end+1:end+2) = [115 sz_z]; % removing data above 2290m at D4; removing all D5 data;
                
            end
            
            tmp(isnan(tmp)) = [];
            tmpz(isnan(tmpz)) = [];
            
            % just in case...
            if isempty(tmp)
                error('ODI 918!')
            end
            
            % **************** EDIT ABOVE THIS LINE ***********
            
            
            % loop through listed mooring
            for ii = 1:length(tmp)
                
                % found mooring
                tmp1 = find(osnap_mr.lon== tmp(ii));
                tmp2 = tmpz(ii);
                
                if isempty(tmp1)
                    error('JDA 910!')
                end
                
                % nan the data
                osnap_mr.t(tmp1,1:tmp2) = nan;
                osnap_mr.s(tmp1,1:tmp2) = nan;
                osnap_mr.d(tmp1,1:tmp2) = nan;
                osnap_mr.u(tmp1,1:tmp2) = nan;
                osnap_mr.v(tmp1,1:tmp2) = nan;

            end
            
            clear tmp tmpz ii tmp1 tmp2
            
        end
% =========================================================================




        %% *************************
        %% *** SSH and WIND data ***
        %% *************************

        % For this time interval
        ipt1 = find(time_adt<= time(i_t),1,'last'); 
        ipt2 = find(abs(time_wind-(time(i_t))) == min(abs(time_wind-(time(i_t))))); 

        
%         % LOG
%         if cnt_mc == 1
%             disp([' Found time_adt= ',datestr(time_adt(ipt1),'yyyy-mm-dd HH:MM:SS')]);
%             disp([' Found time_wind= ',datestr(time_wind(ipt2),'yyyy-mm-dd HH:MM:SS')]);
%         end
        
% =========================================================================
% FLAG_MC =================================================================
% =======================================================================
% IF FLAG_MC
        if flag_MC == 1 && flag_data_interval > 1
            % randomly draw values
        
            if flag_MC_SE == 1
                % draw from mean +/- SE 
                ssh_t  = ssh(:,ipt1,1)  + ssh(:,ipt1,3).*randn(sz_n,1);
                taux_t = taux(:,ipt2,1) + taux(:,ipt2,3).*randn(sz_n,1);
                tauy_t = tauy(:,ipt2,1) + tauy(:,ipt2,3).*randn(sz_n,1);
                u10_t  = u10(:,ipt2,1)  + u10(:,ipt2,3).*randn(sz_n,1);
                v10_t  = v10(:,ipt2,1)  + v10(:,ipt2,3).*randn(sz_n,1);

            elseif flag_MC_SE == 0
                % draw from mean +/- SD
                ssh_t  = ssh(:,ipt1,1)  + ssh(:,ipt1,2).*randn(sz_n,1);
                taux_t = taux(:,ipt2,1) + taux(:,ipt2,2).*randn(sz_n,1);
                tauy_t = tauy(:,ipt2,1) + tauy(:,ipt2,2).*randn(sz_n,1);
                u10_t  = u10(:,ipt2,1)  + u10(:,ipt2,2).*randn(sz_n,1);
                v10_t  = v10(:,ipt2,1)  + v10(:,ipt2,2).*randn(sz_n,1);

            else
                error('ADC 910!');
                
            end



        else
            
            % Temporally means or daily values
            ssh_t = ssh(:,ipt1,1);
            taux_t = taux(:,ipt2,1);
            tauy_t = tauy(:,ipt2,1);
            u10_t = u10(:,ipt2,1);
            v10_t = v10(:,ipt2,1);
                

% =========================================================================
% FLAG_TEST_AGV ===========================================================
% =====================================================================
            if flag_test_AGV== 1
                ipt3 = find(floor(time_agv) == floor(time(i_t))); % note: AGV's time is at 00:00:00 of each day
                agv_t = agv(:,ipt3,1);
            end
% =====================================================================
% =========================================================================



% END IF FLAG_MC   
        end
% =========================================================================
% =========================================================================

        

        %% **********************
        %% *** Load to t_grid ***
        %% **********************
        % Initialize t-grid
        t_t = nan(sz_n,sz_z);
        s_t = nan(sz_n,sz_z);
        u_t = nan(sz_n,sz_z);
        v_t = nan(sz_n,sz_z);
        d_t = nan(sz_n,sz_z);
        pd_t = nan(sz_n,sz_z);
        pt_t = nan(sz_n,sz_z);


        % Retrive gridded observations
        t_o = osnap_mr.t;
        s_o = osnap_mr.s;
        d_o = osnap_mr.d;
        u_o = osnap_mr.u;
        v_o = osnap_mr.v;
        lon_o = osnap_mr.lon;
        lat_o = osnap_mr.lat;
        sz_o = length(lon_o);

        % variable to store the geo-location of OSNAP instruments
        ip_obs = nan(sz_o,1); 

        
        % Fill the t-grid with observed values
% LOOP MEASURED PROFILES
        for i_o = 1:sz_o        
            
            % Found mooring indices
            ip_obs(i_o) = find(lon==lon_o(i_o)); 

            % Load observations
            t_t(ip_obs(i_o),:) = t_o(i_o,:);  
            s_t(ip_obs(i_o),:) = s_o(i_o,:);
            d_t(ip_obs(i_o),:) = d_o(i_o,:);
            u_t(ip_obs(i_o),:) = u_o(i_o,:);
            v_t(ip_obs(i_o),:) = v_o(i_o,:);  

            % Calculate dens/pden/ptmp from each obs. T/S profile
            pres = gsw_p_from_z(-depth,lat_o(i_o))';
            [SA,~,~] = gsw_SA_Sstar_from_SP(s_o(i_o,:),pres,lon_o(i_o),lat_o(i_o)); % absolute salinity             
            pd_t(ip_obs(i_o),:) = gsw_pot_rho_t_exact(SA,t_o(i_o,:),pres,0); % potential density 
            pt_t(ip_obs(i_o),:) = gsw_pt_from_t(SA,t_o(i_o,:),pres,0); % potential temperature      

            
% END LOOP MEASURED PROFILES                
        end
        
        
        
        
        
        %% *********************
        %% *** Update t_grid ***
        %% *********************
        % //////////////////////////////////////////////////////////
        % // Update observations ///////////////////////////////////
        % //////////////////////////////////////////////////////////
        
% Disabled in July 2019
%         % *******************************
%         % OM6 (C3) 2016-2018 CM data loss
%         % *******************************
%         % Added on 3/31/2019
%         %
%         % Copy velocities from OM7 (K7) to OM6
%         % when there are no data returned from the ADCP on OM6;
%         if strcmp(flag_section,'west') || strcmp(flag_section,'all')
%             ipx = 5; % mr index for OM6
%             rg = mask_land(ip_obs(ipx),:)==1;
%             u_t(ip_obs(ipx),rg) = u_t(ip_obs(ipx+1),rg); % copy U 
%             v_t(ip_obs(ipx),rg) = v_t(ip_obs(ipx+1),rg); % copy V
%         end
        
        
        % *******************************
        % OM5 (C2) 2016-2018 ADCP loss
        % *******************************
        % Added on 3/31/2019;
        % Updated on 8/6/2019, only copy to where nans appear;
        %
        % Copy velocities from OM6 to OM5 (high correlation),
        % when there are no data returned from the ADCP on OM3;
        if strcmp(flag_section,'west') || strcmp(flag_section,'all')
            ipx = 4; % mr index for OM5
            rg = mask_land(ip_obs(ipx),:)==1 & isnan(v_t(ip_obs(ipx),:));
            u_t(ip_obs(ipx),rg) = u_t(ip_obs(ipx+1),rg); % copy U 
            v_t(ip_obs(ipx),rg) = v_t(ip_obs(ipx+1),rg); % copy V
        end
        
        
        % ****************************************
        % OM3 (C1) 2015-2016 ADCP reduced range
        % ****************************************
        % May not needed anymore, data are filled to the bottom;
        % Updated on 8/6/2019, only copy to where nans exist;
        %
        % Copy velocities from OM5 to OM3, because
        % OM3 ADCP only had data returned for the upper layers;
        %
        if strcmp(flag_section,'west') || strcmp(flag_section,'all')
            ipx = 2; % mr index for OM2
            rg = mask_land(ip_obs(ipx),:)==1 & isnan(v_t(ip_obs(ipx),:));
            u_t(ip_obs(ipx),rg) = u_t(ip_obs(ipx+2),rg); % copy U 
            v_t(ip_obs(ipx),rg) = v_t(ip_obs(ipx+2),rg); % copy V
        end
        
        
        % ***************
        % OM2 (C1)
        % ***************
        % Copy velocities from OM3 to OM2 because
        % OM2 has only T/S records and it is west of OM3 
        % so cannot be interpolated from OM3;
        if strcmp(flag_section,'west') || strcmp(flag_section,'all')
            ipx = 1; % mr index for OM2
            rg = mask_land(ip_obs(ipx),:)==1;
            u_t(ip_obs(ipx),rg) = u_t(ip_obs(ipx+1),rg); % copy U 
            v_t(ip_obs(ipx),rg) = v_t(ip_obs(ipx+1),rg); % copy V
        end
        
        
        % ***************
        % OM24 (CF1)
        % ***************
        % Copy velocities from CF2 to CF1 
        % when there is no data from the ADCP on CF1;
        if strcmp(flag_section,'east') 
            ipx = 1; % mr index for OM24
            rg = mask_land(ip_obs(ipx),:)==1;
            u_t(ip_obs(ipx),rg) = u_t(ip_obs(ipx+1),rg); % copy U 
            v_t(ip_obs(ipx),rg) = v_t(ip_obs(ipx+1),rg); % copy V
        
        elseif strcmp(flag_section,'all')
            ipx = 23; % mr index for OM24
            rg = mask_land(ip_obs(ipx),:)==1;
            u_t(ip_obs(ipx),rg) = u_t(ip_obs(ipx+1),rg); % copy U 
            v_t(ip_obs(ipx),rg) = v_t(ip_obs(ipx+1),rg); % copy V
        end
        
        
        
%         % ***************
%         % RTADCP
%         % ***************
%         % When there are no data returned from RTADCP:
%         % Copy velocity from EB1 to RTADCP...
%         %   This is the same as copying EB1 data
%         %       into the wedge east of the mooring (EB1)
%         %       when no ADCP data are available;
%         %   This has been used as 'default' before 
%         %       10/5/2018; afterwards, this doesn't change anything 
%         %       because there are always data 'returned' i.e.,
%         %       using the time-mean RTADCP data;
%         if strcmp(flag_section,'east') || strcmp(flag_section,'all')
%             % Found mr index for RTADCP
%             switch flag_section
%                 case 'east'
%                     ipx = 31;
% 
%                 case 'all'
%                     ipx = 53;
%             end
% 
%             rg = isnan(v_t(ip_obs(ipx),:)) & mask_land(ip_obs(ipx),:)==1;
%             
%             if any(rg~=0) && flag_mean_RTADCP== 1
%                 % if flag_mean_RTADCP= 1, there must always be data
%                 % returned;
%                 % need to check fun_osnap_mooring_vinterp_random_v*
%                 error('DSK 291!');
%             end
%             
%             u_t(ip_obs(ipx),rg) = u_t(ip_obs(ipx-1),rg); % copy U from EB1
%             v_t(ip_obs(ipx),rg) = v_t(ip_obs(ipx-1),rg); % copy V from EB1
%         end 
        % //////////////////////////////////////////////////////////
        % //////////////////////////////////////////////////////////

        
% % =========================================================================
% % FLAG_TEST_RTADCP ========================================================
% % =====================================================================
%         if flag_test_RTADCP== 1
%         % Copying the data from EB1 to overwrite RTADCP    
%             if strcmp(flag_section,'east') || strcmp(flag_section,'all')
%                 switch flag_section
%                     case 'east'
%                         ipx = 34;
% 
%                     case 'all'
%                         ipx = 56;
%                 end 
%                 
%                 % copying u and v from EB1
%                 rg = mask_land(ip_obs(ipx),:)==1;
%                 u_t(ip_obs(ipx),rg) = u_t(ip_obs(ipx-1),rg);
%                 v_t(ip_obs(ipx),rg) = v_t(ip_obs(ipx-1),rg);
%             end
%             
%             
%         elseif flag_test_RTADCP== 2
%         % Using time-mean ADCP data ...    
%             
%             % This test requires modification 
%             % in fun_osnap_mooring_vinterp_random_v2.m!!!
%             % USE "if om== 58" INSTEAD OF
%             % "if om== 58 && all(isnan(vp_om))"!!!
%             
%             if i_t == i_t_begin
%                 disp('....... Make sure for flag_test_RTADCP= 2, modification''s been made to fun_osnap_mooring_vinterp_random_v2.m!');
%             end
%             
%         end
% % =========================================================================
% % =========================================================================   
        
     

        % ---------------------------------------------
        %% (2.2) CREATE MASKS
        % ---------------------------------------------

        %% ******************
        %% *** Mode masks ***
        %% ******************
        %   0- no calculation will be done; 1- mode specified.
        mask_mode_d = zeros(sz_n,1);
        mask_mode_g = zeros(sz_n,1);
        
% FLAG_SECTION            
        switch flag_section

            case 'west'
            % *** WEST ***
            % Same for all input
%
%	[Table 1]  24 waypoints for OSNAP West. 
%	------------------------------------------------------------
%	WEST#       | Mooring ID
%	------------+-----------------------------------------------
%	1-5         |   OM2-6 (C1,C2,C3)
%   6-9         |   OM7-8 (K7-8), OM9 (DSOW1), OM10 (K9), 
%	10-12       |   OM11 (DSOW2), OM12 (K10), OM13 (DSOW5)
%   13          |   OM14 (DSOW3)
%	14,15,16-22	|   OM15 (LS8), OM16(DSOW4), OM17-23(LS7-1)
%   23-24       |   OM64-65 (LSB, LSA)
%	------------------------------------------------------------
%   * The # is to use with ip_obs and *_o variables.

% LOOP T-GRID
                for i_n = 1:sz_n
                    % Segments with directly measured velocities
                    if i_n>=ip_obs(1) && i_n<=ip_obs(12) ...      % OM2-DSOW5
                            || i_n>=ip_obs(13) && i_n<=ip_obs(24) % DSOW3-LSA
                        mask_mode_d(i_n) = 1;
                    end

                    % Segments with geostrophic velocities
                    if i_n>=ip_obs(11) && i_n<=ip_obs(14)	% K10-LS8
                        mask_mode_g(i_n) = 1;
                    end
% END LOOP T-GRID
                end


            case 'east'
            % *** EAST ***
%
%       [Table 2]  35 waypoints for OSNAP East.
%       -------------------------------------------------------
%       EAST#       |   Mooring ID
%       ---------------------------------------------------------------------
%       1-7,8	    |	OM24-30 (CF1-7), OM31 (NOCM1)
%       9,10        |   OM32,OM33 (NOCM2,NOCM3) 
%       11,12       |   OM34 (FLMA), OM39 (FLMB)
%       13,14       |   OM40 (NOCM4), OM41(NOCM5)
%       15,16       |   OM42-43(IC0, IC1)
%       17-21       |   OM59 (IRW), OM44 (IC2), OM60 (IRM), OM45 (IC3), OM61 (IRE)
%       22-23       |   OM46 (IC4), OM47 (UMM1)
%       24-26       |   OM48-OM50 (UMD1,UMD2,UMD3) 
%       27-29       |   OM51 (UMM2), OM52 (UMD4), OM62 (UMD5)
%       30-32       |   OM53 (UMM3), OM54 (UMM4), OM63 (IB5)
%       33-34       |   OM55 (RTWB1), OM56 (RTWB2)
%       35-36       |   OM57 (RTEB1), OM58 (RTADCP)
%	--------------------------------------------------------
%       * The # is to use with ip_obs and *_o variables.
%
% LOOP T-GRID
                for i_n = 1:sz_n
                    % Segments with directly measured velocities
                    if i_n>=ip_obs(1) && i_n<=ip_obs(14) ...          	% CF1-(NOC)M5
                            || i_n>=ip_obs(15) && i_n<=ip_obs(30) ...	% IC0-(UM)M3
                            || i_n>=ip_obs(33) && i_n<=ip_obs(34) ...   % RTWB1-RTWB2
                            || i_n>=ip_obs(35) && i_n<=ip_obs(36)       % RTEB1-RTADCP
                        mask_mode_d(i_n) = 1;
                    end

                    % Segments with geostrophic velocities
                    if i_n>=ip_obs(8) && i_n<=ip_obs(16) ...            % (NOC)M1-IC1
                            || i_n>=ip_obs(23) && lon(i_n)<= -14.7 ...	% (UM)M1- SAMS glider eastern endpoint
                            || i_n>=ip_obs(33) && i_n<= ip_obs(35)      % RTWB1-RTEB1
                        mask_mode_g(i_n) = 1;
                    end
% END LOOP T-GRID
                end
                        
            case 'all'
            % *** All ***   
%
%	[Table 3] Total of 60 moorings for 2014-2020:
%	---------------------------------------------------------------------
% 	All#        |   Mooring ID
%	------------+--------------------------------------------------------
%	1-5         |   OM2-6 (C1,C2,C3)
%   6-9         |   OM7-8 (K7,K8), OM9 (DSOW1), OM10 (K9), 
%	10-12       |   OM11 (DSOW2), OM12 (K10), OM13 (DSOW5)
%   13          |   OM14 (DSOW3)
%	14,15,16-22	|   OM15 (LS8), OM16(DSOW4), OM17-23(LS7-1)
%   23-24       |   OM64-65 (LSB, LSA)
%	---------------------------------------------------------------------
% 	25-31,32	|	OM24-30 (CF1-7), OM31 (NOCM1)
%	33,34       |   OM32,OM33 (NOCM2,NOCM3) 
% 	35,36       |   OM34 (FLMA), OM39 (FLMB)
% 	37,38       |   OM40 (NOCM4), OM41(NOCM5)
%   39,40       |   OM42-43(IC0, IC1)
%   41-45       |   OM59 (IRW), OM44 (IC2), OM60 (IRM), OM45 (IC3), OM61 (IRE)
%   46-47       |   OM46 (IC4), OM47 (UMM1)
% 	48-50       |   OM48-OM50 (UMD1,UMD2,UMD3) 
% 	51-53       |   OM51 (UMM2), OM52 (UMD4), OM62 (UMD5)
%  	54-56       |   OM53 (UMM3), OM54 (UMM4), OM63 (IB5)
%	57-58       |   OM55 (RTWB1), OM56 (RTWB2)
%	59-60       |   OM57 (RTEB1), OM58 (RTADCP)
%	--------------------------------------------------------
%   * The # is to use with ip_obs and *_o variables.


% LOOP T-GRID
                for i_n = 1:sz_n
                    % Segments with directly measured velocities
                    if i_n>=ip_obs(1) && i_n<=ip_obs(12) ...            % OM2-DSOW5
                            || i_n>=ip_obs(13) && i_n<=ip_obs(24) ...   % DSOW3-LSA
                            || i_n>=ip_obs(25) && i_n<=ip_obs(38) ...	% CF1-(NOC)M5
                            || i_n>=ip_obs(39) && i_n<=ip_obs(54) ...	% IC0-(UM)M3
                            || i_n>=ip_obs(57) && i_n<=ip_obs(58) ...   % RTWB1-RTWB2
                            || i_n>=ip_obs(59) && i_n<=ip_obs(60)       % RTEB1-ADCP2
                        mask_mode_d(i_n) = 1;
            %                 % log
            %                 disp(['+ Mode mask for direct measurements: ',num2str(i_n)]);
                    end

                    % Segments with dynamic topography profiles (geostrophy)
                    if i_n>=ip_obs(11) && i_n<=ip_obs(14) ...           % K10-LS8
                            || i_n>=ip_obs(32) && i_n<=ip_obs(40) ...	% (NOC)M1-IC1
                            || i_n>=ip_obs(47) && lon(i_n)<= -14.7 ...	% (UM)M1- SAMS glider eastern endpoint
                            || i_n>=ip_obs(57) && i_n<= ip_obs(59)      % RTWB1-RTEB1
                        mask_mode_g(i_n) = 1;
            %                 % log
            %                 disp(['+ Mode mask for geostrophy: ',num2str(i_n)]); 
                    end
% END LOOP T-GRID 
                end



% END FLAG_SECTION
        end 

        
        
        %% ************************
        %% *** Update mode mask ***
        %% ************************

% % ========================================================================= 
% % FLAG_TEST_NORR (2/2) ====================================================
% % =========================================================================  
% % Modifying mode masks...
%         if flag_test_noRR~= 0
% 
%             % geostrophy between NOCM5-RTADCP
%             tmp = [-037.7992 +059.5771 3214;... % NOC-M5 (OM41)
%                    -014.7000 +058.0000 nan]; % SAMS Glider East Waypoint
%                
%                
%             % find data indices
%             ipx = find(lon>= tmp(1,1) & lon<= tmp(2,1));
%             
%             
%             % modify the mask
%             mask_mode_g(ipx) = 1;
%             mask_mode_d(ipx) = 0;
%             
%             
%             clear tmp tmp1
% 
%         end
% % =========================================================================   

% =========================================================================
% FLAG_TEST_M2M3_GEO (2/4) ================================================
% =========================================================================
% Modify mask_mode_d to keep only 
% geostrophic velocities between M2 and M3;
        if flag_test_M2M3_geo== 1 && ~strcmp(flag_section,'west')
            
            % between UMM2 and UMM3
            tmp = [-28.0200	58.0362;... % UM-M2
                   -24.4220	58.0130];   % UM-M3
               
               
            % find data indices
            ipx = find(lon> tmp(1,1) & lon<= tmp(2,1));
            
            
            % modify the mask
            mask_mode_d(ipx) = 0;
            
        end
% =========================================================================  


% =========================================================================
% FLAG_TEST_D5_NOUSE (2/2) ================================================
% =========================================================================
% If no D5 data are used, prevent interpoalting velocities 
% between D4 and UMM3 by updating mask_mode_d
        if flag_test_d5_nouse== 1 && ~strcmp(flag_section,'west')
            
            % between D4 and UMM3
            tmp = [-26.9673	58.0103;... % UM-D4
                   -24.4220	58.0130];   % UM-M3
               
               
            % find data indices
            ipx = find(lon> tmp(1,1) & lon<= tmp(2,1));
            
            
            % modify the mask
            mask_mode_d(ipx) = 0;
            
        end

% =========================================================================      
        
        

% =========================================================================
% flag_test_Vgeo_only_in_CMarray (1/2) ====================================
% =========================================================================
% Modifying mode masks to calculate 
% the geostrophic velocities in the area
% by disabling the use of any direct 
% current meter observations
        if flag_test_Vgeo_only_in_CMarray == 1
            
            % Found indices to update
            if strcmp(flag_section,'all')
                
                % ipo0 is unique througout this script
                % and will be called in the 2nd
                % part of this FLAG
                ipo0 = [6 12; ... % K7-DSOW5
                    13 19; ... % DSOW3-LS5
                    28 36; ... % CF6-NOCM5
                    37 52];    % IC0-UMM3
                
                % Found indices
                ipx = [ip_obs(ipo0(1,1)):ip_obs(ipo0(1,2)) ...	
                    ip_obs(ipo0(2,1)):ip_obs(ipo0(2,2)) ...       
                    ip_obs(ipo0(3,1)):ip_obs(ipo0(3,2)) ...       
                    ip_obs(ipo0(4,1)):ip_obs(ipo0(4,2))];         
            else
                error('XGI 001!');
            end
            
            
            
            % Update the masks
            mask_mode_g(ipx) = 1;
            mask_mode_d(ipx) = 0;
            
            
            
        end
% =========================================================================   
        


        %% ******************
        %% *** Data masks ***
        %% ******************
        % Indicate where observations are available
        %   1- with observations; 0- w/o observations.

        % Initialization
        mask_obs_uv = ones(sz_n,sz_z);
        mask_obs_uv(isnan(v_t)) = 0;

        mask_obs_ts = ones(sz_n,sz_z);
        mask_obs_ts(isnan(t_t)) = 0;
   
        
        % Adjustment #1
        % @ mask_obs_uv
        %   extend all short moorings as tall moorings in the
        %   geostrophy segments to avoid velocity interpolation
        %   over the short moorings from velocities at the bounding 
        %   tall dynamic height moorings;
        for i_p = 1:length(ip_obs)
            i_n = ip_obs(i_p);
            if mask_mode_g(i_n)==1 && mask_mode_d(i_n)==1
                % where short mooring are availble in the
                % geostrophy segments;
                if any(mask_obs_uv(i_n,:)==1)
                    mask_obs_uv(i_n,:) = 1;
                end
                
            elseif mask_mode_g(i_n)==0 && mask_mode_d(i_n)==1
                % boundary moorings segments;
                % do nothing;
                % interp over short moorings is fine;
            end
        end

        
        
        % Adjustment #2
        % @ mask_obs_ts & @ mask_obs_uv
        % All dry points are set to 1. 
        %   this allows horizontal interpolation
        %   between land and the measurements 
        %   in the bottom triangles;
        %   
        mask_obs_ts(mask_land==0) = 1;
        mask_obs_uv(mask_land==0) = 1;
        
        
        
        % Adjustment #3
        % @ mask_obs_ts & @ mask_obs_uv
        %   set to 0 for the dry nodes between 
        %   adjacent moorings.  This is to
        %   allow the horizontal interpolation
        %   above the common maximum depth;
        for i = 1:length(ip_obs)-1
            j = i+1;
            if mask_mode_d(ip_obs(i))==1 && mask_mode_d(ip_obs(j))==1
%                     disp([num2str(ip_obs(i)),' ',num2str(ip_obs(j))]);
                for m= ip_obs(i)+1:ip_obs(j)-1
                    mask_obs_uv(m,:) = 0;
                    mask_obs_ts(m,:) = 0;
                end
            end
        end

        
        
        % *** Ekman flux mask ***
        %   Places where geostrophy is applied
        %  
        mask_ek_i = zeros(sz_n-1,1);

        
        % *** Flux mask ***
        %   where geostrophic velocities are used and 
        %   flux correction will be
        %   needed
        %
        mask_flux_i = zeros(sz_n-1,sz_z-1); % 1-need, 0-no deed
        
        
% =========================================================================
% flag_test_M2M3_geo (3/4) ================================================
% =========================================================================
% Modifying mask_obs_ts to enable the calculation of 
%   velocity shears between M2(extended) and M3 over the dry nodes; 
        if flag_test_M2M3_geo == 1 && ~strcmp(flag_section,'west')

            clear tmp tmp1
            
            tmp = [-28.0200 -24.4220]; % UMM2, UMM3

            % nan
            tmp1 = lon> tmp(1) & lon< tmp(2);
            mask_obs_ts(tmp1,:) = 0;

        end
% =========================================================================
        
        
% % =========================================================================
% % flag_test_no_OOI_data ===================================================
% % =====================================================================
%         if flag_test_no_OOI_data == 1
%             
%             % Found mr index for FLMA and FLMB
%             switch flag_section
%                 case 'east'
%                     ipo1 = 11;
%                     ipo2 = 12;                 
%                     
%                 case 'all'
%                     ipo1 = 33;
%                     ipo2 = 34;
%                     
%             end
% 
% 
%             % *** FLMA ***
%             % Found depth indices for the full water column
%             ipz = mask_land(ip_obs(ipo1),:) == 1;
%             
%             % update data
%             u_t(ip_obs(ipo1),ipz) = nan;
%             v_t(ip_obs(ipo1),ipz) = nan;
%             
%             t_t(ip_obs(ipo1),ipz) = nan;
%             s_t(ip_obs(ipo1),ipz) = nan;
%             d_t(ip_obs(ipo1),ipz) = nan;
%             pd_t(ip_obs(ipo1),ipz) = nan;
%             pt_t(ip_obs(ipo1),ipz) = nan;
%             
%             % update masks
%             mask_obs_ts(ip_obs(ipo1),ipz) = 0;
%             mask_obs_uv(ip_obs(ipo1),ipz) = 0;
% 
%             
%             
%             % *** FLMB ***
%             ipz = mask_land(ip_obs(ipo2),:) == 1; 
%             
%             % update data
%             u_t(ip_obs(ipo2),ipz) = nan;
%             v_t(ip_obs(ipo2),ipz) = nan;
%             
%             t_t(ip_obs(ipo2),ipz) = nan;
%             s_t(ip_obs(ipo2),ipz) = nan;
%             d_t(ip_obs(ipo2),ipz) = nan;
%             pd_t(ip_obs(ipo2),ipz) = nan;
%             pt_t(ip_obs(ipo2),ipz) = nan;
%             
%             % update masks
%             mask_obs_ts(ip_obs(ipo2),ipz) = 0;
%             mask_obs_uv(ip_obs(ipo2),ipz) = 0;
%                     
%         end
% % =========================================================================
% % =========================================================================
        
        
% % ========================================================================= 
% % FLAG_TEST_NORR (extra) ==================================================
% % =========================================================================  
% % Modifying data mask to eliminate the RR...
%         if flag_test_noRR== 2
% 
%             % geostrophy between NOCM5-RTADCP
%             tmp = [-037.7992 +059.5771 3214;... % NOC-M5 (OM41)
%                     -024.4287 +058.0128 2850]; % UM-M3 (OM53)
%                
%                
%             % find data indices
%             ipx = find(lon> tmp(1,1) & lon< tmp(2,1));
%             
%             
%             % modify the mask
%             mask_obs_ts(ipx,:) = 0;
%             
%             
%             clear tmp tmp1
% 
%         end
% % =========================================================================  


%         % PLOT
%         pcolor(lon,depth,s_t'); shading flat;
%         axis ij; colorbar
%         title('(2.2) t\_t');
        

        % ----------------------------------------------------
        %% (2.3) FILL THE REMAINING T-GRID POINTS
        % ----------------------------------------------------
    %     % LOG
    %     display('+ Filling background t/s/dens/pden/ptmp to the entire section ...');   

% =========================================================================
% FLAG_OA =================================================================
% =======================================================================
% Fill the entire section with analysis
%
% FLAG_OA
        if flag_OA== 1
        % Fill analysis prepared separately
        
        	% Find the analysis fields averaged
            % in the same time interval
            ipt = find(time_analysis<= time(i_t)+ flag_data_interval_4calculation/2,1,'last'); 
            
%             if cnt_mc == 1
%                 % LOG
%                 disp([' Found time_analysis = ',datestr(time_analysis(ipt),'yyyy-mm-dd HH:MM:SS')]);
%             end
            
            t_t_analysis = t_analysis(:,:,ipt);
            s_t_analysis = s_analysis(:,:,ipt);
            d_t_analysis = d_analysis(:,:,ipt);
            pd_t_analysis = pd_analysis(:,:,ipt);
            pt_t_analysis = pt_analysis(:,:,ipt);
            
            
            % Fill gaps on t-grid with the analysis
            t_t(isnan(t_t)) = t_t_analysis(isnan(t_t));
            s_t(isnan(s_t)) = s_t_analysis(isnan(s_t));
            d_t(isnan(d_t)) = d_t_analysis(isnan(d_t));
            pd_t(isnan(pd_t)) = pd_t_analysis(isnan(pd_t));
            pt_t(isnan(pt_t)) = pt_t_analysis(isnan(pt_t));

        else
            
            % Fill with updated WOA13 monthly climatology
            t_t(isnan(t_t)) = t_t_c(isnan(t_t));
            s_t(isnan(s_t)) = s_t_c(isnan(s_t));
            pt_t(isnan(pt_t)) = pt_t_c(isnan(pt_t));
            pd_t(isnan(pd_t)) = pd_t_c(isnan(pd_t));
            d_t(isnan(d_t)) = d_t_c(isnan(d_t));

% END FLAG_OA
        end
% =======================================================================
% =========================================================================
       




% =========================================================================
% FLAG_TEST_M2M3_GEO (4/4)=================================================
% =========================================================================
% Modifying mask_land to accomendate M2-extended;
% Otherwise, d_t on land will be discarded when calculating g_v;
        if flag_test_M2M3_geo == 1 && ~strcmp(flag_section,'west')
            mask_land = d_t;
            mask_land(~isnan(d_t)) = 1;
            mask_land(isnan(d_t)) = 0;
        end
% =========================================================================



% =========================================================================
% FLAG_TEST_MODEL_CLIM_CASE1 ==============================================
% =======================================================================
% Test of using all model 
% climatology u/v everywhere ...
        if flag_test_model_clim_case1 == 1
            
            disp('+ All U/V/T/S are filled with FLAME UV only!');
            
            % month
            im = str2double(datestr(time(i_t)+flag_data_interval_4calculation/2,'mm')); 
                
%             
%             % model climatology
%             t_t = nanmean(t_clim(:,:,im),3);
%             s_t = nanmean(s_clim(:,:,im),3);
% 
%             % Derive dens, pden and ptmp
%             % from climatological T/S 
%             pres = gsw_p_from_z(-repmat(depth,sz_n,1),lat'); % pres
%             [SA,~,~] = gsw_SA_Sstar_from_SP(s_t,pres,lon',lat'); % absolute salinity 
%             pt_t = gsw_pt_from_t(SA,t_t,pres,0); % potential temperature  
%             d_t = gsw_rho_t_exact(SA,t_t,pres); % density 
%             pd_t = gsw_pot_rho_t_exact(SA,t_t,pres,0); % potential density
            
            
            % FLAME u/v
            u_t = nanmean(u_clim_model(:,:,im),3);
            v_t = nanmean(v_clim_model(:,:,im),3);
 
            
            % update the masks
            mask_obs_ts = mask_land;
            mask_obs_uv = mask_land;
            mask_mode_d = ones(sz_n,1);
            mask_mode_g = zeros(sz_n,1);
        end
% =======================================================================
% =========================================================================



%         % Remove invalid land point;
%         t_t(mask_land==0) = nan;
%         s_t(mask_land==0) = nan;
%         d_t(mask_land==0) = nan;
%         pd_t(mask_land==0) = nan;
%         pt_t(mask_land==0) = nan;
        

        
%         % PLOT
%         pcolor(lon,depth,s_t'); shading flat;
%         axis ij; colorbar
%         title('(2.3) t\_t');
        

        % :::::::::::::::::::::::::::::::::::::::::: NOTE
        % #1# t_t, s_t, u_t, v_t, d_t, pd_t, pt_t is filled;
        % #2# All masks used for control the calculations are created.


        
        % -----------------------------------------------
        %% (3) Convert t-grid to i-grid or v-grid
        % -----------------------------------------------
        % Update the i-grid (see below) 
        %
        % -------------------------
        %  
        %  i_n i_n i_n+1
        %    |  |  |
        %
        %    t  v  t  v  t   --  i_z
        %       i     i      --  i_z
        %    t  v  t  v  t   --  i_z+1
        %       i     i      --  i_z+1
        %    t  v  t  v  t    
        % -------------------------

        % v-grid (an intermediate grd)
        t_v = nan(sz_n-1,sz_z);
        s_v = nan(sz_n-1,sz_z);
        d_v = nan(sz_n-1,sz_z);
        pd_v = nan(sz_n-1,sz_z);
        pt_v = nan(sz_n-1,sz_z);
        
        v_v = nan(sz_n-1,sz_z);
        g_v = nan(sz_n-1,sz_z);

        % i-grid
        t_i = nan(sz_n-1,sz_z-1);
        pt_i = nan(sz_n-1,sz_z-1);
        s_i = nan(sz_n-1,sz_z-1); 
        d_i = nan(sz_n-1,sz_z-1); 
        pd_i = nan(sz_n-1,sz_z-1);
        v_i = nan(sz_n-1,sz_z-1);



        % ------------------------------------------------------
        %% (3.1) AVERAGE PROPERTY FROM T-GRID TO I-GRID
        % ------------------------------------------------------

% LOOP I-GRID     
        for i_n = 1:sz_n-1
% LOOP I-DEPTH
            for i_k = 1:sz_z-1

                % Domain averages
                % 
                % Using nanmean keeps the 
                % same values from the wet nodes
                % for the nodes between dry 
                % and wet nodes
                t_i(i_n,i_k) = nanmean(nanmean(t_t(i_n:i_n+1,i_k:i_k+1),1),2);
                s_i(i_n,i_k) = nanmean(nanmean(s_t(i_n:i_n+1,i_k:i_k+1),1),2);
                d_i(i_n,i_k) = nanmean(nanmean(d_t(i_n:i_n+1,i_k:i_k+1),1),2);
                pd_i(i_n,i_k) = nanmean(nanmean(pd_t(i_n:i_n+1,i_k:i_k+1),1),2);
                pt_i(i_n,i_k) = nanmean(nanmean(pt_t(i_n:i_n+1,i_k:i_k+1),1),2);

% END LOOP I-DEPTH
            end
            
            % Copy the deepest valid value into the bottom
            %   This is necessary since the bathymetry changes 
            %   between T and V/I grids;
            ipz = find(~isnan(t_i(i_n,:)),1,'last');
            t_i(i_n,ipz:end) = t_i(i_n,ipz);
            s_i(i_n,ipz:end) = s_i(i_n,ipz);
            d_i(i_n,ipz:end) = d_i(i_n,ipz);
            pd_i(i_n,ipz:end) = pd_i(i_n,ipz);
            pt_i(i_n,ipz:end) = pt_i(i_n,ipz);

            
% END LOOP I-GRID        
        end


        % remove land points
        t_i(mask_land_i==0) = nan;
        s_i(mask_land_i==0) = nan;
        d_i(mask_land_i==0) = nan;
        pd_i(mask_land_i==0) = nan;
        pt_i(mask_land_i==0) = nan;


        % :::::::::::::::::::::::::::::::::::::::::::: NOTE
        % #1# t_i, s_i, d_i, pd_i, and pt_i are obtained. 
        % #2# u_t, v_t contain velocity measurements only.




        % ---------------------------------------------------------------
        %% (3.2) INTERPOLATE MOORED VELOCITIES FROM T-GRID TO V-GRID
        % ---------------------------------------------------------------
        %
        
        % Interpolate between moorings;
        % V = 0 below the common depth of two moorings;
        %
    %     display('+ Interpolating u/v ...');
                
%         % Fill u/v on the dry points with zeros
%         u_t(mask_land==0) = 0;
%         v_t(mask_land==0) = 0;


% LOOP T-DEPTH
        for i_k = 1:sz_z

            % counters
            i = 1; % Left node counter
            j = 0; % Right node counter
            o = 0; % Count valid regions (between endpoints)

% SEARCH I ON T-GRID         
            while i <= sz_n-2
                if mask_obs_uv(i,i_k)==1
                    % Found left active node;
    %                 display(['+ Found ACTIVE, WET node i= ',num2str(i),' k= ',num2str(i_k)]);

% SEARCH J ON T-GRID
                    % right active node
                    j = i + 1;
                    while j <= sz_n-1
                        if mask_obs_uv(j,i_k)==1
                            % Found right active node;
    %                       display(['+ Found ACTIVE, WET node j= ',num2str(i),' k= ',num2str(i_k)]);

                            % Segments with direct u/v measurements
                            if (mask_mode_d(i+1) == 1 && mask_mode_d(j-1) == 1) 

                                % Compute velocity normal to section at two endpoints 
                                % * left node *
                                upl = u_t(i,i_k);
                                vpl = v_t(i,i_k);
                                v_l = -upl.*sin(ang_i(i))+vpl.*cos(ang_i(i));
                                
                                % * right node *
                                upr = u_t(j,i_k);
                                vpr = v_t(j,i_k);
                                v_r = -upr.*sin(ang_i(j))+vpr.*cos(ang_i(j)); 
                                                             
 
                                % Fill nodes between the bounding moorings                            
% LOOP M BETWEEN ENDPOINTS I AND J ON V-GRID  
                                for m = i:j-1
                                    % Update distances to 
                                    % the left and right active 
                                    % nodes
                                    dsl = gsw_distance([lon(i) lon_i(m)],[lat(i) lat_i(m)]);
                                    dsr = gsw_distance([lon(j) lon_i(m)],[lat(j) lat_i(m)]);

                                    % weights
                                    wl = dsr/(dsl+dsr);
                                    wr = dsl/(dsl+dsr);

                                    % Linear interpolation
                                    v_v(m,i_k) = (v_l*wl + v_r*wr)*mask_land_v(m,i_k);

% END LOOP M BETWEEN ENDPOINTS I AND J ON V-GRID 
                                end

                            end
                            % Start next search loop
                            i = j;
                            break;
                        else
                            % keep searching for the right active node;
                            j = j + 1;
                        end

% END SEARCH J ON T-GRID                       
                    end
                else
                    % keep searching for the left active node;
                    i = i + 1;
                end

                % no more active nodes
                if j == sz_n
    %                 display('+ End search loop, no more right node found.');
                    break;
                end

% END SEARCH I ON T-GRID  
            end

% END LOOP T-DEPTH         
        end
        %%

        % Fill the bottom triangles by
        % copying values from above
        
% LOOP V-GRID
        for i_n = 1:sz_n-1
            if mask_mode_d(i_n)== 1 || mask_mode_d(i_n+1)== 1
                
                % temporary variable
                tmp = v_v(i_n,:);
                
                % at least two data pints needed
                if sum(~isnan(tmp))>= 2

                    % fill any gaps
                    tmp(isnan(tmp)) = interp1(find(~isnan(tmp)), tmp(~isnan(tmp)), find(isnan(tmp)),'pchip',nan);
                    v_v(i_n,:) = tmp;

                    % copy to the bottom
                    ipz = find(~isnan(tmp),1,'last');
                    v_v(i_n,ipz+1:end) = v_v(i_n,ipz).*mask_land_v(i_n,ipz+1:end);
                    
                end
            end
% END LOOP V-GRID            
        end

        

        % //////////////////////////////////////////////////////////
        % // Fill the wedges in Rockall Trough /////////////////////
        % //////////////////////////////////////////////////////////
        % Directly copy measurements into the wedges
        %
        if strcmp(flag_section,'east') || strcmp(flag_section,'all')     
            
            % *********************
            % West of RTWB1
            % *********************
            % Copy the velocities from RTWB1 into the wedge
            % between 13.0323 W and RTWB1 (The number 13.0323 W is from
            % the report writtend by Chris at SAMS).

            % Find mr# for RTWB1
            
            switch flag_section
                case 'east'
                    ipo = 33;
                case 'all'
                    ipo = 57;
            end

            % Locate the wedge
            ipx = find(lon_i>= -13.0323 & lon_i<= lon_o(ipo));

            % Velocities at WB1
            tmp = -u_t(ip_obs(ipo),:).*sin(ang_i(ip_obs(ipo)))...
                +v_t(ip_obs(ipo),:).*cos(ang_i(ip_obs(ipo)));

            % Copy into the wedge
            v_v(ipx,:) = repmat(tmp,length(ipx),1)...
                .*mask_land_v(ipx,:);

            clear ipx tmp

            
            % ******************************************************
            % Between EB1 and RTADCP and below the common depths
            % ******************************************************
            % Copy the velocities from EB1 into the wedge
            % between EB1 and RTADCP and below the mximum common
            % depth of the two moorings.

            % Find indices for RTEB1 and RTADCP
            switch flag_section
                case 'east'
                    ipo1 = 35;
                    ipo2 = 36;
                case 'all'
                    ipo1 = 59;
                    ipo2 = 60;
            end

            % Locate the wedge
            % between the two moorings
            % and below the common maximum
            % depth
            ipx = find(lon_i>= lon_o(ipo1) & lon_i<= lon_o(ipo2));
            ipz = find(isnan(v_t(ip_obs(ipo2),:))); % depths below RTADCP

            % Velocities at EB1 rotated to the section
            tmp = -u_t(ip_obs(ipo1),:).*sin(ang_i(ip_obs(ipo1)))...
                +v_t(ip_obs(ipo1),:).*cos(ang_i(ip_obs(ipo1)));

            % Copy into the wedge
            v_v(ipx,ipz) = repmat(tmp(ipz),length(ipx),1)...
                .*mask_land_v(ipx,ipz);

            clear ipx ipz tmp
            
            
            % *********************
            % East of RTADCP
            % *********************
            % Copy the velocities from RTADCP into the wedge between
            % the mooring and 9.2W (this eastern limit is from
            % the report by Chris at SAMS).

            % Find mr# for RTADCP1
            switch flag_section
                case 'east'
                    ipo = 36;
                case 'all'
                    ipo = 60;
            end

            % Locate the wedge
            ipx = find(lon_i>= lon_o(ipo) & lon_i<= -9.2);

            % Velocities at RTADCP
            tmp = -u_t(ip_obs(ipo),:).*sin(ang_i(ip_obs(ipo)))...
                +v_t(ip_obs(ipo),:).*cos(ang_i(ip_obs(ipo)));

            % Copy into the wedge
            v_v(ipx,:) = repmat(tmp,length(ipx),1)...
                .*mask_land_v(ipx,:);

            clear ipx tmp
        end
        % //////////////////////////////////////////////////////////
        % //////////////////////////////////////////////////////////

        

% =========================================================================
%% FLAG_TOPO_INTERP =======================================================
% =======================================================================
% In area with the boundary current moorings and for the depths 300
% meters above the seafloor and only in the water column d>1000m;
%
% 7/24/2019: Update to include UMD5, UMM3;

% FLAG_TOPO_INTERP
       if flag_topo_interp== 1

           switch flag_section
               case 'west'
                   rg_v_interp = [1:11 13:21]; % OM2-OM13, OM14-OM23
               case 'east'
                   rg_v_interp = [1:13 15:29 33]; % OM24(CF1)-OM41(NOCM5), OM42(IC0)-OM53(UMM3), RTWB1-RTWB2
               case 'all'
                   rg_v_interp = [1:11 13:21 25:38 39:53 58]; 
           end

% LOOP MOORINGS
            for i= rg_v_interp

    %             % LOG
    %             disp(['- mooring#',num2str(i)]);

                % indices on t-grid
                l= ip_obs(i);
                r= ip_obs(i+1);
                
                % in case no data available from OOI, or RREX, which is
                % defined as the mooring on the west for intepolating;
                % continue;
                if all(isnan(v_t(l,:)))
                    continue;
                end
                
                % in case no data available from OOI, or RREX, which is
                % defined as the mooring on the east for intepolating;
                % skip;
                if all(isnan(v_t(r,:)))
                    
                    % found new mooring on the east;
                    r1 = ip_obs(i+2); % to skip 1 mooring...
                    r2 = ip_obs(i+3); % to skip 2 moorings...
                    
                    if any(~isnan(v_t(r1,:))) && r1<= ip_obs(max(rg_v_interp)+1)
                        r= r1;
                    elseif any(~isnan(v_t(r2,:))) && r2<= ip_obs(max(rg_v_interp)+1)
                        r= r2;
                    else
                        continue; % newly added
                    end
                    
                end
                

                % topography-following interp for 
                % the bottom 300 meters in the water
                % columns
% IF MODE_D SEGMENT
                if all(mask_mode_d(l:r)==1)
                    
% LOOP DEPTH FROM BOTTOM    
                    for k= 1:15 % 15 layers only

                        % search from bottom
                        if k== 1
    %                         % bottom node (mask_land may have been 
    %                         %     modified to accomendate UMM2_extended)
    %                         ibl = find(mask_land(l,:)==1,1,'last');
    %                         ibr = find(mask_land(r,:)==1,1,'last');

                            % bottom node
                            ibl = find(~isnan(v_t(l,:)),1,'last');
                            ibr = find(~isnan(v_t(r,:)),1,'last');
                        else
                            % shift up one node
                            ibl = ibl-1;
                            ibr = ibr-1;
                        end

                        % both moorings are shallower than 1000m or any of
                        % the two is out of surface
                        if depth(ibl)<=1000 && depth(ibr)<=1000% || ibl==0 || ibr==0
                            break;
                        end


                        % u/v at the left and right nodes  
                        ul = u_t(l,ibl);
                        vl = v_t(l,ibl);

                        ur = u_t(r,ibr);
                        vr = v_t(r,ibr);


    %                     % properties at the left and right nodes
    %                     tl =   t_t(l,ibl);
    %                     sl =   s_t(l,ibl);
    %                     ptl = pt_t(l,ibl);
    %                     dl =   d_t(l,ibl);
    %                     pdl = pd_t(l,ibl);
    %                     
    %                     tr =   t_t(r,ibr);
    %                     sr =   s_t(r,ibr);
    %                     ptr = pt_t(r,ibr);
    %                     dr =   d_t(r,ibr);
    %                     pdr = pd_t(r,ibr);



                        % above the shallowest instrument,
                        % go to the next pair of moorings
                        if isnan(ul) || isnan(ur) || isnan(vr) || isnan(vr)
        %                     disp('+ ...above the shallowest instrument!');
                            break;
                        end

        %                 % LOG
        %                 disp(['- ibl= ',num2str(ibl),' ibr= ',num2str(ibr)]);


    % LOOP M ON V-GRID   
                        for m= l:r-1

                            % target node
                            if k== 1
                                ib(m) = find(mask_land_v(m,:)==0,1,'first');
                            end

                            % bottom grid
                            ib(m) = ib(m) - 1;

                            % reached the sea surface
                            if ib(m) <= 0
                                break;
                            end

        %                     % LOG
        %                     disp([' ib= ',num2str(ib(m))]);


                            % horizontal distances to the left and right active 
                            % nodes on T-Grid
                            dsl = gsw_distance([lon(l) lon_i(m)],[lat(l) lat_i(m)]);
                            dsr = gsw_distance([lon(r) lon_i(m)],[lat(r) lat_i(m)]);


                            % absolute distances to the left and right active
                            % nodes
                            dsl = sqrt(dsl^2 + (depth(ib(m))-depth(ibl))^2);
                            dsr = sqrt(dsr^2 + (depth(ib(m))-depth(ibr))^2);

                            % weights
                            wl = dsr/(dsl+dsr);
                            wr = dsl/(dsl+dsr);

                            % velocity normal to the section
                            % at the left and right active nodes
                            vpl = -ul.*sin(ang_i(m))+vl.*cos(ang_i(m));
                            vpr = -ur.*sin(ang_i(m))+vr.*cos(ang_i(m));                      

        %                     disp(['vpl= ',num2str(vpl),' vpr=',num2str(vpr)]);

                            % interpolate/overwrite
                            % *** velocity ***
                            v_v(m,ib(m)) = (vpl*wl + vpr*wr)*mask_land_v(m,ib(m));

%                             % *** property *** TESTED, NOT USED
%                             t_v(m,ib(m)) = (tl*wl + tr*wr)*mask_land_v(m,ib(m));
%                             s_v(m,ib(m)) = (sl*wl + sr*wr)*mask_land_v(m,ib(m));
%                             d_v(m,ib(m)) = (dl*wl + dr*wr)*mask_land_v(m,ib(m));
%                             pt_v(m,ib(m)) = (ptl*wl + ptr*wr)*mask_land_v(m,ib(m));
%                             pd_v(m,ib(m)) = (pdl*wl + pdr*wr)*mask_land_v(m,ib(m));


% END LOOP M ON V-GRID       
                        end
% END LOOP DEPTH FROM BOTTOM        
                    end
% END IF MODE_D SEGMENT                 
                end


% END LOOP MOORINGS    
            end
% END FLAG_TOPO_INTERP        
       end

       clear ib dsl dsr wl wr vpl vpr ul ur vl vr l r m k i
       

        % Fill dry velocity nodes with zeros,
        %   for the interpolating purpose;
        % Not to the property nodes (zero could be 
        %   a meaningless number for properties;
        %
        % Since peropties *_v are not used here, no 
        %   changes need to make (commented out);
        v_v(mask_land_v==0) = 0;
%         t_v(mask_land_v==0) = nan;
%         s_v(mask_land_v==0) = nan;
%         d_v(mask_land_v==0) = nan;
%         pt_v(mask_land_v==0) = nan;
%         pd_v(mask_land_v==0) = nan;

%% =======================================================================
% =========================================================================    




%         %% PLOT
%         figure;
%         pcolor(lon_i,depth,v_v'); shading flat; axis ij;
%         title('(3.2) v\_v');
        
        
        % -----------------------------------------------------------------
        %% (3.3) INTERPOLATE PEROPERTIES AND OBTAIN VERTICAL VELOCITY SHEAR
        % -----------------------------------------------------------------
        % To interpolate properties from the adjacent current meter
        % moorings, and calculate the geostrophic velocity shears from the
        % adjacent dynamic height moorings (including the density profiles
        % in the glider domain);
        %
        % (3.3.1)
        % Search for a pair of active nodes (mask_obs_ts is 1), and
        % fill t/s/rho/sigma/theta in between in the area covered by 
        % current meter moorings. This prevents interpolation
        % for the area between adjacent dynamic height moorings.
        %
        % (3.3.2)
        % Calculate geostrophic velocity shears when mask_mode_g is 1 for both
        % active endpoints.  In the meantime, calulte the surface geostrophic 
        % velocities from ssh.
        %



        % ----------------------------------------------------------
        %% (3.3.1) Update property between current meter moorings
        % ----------------------------------------------------------
        % In the area covered by the boudnary current moorings, interpolate 
        % the property profiles to fill the gap between the moorings;
        %
    %     % LOG
    %     display('+ Interpolating t/s ...');


% LOOP T-DEPTH       
        for i_k = 1:sz_z
    %         display([' + Interpolate t/s on depth k= ',num2str(i_k)]);

            % counters
            i = 1; % Left node counter
            j = 0; % Right node counter
            o = 0; % Count valid regions (between endpoints)

% SEARCH I ON T-GRID           
            while i <= sz_n-2
                if mask_obs_ts(i,i_k)==1 && mask_mode_d(i)== 1
                    % Found the left active node;
%                     display(['+ Found ACTIVE node i= ',num2str(i),' k= ',num2str(i_k)]);
                    
% SEARCH J ON T-GRID
                    % right node counter
                    j = i + 1;
                    while j <= sz_n-1
                        if mask_obs_ts(j,i_k)==1 && mask_mode_d(j)== 1
                            % Found the right active node;
%                             display(['+ Found ACTIVE node j= ',num2str(j),' k= ',num2str(i_k)]);
                            
                            
                            % t/s at the active nodes
                            tpl =   t_t(i,i_k);
                            spl =   s_t(i,i_k);
                            rhopl = d_t(i,i_k);
                            pdpl = pd_t(i,i_k);
                            ptpl = pt_t(i,i_k);
                            
                            tpr =   t_t(j,i_k);
                            spr =   s_t(j,i_k);  
                            rhopr = d_t(j,i_k);
                            pdpr = pd_t(j,i_k);
                            ptpr = pt_t(j,i_k);

                            
                            % Fill inactive nodes between the two endpoints:
                            % (1) dry on at least one side, this is 
                            %     a point in the bottom triangles, don't 
                            %     interpolate, copy the values or keep 
                            %     climtology;
                            % (2) two sides are wet, linear interpolation for
                            % the area between adjacent mooring;
                            %
                            if mask_land(i,i_k)==0 || mask_land(j,i_k)==0 ...
                                    || isnan(tpl) || isnan(tpr)
                                % (1) wet node on one side
                                %   For the bottom triangles, either copy
                                %   values from the other side or do thing
                                %   to keep the preloaded fields.
                                % 
%                                 % ===== [A] copy t/s from the other side ====
%                                 if mask_land(i,i_k)==1 || mask_land(j,i_k)==0
%                                     for m = i:j-1
%                                         t_v(m,i_k) = tpl;
%                                         s_v(m,i_k) = spl;
%                                         d_v(m,i_k) = rhopl;
%                                         pd_v(m,i_k) = pdpl;
%                                         pt_v(m,i_k) = ptpl;
%                                     end
%                                 elseif mask_land(i,i_k)==0 || mask_land(j,i_k)==1
%                                     for m = i:j-1
%                                         t_v(m,i_k) = tpr;
%                                         s_v(m,i_k) = spr;
%                                         d_v(m,i_k) = rhopr;
%                                         pd_v(m,i_k) = pdpr;
%                                         pt_v(m,i_k) = ptpr;
%                                     end     
%                                 end
                                % ===== [B] keep background t/s ==============
                                % Do nothing to upset the pre-filled
                                % climatology fields;
                                % ============================================
                                
                            elseif mask_land(i, i_k)==1 && mask_land(j, i_k)==1
                                %
                                % (2) wet nodes on both sides
                                %

% LOOP M BETWEEN ENDPOINTS I AND J ON T-GRID
                                for m = i:j-1
    %                                 display(['+ Interpolate to INACTIVE node m= ',num2str(m),' k= ',num2str(i_k)]);

                                    % Fill the grid points where v_v
                                    %   are filled with measured values;
                                    % Make sure they are wet nodes, using
                                    %   mask_land_v;
                                    if ~isnan(v_v(m,i_k)) && mask_land_v(m,i_k)== 1
                                        % Update distances to left and right active 
                                        % nodes on the t-grid             
                                        dsl = gsw_distance([lon(i) lon_i(m)],[lat(i) lat_i(m)]);
                                        dsr = gsw_distance([lon(j) lon_i(m)],[lat(j) lat_i(m)]);

                                        % weights
                                        wl = dsr/(dsl+dsr);
                                        wr = dsl/(dsl+dsr); 

                                        % distance-weighted interpolation
                                        % NOTE
                                        %   Don't overwrite
                                        %   topography-following interp.
                                        %   values if any;
                                        if isnan(s_v(m,i_k))
                                            t_v(m,i_k) = (tpl*wl + tpr*wr);
                                            s_v(m,i_k) = (spl*wl + spr*wr);
                                            d_v(m,i_k) = (rhopl*wl + rhopr*wr);
                                            pd_v(m,i_k) = (pdpl*wl + pdpr*wr);
                                            pt_v(m,i_k) = (ptpl*wl + ptpr*wr);
                                        end

                                        % v-grid to i-grid
                                        if i_k > 1
                                            t_i(m,i_k-1) = nanmean(t_v(m,i_k-1:i_k),2)*mask_land_i(m,i_k-1);
                                            s_i(m,i_k-1) = nanmean(s_v(m,i_k-1:i_k),2)*mask_land_i(m,i_k-1);
                                            d_i(m,i_k-1) = nanmean(d_v(m,i_k-1:i_k),2)*mask_land_i(m,i_k-1);
                                            pd_i(m,i_k-1) = nanmean(pd_v(m,i_k-1:i_k),2)*mask_land_i(m,i_k-1);
                                            pt_i(m,i_k-1) = nanmean(pt_v(m,i_k-1:i_k),2)*mask_land_i(m,i_k-1);
                                        end
                                        
                                    else
                                        % do nothing to upset the prefilled
                                        % OA fields;   
                                        
                                    end

% END LOOP M BETWEEN ENDPOINTS I AND J ON T-GRID
                                end 

                            end % Fill inactive nodes

                            % Start next search loop
                            i = j;
                            break;
                        else
                            % keep searching for the right active node;
                            j = j + 1;
                        end

% END SEARCH J ON T-GRID                     
                    end
                else
                    % keep searching for the left active node;
                    i = i + 1;
                end

                % no more active nodes to the right of current left node.
                if j == sz_n
    %                 display('+ End search loop, no more right node found.');
                    break;
                end

% END SEARCH I ON T-GRID             
            end

% END LOOP T-DEPTH           
        end
        
        
        
        
%         % PLOT
%         figure;
%         pcolor(lon_i,depth,s_v'); shading flat; axis ij;
%         title('(3.3.1) s\_v');
%         
%         figure;
%         pcolor(lon_i,depth_i,s_i'); shading flat; axis ij;
%         title('(3.3.1) s\_i');
        

        
        
        %% remove invalid data points
        t_i(mask_land_i==0) = nan;
        s_i(mask_land_i==0) = nan;
        pd_i(mask_land_i==0) = nan;
        pt_i(mask_land_i==0) = nan;
        d_i(mask_land_i==0) = nan;


%         % PLOT
%         figure;
%         pcolor(lon_i,depth_i,s_i'); shading flat; axis ij;
%         title('(3.3.1) Land points nan''ed s\_i');
        


        % ----------------------------------------------
        %% (3.3.2) Compute thermal wind shears
        % ----------------------------------------------
        
% =========================================================================
% FLAG_BARO_CORR (1/4) =========================================================
% ========================================================================
% Initialize variables: surface geostrophic and barotropic velocities
%
% IF FLAG_BARO_CORR            
        if flag_baro_corr==1 

            % Initialization
            if i_t == i_t_begin  
                v_ssh = nan(sz_n-1,1);
                v_baro = nan(sz_n-1,1);
            end               

% END IF FLAG_BARO_CORR
        end
% ========================================================================
% ========================================================================= 


% =========================================================================
% FLAG_FULL_GRID_SHEARS ===================================================
% ========================================================================
% Set the T/S observation mask all to 1 for allowing all t/s 
% in the geostrophy segments to be used for calculating the 
% thermal wind shears.
%
% IF FLAG_FULL_GRID_SHEARS
        if flag_full_grid_shears == 1   
        % Apply full grid geostrophy
        % everywhere across the section
        %
% LOOP V-GRID        
            for i_n = 1:sz_n-1
% LOOP V-GRID DEPTH            
                for i_k = 1:sz_z

                    % Only in the geostrophy segments
                    if isnan(v_v(i_n,i_k))
                        %
                        mask_obs_ts(i_n:i_n+1,i_k)= 1;
                    end

% END LOOP V-GRID DEPTH                
                end
% END LOOP V-GRID            
            end

            % Set 1 on all dry nodes
            mask_obs_ts(mask_land==0) = 1;
            
            
            % Create p_geos for consistency
            p_geos = 1:sz_n-1;
            
            
            
        elseif flag_full_grid_shears == 2
        % [default] Applying full-grid geostrophy only
        % over the SAMS gldier domain;
        
            % *** East/All ***
            if strcmp(flag_section,'east') || strcmp(flag_section,'all')
                % Applying the full-grid 
                % geostrophic method between 
                % UMM4 and SAMS_glider_endpoint;
                
                % Found mr# for UMM4
                if flag_test_IB5==1 && time(i_t)>=datenum(2018,6,1)
                    switch flag_section
                        case 'east'
                            ipo2 = 32;
                        case 'all'
                            ipo2 = 56;
                    end
                    
                elseif flag_test_IB5==0 || (flag_test_IB5==1 && time(i_t)<datenum(2018,6,1))
                    switch flag_section
                        case 'east'
                            ipo2 = 31;
                        case 'all'
                            ipo2 = 55;
                    end
                end
                
                % Determine the western and eastern boundaries
                p2 = ip_obs(ipo2); % (UM)M4 or IB5
                p3 = find(lon<= -14.7,1,'last'); % SAMS Gliders' Eastern Endpoint

                % Always use SAMS glider data
                p_geos = p2:p3;
                    

                % Update the mask on T-Grid
                mask_obs_ts(p_geos,:) = 1;

                
            end% if strcmp

        elseif flag_full_grid_shears == 0
        % Apply full-grid geostrophy
        % only in the OUC and SAMS gliders domains
        %
        
            % *** East/All ***
            if strcmp(flag_section,'east') || strcmp(flag_section,'all')
                % By default, applying the full-grid 
                % geostrophic method in the following
                % areas:
                % (1) M3 to M4, for 2015-06-27 to 2016-06-05;
                % (2) M4 to SAMS_glider_endpoint;
                
                
                % Found mr# for UMM3 and UMM4
                switch flag_section
                    case 'east'
                        ipo1 = 30;
                        ipo2 = 31;
                    case 'all'
                        ipo1 = 52;
                        ipo2 = 53;
                end

                % Determine the western and eastern boundaries
                p1 = ip_obs(ipo1); % (UM)M3
                p2 = ip_obs(ipo2); % (UM)M4
                p3 = find(lon<= -14.7,1,'last'); % SAMS Gliders' Eastern Endpoint


                
                % Found indices for the geostrophic segments ...
                if time(i_t)>= datenum(2015,6,28) && ...
                        time(i_t)< datenum(2016,6,6)
                    % WHOI/OUC glider was in place 
                    % between 27-Jun-2015 and
                    % 05-Jun-2016;
                    p_geos = p1:p3;
                    
                else
                    % Always use SAMS glider data
                    p_geos = p2:p3;
                    
                    
                end
                

                % Update the mask on T-Grid
                mask_obs_ts(p_geos,:) = 1;

                
            end% if strcmp
            
        else
            error('NOW 001!');


% END IF FLAG_FULL_GRID_SHEARS
        end 
% ========================================================================
% =========================================================================    

    
    
% =========================================================================
% flag_OOI_for_tw =====================================================
% ========================================================================
% Update mask_obs_ts to avoid the use of density profiles for
% calculating the thermal wind shears
%
        if flag_OOI_for_tw == 0 && (strcmp(flag_section,'east') || strcmp(flag_section,'all'))
                        
            % *** OOI ***
            % Don't use OOI moorings for the thermal
            % wind calculations
            
            % LOG
            disp('>> OOI FLMA and FLMB are not used as dynamic height moorings!');
                    
            switch flag_section
                case 'east'
                    ipo1 = 11;
                    ipo2 = 12;                 
                    
                case 'all'
                    ipo1 = 33;
                    ipo2 = 34;
                    
            end

            % nan
            ipz = isnan(v_t(ip_obs(ipo1),:));
            mask_obs_ts(ip_obs(ipo1),ipz) = 0;

            ipz = isnan(v_t(ip_obs(ipo2),:));
            mask_obs_ts(ip_obs(ipo2),ipz) = 0;

        end 
% ========================================================================
% ========================================================================= 
    
    
    
% LOOP T-DEPTH   
        for i_k = 1:sz_z
    %         display([' + Interpolate t/s on depth k= ',num2str(i_k)]);

            % counters
            i = 1; % Left node counter
            j = 0; % Right node counter
            o = 0; % Count valid regions (between endpoints)

% SEARCH I ON T-GRID        
            while i <= sz_n-2
                if mask_obs_ts(i,i_k)==1
                    % Found the left active node;
    %                 display(['+ Found ACTIVE node i= ',num2str(i),' k= ',num2str(i_k)]);

% SEARCH J ON T-GRID
                    % right node counter
                    j = i + 1;
                    while j <= sz_n-1
                        if mask_obs_ts(j,i_k)==1
                            % Found the right active node;
    %                         display(['+ Found ACTIVE node j= ',num2str(j),' k= ',num2str(i_k)]);

                            % Density
                            rhopl = d_t(i,i_k);
                            rhopr = d_t(j,i_k);

                            % Determine if I and J are both wet nodes
                            if mask_land(i, i_k)==1 && mask_land(j, i_k)==1

                                % Dynamic topography profiles,
                                % compute geostrophy;
                                if mask_mode_g(i+1) == 1 && mask_mode_g(j-1) == 1
    %                             display(['+ Found geostrophy flag between i= ',num2str(i),' and j= ',num2str(j)]);


                                    % Find separation [m] between two density
                                    % profiles
                                    ds = gsw_distance([lon(i) lon(j)],[lat(i) lat(j)]);

                                    
                                    % Compute vertical geostrophic shear 
                                    % with dv/dz = -gr/(f rho_0)*(drho/dx)
                                    % relative to BOTTOM!
                                    %
% LOOP M BETWEEN ENDPOINTS I AND J ON V-GRID
                                    for m=i:j-1
    %                                     display(['+ Calculate velocity shears for node m= ',num2str(m),' k= ',num2str(i_k)]);
                                        g_v(m,i_k) = -1.0*(gr/(f_i(m)*rho0))*((rhopr-rhopl)/ds)*mask_land_v(m,i_k);

% END LOOP M BETWEEN ENDPOINTS I AND J ON V-GRID    
                                    end

                                end %  END if mask_mode_g
                            end %  END if mask_obs_ts

                            % Start next search loop
                            i = j;
                            break;
                        else
                            % keep searching for the right active node;
                            j = j + 1;
                        end

% END SEARCH J ON T-GRID                    
                    end
                else
                    % keep searching for the left active node;
                    i = i + 1;
                end

                % no more active nodes to the right of current left node.
                if j == sz_n
    %                 display('+ End search loop, no more right node found.');
                    break;
                end

% END SEARCH I ON T-GRID            
            end

% END LOOP T-DEPTH        
        end


%         return
% =========================================================================
% FLAG_BARO_CORR (2/4) ====================================================
% =========================================================================
% Calculating the surface geostrophic velocity
% from the sea surface slope;

% FLAG_BARO_CORR
        if flag_baro_corr == 1

            % Using the sea surface slope between the bounding 
            % density profiles that are from tall moorings or 
            % in the SAMS glider domain;
            ipx = find(mask_obs_ts(:,1)==1)';
            
%             % Using the sea surface height at the locations of 
%             % any pair of density profiles existed throughout
%             % the water column (including deep moorings);
%             mask_tmp = mask_obs_ts;
%             mask_tmp(mask_land==0) = 0;
%             ii= 1;
%             ipx = nan;
%             for i_n= 1:sz_n
%                 if any((mask_tmp(i_n,:))==1)
%                     ipx(ii) = i_n;
%                     ii = ii+1;
%                 end
%             end
%             
%             disp('Test of using Vssh at each pair of density profiles used including that from short moorings.');
            
            % Loop geostrophy segments
            for ix = 1:length(ipx)-1
                % left and right nodes
                i = ipx(ix);
                j = ipx(ix+1);

                % for geostrophic segment
                if mask_mode_g(i)== 1 && mask_mode_g(j)== 1

                    % sea surface slope
                    dh = ssh_t(j) - ssh_t(i);

                    % distance between two end points
                    ds = gsw_distance([lon(i) lon(j)],[lat(i) lat(j)]);

                    % surface geostrophic velocities
                    % derived from the sea surface slope
                    for m = i:j-1
                        v_ssh(m) =  gr/f_i(m)*dh/ds;
%                         disp(num2str(m));
                    end
                    
                    
% =========================================================================
% FLAG_TEST_AGV ===========================================================
% =====================================================================
                    if flag_test_AGV== 1
                        
                        % Integrated velocity between i and j
                        sumv = nansum(agv_t(i:j-1));
                        
                        % Distance weighted velocities
                        % for each grid between i and j
                        wi = dist_i(i:j-1)./sum(dist_i(i:j-1));
                        v_ssh(i:j-1) = wi.*sumv;
                        
                    end
% =====================================================================
% =========================================================================
                    
                end
            end
            
        
% % =========================================================================
% % FLAG_AVISO_SMOOTHING ====================================================
% % =====================================================================           
% % Spatial averaging on distances ~100km
% %   to reduce redundancy (see Gourcuff et al. 2011);  
% %
% % Only in the areas where the full-grid shears are 
% %   calculated, designated by p_geos (e.g., glider domain);
% %
%             if flag_aviso_smoothing == 1
% 
% 
%                 % segments with lengths O~100km
%                 distcum = cumsum(dist_i(p_geos)).*1e-3; % km
%                 ip_seg = nan(floor(distcum(end)/100),1);
%                 for i_s = 1:length(ip_seg)
%                     ip_seg(i_s) = find((distcum-i_s*100)>0,1,'first');
%                 end
%                 ip_seg(end+1) = length(p_geos);
% 
% 
%                 % average in each segment
%                 tmp = v_ssh(p_geos);
% 
%                 for i_s = 1:length(ip_seg)
% 
%                     if i_s == 1
%                         ipn = 1:ip_seg(1);
%                     else
%                         ipn = ip_seg(i_s-1)+1:ip_seg(i_s);
%                     end
% 
%                     % averaging
%                     tmp(ipn) = nanmean(tmp(ipn));
%                 end
% 
% 
%                 % update v_ssh
%                 v_ssh(p_geos) = tmp;
% 
% 
%                 % clean up
%                 clear agv_seg tmp ip_seg i_s i_x distcum 
%                     
%                 
% % END FLAG_AVISO_SMOOTHING
%             end
% % =====================================================================
% % =========================================================================
            
            
            
% END FLAG_BARO_CORR
        end
% =========================================================================
% =========================================================================




% =========================================================================
% FLAG_FULL_GRID_BARO_CORR ================================================
% =========================================================================
% Calculate the surface geostrophic 
%   velocity across each grid cell;

% IF FLAG_FULL_GRID_BARO_CORR
        if flag_full_grid_baro_corr == 1 && flag_baro_corr == 1
            
            % Loop geostrophy segments
            for i = 1:sz_n-1
                % left and right nodes
                j = i+1;

                % for geostrophic segment
                if mask_mode_g(i)== 1 && mask_mode_g(j)== 1

                    % sea surface slope
                    dh = ssh_t(j) - ssh_t(i);

                    % distance between two end points
                    ds = gsw_distance([lon(i) lon(j)],[lat(i) lat(j)]);

                    % surface geostrophic velocities
                    % derived from the sea surface slope
                    for m = i:j-1
                        v_ssh(m) =  gr/f_i(m)*dh/ds;
%                         disp(num2str(m));
                    end
                    
                    
% =========================================================================
% FLAG_TEST_AGV ===========================================================
% =====================================================================
                    if flag_test_AGV== 1
                        
                        % Integrated velocity between i and j
                        sumv = nansum(agv_t(i:j-1));
                        
                        % Distance weighted velocities
                        % for each grid between i and j
                        wi = dist_i(i:j-1)./sum(dist_i(i:j-1));
                        v_ssh(i:j-1) = wi.*sumv;
                        
                    end
% =====================================================================
% =========================================================================
                    
                end
            end
            
            
% =========================================================================
% FLAG_AVISO_SMOOTHING ====================================================
% =====================================================================           
% Spatial averaging on distances ~100km across the full array
%   to reduce redundancy (see Gourcuff et al. 2011);

            if flag_aviso_smoothing== 1 && flag_vbaro_mean== 0
            % only needed when v_baro_mean is not used...    

                % segments with lengths O~100km
                distcum = cumsum(dist_i).*1e-3; % km
                ip_seg = nan(floor(distcum(end)/100),1);
                for i_s = 1:length(ip_seg)
                    ip_seg(i_s) = find((distcum-i_s*100)>0,1,'first');
                end
                ip_seg(end+1) = sz_n-1;


                % average in each segment
                tmp = v_ssh;

                for i_s = 1:length(ip_seg)

                    if i_s == 1
                        ipn = 1:ip_seg(1);
                    else
                        ipn = ip_seg(i_s-1)+1:ip_seg(i_s);
                    end

                    % averaging
                    tmp(ipn) = nanmean(tmp(ipn));
                end

                
                % update v_ssh
                v_ssh = tmp;
                
                
                % clean up
                clear agv_seg tmp ip_seg i_s i_x distcum 
                
% END FLAG_AVISO_SMOOTHING
            end
% =====================================================================
% =========================================================================
        

% END FLAG_FULL_GRID_BARO_CORR        
        end
% =========================================================================
% =========================================================================    



% =========================================================================
% flag_test_glider_shear
% =========================================================================
% Set the vertical velocity shears to be zero 
% below 1500m
        if flag_test_glider_shear== 1 && ~strcmp(flag_section,'west')
            % Applying the full-grid 
                % geostrophic method between 
                % UMM4 and SAMS_glider_endpoint;
                
                % Found mr# for IB4 or IB5
                if flag_test_IB5==1 && time(i_t)>=datenum(2018,6,1)
                    switch flag_section
                        case 'east'
                            ipo2 = 32;
                        case 'all'
                            ipo2 = 56;
                    end
                     % depth range
                     ipz = find(depth>= 2000);
                elseif flag_test_IB5==0 || (flag_test_IB5==1 && time(i_t)<datenum(2018,6,1))
                    switch flag_section
                        case 'east'
                            ipo2 = 31;
                        case 'all'
                            ipo2 = 55;
                    end
                     % depth range
                     ipz = find(depth>= 2000);
                end

                % Determine the western and eastern boundaries
                p2 = ip_obs(ipo2); % (UM)M4/IB5
                p3 = find(lon<= -14.7,1,'last'); % SAMS Gliders' Eastern Endpoint

                % SAMS glider domain
                ipn = p2:p3;
                
                % Update the shears
                g_v(p2,:) = g_v(p2+1,:); % abnormal shears typically observed right next to UMM4
                g_v(ipn,ipz) = 0;
                
        end
% =========================================================================
% =========================================================================



% =========================================================================
% flag_test_SAMSgdr_geov
% =========================================================================
% In the SAMS glider domain, replacing v_ssh with geov derived from glider
% DAC (Loic et al. 2018 JGR);
        if flag_test_SAMSgdr_geov== 1 && ~strcmp(flag_section,'west')
            
            % load data
            if ~exist('meanmap','var')
                
                data_glider_geov = [dirio,'/glider/velocity/meanUKOSNAP_fordan.mat'];
                load(data_glider_geov);
                disp(['+ data_glider_geov= ',data_glider_geov]);
                
            end
            
            
            % data on the original grid
            v_gdr = meanmap.geov(2,:); % at 3m
            lon_gdr = meanmap.lon;
            lat_gdr = meanmap.lat;
            
            
            % interpolating to the I-grid
            tmp = find(lon_i>= lon_gdr(1) & lon_i<= lon_gdr(end));
            
            
%             % [A] interp accord to distance
%             lon_target(1) = lon_gdr(1);
%             lat_target(1) = lat_gdr(1);
%             lon_target(2:1+length(tmp)) = lon_i(tmp);
%             lat_target(2:1+length(tmp)) = lat_i(tmp);
           
%             % tested but not used; the new osnap line doesn't follow 58N
%             % closely because of the location of IB5;
%             dist_target = gsw_distance(lon_target, lat_target);
%             dist_target = cumsum([0 dist_target]);
%             dist_gdr = gsw_distance(lon_gdr, lat_gdr);
%             dist_gdr = cumsum([0 dist_gdr]);
            
            % [B] interp accord to longitude
            v_gdr_int = interp1(lon_gdr, v_gdr, lon_i(tmp), 'pchip');
            
            
            % replacing v_ssh
            v_ssh(tmp) = v_gdr_int;
            
            
            
        end
% =========================================================================
% =========================================================================




        % ----------------------------------------------------------------
        %% (3.4) INTEGRATE AND REFERENCE GEOSTROPHIC SHEAR TO GET VELOCITY
        % ----------------------------------------------------------------
        

% =========================================================================
% FLAG_BARO_CORR (3/4) ====================================================
% =======================================================================
% Before referencing the shears, set up 
% a mask to indicate the grid cells where the 
% SSH-derived barotropic velocities will be
% needed;
        if flag_baro_corr == 1
            
            % grid cells w/o velocity values
            mask_v_baro = v_v;
            mask_v_baro(~isnan(mask_v_baro)) = 0;
            mask_v_baro(isnan(mask_v_baro)) = 1;
            
            % remove invalid land points (may not need)
            mask_v_baro(mask_land_v==0) = 0;
            
        end
% =======================================================================
% =========================================================================


%         % PLOT
%         figure;
%         pcolor(lon_i,depth,mask_v_baro'); shading flat; colorbar; axis ij;
%         title('(3.4) v\_v derived mask\_v\_baro');
        
        
% =========================================================================
% flag_keep_vrefshortmoorings_all =========================================
% =======================================================================
% Updating mask_v_baro so that 
% in the areas with deep
% moorings, the shears referend to velocities 
% from the deep moorings won't be overwritten
% later;
        if flag_keep_vrefshortmoorings_all== 1 && ...
                flag_keep_vrefshortmoorings_selected== 0 && ...
                flag_baro_corr== 1 ...
       
            for i_n = 1:sz_n-1
                if mask_mode_d(i_n)==1 && mask_mode_d(i_n+1)==1
                    mask_v_baro(i_n,:) = 0;
                end
            end
        end
% =======================================================================
% =========================================================================


% =========================================================================
% flag_keep_vrefshortmoorings_selected ====================================
% =======================================================================
% Updating mask_v_baro to 0 to keep the 
% absolute geostrophic velocity referenced to the 
% deep moorings - avoiding the barotropic correction;
%
        if flag_keep_vrefshortmoorings_selected == 1 && flag_baro_corr == 1
            
            % Find mr# to specifiy the western and
            % eastern boundaries of an area for
            % keeping the deep mooring reference and then
            % update mask_v_baro to zero in this area
 
%             % *** K10 to DSOW5 ***
%             ipo1= -50.2548; 
%             ipo2= -49.7807; 
% 
%             % Set zeros in the mask
%             mask_v_baro(lon_i >= ipo1 & lon_i <= ipo2,:) = 0;


            % *** DSOW3 to LS8 ***
            ipo1= -47.5645; 
            ipo2= -47.3337;

            % Set zeros in the mask
            mask_v_baro(lon_i>= ipo1 & lon_i<= ipo2,:) = 0;


            % *** NOCM1 to NOCM5 ***
            ipo1= -41.1118; 
            ipo2= -37.7995; 
            if flag_test_remove_NOCM5==1 || time(i_t)>=datenum(2018,6,1)
            ipo2= -38.5665;
            end

            % Set zeros in the mask
            mask_v_baro(lon_i>= ipo1 & lon_i<= ipo2,:) = 0;


            % *** IC0 to IC1 ***
            ipo1= -35.1252; 
            ipo2= -33.6867; 

            % Set zeros in the mask
           mask_v_baro(lon_i>= ipo1 & lon_i<= ipo2,:) = 0;


            % *** UMM1 to UMM2 ***
            ipo1= -30.5293; % UMM1
            ipo2= -28.0200; % UMM2

            % Set zeros in the mask
            mask_v_baro(lon_i>= ipo1 & lon_i<= ipo2,:) = 0;


%             % *** UMM2 to UMM3 *** 
%             % disabled on 9/24/2019
%             % enabled on 7/22/2019
%             ipo1= -28.0200;
%             ipo2= -24.4220;
% 
%             % Set zeros in the mask
%             mask_v_baro(lon_i>= ipo1 & lon_i<= ipo2,:) = 0;



            % *** RTWB1 to RTWB2 ***
            ipo1= -12.7108; 
            ipo2= -12.3122;

            % Set zeros in the mask
            mask_v_baro(lon_i>= ipo1 & lon_i<= ipo2,:) = 0;

        end
% =======================================================================
% =========================================================================



%         % PLOT
%         figure;
%         pcolor(lon_i,depth,mask_v_baro'); shading flat; colorbar; axis ij;
%         title('(3.4) Updated mask\_v\_baro');

        

        % ----------------------------------------------------------
        %% (3.4.1) Reference the shears to the bottom
        % ----------------------------------------------------------
        % The refernece velocity is either zero at the sea floor 
        % or the velocity at the top of deep moorings;
        
    %     % LOG
    %     display('+ Integrating to get the absolute geostrophic velocities...');


% LOOP V-GRID
        for i_n = 1:sz_n-1

% IN G SEGMENTS        
             % For the geostrophy segments only
            if mask_mode_g(i_n)==1 && mask_mode_g(i_n+1)==1

    %             % LOG
    %             display(num2str(i_n));

                % dv/dz profile on the i-grid depth levels
                g_prof = 0.5*(g_v(i_n,2:end)+g_v(i_n,1:end-1));
                
                % fill dv/dz to be zero in the bottom triangles
                if any(~isnan(g_prof))
                    g_prof(isnan(g_prof)) = 0;
                else
                    continue;
                end

                % Search for the reference velocity
                i_ref = find(~isnan(v_v(i_n,:)),1,'first');     % bottom velocity
                i_seafloor = find(mask_land_v(i_n,:)==0,1,'first');  % the seafloor

    %             % LOG
                %	display(['i_ref= ',num2str(i_ref),' i_seafloor= ',num2str(i_seafloor),' v_v= ',num2str(v_v(i_n,i_ref))]);      


                % Update Ekman flux flag to 1 in the geostrophy segments
                %   this adds Ekman velocity onto the Ekman layers for all
                %   geostrophic sections;
                % 
                mask_ek_i(i_n) = 1;


                % Integrate the velocity shears with depth to obtain 
                %   the absolute geostrophic velocity profile- upward 
                %   to the seasurface and downward to the seafloor;
                % Do not overwrite the directly measured velocities;

                % upward integration
                for i_k = i_ref-1:-1:1
                    if isnan(v_v(i_n,i_k)) && mask_land_v(i_n,i_k)== 1

    %                     % LOG
    %                     display([' for i_k= ',num2str(i_k),' upward interp...']);
    
    
% =========================================================================
% FLAG_VREF_SHORTMOORINGS =================================================
% =======================================================================
% If this flag is zero, force to reference 
% the shears to a bottom level of 
% no motion across the section, i.e., 
% at the sea floor AND 
% at the top of deep moorings;
                        if flag_vref_shortmoorings == 0 && i_k == i_ref-1
                            
%                             % LOG
%                             display([' for i_n= ',num2str(i_n),' vref= 0!']);
                            
                            % assuming a bottom level of no motion
                            v_v(i_n,i_k) = g_prof(i_k)*(depth(i_k+1)-depth(i_k));   
                            
                            
                            % update flux flag here
                            mask_flux_i(i_n,i_k) = 1;
                            
                            
                            % next level
                            continue;
                            
                        end
% =======================================================================
% =========================================================================



                        % calculate geostrophic velocity
                        v_v(i_n,i_k) = v_v(i_n,i_k+1) + g_prof(i_k)*(depth(i_k+1)-depth(i_k));   

                        % update flux flag here
                        mask_flux_i(i_n,i_k) = 1;
                        
                        
                    else
                        % do nothing to keep observations;
                    end
                end

                
%                 % downward integration
%                 for i_k = i_ref+1:i_seafloor-1
%                     if isnan(v_v(i_n,i_k))
% 
%     %                     % LOG
%     %                     display([' for i_k= ',num2str(i_k),' DOWNWARD interp...']);
% 
%                         % get absolute geostrophic velocity
%                         v_v(i_n,i_k) = v_v(i_n,i_k-1) - g_prof(i_k-1)*(depth(i_k)-depth(i_k-1)); 
% 
%                         % update flux flag here
%                         mask_flux_i(i_n,i_k) = 1;
%                     else
%                         % do nothing to keep observations;
%                     end
%                 end
                
% END IN G SEGMENTS
            end

% END LOOP V-GRID        
        end  

           
%         % PLOT
%         figure;
%         pcolor(lon_i,depth_i,mask_flux_i'); shading flat;
%         axis ij; colorbar;
%         title('(3.4.1) mask\_flux\_i');
        

        % ----------------------------------------------------------
        %% (3.4.2) Apply the barotropic correction
        % ----------------------------------------------------------
        % Correct all the referenced geostrophic profiles by adding a
        % bartropic component obtained by differencing the sea-surface
        % slope derived surface geostrophic velocity and the referenced
        % geostrophic velocity at the sea surface - even if with a BLONM;
        %

%         return
% =========================================================================
% flag_baro_corr (4/4) =========================================================
% =======================================================================
% Calculate the barotropic velocity and 
% apply to all geostropic segments
%
% flag_baro_corr
        if flag_baro_corr== 1
            % LOG
        %     display('+ Incorporating the barotropic velocity ...');



                    % Add barotropic velocity 
                    % obtained by the difference between 
                    % the surface geostrophic velocity
                    % derived from the sea surface height
                    % and the one derived from the thermal 
                    % wind referenced to a velocity at bottom
                    % 

% LOOP I-GRID
                    for i_n= 1:sz_n-1
                        % geostrophy segments only
                        if mask_mode_g(i_n)==1 && mask_mode_g(i_n+1)==1
                            % The surface geo. velo from the sea-surface
                            % slope was calcualted when getting the
                            % velocity shears (3.3.2).


                            % Calcualte the barotropic velocity
                            if isnan(v_ssh(i_n))
                                % no data, so
                                % no barotropic velocity to be applied
                %                 disp([' i_n= ',num2str(i_n)]);
                                v_baro(i_n) = 0;
                            else
                                % barotropic velocity
                                v_baro(i_n) = v_ssh(i_n) - v_v(i_n,1);
                            end

% continue

% =========================================================================
% flag_vbaro_mean (2/2) =========================================
% =======================================================================
% Use the pre-prepared mean barotropic velocties
% to replace the calculated time-varying v_baro

                            if flag_vbaro_mean == 1
% IF FLAG_MC
                                if flag_MC == 1
                                    
                                    if flag_MC_SE == 1
                                        % drawing from mean w/ SE
                                        v_baro(i_n) = v_baro_mean(i_n,1) + v_baro_mean(i_n,3).*randn(1);
                                        
                                    elseif flag_MC_SE == 0
                                        % drawing from mean w/ SD
                                        v_baro(i_n) = v_baro_mean(i_n,1) + v_baro_mean(i_n,2).*randn(1);
                                    end
                                    
                                elseif flag_MC == 0
                                    v_baro(i_n) = v_baro_mean(i_n,1);
                                    
% END IF FLAG_MC          
                                end
                                
                            end
% =======================================================================
% =========================================================================  



                            % Add the barotropic velocity 
                            % to the referenced geostrophic 
                            % profiles
                            rg = mask_v_baro(i_n,:)==1;
                            v_v(i_n,rg) = v_v(i_n,rg) + v_baro(i_n);

                        end % END CHECK MODE G

% END LOOP I-GRID
                    end

                    
                % remove invalid data points;
                v_v(mask_land_v==0) = 0;


% END flag_baro_corr    
        end
        
% =======================================================================
% =========================================================================




        % ////////////////////////////////////////////////////////////
        % // Changes to velocity field ///////////////////////////////
        % ////////////////////////////////////////////////////////////
        % Above the Rockall Plateau,
        % between SAMS glider's eastern 
        % endpoint and west of the RT wedges
        if strcmp(flag_section,'east') || strcmp(flag_section,'all')
            
            % Find mr# for RTWB1
            ipo1 = -12.7108;
            ipx = find(lon>=-14.7 & lon<= ipo1);
            
            
            % Set velocities to zero
            for i_n = ipx
                ipz = find(isnan(v_v(i_n,:)));
                v_v(i_n,ipz) = 0;
                
                % Update the flux mask
                mask_flux_i(i_n,ipz) = 1;
            end
        
        end
        % ////////////////////////////////////////////////////////////
        % ////////////////////////////////////////////////////////////
        
        

        % ------------------------
        %% (3.5) FILL BOUNDARIES
        % ------------------------
%     % LOG
%     display('+ Filling the inshore unmeasured areas with clim t/s/u/v...');

        % Fining boundaries
        switch flag_section
            case {'west', 'east'}
                ip_b = find(lon_i< lon_o(1) | lon_i> lon_o(end));
            case 'all' 
                % west of OM2, between OM65 and OM24, and east of OM58
                % updated on 20220125 
                if time(i_t)<datenum(2018,9,1) || flag_test_remove_LSAB==1
                    ip_b = find(lon_i< lon_o(1) |...
                    (lon_i> lon_o(22) & lon_i< lon_o(25)) |...
                    lon_i> lon_o(end));
                elseif flag_test_remove_LSAB==0&&time(i_t)>=datenum(2018,9,1)
                    ip_b = find(lon_i< lon_o(1) |...
                        (lon_i> lon_o(24) & lon_i< lon_o(25)) |...
                        lon_i> lon_o(end));
                end
        end
        
        
        % By default, keep original angles to
        % rotate the velocities
        ang_i_tmp = ang_i;
        
        
% =========================================================================
% FLAG_OPA_CLIM (2/2) =====================================================
% =======================================================================
% Fill the inshore unmeasured LC componenet
% with OPA model velocity climatology;
        if flag_opa_clim== 1
            
            % Loading LC climatology provided by Brad
            if strcmp(flag_section,'west') || strcmp(flag_section,'all')
                
                % found the month
                im = str2double(datestr(time(i_t)+flag_data_interval_4calculation/2,'mm'));

                % replace FLAME climatology above the Labrador shelf
                u_t_c(1:21,:) = u_LC_model_clim(:,:,im);
                v_t_c(1:21,:) = v_LC_model_clim(:,:,im);

                % fill uv dry nodes with zeros
                u_t_c(isnan(u_t_c)) = 0;
                v_t_c(isnan(v_t_c)) = 0;
                
                u_t_c(mask_land==0) = 0;
                v_t_c(mask_land==0) = 0;

            end
        end
% =======================================================================
% =========================================================================




% =========================================================================
% FLAG_TEST_LC_UV_CLIM ====================================================
% =========================================================================
% Fill the inshore unmeasured LC componenet
% with  model velocity climatology;
        
        if flag_test_LC_UV_clim ~=0 
            
            % Loading LC UV climatology 
            switch flag_test_LC_UV_clim
            
                case 1
                    data_LC_UV_clim = [dirio,'data/climatology/0115_NEMOnew_2014-2017clim_mm_OSNAP_West.mat'];
                    load(data_LC_UV_clim);
                    u_LC_model_clim = nemo_clim.u(1:21,:,:);
                    v_LC_model_clim = nemo_clim.v(1:21,:,:);
                    
                case 2
                    data_LC_UV_clim = [dirio,'data/climatology/0115_GLORYS12v1_2014-2018clim_mm_LabSeaShelf.mat'];
                    load(data_LC_UV_clim);
                    u_LC_model_clim = glorys_clim.u(1:21,:,:);
                    v_LC_model_clim = glorys_clim.v(1:21,:,:);
                    
                case 3 
                    data_LC_UV_clim = [dirio,'data/climatology/0206_VITALS_2014-2018clim_mm_LabSeaShelf.mat'];
                    load(data_LC_UV_clim);
                    u_LC_model_clim = vitals_clim.u(1:21,:,:);
                    v_LC_model_clim = vitals_clim.v(1:21,:,:);
                    
                case 4
                    data_LC_UV_clim = [dirio,'data/climatology/0206_Multimodel_clim_mm_LabSeaShelf.mat'];
                    load(data_LC_UV_clim);
                    u_LC_model_clim = modelavg_clim.u(1:21,:,:);
                    v_LC_model_clim = modelavg_clim.v(1:21,:,:);
                    
                otherwise
                    error('IQU 187!');
                    
            end
            
            % display
            if i_t == i_t_begin && cnt_mc == 1
                disp(['+ data_sect_UV for Labrador Current= ', data_LC_UV_clim]);
            end
            
            
            % Filling boundary grids over the Labrador shelf
            if strcmp(flag_section,'west') || strcmp(flag_section,'all')
                
                % found the month
                im = str2double(datestr(time(i_t)+flag_data_interval_4calculation/2,'mm'));

                % replace climatology above the Labrador shelf
                u_t_c(1:21,:) = u_LC_model_clim(1:21,:,im);
                v_t_c(1:21,:) = v_LC_model_clim(1:21,:,im);

                % fill uv dry nodes with zeros
                u_t_c(isnan(u_t_c)) = 0;
                v_t_c(isnan(v_t_c)) = 0;
                
                u_t_c(mask_land==0) = 0;
                v_t_c(mask_land==0) = 0;

            end
        end
% =========================================================================
% =========================================================================







% =========================================================================
% FLAG_LC_UV (2/2) ========================================================
% =======================================================================
% Fill the inshore unmeasured LC componenet
% with model time-varying output

        if flag_LC_UV == 1
        % with NEMO velocity output;
        % NOTE
        %   NEMO output only available till 12/26/2017
            
            if strcmp(flag_section,'west') || strcmp(flag_section,'all')
                
                % found the time period and average
                if time(i_t)< datenum(2018,1,1)
                    ipt = find(time_LC_uv>= time(i_t) & time_LC_uv< time(i_t) + flag_data_interval_4calculation);
                else
                    % recycle the 2017 model data for 2018
                    ipt = find((time_LC_uv+365)>= time(i_t) & (time_LC_uv+365)< time(i_t) + flag_data_interval_4calculation);
                end

%                 % LOG
%                 if cnt_mc== 1
%                     disp(['Found time_LC_uv= ',datestr(time_LC_uv(ipt(1))), ' to ',datestr(time_LC_uv(ipt(end)))]);
%                 end
                
                % replace FLAME climatology above the Labrador shelf
                % average 
                u_tmp = nanmean(u_LC_model(1:21,:,ipt),3);
                v_tmp = nanmean(v_LC_model(1:21,:,ipt),3);

                % overwrite background
                u_t_c(1:21,:) = u_tmp;
                v_t_c(1:21,:) = v_tmp;

                % fill uv dry nodes with zeros, for the subsequent
                % averaging purpose...
                u_t_c(isnan(u_t_c)) = 0;
                v_t_c(isnan(v_t_c)) = 0;
                
                u_t_c(mask_land==0) = 0;
                v_t_c(mask_land==0) = 0;

            end
            
            
        elseif flag_LC_UV == 2
        % Fill monthly output from GloSea5
        
            if strcmp(flag_section,'west') || strcmp(flag_section,'all')
                
                % found the time period and average
                ipt = find(time_LC_uv<= time(i_t)+flag_data_interval_4calculation/2,1,'last');
                
                
                % replace FLAME climatology above the Labrador shelf
                % average 
                u_tmp = u_LC_model(1:21,:,ipt);
                v_tmp = v_LC_model(1:21,:,ipt);

                % overwrite background
                u_t_c(1:21,:) = u_tmp;
                v_t_c(1:21,:) = v_tmp;

                % fill uv dry nodes with zeros, for the subsequent
                % averaging purpose...
                u_t_c(isnan(u_t_c)) = 0;
                v_t_c(isnan(v_t_c)) = 0;
                
                u_t_c(mask_land==0) = 0;
                v_t_c(mask_land==0) = 0;
               
                
                % reset new angles above the shelf
                % that is, velocities need not to be rotated!
                ang_i_tmp = ang_i;
                ang_i_tmp(1:20) = 0;
            end
        
        end
% =======================================================================
% =========================================================================


% =========================================================================
% FLAG_LC_TS (2/2) ========================================================
% =======================================================================
% Fill the inshore unmeasured LC componenet TS
       if flag_LC_TS == 1
            
            if strcmp(flag_section,'west') || strcmp(flag_section,'all')
                
                % found the time period and average
                ipt = find(time_LC_ts<= time(i_t)+flag_data_interval_4calculation/2,1,'last');
                
%                 % LOG
%                 if cnt_mc== 1
%                     disp(['Found time_LC_ts= ',datestr(time_LC_ts(ipt))]);
%                 end
                
                % overwrite background
                t_t_c(1:21,:) = t_LC_rea(:,:,ipt);
                s_t_c(1:21,:) = s_LC_rea(:,:,ipt);

            end
        end
% =======================================================================
% =========================================================================

% return
% LOOP BOUNDARY GRID
        for i_x = 1:length(ip_b)
            
            % loop
            i_n = ip_b(i_x);
  

%             display(['+ Working on boundary point i_n= ',num2str(i_n)]);

            % *** u/v ***
            % Don't overwrite any calculated velocities, like in the
            % wedge east of RTADCP2
            
            % U/V on the adjacent nodes
            % Left node u/v
            upl = u_t_c(i_n,:);
            vpl = v_t_c(i_n,:);
            v_l = -upl.*sin(ang_i_tmp(i_n))+vpl.*cos(ang_i_tmp(i_n));
            % Right node u/v
            upr = u_t_c(i_n+1,:);
            vpr = v_t_c(i_n+1,:);
            v_r = -upr.*sin(ang_i_tmp(i_n))+vpr.*cos(ang_i_tmp(i_n));

            
            % Check if velocities were calculated
            if any(isnan(v_v(i_n,:)))
                % Update v_v
                v_v(i_n,:) = 0.5*(v_l + v_r).*mask_land_v(i_n,:);
                
                % update flux flag here
                mask_flux_i(i_n,:) = 1.*mask_land_i(i_n,:);
            end


            
            % *** property ***
            % Overwrite pre-loaded analysis fields with the model
            % climatology
            
            % Update t_i
            t_v(i_n,:) = nanmean(t_t_c(i_n:i_n+1,:),1);
            for i_z = 1:sz_z-1
                t_i(i_n,i_z) = nanmean(t_v(i_n,i_z:i_z+1),2);
            end
            
            % Update s_i
            s_v(i_n,:) = nanmean(s_t_c(i_n:i_n+1,:),1);
            for i_z = 1:sz_z-1
                s_i(i_n,i_z) = nanmean(s_v(i_n,i_z:i_z+1),2);
            end

            % Update d_i
            d_v(i_n,:) = nanmean(d_t_c(i_n:i_n+1,:),1);
            for i_z = 1:sz_z-1
                d_i(i_n,i_z) = nanmean(d_v(i_n,i_z:i_z+1),2);
            end

            % Update pd_i
            pd_v(i_n,:) = nanmean(pd_t_c(i_n:i_n+1,:),1);
            for i_z = 1:sz_z-1
                pd_i(i_n,i_z) = nanmean(pd_v(i_n,i_z:i_z+1),2);
            end

            % Update pt_i
            pt_v(i_n,:) = nanmean(pt_t_c(i_n:i_n+1,:),1);
            for i_z = 1:sz_z-1
                pt_i(i_n,i_z) = nanmean(pt_v(i_n,i_z:i_z+1),2);
            end
        
            
% END LOOP BOUNDARY GRID          
        end

        
        
% ===============================================================
% TEMPORARY - Fill model climatology u/v
% ===============================================================
% % If not all OSNAP data are available, using
% % FLAME model climatology instead;
% 
% % LOOP V/I-GRID
%         for i_n = 1:sz_n-1
%             
%             if isnan(v_v(i_n,1))
%                 
%                 % Update flux masks
%                 mask_ek_i(i_n) = 0;
%                 
%                 % LOG
%                 if cnt_mc== 1
%                     disp(['NOTE u/v_clim_model used at i_n= ',num2str(i_n),', lon_i= ',num2str(lon(i_n))]);
%                 end
%                 
% 
%                 % U/V on the adjacent nodes
%                 % Left node u/v
%                 upl = u_t_c(i_n,:);
%                 vpl = v_t_c(i_n,:);
%                 v_l = -upl.*sin(ang_i(i_n))+vpl.*cos(ang_i(i_n));
%                 % Right node u/v
%                 upr = u_t_c(i_n+1,:);
%                 vpr = v_t_c(i_n+1,:);
%                 v_r = -upr.*sin(ang_i(i_n))+vpr.*cos(ang_i(i_n));
% 
%                 % Interpolate to V-grid
%                 v_tmp = 0.5*(v_l + v_r);
%                 
%                 % Fill nans
%                 ipz = find(isnan(v_v(i_n,:)));
%                 v_v(i_n,ipz) = v_tmp(ipz).*mask_land_v(i_n,ipz);
%             end
% % END LOOP V/I-GRID            
%         end
        
        
% Modified on 4/1/2019: Don't fill with model uv clim; simply skip this
% time interval;        
        if any(isnan(v_v(:)))
            disp('NOTE   OSNAP data are not available across the whole section -- VOID this time step (v_v=nan);');
            v_v = nan(sz_n-1,sz_z);
        end
% ===============================================================
% ===============================================================

%         save ./temp/v_v2 v_v
%         return

% =========================================================================
% flag_test_model_clim_case2 ================================================
% =======================================================================
% Use property climatololgy in the 
% areas with direct measurements
        if flag_test_model_clim_case2 == 1
            
        % LOOP V/I-GRID
                for i_n = 1:sz_n-1

                    if mask_mode_d(i_n)== 0 || mask_mode_d(i_n+1)== 0
                        % *** t/s/d/pd/pt ***
                        t_v(i_n,:) = nanmean(t_t_c(i_n:i_n+1,:),1);
                        for i_z = 1:sz_z-1
                            t_i(i_n,i_z) = nanmean(t_v(i_n,i_z:i_z+1),2);
                        end

                        d_v(i_n,:) = nanmean(d_t_c(i_n:i_n+1,:),1);
                        for i_z = 1:sz_z-1
                            d_i(i_n,i_z) = nanmean(d_v(i_n,i_z:i_z+1),2);
                        end

                        pd_v(i_n,:) = nanmean(pd_t_c(i_n:i_n+1,:),1);
                        for i_z = 1:sz_z-1
                            pd_i(i_n,i_z) = nanmean(pd_v(i_n,i_z:i_z+1),2);
                        end

                        pt_v(i_n,:) = nanmean(pt_t_c(i_n:i_n+1,:),1);
                        for i_z = 1:sz_z-1
                            pt_i(i_n,i_z) = nanmean(pt_v(i_n,i_z:i_z+1),2);
                        end

                    end

        % END LOOP V/I-GRID            
                end
        end
% =======================================================================
% =========================================================================



% =========================================================================
% flag_test_Vgeo_only_in_CMarray (2/2) ====================================
% =======================================================================
% Apply topography-following interpolation with the geostrophic velocities
% if Vgeo happens to be over the slopes;
%
% In area with the boundary current moorings and for the depths 300
% meters above the seafloor and only in the water column d>1000m;
%
% NOTE on the grid numbering
%   T-grid# [m:n] - original T-grid cells
%   V-grid# [m:n-1] - corresponding V-grid cells between those T-grid cells
%   T-grid# [m+1:n-1] - T-grid cells interpolated from above V-grid cells
%   V-gird# [m+1:n-2] - V-grid cells interpolated from above T-grid cells
%

% IF FLAG
        if flag_test_Vgeo_only_in_CMarray== 1

            % Checking the total number of segments
            if exist('ipo0','var')
                sz_tmp = size(ipo0);
            else
                error('IPO 918!');
            end
                
           
% LOOP TEST SEGMENT 
            for kk = 1:sz_tmp(1)
               
                % Indices
                rg_v_interp = ip_obs(ipo0(kk,1)):ip_obs(ipo0(kk,2));
                
               
                % ====================================
                % Interpolating from V-Grid to T-Grid
                % ====================================
                % Temporary velocity field on T-grid
                v_t_tmp = zeros(size(v_t));
            
% LOOP V-GRID
                for ii= rg_v_interp(1:end-2)

            %             % LOG
            %             disp(['- mooring#',num2str(i)]);

                    % indices on t-grid
                    l= ii;
                    r= ii+1;


                    % topography-following interp for 
                    % the bottom 300 meters in the water
                    % columns
% LOOP DEPTH FROM BOTTOM    
                    for k= 1:15 % 15 layers only

                        % search from bottom
                        if k== 1
                            % bottom node
                            ibl = find(mask_land_v(l,:)==1,1,'last');
                            ibr = find(mask_land_v(r,:)==1,1,'last');
                        else
                            % shift up one node
                            ibl = ibl-1;
                            ibr = ibr-1;
                        end

                        % both profiles are shallower than 1000m
                        if depth(ibl)<=1000 && depth(ibr)<=1000
                            break;
                        end


                        % cross-sectional velocity 
                        % at the left and right nodes  
                        vpl = v_v(l,ibl);
                        vpr = v_v(r,ibr);


            % LOOP M ON T-GRID   
                        for m= l+1:r % T-grid centers are at the left faces of the V-grid cells

                            % target node on T-Grid
                            if k== 1
                                ib(m) = find(mask_land(m,:)==0,1,'first');
                            end

                            % bottom grid
                            ib(m) = ib(m) - 1;

                            % reached the sea surface
                            if ib(m) <= 0
                                break;
                            end


                            % horizontal distances to the left and right active 
                            % nodes on V-Grid
                            dsl = gsw_distance([lon_i(l) lon(m)],[lat_i(l) lat(m)]);
                            dsr = gsw_distance([lon_i(r) lon(m)],[lat_i(r) lat(m)]);


                            % absolute distances to the left and right active
                            % nodes
                            dsl = sqrt(dsl^2 + (depth(ib(m))-depth(ibl))^2);
                            dsr = sqrt(dsr^2 + (depth(ib(m))-depth(ibr))^2);

                            % weights
                            wl = dsr/(dsl+dsr);
                            wr = dsl/(dsl+dsr);


                            % interpolate/overwrite
                            v_t_tmp(m,ib(m)) = (vpl*wl + vpr*wr) * mask_land(m,ib(m));

% END LOOP M ON T-GRID       
                        end
% END LOOP DEPTH FROM BOTTOM        
                    end

% END LOOP V-GRID 
                end



                % ========================================
                % Interpolating from T-Grid back to V-Grid
                % ========================================
                % Overwrite values on V-Grid by applying
                % another round of topography-following 
                % interpolation ...

% LOOP T-GRID
                for ii= rg_v_interp(2:end-2)

                    % indices on t-grid
                    l= ii;
                    r= ii+1;


                    % topography-following interp for 
                    % the bottom 300 meters in the water
                    % columns
% LOOP DEPTH FROM BOTTOM    
                    for k= 1:15 % 15 layers only

                        % search from bottom
                        if k== 1
                            % bottom node
                            ibl = find(mask_land(l,:)==1,1,'last');
                            ibr = find(mask_land(r,:)==1,1,'last');
                        else
                            % shift up one node
                            ibl = ibl-1;
                            ibr = ibr-1;
                        end

                        % both profiles are shallower than 1000m
                        if depth(ibl)<=1000 && depth(ibr)<=1000
                            break;
                        end


                        % cross-sectional velocity at the left and right nodes  
                        vpl = v_t_tmp(l,ibl);
                        vpr = v_t_tmp(r,ibr);


% LOOP M ON V-GRID   
                        for m= l:r-1

                            % target node on V-Grid
                            if k== 1
                                ib(m) = find(mask_land_v(m,:)==0,1,'first');
                            end

                            % bottom grid
                            ib(m) = ib(m) - 1;

                            % reached the sea surface
                            if ib(m) <= 0
                                break;
                            end


                            % horizontal distances to the left and right active 
                            % nodes on V-Grid
                            dsl = gsw_distance([lon(l) lon_i(m)],[lat(l) lat_i(m)]);
                            dsr = gsw_distance([lon(r) lon_i(m)],[lat(r) lat_i(m)]);


                            % absolute distances to the left and right active
                            % nodes
                            dsl = sqrt(dsl^2 + (depth(ib(m))-depth(ibl))^2);
                            dsr = sqrt(dsr^2 + (depth(ib(m))-depth(ibr))^2);


                            % weights
                            wl = dsr/(dsl+dsr);
                            wr = dsl/(dsl+dsr);


                            % interpolate/overwrite
                            v_v(m,ib(m)) = (vpl*wl + vpr*wr) * mask_land_v(m,ib(m));

% END LOOP M ON V-GRID       
                        end
% END LOOP DEPTH FROM BOTTOM        
                    end

% END LOOP T-GRID 
                end


                % Cleaning up ...
                clear v_t_tmp
                clear ib dsl dsr wl wr vpl vpr ul ur vl vr l r m k i upl upr
          
               
% END TEST SEGMENT
            end
          
           
% END FLAG        
        end

% =======================================================================
% =========================================================================



        % ----------------------------------
        %% (3.6) FINISH DATA ON THE I-GRID
        % ----------------------------------
        % Fill land points with zero
        v_v(mask_land_v==0) = 0;
        
        
        % Convert velocity from V-grid to I-grid
        v_i = 0.5*(v_v(:,1:end-1) + v_v(:,2:end)).*mask_land_i;

        
        % seawater heat capacity at 1 standard atmosphere pressure
        pres = gsw_p_from_z(-repmat(depth_i,sz_n-1,1),lat_i);
        pres(mask_land_i==0) = nan;
        [SA,~,~] = gsw_SA_Sstar_from_SP(s_i,pres,lon_i,lat_i); % absolute salinity
        swcp = gsw_cp_t_exact(SA,t_i,0); % isobaric heat capacity 
        
        
        
% =========================================================================
% flag_test_pden ==========================================================
% =======================================================================
        if flag_test_pden == 1
            gamma_n = eos80_legacy_gamma_n(s_i,t_i,pres,mean(lon_i),mean(lat_i)) + 1000;
%             gamma_n(gamma_n< 1020 | gamma_n> 1030) = nan; % create 'gaps'!
            gamma_n(mask_land_i==0) = nan;
            pd_i = gamma_n;
            
        elseif flag_test_pden == 2
            pd_i = gsw_pot_rho_t_exact(SA,t_i,pres,2000); % potential density

        end
% =========================================================================
% =========================================================================

        
        % Make sure no invalid data points...
        v_i(mask_land_i==0) = nan;
        t_i(mask_land_i==0) = nan;
        s_i(mask_land_i==0) = nan;
        d_i(mask_land_i==0) = nan;        
        pd_i(mask_land_i==0) = nan;
        pt_i(mask_land_i==0) = nan;

       
        % Just in case there are nans at bottom pden grids
        [xi,yi] = find(isnan(pd_i) & ~isnan(v_i));
        if ~isempty(xi)
            
%             % LOG
%             disp(['- Found',num2str(length(xi)),' empty density bins...']);
            
            for jj= 1:length(xi)
                pd_i(xi(jj),yi(jj)) = pd_i(xi(jj),yi(jj)-1); % copy from above
                d_i(xi(jj),yi(jj)) = d_i(xi(jj),yi(jj)-1); % copy from above
                t_i(xi(jj),yi(jj)) = t_i(xi(jj),yi(jj)-1); % copy from above
                s_i(xi(jj),yi(jj)) = s_i(xi(jj),yi(jj)-1); % copy from above
                pt_i(xi(jj),yi(jj)) = pt_i(xi(jj),yi(jj)-1); % copy from above
            end
        end
        
        
        % :::::::::::::::::::::::::::::::::::::::::::::::::: NOTE
        % v_i field without Ekman velocity finished.
        % <property>_i ready for the flux calculations.


        
        
        
        % ----------------------------
        %% (3.7) EKMAN FLUX CORRECTION
        % ----------------------------
    %     display('+ Computing Ekman flux and velocity ...');

        % Flux mask!
        % 
        % construction is done in Step (3.5)

        % Ekman transport
        ue = (tauy_t(1:end-1)+tauy_t(2:end))*0.5./(f_i*rho0);
        ve = -(taux_t(1:end-1)+taux_t(2:end))*0.5./(f_i*rho0);
        vflux_ek_i = (-ue.*sin(ang_i)+ve.*cos(ang_i)).*dist_i;
        
        
        % Ekman layer depth  
% =========================================================================
% FLAG_VARY_DEK ===========================================================
% =======================================================================
% Determine the way of applying the Ekman layer depth ...
        switch flag_vary_dek 
            case 0
                d_ek = 50*ones(sz_n-1,1); % constant
                
            case 1
                U_wind = sqrt((0.5*(u10_t(1:end-1)+u10_t(2:end))).^2 ...
                    + (0.5*(v10_t(1:end-1)+v10_t(2:end))).^2);
                d_ek = 7.6.*U_wind./sqrt(sind(lat_i));
                
            otherwise
                error('IYT 875!');
        end
% =======================================================================
% =========================================================================


        % Ekman velocity
        v_e_i = vflux_ek_i./d_ek./dist_i;

        
% LOOP I-GRID 
        for i_n = 1:sz_n-1

            % Found the indices
            z_ek = find(abs(depth_i-d_ek(i_n))<10,1,'first');
            
            % Just in case...
            if isempty(z_ek)
                z_ek = 1;
            end
            
% LOOP EKMAN DEPTH
            for i_k = 1:z_ek   
                % Add this Ekman velocity to 
                % the velocity field. So Ekman flux correction is directly 
                % coverted to velocity correction to the velocity field 
                % in the Ekman layer.
                %
                if mask_ek_i(i_n)==1
                    % Incorporate Ekman velocity
                    v_i(i_n,i_k) = v_i(i_n,i_k) + v_e_i(i_n);
                end

% END LOOP EKMAN DEPTH
            end

% END LOOP I-GRID
        end

        % Ekman flux for the entire section
        vflux_ek = nansum(vflux_ek_i); % [m^3/s]
        
        
        
% =========================================================================
% flag_test_direct_velo_only ==============================================
% =======================================================================
    if flag_test_direct_velo_only~= 0
        
        % load data
%         data_sect_mean_velo = [dirio,'/20180625_test5_v_mean.mat']; % used in the 2019 Science paper
%         data_sect_mean_velo = [dirio,'/20181009_test1_v_mean.mat'];
        data_sect_mean_velo = [dirio,'/20181009_test1_v_mean_first11mon.mat'];
        
        stc_tmp = load(data_sect_mean_velo,'v_mean');
        v_mean = stc_tmp.v_mean;
        
        % LOG
        if i_t == i_t_begin
            disp(['+ data_sect_mean_velo= ',data_sect_mean_velo]);
        end
        
        
        switch flag_test_direct_velo_only
            case 1
                % fill the areas without direct measurements;
                v_i(mask_flux_i== 1) = v_mean.mean(mask_flux_i== 1);
        
            case 2
                % find indices between UMM4 and 14.7W
                ipx = find(lon_i>=-021.1435 & lon_i<= -14.7);
            
                % fill with the time-mean velocity
                v_i(ipx,:) = v_mean.mean(ipx,:);
            
            case 3
                % fill the areas that have direct measurements;
                % in opposition to case 1;
                v_i(mask_flux_i== 0) = v_mean.mean(mask_flux_i== 0);
            
            otherwise
                error('XUG 010!');
        end
        
    end
% ========================================================================
% =========================================================================

        

        % :::::::::::::::::::::::::::::::::::::::::::::: NOTE
        % v_i field is ready for computing the volume flux.
        % 


                
        
        % ---------------------------------
        %% (4) Flux corrections
        % ---------------------------------

        
% =========================================================================
% FLAG_VCOMP_EVERYWHERE ===================================================
% =======================================================================
% Determine where to apply Vcomp ...
            if flag_vcomp_everywhere == 1
                % to the entire section
                mask_flux_i = mask_land_i;
                
            elseif flag_vcomp_everywhere== 0
                % do nothing;

            end
% =======================================================================
% =========================================================================


% =========================================================================
% flag_vcomp_selected =====================================================
% =======================================================================
% Apply Vcomp only to a specific part of the 
% section by modifying mask_flux_i;
% Always apply Vcomp to the west section;
%
        if flag_vcomp_selected== 0
            % do thing;
            
            
        elseif flag_vcomp_selected== 1
            
            % NO flux correction is applied 
            % to the geostrophic profiles that have been
            % referenced to the deep moorings;
            ipx = mask_v_baro(:,1)== 0;
            mask_flux_i(ipx,:) = 0;
            
            
        elseif flag_vcomp_selected== 2

            % Apply flux correction 
            % to the area west of the RR (UMM1) 
            % where diret velocity 
            % measurements are not available;
            ipx = lon_i>= -30.5293; % UMM1
            mask_flux_i(ipx,:) = 0;
            
        else 
            error('ISY 819!');
        end
        
        
%         % PLOT
%         figure;
%         pcolor(lon_i,depth_i,mask_flux_i'); shading flat;
%         axis ij; colorbar;
%         title('(3.4.1) mask\_flux\_i');
% =========================================================================
% =========================================================================
        



        % -----------------------------------
        %% (4.1) COMUPTE FLUXES - DEPTH SPACE
        % -----------------------------------
        % (a) the net volume transport constraints
        %

% =========================================================================
% FLAG_VCOMP_CONSTRAINT ===================================================
% =======================================================================
% Net throughflow at the section

        % observed long-term mean transports with uncertainty (reported
        % values)
        vnet_DS = [-1.6 0.20];   % [Sv]  Curry et al. 2014, JPO
        vnet_BS = [-1.0 0.05];    % [Sv]  Woodgate 2018, PiO
        
        
        
        % varying the numbers ...
        if flag_MC == 0
            
            vnet_DS = vnet_DS(1);
            vnet_BS = vnet_BS(1);
            
        elseif flag_MC == 1
            
            vnet_DS = vnet_DS(1) + vnet_DS(2)*randn(1);
            vnet_BS = vnet_BS(1) + vnet_BS(2)*randn(1);
            
        end
        
            
% IF FLAG_VCOMP_CONSTRAINT
        if flag_vcomp_constraint == 0
            % Compensate the whole section:
            vflux_net = 0;
        
        elseif flag_vcomp_constraint == 1
            % Compensate at W and E seperately:
            % -1.6 Sv at W and +0.6 Sv at E;
            
            switch flag_section
                case 'west'
                    vflux_net = vnet_DS; 
                    
                case 'east'
                    vflux_net = vnet_BS - vnet_DS;
                    
                case 'all'
                    vflux_net_w = vnet_DS;
                    vflux_net_e = vnet_BS - vnet_DS;
                    
            end
            
        elseif flag_vcomp_constraint == 2 
            % Compensate at W and E seperately:
            % A zero-net-flow constraint across the full array;

            switch flag_section
                case 'west'
                    vflux_net = vnet_DS;
                    
                case 'east'
                    vflux_net = -vnet_DS; 
                    
                case 'all'
                    vflux_net_w = vnet_DS;
                    vflux_net_e = -vnet_DS; 
                    
            end
            
        else
            error('NEW 002!');
                
            
% END FLAG_VCOMP_CONSTRAINT
        end
% =======================================================================
% ========================================================================= 
    
% return

        % (b) Calculate and apply the compensation transport/velocity
        if exist('vflux_net','var')
        % *** One vflux_comp for the whole section ***    
        
            % total and compensation transport
            vflux = v_i.*area_i;  % [m^3/s]
            vflux_total = nansum(nansum(vflux,1),2);  % [m^3/s]
            vflux_comp = vflux_net*1e6 - vflux_total; % [m^3/s]

            % apply to the area where mask_flux_i== 1
            area_sum = nansum(nansum(area_i(mask_flux_i==1),1),2);
            velo_comp = vflux_comp/area_sum;
            v_i(mask_flux_i==1) = v_i(mask_flux_i==1) + velo_comp;

      
        elseif exist('vflux_net_w','var') && exist('vflux_net_e','var')
        % *** Different vflux_comp for W and E ***   
        
            % flux masks for West and East, seperately
            ipx_w = lon_i<= -45;
            area_w = area_i(ipx_w,:);
            mask_flux_w = mask_flux_i(ipx_w,:);

            ipx_e = lon_i>= -44;
            area_e = area_i(ipx_e,:);
            mask_flux_e = mask_flux_i(ipx_e,:);

            
            % total and compensation transport for West and East,
            % seperately
            vflux = v_i.*area_i;  % [m^3/s]
            vflux_total_w = nansum(nansum(vflux(ipx_w,:),1),2);  % [m^3/s]
            vflux_total_e = nansum(nansum(vflux(ipx_e,:),1),2);  % [m^3/s]
        
            vflux_comp_w = vflux_net_w*1e6 - vflux_total_w; % [m^3/s]
            vflux_comp_e = vflux_net_e*1e6 - vflux_total_e; % [m^3/s]
            
            
            % apply to the area where mask_flux_i== 1
            % but for west and east seperately
            velo_comp_w = vflux_comp_w/(sum(sum(area_w(mask_flux_w==1))));
            v_i(ipx_w,:) = v_i(ipx_w,:) + velo_comp_w.*mask_flux_w;

            velo_comp_e = vflux_comp_e/(sum(sum(area_e(mask_flux_e==1))));
            v_i(ipx_e,:) = v_i(ipx_e,:) + velo_comp_e.*mask_flux_e;
                
% =======================================================================
% =========================================================================
        end
        
        

        % (e) Compensated volume, temp(or heat) and salt flux at each grid point
        vflux = v_i.*area_i*1e-6; % volume transport [Sv]
        hflux = vflux.*d_i.*swcp.*pt_i*1e-9; % temperature transport [PW]
        sflux = vflux.*s_i; % salt transport [Sv]


        % ------------------------------------
        %% (4.2) COMUPTE FLUX - DENSIGY SPACE
        % ------------------------------------
        
        % convert from depth space to to density space
        vflux_d = zspace_to_sigmaspace(vflux,pd_i,dens_bins,'sum'); % changed zspace_to_sigmaspace to zspace_to_sigmaspace_ts
	      
        
%         % LOG
%         disp(['+ cnt_mc= ',num2str(cnt_mc)]);
% %         disp(['+ Total Ekman transport = ',num2str(vflux_ek*1e-6),' Sv']);
% %         disp(['+ Net volume transport constraint = ',num2str(vflux_net),' Sv']);
% %         disp(['+ Net volume transport before applying the constraint = ',num2str(vflux_total*1e-6),' Sv']);
%         disp(['+ Compensation transport required = ',num2str(vflux_comp*1e-6),' Sv']);
% %         disp(['+ Total area = ', num2str(nansum(nansum(area_i,1),2)),' m^2']);
% %         disp(['+ Area to apply flux correction = ',num2str(area_sum),' m^2']);
% %         disp(['+ Compensation velocity = ',num2str(velo_comp),' m/s']);
%         disp(['+ Net, CORRECTED, volume transport = ',num2str(nansum(nansum(vflux))),' Sv']);
%         disp(['+ Net, CORRECTED, volume transport, in density space= ',num2str(nansum(nansum(vflux_d))),' Sv']);

        


        % ---------------------------------
        %% (5) Calculate fluxes
        % ---------------------------------
            
        % ---------------------------------------
        %% (5.1) FLUX ESTIMATE FOR THIS ITERATION
        % ---------------------------------------
        
        % MOC [Sv]
% =========================================================================
% FLAG_MOC_DEF ============================================================
% =========================================================================  
% ========================================================================= 
        % Section-width integrated transport
        vflux_prof_d = squeeze(nansum(vflux_d,1)); % [sigma]
        vflux_prof = squeeze(nansum(vflux,1)); % [depth]
        
        
        % Section-width integrated transport, West and East
        if strcmp(flag_section,'all')
            ipx_w = lon_i<= -45;
            vflux_prof_d_w = squeeze(nansum(vflux_d(ipx_w,:),1)); % [sigma]
            vflux_prof_w = squeeze(nansum(vflux(ipx_w,:),1)); % [depth]
            
            ipx_e = lon_i>= -44;
            vflux_prof_d_e = squeeze(nansum(vflux_d(ipx_e,:),1)); % [sigma]
            vflux_prof_e = squeeze(nansum(vflux(ipx_e,:),1)); % [depth]
        end
                
        
        
        % Deriving MOC, MHT, MFT, MST, Tek... 
        if any(flag_moc_def == 0)
            
            % [sigma] maximum of the streamfunction   
            MOC_n(cnt_mc) = max(cumsum(vflux_prof_d));   % [Sv]
            
            if strcmp(flag_section,'all') && flag_test_AllWE== 1
                MOC_n_w(cnt_mc) = max(cumsum(vflux_prof_d_w));   % [Sv]
                MOC_n_e(cnt_mc) = max(cumsum(vflux_prof_d_e));   % [Sv]
            end
            
% =========================================================================
% FLAG_TEST_PD_MOCWEST ====================================================
% =========================================================================  
            if flag_test_pd_MOCwest== 1 && strcmp(flag_section,'west')
                
                % [sigma_modified] maximum of the streamfunction  
                
                clear tmp tmp1
                
                % pd range for MOC
                tmp = dens_bins>= 1027.0 & dens_bins< 1027.8;
                
                % seraching for MOC in the pd range
                tmp1 = cumsum(vflux_prof_d);
                MOC_n(cnt_mc) = max(tmp1(tmp));
                
            end
% =========================================================================
            
        end
        
        
        if any(flag_moc_def == 1)
            
            % [sigma] Sum of all northward transport
            MOC1_n(cnt_mc) = sum(vflux_prof_d(vflux_prof_d>0));  % [Sv]
            
            if strcmp(flag_section,'all') && flag_test_AllWE== 1
                MOC1_n_w(cnt_mc) = sum(vflux_prof_d_w(vflux_prof_d_w>0));   % [Sv]
                MOC1_n_e(cnt_mc) = sum(vflux_prof_d_e(vflux_prof_d_e>0));   % [Sv]
            end
            
        end
        
        
        if any(flag_moc_def == 2)
            
            % [depth] maximum of the streamfunction   
            MOC2_n(cnt_mc) = max(cumsum(vflux_prof));   % [Sv]
            
            if strcmp(flag_section,'all') && flag_test_AllWE== 1
                MOC2_n_w(cnt_mc) = max(cumsum(vflux_prof_w));   % [Sv]
                MOC2_n_e(cnt_mc) = max(cumsum(vflux_prof_e));   % [Sv]
            end
            
        end
        
        
        if any(flag_moc_def == 3)
            
            % [depth] Sum of all northward transport
            MOC3_n(cnt_mc) = sum(vflux_prof(vflux_prof>0));  % [Sv]
            
            if strcmp(flag_section,'all') && flag_test_AllWE== 1
                MOC3_n_w(cnt_mc) = sum(vflux_prof_w(vflux_prof_w>0));   % [Sv]
                MOC3_n_e(cnt_mc) = sum(vflux_prof_e(vflux_prof_e>0));   % [Sv]
            end
            
        end
% =========================================================================
% =========================================================================
        

        % MHT [PW]
        MHT_n(cnt_mc) = nansum(hflux(:));  
        
        if strcmp(flag_section,'all') && flag_test_AllWE== 1
            MHT_n_w(cnt_mc) = nansum(nansum(hflux(ipx_w,:)));   % [Sv]
            MHT_n_e(cnt_mc) = nansum(nansum(hflux(ipx_e,:)));   % [Sv]
        end

        
        
        % MFT [Sv]
        fflux = - (sflux - vflux.*s0)./s0;  
        MFT_n(cnt_mc) = nansum(fflux(:)); 
        
        if strcmp(flag_section,'all') && flag_test_AllWE== 1
            MFT_n_w(cnt_mc) = nansum(nansum(fflux(ipx_w,:)));   % [Sv]
            MFT_n_e(cnt_mc) = nansum(nansum(fflux(ipx_e,:)));   % [Sv]
        end
        
        
        
        % MST [Sv]
        MST_n(cnt_mc) = nansum(sflux(:));
        
        if strcmp(flag_section,'all') && flag_test_AllWE== 1
            MST_n_w(cnt_mc) = nansum(nansum(sflux(ipx_w,:)));   % [Sv]
            MST_n_e(cnt_mc) = nansum(nansum(sflux(ipx_e,:)));   % [Sv]
        end
        
        
        
        % Ekman transport [Sv]
        Tek_n(cnt_mc) = vflux_ek*1e-6;
        
        if strcmp(flag_section,'all') && flag_test_AllWE== 1
            Tek_n_w(cnt_mc) = nansum(vflux_ek_i(ipx_w,:))*1e-6;   % [Sv]
            Tek_n_e(cnt_mc) = nansum(vflux_ek_i(ipx_e,:))*1e-6;   % [Sv]
        end
        
        
        
        % Compensation transport [Sv]
        if exist('vflux_net','var')
            Text_n(cnt_mc) = vflux_comp.*1e-6;

        elseif exist('vflux_net_w','var') && exist('vflux_net_e','var')
            Text_n(cnt_mc) = (vflux_comp_w+vflux_comp_e).*1e-6;
            Text_n_w(cnt_mc) = (vflux_comp_w).*1e-6;
            Text_n_e(cnt_mc) = (vflux_comp_e).*1e-6;
        end
            
        
        
        % Amend transport profiles
        Tprof_d_n(:,cnt_mc) = vflux_prof_d;
        Tprof_n(:,cnt_mc) = vflux_prof;
        
        if strcmp(flag_section,'all') && flag_test_AllWE== 1
            Tprof_d_n_w(:,cnt_mc) = vflux_prof_d_w;
            Tprof_n_w(:,cnt_mc) = vflux_prof_w;
            
            Tprof_d_n_e(:,cnt_mc) = vflux_prof_d_e;
            Tprof_n_e(:,cnt_mc) = vflux_prof_e;
        end

        
        clear vflux_prof_d*
        
        % -------------------------------------
        %% (5.2) DIFFERENCE IN RUNNING AVERAGES
        % -------------------------------------


% =========================================================================
% FLAG_MC =================================================================
% =======================================================================
        if flag_MC==0
            
            delta_MOC = 0;
            delta_MOC1 = 0;
            delta_MOC2 = 0;
            delta_MOC3 = 0;
            delta_MHT = 0;
            delta_MFT = 0;
            
            disp(['+ cnt_mc= ',num2str(cnt_mc),', will quit MC ...']);
            
        elseif flag_MC == 1
            
            if cnt_mc > 1
                
                % difference in the successive running averages
                delta_MOC = mean(MOC_n(1:cnt_mc)) - mean(MOC_n(1:cnt_mc-1)); 
                delta_MHT = mean(MHT_n(1:cnt_mc)) - mean(MHT_n(1:cnt_mc-1));
                delta_MFT = mean(MFT_n(1:cnt_mc)) - mean(MFT_n(1:cnt_mc-1));
                
                if isnan(MOC1_n)
                    delta_MOC1 = 0;
                else
                    delta_MOC1 = mean(MOC1_n(1:cnt_mc)) - mean(MOC1_n(1:cnt_mc-1));
                end
                
                if isnan(MOC2_n)
                    delta_MOC2 = 0;
                else
                    delta_MOC2 = mean(MOC2_n(1:cnt_mc)) - mean(MOC2_n(1:cnt_mc-1));
                end
                
                if isnan(MOC3_n)
                    delta_MOC3 = 0;
                else
                    delta_MOC3 = mean(MOC3_n(1:cnt_mc)) - mean(MOC3_n(1:cnt_mc-1));
                end
                
            else
                % do nothing;
                
            end
            
        end
% =======================================================================
% =========================================================================        
        

%         % LOG
% %         disp(['+   cnt_mc= ',num2str(cnt_mc)]);
%         disp(['+ MOC= ',num2str(MOC_n(cnt_mc))]);
% %         disp(['+     delta_MOC= ',num2str(delta_MOC),' Sv']);
% %         disp('+');



% =========================================================================
% flag_HT_FWT_decomposition ===============================================
% =======================================================================
% MHT and MFT decompsition in both density and depth space
        if flag_HT_FWT_decomposition == 1

            % decomposing...
            [~,mht_z,mht_d,~,mft_z,mft_d] = ...
                fun_HT_FWT_decomposition(depth,dens_bins,area_i,...
                v_i,pd_i,pt_i,s_i,vflux,vflux_d,s0,0);

            
            % saving the components from this iteration ...
            % *** Heat ***
            MHT_n_dia(cnt_mc) = mht_d.diapycnal;
            MHT_n_iso(cnt_mc) = mht_d.isopycnal;
            MHT_n_net(cnt_mc) = mht_d.net;
            
            MHT_n_over(cnt_mc) = mht_z.overturning;
            MHT_n_gyre(cnt_mc) = mht_z.gyre;
            MHT_n_netz(cnt_mc) = mht_z.net;
            
            
            % *** Freshwater ***
            MFT_n_dia(cnt_mc) = mft_d.diapycnal;
            MFT_n_iso(cnt_mc) = mft_d.isopycnal;
            MFT_n_net(cnt_mc) = mft_d.net;
            
            MFT_n_over(cnt_mc) = mft_z.overturning;
            MFT_n_gyre(cnt_mc) = mft_z.gyre;
            MFT_n_netz(cnt_mc) = mft_z.net;
            
            
            clear mht_z mht_d mft_z mft_d
            
            
            
            if strcmp(flag_section,'all') && flag_test_AllWE== 1
                % west
                
                % decomposing...
                [~,mht_z,mht_d,~,mft_z,mft_d] = ...
                    fun_HT_FWT_decomposition(depth,dens_bins,area_i(ipx_w,:),...
                    v_i(ipx_w,:),pd_i(ipx_w,:),pt_i(ipx_w,:),s_i(ipx_w,:),vflux(ipx_w,:),vflux_d(ipx_w,:),s0,0);


                % saving the components from this iteration ...
                % *** Heat ***
                MHT_n_dia_w(cnt_mc) = mht_d.diapycnal;
                MHT_n_iso_w(cnt_mc) = mht_d.isopycnal;
                MHT_n_net_w(cnt_mc) = mht_d.net;

                MHT_n_over_w(cnt_mc) = mht_z.overturning;
                MHT_n_gyre_w(cnt_mc) = mht_z.gyre;
                MHT_n_netz_w(cnt_mc) = mht_z.net;


                % *** Freshwater ***
                MFT_n_dia_w(cnt_mc) = mft_d.diapycnal;
                MFT_n_iso_w(cnt_mc) = mft_d.isopycnal;
                MFT_n_net_w(cnt_mc) = mft_d.net;

                MFT_n_over_w(cnt_mc) = mft_z.overturning;
                MFT_n_gyre_w(cnt_mc) = mft_z.gyre;
                MFT_n_netz_w(cnt_mc) = mft_z.net;
            

            
                clear mht_z mht_d mft_z mft_d
            
            
                % east
                
                % decomposing...
                [~,mht_z,mht_d,~,mft_z,mft_d] = ...
                    fun_HT_FWT_decomposition(depth,dens_bins,area_i(ipx_e,:),...
                    v_i(ipx_e,:),pd_i(ipx_e,:),pt_i(ipx_e,:),s_i(ipx_e,:),vflux(ipx_e,:),vflux_d(ipx_e,:),s0,0);


                % saving the components from this iteration ...
                % *** Heat ***
                MHT_n_dia_e(cnt_mc) = mht_d.diapycnal;
                MHT_n_iso_e(cnt_mc) = mht_d.isopycnal;
                MHT_n_net_e(cnt_mc) = mht_d.net;

                MHT_n_over_e(cnt_mc) = mht_z.overturning;
                MHT_n_gyre_e(cnt_mc) = mht_z.gyre;
                MHT_n_netz_e(cnt_mc) = mht_z.net;


                % *** Freshwater ***
                MFT_n_dia_e(cnt_mc) = mft_d.diapycnal;
                MFT_n_iso_e(cnt_mc) = mft_d.isopycnal;
                MFT_n_net_e(cnt_mc) = mft_d.net;

                MFT_n_over_e(cnt_mc) = mft_z.overturning;
                MFT_n_gyre_e(cnt_mc) = mft_z.gyre;
                MFT_n_netz_e(cnt_mc) = mft_z.net;
                
                clear mht_z mht_d mft_z mft_d
                
            end
            
        end
% =========================================================================
% =========================================================================





        % --------------------------------
        %% (6) Accumulate data fields
        % --------------------------------
        if cnt_mc == 1
            
            % Initialize
            if flag_baro_corr==1
                v_ssh_accu = v_ssh;
                v_baro_accu = v_baro;
            end
            
            velo_accu = v_i;
            temp_accu = t_i;
            ptmp_accu = pt_i;
            salt_accu = s_i;
            dens_accu = d_i; 
            pden_accu = pd_i;
         
            vflux_accu = vflux;
            hflux_accu = hflux;
            sflux_accu = sflux;
            
%             vflux_ek_accu = vflux_ek*1e-6; % [Sv]
%             
%             if exist('vflux_net','var')
%                 vflux_ext_accu = vflux_comp.*1e-6; % Compensation transport [Sv]
%                 
%             elseif exist('vflux_net_w','var') && exist('vflux_net_e','var')
%                 vflux_ext_accu = [vflux_comp_w vflux_comp_e].*1e-6;
%             end



        else
            
            % Accumulate
            if flag_baro_corr==1
                v_ssh_accu = v_ssh_accu + v_ssh;
                v_baro_accu = v_baro_accu + v_baro;
            end
    
            velo_accu = velo_accu + v_i;
            temp_accu = temp_accu + t_i;
            ptmp_accu = ptmp_accu + pt_i;
            salt_accu = salt_accu + s_i;
            dens_accu = dens_accu + d_i;
            pden_accu = pden_accu + pd_i;     
            
            vflux_accu = vflux_accu + vflux;
            hflux_accu = hflux_accu + hflux;
            sflux_accu = sflux_accu + sflux;
            
%             vflux_ek_accu = vflux_ek_accu + vflux_ek*1e-6; % [Sv]
%             
%             if exist('vflux_net','var')
%                 vflux_ext_accu = vflux_ext_accu + vflux_comp.*1e-6; % Compensation transport [Sv]
%                 
%             elseif exist('vflux_net_w','var') && exist('vflux_net_e','var')
%                 vflux_ext_accu = vflux_ext_accu + [vflux_comp_w vflux_comp_e].*1e-6; % Compensation transport [Sv]
%             end

        end 

        
%         % LOG for checking purpose ...
%         disp(['+    cnt_mc= ',num2str(cnt_mc)]);
%         disp(['+    CORRECTED, vflux = ',num2str(nansum(nansum(vflux,1),2)),' Sv']);
%         disp(['+    CORRECTED, vflux_accu./cnt_mc = ',num2str(nansum(nansum(vflux_accu./cnt_mc,1),2)),' Sv']);
        

        % Double check the net throughflow at this iteration ...
        clear vnet_tmp
        if exist('vflux_net','var')
            vnet_tmp = vflux_net;
        elseif exist('vflux_net_w','var') && exist('vflux_net_e','var')
            vnet_tmp = vflux_net_w + vflux_net_e;
        end
        
        if abs(nansum(nansum(vflux,1),2) - vnet_tmp) > 1e-5
            error('SOS 701!');
        end
        
        
        % Update the counter for #iteration
        cnt_mc = cnt_mc + 1;
                   
  
% END LOOP MONTE CARLO SIMULATION   
    end
    

    
    % ------------------------------------
    %% (7) Outputs for this time interval
    % ------------------------------------
	% Mean fluxes with uncertainty [Sv]
    
    % Transport profile
    Tprof_mc(i_per,:,:) = [nanmean(Tprof_n,2) nanstd(Tprof_n,0,2)];
    Tprof_d_mc(i_per,:,:) = [nanmean(Tprof_d_n,2) nanstd(Tprof_d_n,0,2)];

    if strcmp(flag_section,'all') && flag_test_AllWE== 1
        % *** West ***
        Tprof_mc_w(i_per,:,:) = [nanmean(Tprof_n_w,2) nanstd(Tprof_n_w,0,2)];
        Tprof_d_mc_w(i_per,:,:) = [nanmean(Tprof_d_n_w,2) nanstd(Tprof_d_n_w,0,2)];

        % *** East ***
        Tprof_mc_e(i_per,:,:) = [nanmean(Tprof_n_e,2) nanstd(Tprof_n_e,0,2)];
        Tprof_d_mc_e(i_per,:,:) = [nanmean(Tprof_d_n_e,2) nanstd(Tprof_d_n_e,0,2)];

    end
    
    
    
    
    % Ekman transport 
    Tek_n(isnan(Tek_n)) = []; 
    Tek_mc(i_per,:) = [mean(Tek_n) std(Tek_n)];
    
    if strcmp(flag_section,'all') && flag_test_AllWE== 1
        Tek_n_w(isnan(Tek_n_w)) = [];  
        Tek_mc_w(i_per,:) = [mean(Tek_n_w) std(Tek_n_w)];
        
        Tek_n_e(isnan(Tek_n_e)) = [];  
        Tek_mc_e(i_per,:) = [mean(Tek_n_e) std(Tek_n_e)];
    end
    
    
    
    
    
    % Compensation transport
    if exist('vflux_net','var')
        Text_n(isnan(Text_n)) = []; 
        Text_mc(i_per,:) = [mean(Text_n) std(Text_n)];
        
        clear Text_n_w Text_n_e Text_mc_w Text_mc_e
        
    elseif exist('vflux_net_w','var') && exist('vflux_net_e','var')
        Text_n(isnan(Text_n)) = []; 
        Text_mc(i_per,:) = [mean(Text_n) std(Text_n)];
    
        Text_n_w(isnan(Text_n_w)) = []; 
        Text_mc_w(i_per,:) = [mean(Text_n_w) std(Text_n_w)];
        
        Text_n_e(isnan(Text_n_e)) = []; 
        Text_mc_e(i_per,:) = [mean(Text_n_e) std(Text_n_e)];
        
    end
        
    
    
    
    
    % Fluxes
    MOC_n(isnan(MOC_n)) = [];
    MOC_mc(i_per,:) = [mean(MOC_n) std(MOC_n)];
    
    
    MOC1_n(isnan(MOC1_n)) = [];
    MOC1_mc(i_per,:) = [mean(MOC1_n) std(MOC1_n)];
    
    MOC2_n(isnan(MOC2_n)) = [];
    MOC2_mc(i_per,:) = [mean(MOC2_n) std(MOC2_n)];
    
    MOC3_n(isnan(MOC3_n)) = [];
    MOC3_mc(i_per,:) = [mean(MOC3_n) std(MOC3_n)];
    
    MHT_n(isnan(MHT_n)) = [];
    MHT_mc(i_per,:) = [mean(MHT_n) std(MHT_n)];
    
    MFT_n(isnan(MFT_n)) = [];
    MFT_mc(i_per,:) = [mean(MFT_n) std(MFT_n)];
    
    MST_n(isnan(MST_n)) = [];
    MST_mc(i_per,:) = [mean(MST_n) std(MST_n)];
    
    
    
    if strcmp(flag_section,'all') && flag_test_AllWE== 1
        % *** West ***
        % Fluxes
        MOC_n_w(isnan(MOC_n_w)) = [];
        MOC_mc_w(i_per,:) = [mean(MOC_n_w) std(MOC_n_w)];

        MOC1_n_w(isnan(MOC1_n_w)) = [];
        MOC1_mc_w(i_per,:) = [mean(MOC1_n_w) std(MOC1_n_w)];

        MOC2_n_w(isnan(MOC2_n_w)) = [];
        MOC2_mc_w(i_per,:) = [mean(MOC2_n_w) std(MOC2_n_w)];

        MOC3_n_w(isnan(MOC3_n_w)) = [];
        MOC3_mc_w(i_per,:) = [mean(MOC3_n_w) std(MOC3_n_w)];

        MHT_n_w(isnan(MHT_n_w)) = [];
        MHT_mc_w(i_per,:) = [mean(MHT_n_w) std(MHT_n_w)];

        MFT_n_w(isnan(MFT_n_w)) = [];
        MFT_mc_w(i_per,:) = [mean(MFT_n_w) std(MFT_n_w)];

        MST_n_w(isnan(MST_n_w)) = [];
        MST_mc_w(i_per,:) = [mean(MST_n_w) std(MST_n_w)];
    
    
        % *** East ***
        % Fluxes
        MOC_n_e(isnan(MOC_n_e)) = [];
        MOC_mc_e(i_per,:) = [mean(MOC_n_e) std(MOC_n_e)];

        MOC1_n_e(isnan(MOC1_n_e)) = [];
        MOC1_mc_e(i_per,:) = [mean(MOC1_n_e) std(MOC1_n_e)];

        MOC2_n_e(isnan(MOC2_n_e)) = [];
        MOC2_mc_e(i_per,:) = [mean(MOC2_n_e) std(MOC2_n_e)];

        MOC3_n_e(isnan(MOC3_n_e)) = [];
        MOC3_mc_e(i_per,:) = [mean(MOC3_n_e) std(MOC3_n_e)];

        MHT_n_e(isnan(MHT_n_e)) = [];
        MHT_mc_e(i_per,:) = [mean(MHT_n_e) std(MHT_n_e)];

        MFT_n_e(isnan(MFT_n_e)) = [];
        MFT_mc_e(i_per,:) = [mean(MFT_n_e) std(MFT_n_e)];

        MST_n_e(isnan(MST_n_e)) = [];
        MST_mc_e(i_per,:) = [mean(MST_n_e) std(MST_n_e)];
    end
    
        
    
% =========================================================================
% flag_HT_FWT_decomposition ===============================================
% =======================================================================
% MHT and MFT decompsition in both density and depth space
    if flag_HT_FWT_decomposition == 1

        % *** Heat ***
        MHT_n_dia(isnan(MHT_n_dia)) = [];
        MHT_mc_dia(i_per,:) = [mean(MHT_n_dia) std(MHT_n_dia)];

        MHT_n_iso(isnan(MHT_n_iso)) = [];
        MHT_mc_iso(i_per,:) = [mean(MHT_n_iso) std(MHT_n_iso)];
        
        MHT_n_net(isnan(MHT_n_net)) = [];
        MHT_mc_net(i_per,:) = [mean(MHT_n_net) std(MHT_n_net)];

        MHT_n_over(isnan(MHT_n_over)) = [];
        MHT_mc_over(i_per,:) = [mean(MHT_n_over) std(MHT_n_over)];

        MHT_n_gyre(isnan(MHT_n_gyre)) = [];
        MHT_mc_gyre(i_per,:) = [mean(MHT_n_gyre) std(MHT_n_gyre)];

        MHT_n_netz(isnan(MHT_n_netz)) = [];
        MHT_mc_netz(i_per,:) = [mean(MHT_n_netz) std(MHT_n_netz)];


        % *** Freshwater ***
        MFT_n_dia(isnan(MFT_n_dia)) = [];
        MFT_mc_dia(i_per,:) = [mean(MFT_n_dia) std(MFT_n_dia)];

        MFT_n_iso(isnan(MFT_n_iso)) = [];
        MFT_mc_iso(i_per,:) = [mean(MFT_n_iso) std(MFT_n_iso)];
        
        MFT_n_net(isnan(MFT_n_net)) = [];
        MFT_mc_net(i_per,:) = [mean(MFT_n_net) std(MFT_n_net)];

        MFT_n_over(isnan(MFT_n_over)) = [];
        MFT_mc_over(i_per,:) = [mean(MFT_n_over) std(MFT_n_over)];

        MFT_n_gyre(isnan(MFT_n_gyre)) = [];
        MFT_mc_gyre(i_per,:) = [mean(MFT_n_gyre) std(MFT_n_gyre)];

        MFT_n_netz(isnan(MFT_n_netz)) = [];
        MFT_mc_netz(i_per,:) = [mean(MFT_n_netz) std(MFT_n_netz)];



        if strcmp(flag_section,'all') && flag_test_AllWE== 1
            
            % *** West ***
                % *** Heat ***
            MHT_n_dia_w(isnan(MHT_n_dia_w)) = [];
            MHT_mc_dia_w(i_per,:) = [mean(MHT_n_dia_w) std(MHT_n_dia_w)];

            MHT_n_iso_w(isnan(MHT_n_iso_w)) = [];
            MHT_mc_iso_w(i_per,:) = [mean(MHT_n_iso_w) std(MHT_n_iso_w)];
            
            MHT_n_net_w(isnan(MHT_n_net_w)) = [];
            MHT_mc_net_w(i_per,:) = [mean(MHT_n_net_w) std(MHT_n_net_w)];

            MHT_n_over_w(isnan(MHT_n_over_w)) = [];
            MHT_mc_over_w(i_per,:) = [mean(MHT_n_over_w) std(MHT_n_over_w)];

            MHT_n_gyre_w(isnan(MHT_n_gyre_w)) = [];
            MHT_mc_gyre_w(i_per,:) = [mean(MHT_n_gyre_w) std(MHT_n_gyre_w)];

            MHT_n_netz_w(isnan(MHT_n_netz_w)) = [];
            MHT_mc_netz_w(i_per,:) = [mean(MHT_n_netz_w) std(MHT_n_netz_w)];


                % *** Freshwater ***
            MFT_n_dia_w(isnan(MFT_n_dia_w)) = [];
            MFT_mc_dia_w(i_per,:) = [mean(MFT_n_dia_w) std(MFT_n_dia_w)];

            MFT_n_iso_w(isnan(MFT_n_iso_w)) = [];
            MFT_mc_iso_w(i_per,:) = [mean(MFT_n_iso_w) std(MFT_n_iso_w)];
            
            MFT_n_net_w(isnan(MFT_n_net_w)) = [];
            MFT_mc_net_w(i_per,:) = [mean(MFT_n_net_w) std(MFT_n_net_w)];

            MFT_n_over_w(isnan(MFT_n_over_w)) = [];
            MFT_mc_over_w(i_per,:) = [mean(MFT_n_over_w) std(MFT_n_over_w)];

            MFT_n_gyre_w(isnan(MFT_n_gyre_w)) = [];
            MFT_mc_gyre_w(i_per,:) = [mean(MFT_n_gyre_w) std(MFT_n_gyre_w)];

            MFT_n_netz_w(isnan(MFT_n_netz_w)) = [];
            MFT_mc_netz_w(i_per,:) = [mean(MFT_n_netz_w) std(MFT_n_netz_w)];
        
        
            
            % *** East ***
                % *** Heat ***
            MHT_n_dia_e(isnan(MHT_n_dia_e)) = [];
            MHT_mc_dia_e(i_per,:) = [mean(MHT_n_dia_e) std(MHT_n_dia_e)];

            MHT_n_iso_e(isnan(MHT_n_iso_e)) = [];
            MHT_mc_iso_e(i_per,:) = [mean(MHT_n_iso_e) std(MHT_n_iso_e)];
            
            MHT_n_net_e(isnan(MHT_n_net_e)) = [];
            MHT_mc_net_e(i_per,:) = [mean(MHT_n_net_e) std(MHT_n_net_e)];

            MHT_n_over_e(isnan(MHT_n_over_e)) = [];
            MHT_mc_over_e(i_per,:) = [mean(MHT_n_over_e) std(MHT_n_over_e)];

            MHT_n_gyre_e(isnan(MHT_n_gyre_e)) = [];
            MHT_mc_gyre_e(i_per,:) = [mean(MHT_n_gyre_e) std(MHT_n_gyre_e)];

            MHT_n_netz_e(isnan(MHT_n_netz_e)) = [];
            MHT_mc_netz_e(i_per,:) = [mean(MHT_n_netz_e) std(MHT_n_netz_e)];


                % *** Freshwater ***
            MFT_n_dia_e(isnan(MFT_n_dia_e)) = [];
            MFT_mc_dia_e(i_per,:) = [mean(MFT_n_dia_e) std(MFT_n_dia_e)];

            MFT_n_iso_e(isnan(MFT_n_iso_e)) = [];
            MFT_mc_iso_e(i_per,:) = [mean(MFT_n_iso_e) std(MFT_n_iso_e)];
            
            MFT_n_net_e(isnan(MFT_n_net_e)) = [];
            MFT_mc_net_e(i_per,:) = [mean(MFT_n_net_e) std(MFT_n_net_e)];

            MFT_n_over_e(isnan(MFT_n_over_e)) = [];
            MFT_mc_over_e(i_per,:) = [mean(MFT_n_over_e) std(MFT_n_over_e)];

            MFT_n_gyre_e(isnan(MFT_n_gyre_e)) = [];
            MFT_mc_gyre_e(i_per,:) = [mean(MFT_n_gyre_e) std(MFT_n_gyre_e)];

            MFT_n_netz_e(isnan(MFT_n_netz_e)) = [];
            MFT_mc_netz_e(i_per,:) = [mean(MFT_n_netz_e) std(MFT_n_netz_e)];
            
            
        end


    end
% =========================================================================
% =========================================================================



    
% ////////////////////////////////////////////////////
% // For checking purpose only ///////////////////////
% ////////////////////////////////////////////////////
    if i_per == 1
        % initializing structure array 
        % with one field filled with
        % nans ...
        MC_runs = struct('MOC',nan(sz_t_max,1));
        
        if strcmp(flag_section,'all') && flag_test_AllWE== 1
            MC_runs_w = struct('MOC',nan(sz_t_max,1));
            MC_runs_e = struct('MOC',nan(sz_t_max,1));
        end
    end

    % MOC, MHT and MFT
    MC_runs(i_per).MOC = MOC_n;
    MC_runs(i_per).MOC1 = MOC1_n; 
    MC_runs(i_per).MOC2 = MOC2_n; 
    MC_runs(i_per).MOC3 = MOC3_n; 
    MC_runs(i_per).MHT = MHT_n; 
    MC_runs(i_per).MFT = MFT_n;
    MC_runs(i_per).MST = MST_n;
    
    % Ekman transport
    MC_runs(i_per).Tek = Tek_n;
    
    % Compensation transport
    MC_runs(i_per).Text = Text_n;
    if  exist('vflux_net_w','var') && exist('vflux_net_e','var')
        MC_runs(i_per).Text_w = Text_n_w;
        MC_runs(i_per).Text_e = Text_n_e;
    end
    
    
    if strcmp(flag_section,'all') && flag_test_AllWE== 1
        % *** West ***
            % MOC, MHT and MFT
        MC_runs_w(i_per).MOC = MOC_n_w;
        MC_runs_w(i_per).MOC1 = MOC1_n_w; 
        MC_runs_w(i_per).MOC2 = MOC2_n_w; 
        MC_runs_w(i_per).MOC3 = MOC3_n_w; 
        MC_runs_w(i_per).MHT = MHT_n_w; 
        MC_runs_w(i_per).MFT = MFT_n_w;
        MC_runs_w(i_per).MST = MST_n_w;

            % Ekman transport
        MC_runs_w(i_per).Tek = Tek_n_w;
        
        
        % *** East ***
            % MOC, MHT and MFT
        MC_runs_e(i_per).MOC = MOC_n_e;
        MC_runs_e(i_per).MOC1 = MOC1_n_e; 
        MC_runs_e(i_per).MOC2 = MOC2_n_e; 
        MC_runs_e(i_per).MOC3 = MOC3_n_e; 
        MC_runs_e(i_per).MHT = MHT_n_e; 
        MC_runs_e(i_per).MFT = MFT_n_e;
        MC_runs_e(i_per).MST = MST_n_e;

            % Ekman transport
        MC_runs_e(i_per).Tek = Tek_n_e;
    end
    
    
% % =========================================================================
% % flag_HT_FWT_decomposition ===============================================
% % =======================================================================
% % MHT and MFT components
%     if flag_HT_FWT_decomposition == 1
% 
%         % *** Heat ***
%         MC_runs(i_per).MHT_dia = MHT_n_dia;
%         MC_runs(i_per).MHT_iso = MHT_n_iso;
%         MC_runs(i_per).MHT_over = MHT_n_over;
%         MC_runs(i_per).MHT_gyre = MHT_n_gyre;
%         MC_runs(i_per).MHT_net = MHT_n_net;
%         
%         % *** Freshwater ***
%         MC_runs(i_per).MFT_dia = MFT_n_dia;
%         MC_runs(i_per).MFT_iso = MFT_n_iso;
%         MC_runs(i_per).MFT_over = MFT_n_over;
%         MC_runs(i_per).MFT_gyre = MFT_n_gyre;
%         MC_runs(i_per).MFT_net = MFT_n_net;
% %         MC_runs(i_per).MFT_wBS = MFT_n_wBS;
% %         MC_runs(i_per).MFT_div = MFT_n_div;
% %         MC_runs(i_per).MFT_net_div = MFT_n_net_div;
%     end
% % =========================================================================
% % =========================================================================


    clear MOC1_n* MOC2_n* MOC3_n* MOC_n* MHT_n* MFT_n* MST_n* Tek_n* Text_n*
    clear MHT_n_dia* MHT_n_iso* MHT_n_over* MHT_n_gyre* MHT_n_net*
    clear MFT_n_dia* MFT_n_iso* MFT_n_over* MFT_n_gyre* MFT_n_net* MFT_n_wBS* MFT_n_div* MFT_n_net_div*
    

% ////////////////////////////////////////////////////
% // For checking purpose only ///////////////////////
% ////////////////////////////////////////////////////


    
    
    % Mean fields averaged from all iterations
	if flag_baro_corr==1
        v_ssh_accu = v_ssh_accu./(cnt_mc-1);
        v_baro_accu = v_baro_accu./(cnt_mc-1);
	end
%     
    if cnt_mc~=1
        velo_accu = velo_accu./(cnt_mc-1);
        temp_accu = temp_accu./(cnt_mc-1);
        ptmp_accu = ptmp_accu./(cnt_mc-1);
        salt_accu = salt_accu./(cnt_mc-1);
%         dens_accu = dens_accu./(cnt_mc-1);
        pden_accu = pden_accu./(cnt_mc-1);
%         
%         vflux_accu = vflux_accu./(cnt_mc-1);  % [Sv]
%         hflux_accu = hflux_accu./(cnt_mc-1);  % [PW]
%         sflux_accu = sflux_accu./(cnt_mc-1);  % [Sv]
    end
%     vflux_ek_accu = vflux_ek_accu./(cnt_mc-1);  % [Sv]
%     vflux_ext_accu = vflux_ext_accu./(cnt_mc-1);  % Compensation transport [Sv]
    
            
% 	% *** Assemble
    per(i_per) = time(i_t); % time
%     
    if flag_baro_corr==1
        v_ssh_per(:,i_per) = v_ssh_accu;
        v_baro_per(:,i_per) = v_baro_accu;
    end
%     
    velo_per(:,:,i_per) = velo_accu;
    temp_per(:,:,i_per) = temp_accu;
    ptmp_per(:,:,i_per) = ptmp_accu;
    salt_per(:,:,i_per) = salt_accu;
%     dens_per(:,:,i_per) = dens_accu;
    pden_per(:,:,i_per) = pden_accu;
%     
    vflux_per(:,:,i_per) = vflux_accu;
    hflux_per(:,:,i_per) = hflux_accu;
    sflux_per(:,:,i_per) = sflux_accu;
%     
% %     vflux_ek_per(i_per) = vflux_ek_accu;
% %     vflux_ext_per(i_per,:) = vflux_ext_accu;
% 
% 
%     
% 	*** Conversion of fluxes to density space
    vflux_d_accu = zspace_to_sigmaspace(vflux_accu,pden_accu,dens_bins,'sum');
  hflux_d_accu = zspace_to_sigmaspace(hflux_accu,pden_accu,dens_bins,'sum');   
	sflux_d_accu = zspace_to_sigmaspace(sflux_accu,pden_accu,dens_bins,'sum');   
%         
    vflux_d_per(:,:,i_per) = vflux_d_accu;
    hflux_d_per(:,:,i_per) = hflux_d_accu;
    sflux_d_per(:,:,i_per) = sflux_d_accu;



    

    % LOG
    disp(['+ i_per= ',num2str(i_per)]);
    disp(['+ Total MC iterations= ',num2str(cnt_mc-1)]);
    disp(['+ Ekman transport = ',num2str(Tek_mc(i_per,1)),' Sv']);
    disp(['+ Total compensation transport required = ',num2str(Text_mc(i_per,1)),' Sv']);
    disp(['+ Net, CORRECTED, volume transport = ',num2str(nansum(nansum(vflux_accu))),' Sv']);
%     disp(['+ Net, CORRECTED, volume transport, in density space= ',num2str(nansum(nansum(vflux_d_accu))),' Sv']);
    disp(['+ MOC ave= ',num2str(MOC_mc(i_per,1)),' std= ',num2str(MOC_mc(i_per,2))]);
    disp(['+ MHT ave= ',num2str(MHT_mc(i_per,1)),' std= ',num2str(MHT_mc(i_per,2))]);
    disp(['+ MFT ave= ',num2str(MFT_mc(i_per,1)),' std= ',num2str(MFT_mc(i_per,2))]);
    
    
    if strcmp(flag_section,'all') && flag_test_AllWE== 1
        disp(['+ MOC_e ave= ',num2str(MOC_mc_e(i_per,1)),' std= ',num2str(MOC_mc_e(i_per,2))]);
        disp(['+ MOC_w ave= ',num2str(MOC_mc_w(i_per,1)),' std= ',num2str(MOC_mc_w(i_per,2))]);
    end
    
    disp('+');
    

    clear velo_accu temp_accu ptmp_accu salt_accu dens_accu pden_accu
    clear vflux_ek_accu vflux_ext_accu 
    clear vflux_accu hflux_accu sflux_accu
    clear vflux_d_accu hflux_d_accu sflux_d_accu
    
    
    % Update the counter
    i_per = i_per + 1;
    

% END LOOP TIMESTEP    
end




%------------------+
%% SAVE AND CLEAN  |
%------------------+
disp('+');
disp('+');
disp('+------------------+');
disp('| Save and clean   |');
disp('+------------------+');

% *******************
% *** Coordinates ***
% *******************
ose.time = per;
ose.lon = lon_i;
ose.lat = lat_i;
ose.depth = depth_i;
ose.depth_d = dens_bins; % potential density bins
% ose.area = area_i;
% ose.mask_flux = mask_flux_i;  % grid cell to apply flux correction
% 

% *********************************
% *** Along-section data fields ***
% % *******************************
ose.velo = velo_per;    % velocity 
ose.temp = temp_per;	% in-situ temperature 
ose.ptmp = ptmp_per;    % potential temperature
ose.salt = salt_per;	% practical salinity
% ose.dens = dens_per;	% in-situ denisty
ose.pden = pden_per;    % potential density 
ose.vflux = vflux_per;  % volume flux [Sv]
ose.hflux = hflux_per;  % heat flux [PW]
ose.sflux = sflux_per;  % salt flux [Sv]
ose.vflux_d = vflux_d_per;
ose.hflux_d = hflux_d_per;
ose.sflux_d = sflux_d_per;
% 
% % ose.vflux_ek = vflux_ek_per; 
% % ose.vflux_ext = vflux_ext_per; 
% 
if flag_baro_corr==1
    ose.v_baro = v_baro_per; % barotropic velocity
    ose.v_ssh = v_ssh_per;
end


% ************************
% *** Flux time series ***
% ************************
% % Transport profils
% ose.Tprof = Tprof_mc;
% ose.Tprof_d = Tprof_d_mc;

% MOC, MHT and MFT 
ose.MOC = MOC_mc;   % maxMOCsigma
ose.MOC1 = MOC1_mc; % northMOCsigma
ose.MOC2 = MOC2_mc; % maxMOCz
ose.MOC3 = MOC3_mc; % northMOCz

ose.MHT = MHT_mc;
ose.MFT = MFT_mc;
ose.MST = MST_mc;

% Ekman transport 
ose.Tek = Tek_mc; 

% Compensation transport 
if exist('vflux_net','var')
    ose.Text = Text_mc;

elseif exist('vflux_net_w','var') && exist('vflux_net_e','var')
    ose.Text = Text_mc;
    ose.Text_w = Text_mc_w;
    ose.Text_e = Text_mc_e;
end
    
if strcmp(flag_section,'all') && flag_test_AllWE== 1
    % *** West ***
        % Transport profils
    ose.Tprof_w = Tprof_mc_w;
    ose.Tprof_d_w = Tprof_d_mc_w;

        % MOC, MHT and MFT 
    ose.MOC_w = MOC_mc_w;   % maxMOCsigma
    ose.MOC1_w = MOC1_mc_w; % northMOCsigma
    ose.MOC2_w = MOC2_mc_w; % maxMOCz
    ose.MOC3_w = MOC3_mc_w; % northMOCz

    ose.MHT_w = MHT_mc_w;
    ose.MFT_w = MFT_mc_w;
    ose.MST_w = MST_mc_w;

        % Ekman transport 
    ose.Tek_w = Tek_mc_w; 

    
    % *** East ***
        % Transport profils
    ose.Tprof_e = Tprof_mc_e;
    ose.Tprof_d_e = Tprof_d_mc_e;

        % MOC, MHT and MFT 
    ose.MOC_e = MOC_mc_e;   % maxMOCsigma
    ose.MOC1_e = MOC1_mc_e; % northMOCsigma
    ose.MOC2_e = MOC2_mc_e; % maxMOCz
    ose.MOC3_e = MOC3_mc_e; % northMOCz

    ose.MHT_e = MHT_mc_e;
    ose.MFT_e = MFT_mc_e;
    ose.MST_e = MST_mc_e;

        % Ekman transport 
    ose.Tek_e = Tek_mc_e; 

end


% =========================================================================
% flag_HT_FWT_decomposition ===============================================
% =======================================================================
% MHT and MFT components
if flag_HT_FWT_decomposition == 1

    % *** Heat ***
    ose.MHT_dia = MHT_mc_dia;
    ose.MHT_iso = MHT_mc_iso;
    ose.MHT_net = MHT_mc_net;
    ose.MHT_over = MHT_mc_over;
    ose.MHT_gyre = MHT_mc_gyre ;
    ose.MHT_netz = MHT_mc_netz;
    
    % *** FW ***
    ose.MFT_dia = MFT_mc_dia;
    ose.MFT_iso = MFT_mc_iso;
    ose.MFT_netz = MFT_mc_netz;
    ose.MFT_over = MFT_mc_over;
    ose.MFT_gyre = MFT_mc_gyre;
    ose.MFT_net = MFT_mc_net;



    if strcmp(flag_section,'all') && flag_test_AllWE== 1
        % *** West ***
            % *** Heat ***
        ose.MHT_dia_w = MHT_mc_dia_w;
        ose.MHT_iso_w = MHT_mc_iso_w;
        ose.MHT_net_w = MHT_mc_net_w;
        ose.MHT_over_w = MHT_mc_over_w;
        ose.MHT_gyre_w = MHT_mc_gyre_w;
        ose.MHT_netz_w = MHT_mc_netz_w;

            % *** FW ***
        ose.MFT_dia_w = MFT_mc_dia_w;
        ose.MFT_iso_w = MFT_mc_iso_w;
        ose.MFT_net_w = MFT_mc_net_w;
        ose.MFT_over_w = MFT_mc_over_w;
        ose.MFT_gyre_w = MFT_mc_gyre_w;
        ose.MFT_netz_w = MFT_mc_netz_w;
        
        
        % *** East ***
            % *** Heat ***
        ose.MHT_dia_e = MHT_mc_dia_e;
        ose.MHT_iso_e = MHT_mc_iso_e;
        ose.MHT_net_e = MHT_mc_net_e;
        ose.MHT_over_e = MHT_mc_over_e;
        ose.MHT_gyre_e = MHT_mc_gyre_e;
        ose.MHT_netz_e = MHT_mc_netz_e;


            % *** FW ***
        ose.MFT_dia_e = MFT_mc_dia_e;
        ose.MFT_iso_e = MFT_mc_iso_e;
        ose.MFT_net_e = MFT_mc_net_e;
        ose.MFT_over_e = MFT_mc_over_e;
        ose.MFT_gyre_e = MFT_mc_gyre_e;
        ose.MFT_netz_e = MFT_mc_netz_e;
        
    end
    
end
% =========================================================================
% =========================================================================


% *********************************
% *** Fluxes from MC iterations ***
% *********************************
ose.MC_runs = MC_runs; % for testing purpose

if strcmp(flag_section,'all') && flag_test_AllWE== 1
    ose.MC_runs_w = MC_runs_w;
    ose.MC_runs_e = MC_runs_e;
end


% ****************
% *** Metadata ***
% ****************
% Add sufficient notations to ensure reproducibility

% MC stopping criteria
ose.MC_STOP = [epsilon_MOC epsilon_MHT epsilon_MFT];

% OSNAP-Bering Strait boundary mean salinity
ose.OSNAP_BS_MEAN_SALINITY = s0;

% all input directories
if flag_regridded_profiles == 1
    ose.INPUT_FILES.moorings_regridded_profiles = data_mr_osnap;
elseif flag_regridded_profiles == 0
    ose.INPUT_FILES.moorings_raw = data_osnap_mr_daily;
end

ose.INPUT_FILES.adt = data_sect_ADT;
ose.INPUT_FILES.sst = data_sect_SST;
ose.INPUT_FILES.wind = data_sect_WIND;
if flag_test_AGV== 1
    ose.INPUT_FILES.agv = data_sect_AGV;
end
ose.INPUT_FILES.model_climatology = data_sect_model_clim;

if flag_OA == 1
    ose.INPUT_FILES.OA = data_sect_OA;
end

if flag_vbaro_mean== 1
    ose.INPUT_FILES.mean_baro = data_sect_mean_barotropic;
end

% matlab version
ose.MATLAB_FILE = 'Compute_fluxes_osnap_master_2021.m';
ose.MATLAB_VERSION = version;

% creation date
ose.CREATION_DATE = datestr(now,'yyyy-mm-dd HH:MM:SS');
ose.notation = 'Canadian shelf ADCP data complete on Aug 31 2021';



% SAVE
% output since Nov 15, 2021 should include IB5, new barotropic velocity
% (v20211111) and corrected wind (v20211111)
% 20211130 with corrected 53N Arrary also the OA is updated
file_out = [datestr(today,'yyyymmdd'),'_test_Fluxes_OSNAP_All_nMC_mm_fullOA_LSAB_updtM1235_updt53N_20140601-20200630.mat'];
disp(['+ OUTPUT TO ',dirout,file_out]);
save([dirout,file_out],'ose','-v7.3');



% CLEANING UP ...
delete TemporaryFile*mat


% return
% *********************************
%% calculate and output v_baro_mean
% *********************************
% This part can be used seperately
if flag_baro_corr== 1 && flag_vbaro_mean~= 1
    
    time = ose.time;
    ipt = find(time<= datenum(2020,7,26,0,0,0));  % GEOMAR mooring recovered;
%     ipt = find(time<= datenum(2016,3,28,0,0,0)); % End of the first 21month;
    time = time(ipt);
    vbaro = ose.v_baro(:,ipt);
    
    % output
    vmean = nanmean(vbaro,2); % mean
    vstd = nanstd(vbaro,0,2); % stdev
    
    % stderr
    vste = nan(size(vmean));
    
    for i_n = 1:length(ose.lon)
        [tau_i, dof_eff] = fun_DOFeff(vbaro(i_n,:),time);
        vste(i_n) = vstd(i_n)/sqrt(dof_eff);
    end
    
    % structure
    v_baro_mean.mean = vmean;
    v_baro_mean.sd = vstd;
    v_baro_mean.se = vste;
    v_baro_mean.lon = ose.lon;
    v_baro_mean.lat = ose.lat;
    v_baro_mean.TIME_PERIOD = [time(1) time(end)]; % added to show the time period over which the mean vbaro was derived
    v_baro_mean.CREATION_DATE = ose.CREATION_DATE;
    v_baro_mean.DATA_FILE = file_out;
    
    % save as
    disp(['+ OUTPUT v_baro_mean: ',dirio,file_out(1:16),'_v_baro_mean.mat']);
    save([dirout,file_out(1:13),'_v_baro_mean_M5corrected.mat'],'v_baro_mean');
%     disp(['+ OUTPUT v_baro_mean: ',dirio,'/',file_out(1:16),'v_baro_mean_first21mon.mat']);    
%     save([dirio,'/',file_out(1:16),'v_baro_mean_first21mon.mat'],'v_baro_mean');

end




disp('+' );
disp(['+ Now out Compute_fluxes_osnap_2021 :',datestr(now)]);
disp('+' );
disp('+' );

disp('+ now run compare_moc,+')
% compare_moc


