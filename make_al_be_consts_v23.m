function out = make_al_be_consts_v23()

% This function creates and saves a structure with relevant constants and
% external data for the Al-Be exposure age and erosion rate calculators.  
%
% Syntax: make_al_be_consts_v23
% (no arguments)
%
% See the text documentation or the code itself for what the constants
% actually are. 
%
% This update June 2016. Changes include new muon interaction
% cross-sections and new default production rates. 
%
% Written by Greg Balco -- Berkeley Geochronology Center
% balcs@u.washington.edu -- balcs@bgc.org
% February, 2008
% Part of the online exposure age calculators at: 
%      http://hess.ess.washington.edu/
%
% Copyright 2001-2007, University of Washington
% 2008-, Greg Balco
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This software is under development and is not licensed for use or distribution.  
% Insert license info here. 

al_be_consts.version = '2.3';
al_be_consts.prepdate = fix(clock);

% Be-10 decay constant

al_be_consts.l10 = -log(0.5)./1.387e6; % Chmeleff/Korschinek value

dldt = -log(0.5).*(1.387e6^-2);
al_be_consts.dell10 = sqrt((dldt.*0.012e6)^2); % Chmeleff/Korschinek value

% Al-26 decay constant -- value compatible with Nishiizumi standards
% lambda = 9.83e-7 --> t(1/2) = 0.705e6 yr
% See Nishiizumi (2004) for details.

al_be_consts.l26 = 9.83e-7;

al_be_consts.dell26 = 2.5e-8;

% Effective attenuation length for spallation in rock
% Commonly accepted value: 160 g/cm2
% For discussion see Gosse and Phillips (2000)

al_be_consts.Lsp = 160;

% Fsp - fraction of total production by spallation rather than muons
% For use with Lal/Stone scaling scheme in exposure age calculation
% For details, see Stone (2000)
% This aspect of Stone(2000) is de-emphasized in version 2. These constants
% are only in use for historical comparisons and quick initial guesses for 
% the exposure age and erosion rate solvers. 

% Note that they are probably incorrect WRT Be-10 restandardization. Don't
% use these for anything important. 
al_be_consts.Fsp10 = 0.978;
al_be_consts.Fsp26 = 0.974;

% -------------- standardization ---------------------

% Be-10 standards conversion lookup table

al_be_consts.be_stds_names = strvcat('07KNSTD','KNSTD','NIST_Certified','LLNL31000','LLNL10000',...
    'LLNL3000','LLNL1000','LLNL300','NIST_30000','NIST_30200','NIST_30300','NIST_30600','NIST_27900',...
    'S555','S2007','BEST433','BEST433N','S555N','S2007N','STD11');
al_be_consts.be_stds_cfs = [1.0000 0.9042 1.0425 0.8761 0.9042 0.8644 0.9313 0.8562 0.9313 0.9251 0.9221 0.9130 1 0.9124 0.9124 0.9124 1 1 1 1]';

% Same for Al-26. A zero placeholder is also allowed. 

al_be_consts.al_stds_names = strvcat('KNSTD','ZAL94','SMAL11','0','ZAL94N','ASTER','Z92-0222');
al_be_consts.al_stds_cfs = [1.0000 0.9134 1.021 1.0000 1 1.021 1]';

% ----------- end standardization ---------------------

% ----------- reference production rates ----------------

% Default reference production rates at SLHL for spallation. These are 
% now calibrated to the 'primary' data set of Borchers and others (2016). 
%
% The uncertainties are just the scatter for the secondary data set of
% Borchers. If we were to use the scatter of the primary data set, it would
% be less for St and Lm (5%) and more for the others (15%). 

% Be-10 production rates. 

al_be_consts.P10_ref_St = 4.10;
al_be_consts.delP10_ref_St = 0.35;

al_be_consts.P10_ref_Lm = 3.99; 
al_be_consts.delP10_ref_Lm = 0.37;

al_be_consts.P10_ref_De = 3.93;
al_be_consts.delP10_ref_De = 0.47;

al_be_consts.P10_ref_Du = 3.95;
al_be_consts.delP10_ref_Du = 0.46;

al_be_consts.P10_ref_Li = 4.25;
al_be_consts.delP10_ref_Li = 0.48;

% Al-26 production rates are derived from Be-10 production rates 

R2610 = 6.1.*1.106; % Update assumed production ratio

al_be_consts.P26_ref_St = al_be_consts.P10_ref_St*R2610;
al_be_consts.delP26_ref_St = al_be_consts.delP10_ref_St*R2610;
al_be_consts.P26_ref_Lm = al_be_consts.P10_ref_Lm*R2610; 
al_be_consts.delP26_ref_Lm = al_be_consts.delP10_ref_Lm*R2610;
al_be_consts.P26_ref_De = al_be_consts.P10_ref_De*R2610;
al_be_consts.delP26_ref_De = al_be_consts.delP10_ref_De*R2610;
al_be_consts.P26_ref_Du = al_be_consts.P10_ref_Du*R2610;
al_be_consts.delP26_ref_Du = al_be_consts.delP10_ref_Du*R2610;
al_be_consts.P26_ref_Li = al_be_consts.P10_ref_Li*R2610;
al_be_consts.delP26_ref_Li = al_be_consts.delP10_ref_Li*R2610;


% ---------- done with reference production rates --------------

% ---------- muon parameters -----------------------------------

% Muon interaction cross-sections. All follow Heisinger (2002a,b).
% Note that the energy-dependence-of-muon-interaction-cross-section
% exponent alpha is treated as model-dependent -- it's internal to 
% P_mu_total.m and can't be passed.  

al_be_consts.Natoms10 = 2.006e22;
al_be_consts.Natoms26 = 1.003e22;

% Be-10 and Al-26 muon interaction cross-sections. These are derived from
% fitting P_mu_total.m to Beacon Heights data from Borchers and others
% (2016). 

% Be-10 

al_be_consts.k_neg10 = 0.704 * 0.1828 * 0.00157;
al_be_consts.sigma190_10 = 37.8e-30;

% Al-26 

al_be_consts.k_neg26 = 0.296 * 0.6559 * 0.0118;
al_be_consts.sigma190_26 = 521e-30;

% Uncertainties are nominal. Not used here. 
al_be_consts.delsigma190_10 = al_be_consts.sigma190_10.*0.1;
al_be_consts.delk_neg10 = al_be_consts.k_neg10.*0.1;
al_be_consts.delsigma190_26 = al_be_consts.sigma190_26.*0.1;
al_be_consts.delk_neg26 = al_be_consts.k_neg26.*0.1;

% --------- done with muon parameters ---------------------------

% ----------- paleomagnetic data --------------------------------

% Paleomagnetic records for use in time-dependent production rate schemes
% Derived from Nat Lifton's compilation of paleomagnetic data from
% various sources. See Lifton et al. (2006) and Pigati and Lifton (2005).

% Load the magnetic field data
load PMag_Mar07
% Relative dipole moment and time vector
al_be_consts.M = MM0; 
al_be_consts.t_M = t_M; 
% These start at 7500 yr -- time slices are 7500,8500,9500,10500,11500
% in order to use data from Yang et al; subsequent time slices are 
% 12000:1000:800000 for SINT800 data; final two time points are 801000
% and 1e7. Note that Inf is no longer allowed by interp1. 

al_be_consts.t_M(end) = 1e7;

% Cutoff rigidity blocks for past 6900 yr. 
% TTRc and IHRC are lon x lat x time blocks of Rc values for the past 
% 6900 years.
% Both are derived by Nat Lifton from the magnetic field reconstructions of
% Korte and Constable. 
% TTRC has cutoff rigidity obtained by trajectory tracing -- these are for
% the Lifton and Desilets scaling factors. IHRc has cutoff rigidity
% obtained by finding magnetic inclination and horizontal field strength
% from the field model, then applying Equation 2 of Dunai(2001). 
al_be_consts.TTRc = TTRc; % data block
al_be_consts.IHRc = IHRc; % data block
al_be_consts.lat_Rc = lat_Rc; % lat and lon indices for Rc data block
al_be_consts.lon_Rc = lon_Rc;
al_be_consts.t_Rc = t_Rc; % time vector for Rc data block

% Effective pole positions and field strengths inferred from K and C field
% reconstructions for last 7000 yr. These are used in the
% paleomagnetically-corrected implementation of the Lal SF. They are for
% the same times as the RC slices in the data block above. Again,
% generated by Nat Lifton -- hence KCL = Korte-Constable-Lifton. 

al_be_consts.MM0_KCL = MM0_KCL;
al_be_consts.lat_pp_KCL = lat_pp_KCL;
al_be_consts.lon_pp_KCL = lon_pp_KCL;

% Solar variability from Lifton et al. 2005
% Has been averaged and resampled to the same time slices as everything
% else. 

al_be_consts.S = S; 
al_be_consts.SInf = 0.95; % Long-term mean S value;

% --------------- end paleomagnetic data ------------------------

% Finish up

save al_be_consts_v23 al_be_consts

disp(['Al-Be constants version ' al_be_consts.version]);
disp('Saved'); 


