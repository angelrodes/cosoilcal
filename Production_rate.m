function out = Production_rate(site_lat,site_lon,site_elv,shielding,nuclide)

%   Calculates surface production rates for cosmogenic isotopes,
%   apparent attenuation lengths from surface to a few meters of fast
%   and slow muons, and corresponding pressure for a given latitude,
%   longitude and elevation. It uses the following CRUNUS 2.3 code: 
%   NCEPatm_2.m to calculate pressure,
%   antatm.m to calculate pressure in Antarctica,
%   al_be_consts_v23.mat for constants, 
%   stone2000 for spallation production rates,
%   P_mu_total.m for muon production rates.
%   Calculated 10Be production rates are scaled for other isotopes
%   based on published ratios.

%       Inputs:
%       site_lat: latitude (DD). Southern hemisphere is negative.
%       site_lon: longitude (DD). Western hemisphere is negative.
%       site_elv: elevation (m).
%       shielding: topographic sheilding factor between 0 and 1.
%       Nuclide: mass of the cosmogenic nuclide

%       Outputs:
%       Surface spallation production rate (at/g/a) P(1)
%       Surface fast muons production rate (at/g/a) P(2)
%       Surface negative muons production rate (at/g/a) P(3)
%       Decay constant (l)
%       Apparent spallation attenuation length (g/cm2) L(1)
%       Apparent fast muons attenuation length (g/cm2) L(2) bsed on a large
%        interpolation from data generated in CRONUS 2.3
%        Lfm = 900 + 1310*exp(-0.0004048*elv)
%       Apparent slow muons attenuation length (g/cm2) L(3) bsed on a large
%        interpolation from data generated in CRONUS 2.3
%        Lsm = 500 + 823*exp(-0.0005567*elv)
%       Pressure (hPa) (h)

% Written by Ángel Rodés -- SUERC
% angelrodes@gmail.com
% March, 2018

if site_lat<(-55) 
    h=antatm(site_elv);
else
    h=NCEPatm_2(site_lat,site_lon,site_elv);
end
out.h=h; %Pressure

% run make_al_be_consts_v23.m if error
make_al_be_consts_v23
load al_be_consts_v23.mat
consts=al_be_consts;


% 10Be production rate for spallation
P_ref_St = consts.P10_ref_St;
% constants structure for muon production rate
mconsts.Natoms = consts.Natoms10;
mconsts.sigma190 = consts.sigma190_10;
mconsts.delsigma190 = consts.delsigma190_10; % not used
mconsts.k_neg = consts.k_neg10;
mconsts.delk_neg = consts.delk_neg10; % not used
P_St = P_ref_St * stone2000(site_lat,h,1)*shielding;
dflag='yes';

% 10Be muon calculations
surface=P_mu_total(0,h,mconsts,dflag);

% Attenuation lengths (spallation)
out.L(1)=consts.Lsp;
% Attenuation lengths )muons= interpolated from global data using P_mu_total
out.L(2)=900 + 1310*exp(-0.0004048*site_elv); % fast muons
out.L(3)=500 + 823*exp(-0.0005567*site_elv); % negative muons

% Production rates
if nuclide==10 % from CRONUS 2.3 (Balco, 2008)
    out.P(1)=P_St;
    out.P(2)=surface.P_fast;
    out.P(3)=surface.P_neg;
    out.l=consts.l10;
    out.label='Be-10(SiO2)';
    out.label2='^{10}Be';
elseif nuclide==26 % P26Al(SiO2)/P10Be(SiO2) ratios from CRONUS 2.3 (Balco, 2008)
    out.P(1)=P_St*6.7466;
    out.P(2)=surface.P_fast*6.8915;
    out.P(3)=surface.P_neg*11.3387;
    out.l=consts.l26;
    out.label='Al-26(SiO2)';
    out.label2='^{26}Al';
elseif nuclide==21 % P21Ne(SiO2)/P10Be(SiO2) ratios from Balco & Shuster (2009)
    out.P(1)=P_St*4.10;
    out.P(2)=surface.P_fast*4.10;
    out.P(3)=surface.P_neg*0;
    out.l=0;
    out.label='Ne-21(SiO2)';
    out.label2='^{21}Ne';
elseif nuclide==3 % P3He(px)/P10Be(qtz) ratios from Blard et al. (2013)
    out.P(1)=P_St*33.3;
    out.P(2)=surface.P_fast*33.3;
    out.P(3)=surface.P_neg*33.3;
    out.l=0;
    out.label='He-3(px)';
    out.label2='^{3}He';
elseif nuclide==36 % P36Cl(CaCO3)/P10Be(SiO2) ratios from Heisinger & Nolte (2000)
    out.P(1)=P_St*4.98;
    out.P(2)=surface.P_fast*3.15;
    out.P(3)=surface.P_neg*22.5;
    out.l=log(2)/301e3;
    out.label='Cl-36(CaCO3)';
    out.label2='^{36}Cl';
elseif nuclide==14 % P14C(SiO2)/P10Be(SiO2) ratios from Heisinger & Nolte (2000)
    out.P(1)=P_St*2.71;
    out.P(2)=surface.P_fast*4.70;
    out.P(3)=surface.P_neg*30.8;
    out.l=log(2)/5730;
    out.label='C-14(SiO2)';
    out.label2='^{14}C';
else
    out.P=[0,0,0];
    out.l=0;
    out.label='Unknown-isotope';
    out.label2='Unknown-isotope';
end


