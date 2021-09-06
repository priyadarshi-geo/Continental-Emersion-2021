clear; clc; clf; format shorte;
% Date: July, 2021
% Author: Priydarshi Chowdhury (Monash University)
% >> This code is a part of PNAS manuscript (Chowdhury et al.)
% >> Contact: pchowdhury59@gmail.com / priyadarshi.chowdhury@monash.edu

%% User defined variables
% variables that are kept constant for all simulations
rho_a           = 3330;  % kg/m3; asthenosphere density
rho_oc          = 3000;  % kg/m3; oceanic crust density
rho_olm         = 3320;  % kg/m3; oceanic litho. mantle density
rho_cc_ini      = 2980;  % in kg/m3; initial cratonic crust density         [see SI Appendix Methods for the choice of initial & final
rho_cc_final    = 2850;  % in kg/m3; final cratonic crust density            cratonic crust density at 3.5 and 3.1 Ga respec.]
CL_thick_ini    = 80;    % in km; initial cratonic litho. thickness         [at 3.5 Ga; CLM thick. is calc. using iniial cratonic crust thick and this litho. thick.]
OL_thick        = 80;    % in km; oceanic litho. thickness                  [OLM thick. is calc. using oceanic crust thick and this litho. thick.; kept constant over time]
present_sl      = 3.68;  % km; present-day sea-level                        [see SI Appendix Methods; taken from Eatkins & Sharman, 2010, NOAA Data]
archean_sl1     = 4.41;  % km; Archean sea-level w. gradual cont. growth    [see SI Appendix Methods; taken from Rosas & Korenaga, 2021, Nat. Geosci.]
archean_sl2     = 4.68;  % km; Archean sea-level w. rapid cont. growth      [see SI Appendix Methods; taken from Rosas & Korenaga, 2021, Nat. Geosci.]

% variables for which sensitivity of craton elevation is analyzed
rho_clm         = 3305;  % in kg/m3; CLM density                            [varied between: 3330, 3320, 3310, 3305, 3300, 3290, 3280]
h_oc            = 20;    % in km; oceanic crust thickness; 20/25/30 km      [varied between: 20, 25, 30]
clm_thick_max   = 133;   % in km; maximum CLM thickness allowed             [see SI Appendix Methods; varied between 133-160 km; 133 km when growth rate is 0.4x4kt & 160 km when growth rate is 4kt]
clm_growthfac   = 0.4;   % multiplied to the conductive cooling eq.;        [see SI Appendix Methods; varied between 0.4 and 1]
% NOTE: Initial thickness of the whole oceanic (OL) & craton (CL) lithos., 
% it is used to calculate the thickness of the oceanic (OLM) & continental
% (CLM) lithospheric mantle dynamcially within the code once the initial 
% oceanic & continental crustal thicknesses are known
h_olm           = OL_thick - h_oc;          % km; OLM thickness

%% variables for which parameterized time-evolution is known 
% i.e., cratonic crustal thickness and density, & CLM thickness 
time = 3470:-10:3100;     % in Ma; age duration of our interest

% Forming the crustal thickness vs. age array
% change the last term to modulate the curve
% 5424--> fit to the mean thickness trend 
% 5427--> fit to trace element thickness trend
% 5421--> fit to La/Yb thickness trend
h_cc = -4.995e-4.*time.*time + 3.306.*time - 5424;
max_h_cc = max(h_cc);   % to create a smooth profile
for i=1:length(time)
    if time(i)<3300, h_cc(i)= max_h_cc; end
end

% Forming the crustal density vs. age array (a poplynomial function)
% NOTE: calculation has to be done with time in Ga to get accurate values; 
% 'time' in Ma; time values are needed to be factored by 1000
rho_cc = 596.18.*time/1e3.*time/1e3-3565.58.*time/1e3+8173.967;
if rho_cc < rho_cc_final, rho_cc = rho_cc_final; end


% Forming CLM thickness vs age array
h_clm_ini       = CL_thick_ini - h_cc(1);   % km; initial CLM thickness
for i=1:length(time)
    if i==1, h_clm(i) = h_clm_ini; 
    else
        h_clm(i) = sqrt(4 * 1e-6 * (abs(time(i)-time(1))*1e6*365*24*3600));
        h_clm(i) = clm_growthfac * h_clm(i)/1e3 + h_clm_ini;
    end
    if h_clm(i) > clm_thick_max, h_clm(i) = clm_thick_max; end
end

%% Calculating elevation for a fixed CLM-AM density contrast; see L.#20
% calculating crustal elevation directly relative to the seafloor
for i=1:length(time)
    elev(i) = (h_cc(i)*(rho_a-rho_cc(i)) +...
                        h_clm(i)*(rho_a-rho_clm) +...
                            h_oc*(rho_oc-rho_a) +...
                                h_olm*(rho_olm-rho_a))/rho_a;
end

%% calcualting elevation for different CLM-AM density contrast (0-50 kg/m3)
drho_LAB = 0:10:50; % in kg/m3
% calculating crustal elevation directly relative to the seafloor
for k=1:length(drho_LAB)
for i=1:length(time)
    elev2(i,k) = (h_cc(i)*(rho_a-rho_cc(i)) +...
                            h_clm(i)*(drho_LAB(k)) + ...
                                h_oc*(rho_oc-rho_a) +...
                                    h_olm*(rho_olm-rho_a))/rho_a;
end
end
%display ([time',dH_cc',elev',elev_direct',dH_cc2(:,1),elev2(:,1)]);

%% Calculations for plotting
% LaYb based crustal thickness vs age relation as estimated from TTGs;
tt          = 3470:-10:3240;
LaYb_h_cc   = -4.4304e-4.*tt.*tt + 2.9245.*tt - 4776.3;

% crustal thickness vs age relation as obtained from TTG REE patterns
TTG_melt_age= [3470,3400,3399.5,3350,3349.5,3240];
TTG_Press   = [0.93,0.93,1.1,1.1,1.25,1.25];
% converting pressure to depth
for j=1:length(TTG_melt_age)
    % calculating density at time-instants given in TTG_melt_age
    rho_cc1 = 596.18*TTG_melt_age(j)/1e3*TTG_melt_age(j)/1e3...
        - 3565.58*TTG_melt_age(j)/1e3 + 8173.967;
    if rho_cc1 < rho_cc_final, rho_cc1 = rho_cc_final; end
    
    TTG_depth(j) = ((TTG_Press(j)*1e9)/(rho_cc1*9.81))/1e3;
end

% preparing sea level arrays over time for plotting
sealevel_p  = zeros(length(time),1) + present_sl;
sealevel_a1 = zeros(length(time),1) + archean_sl1; 
sealevel_a2 = zeros(length(time),1) + archean_sl2; 

%% Plotting results
figure(1)
% cratonic crustal thickness vs. age
subplot(2,2,1)
plot(time./1e3,h_cc,'-k'); hold on;                 % user defined thick
plot(tt./1e3,LaYb_h_cc,'-g'); hold on;              % LaYb based thick
plot(TTG_melt_age./1e3,TTG_depth,'-b'); hold on;    % TTG REE based thick
title('Cratonic crustal thickness vs Age');
xlim([3.0 3.6]); ylim([20 60]);
hold off; set(gca,'TickDir','out'); grid on;
xlabel('Age (Ga)'); ylabel('Cratonic crustal thickness (km)'); 
% set(gca,'xtick',[],'ytick',[]);
set(gca,'layer','top','linewidth',0.5);

% cratonic crustal density vs. age
subplot(2,2,2)
plot(time./1e3,rho_cc,'-k'); hold on;
plot([3.47,3.30,3.10],[2980,2900,2850],'ok'); hold on;
title('Continental crust density vs Age');
xlim([3.0 3.6]); ylim([2800 3000]);
hold off; set(gca,'TickDir','out'); grid on;
xlabel('Age (Ga)'); ylabel('Continental crust density (kg/m3)'); 
% set(gca,'xtick',[],'ytick',[]);
set(gca,'layer','top','linewidth',0.5);

% CLM thickness vs age
subplot(2,2,3)
plot(time./1e3,h_clm,'-k'); hold on;
title('CLM thickness vs Age');
xlim([3.0 3.6]); ylim([20 200]);
hold off; set(gca,'TickDir','out'); grid on;
xlabel('Age (Ga)'); ylabel('CLM thickness (km)'); 
% set(gca,'xtick',[],'ytick',[]);
set(gca,'layer','top','linewidth',0.5);

% CLM thickness vs age
subplot(2,2,4)
plot(time./1e3,elev,'-k'); hold on;
title('Craton elevation vs Age');
xlim([3.0 3.6]); ylim([0 6]);
hold off; set(gca,'TickDir','out'); grid on;
xlabel('Age (Ga)'); ylabel('Craton elevation from seafloor (km)'); 
% set(gca,'xtick',[],'ytick',[]);
set(gca,'layer','top','linewidth',0.5);


figure(2)
for k=1:length(drho_LAB)
    plot(time./1e3,elev2(:,k),'-','Color',[0.65 0.65 0.65]); hold on;
end
plot(time./1e3,elev,'-r'); hold on;
plot(time./1e3,sealevel_p,'--b'); hold on;
plot(time./1e3,sealevel_a1,'-b'); hold on;
plot(time./1e3,sealevel_a2,'-b'); hold on;

max_hpresent  = round(max(elev)-present_sl,2);
max_harchean1 = round(max(elev)-archean_sl1,2);
max_harchean2 = round(max(elev)-archean_sl2,2);

title({['Continental elevation (relative to seafloor) vs Age']
    ['For Oceanic crust thickness = ',num2str(h_oc),' km']
    []
    [' Max. elev. above modern sealevel (',num2str(present_sl),' km) = ',...
                            num2str(max_hpresent*1e3),' m']
    [' Max. elev. above Archean sealevel-1 (',num2str(archean_sl1),' km) = ',...
                            num2str(max_harchean1*1e3),' m']
    [' Max. elev. above Archean sealevel-2 (',num2str(archean_sl2),' km) = ',...
                            num2str(max_harchean2*1e3),' m']});
xlim([3.0 3.5]); ylim([0 8]);
hold off; set(gca,'TickDir','out'); grid on;
xlabel('Age (Ga)'); 
ylabel('Continental elevation (km)'); 
% set(gca,'xtick',[],'ytick',[]);
set(gca,'layer','top','linewidth',0.5);
%%