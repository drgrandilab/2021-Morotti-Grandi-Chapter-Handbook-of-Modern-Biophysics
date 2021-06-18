%% Morotti et al. Human Atrial Myocyte model - main file

% This file loads initial conditions, defines simulation parameters, calls
% the ode solver, and plots simulation results.

% Please, cite the following paper when using this model:
% Morotti S, McCulloch AD, Bers DM, Edwards AG, Grandi E. (2015).
% Atrial-selective targeting of arrhythmogenic phase-3 early
% afterdepolarizations in human myocytes. J Mol Cell Cardiol.
% doi: 10.1016/j.yjmcc.2015.07.030

clear
close all
clc

%% Initial conditions
load yf_ham_ina_ran_1Hz % steady-state 1-Hz, w/out Acetylcholine
% load yf_ham_ina_ran_ACh0p1_1Hz % steady-state 1-Hz, w/ Acetylcholine

y0 = yfinal;

%% Simulation parameters
% Protocol parameters
prot_index = 1;     % (1-5) DEFINE STIMULATION PROTOCOL INDEX HERE
% INDEX:            Protocol description:
% 1) 'pace_cc'      current clamp (current injection at rate 'prot_rate')
% 2) 'pace_cc_ead'  current clamp (20-s at 10-Hz, then 10-s at 1_Hz)
% 3) 'vc_tri'       train of triangular voltage pulses (rate 'prot_rate')
% 4) 'vc_step'      train of square voltage pulses (rate 'prot_rate', duration 'prot_interval', voltage 'prot_vm')
% 5) 'caffeine'     caffeine protocol (administration at 5 s)
% 6) 'pace_cc_dad'  current clamp (20-s at 3-Hz, then pause)

prot_rate = 1;      % (Hz) used in prot 1, 3 and 4
prot_interval = 25; % (ms) used in prot 4
prot_vm = -10;      % (mV) used in prot 4

% Ranolazine parameters
ran_flag = 0; % (boolean - 0 for no drug, 1 for drug)
if ran_flag == 1
    drug_index = 1; drug_conc = 10 * (1E-6); % (M) DEFINE [RAN] HERE
else
    drug_index = 0; drug_conc = 0; % Drug Free
end

% Other experimental conditions
exp_Temp = 310; % temperature (300 or 310 K)
exp_Nao = 140; % extracellular [Na] (mM)

% Isoproterenol administration
exp_ISO = 0; % (boolean - 0 for no ISO, 1 for ISO)

% Acetylcholine administration
exp_Ach = 0; % (boolean - 0.1 uM if exp_Ach = 1)

% Parameter array for passing nondefault conditions
prot_par = [prot_index prot_rate prot_interval prot_vm];    % 1 2 3 4
drug_par = [drug_index drug_conc];                          % 5 6
exp_par = [exp_Temp exp_Nao exp_ISO exp_Ach];               % 7 8 9 10
p = [prot_par drug_par exp_par]; 

% Sensitivity analysis parameters
p_SA = ones(1,19);

% Simulation duration
duration = 3 * 1e3; % (ms) DEFINE SIMULATION DURATION HERE

%% Single Run Simulation
tic
tspan = [0; duration];
options = odeset('RelTol',1e-6,'MaxStep',1,'Stats','off'); 
[t,y] = ode15s(@morotti_et_al_ham_ina_ran_model_SA,tspan,y0,options,p,p_SA);
yfinal = y(end,:);
toc

%% Final conditions
% Control
% save yf_ham_ina_ran_1Hz yfinal % steady-state 1-Hz, w/out Acetylcholine

% Acetylcholine 0.1 uM
% save yf_ham_ina_ran_ACh0p1_1Hz yfinal

%% Plot 
figure
set(gcf,'color','w','Position',[100 100 600 800])
subplot(4,1,1); hold on, plot(t/1000,y(:,39))
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('E_m (mV)')
subplot(4,1,2); hold on, plot(t/1000,y(:,38)*1000)
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('[Ca^2^+]_i (\muM)')
subplot(4,1,3); hold on, plot(t/1000,y(:,31))
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('Free [Ca^2^+]_S_R (mM)')
subplot(4,1,4); hold on, plot(t/1000,y(:,34))
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('[Na^+]_i (mM)')
xlabel('Time (s)')

Vm = y(:,39); % (mV)
dVm_calc = (Vm(2:end)-Vm(1:end-1))./(t(2:end)-t(1:end-1));
dVm = [dVm_calc; dVm_calc(end)]; % (mV/ms)
%figure,plot(t,dVm)

%% Plot Currents
currents = morotti_et_al_ham_ina_ran_calcCurrents(t,y,p,p_SA);

figure
set(gcf,'color','w','Position',[100 100 600 800])
subplot(4,1,1); hold on, plot(t/1000,y(:,39),'k')
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('E_m (mV)')
subplot(4,1,2); hold on, plot(t/1000,y(:,38)*1000,'b')
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('[Ca^2^+]_i (\muM)')
subplot(4,1,3); hold on, plot(t/1000,currents(:,3),'b',t/1000,currents(:,4),'g',t/1000,currents(:,1),'r')
set(gca,'box','off','tickdir','out','fontsize',12)
legend('I_C_a','I_N_C_X','I_N_a'), ylabel('Current (A/F)')
subplot(4,1,4); hold on, plot(t/1000,y(:,66),'k',t/1000,y(:,70),'k--')
set(gca,'box','off','tickdir','out','fontsize',12)
legend('O','LO'), ylabel('P_O (%)')
xlabel('Time (s)')
