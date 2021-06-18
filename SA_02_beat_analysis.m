% This file analyzes the properties of each model in the population.
% It creates the matrix 'all_outputs' and the arrays 'output_names' and
% 'output_units'.

clear
close all
clc

%% Loading initial conditions
load SA_ICs_matrix_1000_s0p1

% all_ICs
% columns: N state variables
% rows: N trials

%% Parameters
% load matrix all_parameters (columns: N parameters, rows: N trials)
load SA_par_matrix_1000_s0p1 % sigma 0.1

[N_trials, N_par] = size(all_parameters);

%% Simulation parameters
% Protocol parameters
prot_index = 1;     % (1-5) DEFINE STIMULATION PROTOCOL INDEX HERE

prot_rate = 1;      % (Hz) used in prot 1, 3 and 4
prot_interval = 0; % (ms) used in prot 4
prot_vm = 0;      % (mV) used in prot 4

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
duration = 3e3; % (ms)
tspan = [0; duration];
options = odeset('RelTol',1e-6,'MaxStep',1,'Stats','off'); 

%% Outputs definition
% Outputs for logistic regression (1-5):
% 1) flag_abnormal_rep 2) flag_ead 3) flag_normal_AP 4) run_Vm_min 5) run_delta_Na
% Outputs for linear regression (6-22) - only for normal beat:
% 1) dVm_max 2) Vm_max 3) -Vm_min 4) AP_amp 5) APD90 6) APD70 7) APD50 8) APD30
% 9) Ca_max 10) Ca_min 11) CaT_amp 12) CaT_rise 13) CaT_decay_50 14) CaT_decay_63
% 15) Na_min 16) CaSR_max 17) CaSR_min

output_names = {'abn rep', 'ead', 'normal', 'min Em', 'delta [Na]',...
    'UV', 'AP peak', '|Em rest|', 'AP amp', 'APD90',...
    'APD70', 'APD50', 'APD30', 'syst [Ca]', 'diast [Ca]',...
    'CaT amp', 'CaT ttp', 'CaT t50', 'CaT tau', 'diast [Na]',...
    'CaSR max', 'CaSR min'};

output_units = {'-', '-', '-', 'mV', 'mM',...
    'mV/ms', 'mV', 'mV', 'mV', 'ms',...
    'ms', 'ms', 'ms', 'mM', 'mM',...
    'mM', 'ms', 'ms', 'ms', 'mM',...
    'mM', 'mM'};

N_outputs = length(output_names); % number of outputs of beat analysis

all_outputs = zeros(N_trials,N_outputs);

%% Run cycle
tic
parfor ii = 1:N_trials
%for ii = 1:30%:50
    X = sprintf('Run %d on %d',ii,N_trials); disp(X)
    y0 = all_ICs(ii,:);
    %y0 = real(all_ICs(ii,:));
    p_SA = all_parameters(ii,:); % 19 parameters

    [t,y] = ode15s(@morotti_et_al_ham_ina_ran_model_SA,tspan,y0,options,p,p_SA);
    
    time = t; % (ms)
    Vm = y(:,39); % (mV)
    Ca = y(:,38); % (mM)
    CaSR = y(:,31); % (mM)
    Na = y(:,34); % (mM)
    dVm_calc = (Vm(2:end)-Vm(1:end-1))./(t(2:end)-t(1:end-1));
    dVm = [dVm_calc; dVm_calc(end)]; % (mV/ms)
    period = 1000/prot_rate; % ms
    AP_index = 1; % w/ 1 first AP, otherwise last AP
    
    %ii
    outputs = function_beat_analysis_EAD(time,Vm,Ca,CaSR,Na,dVm,period,AP_index)%;
    all_outputs(ii,:) = outputs;
    
    figure%(100)
    subplot(3,1,1),hold on,plot(time,Vm)
    subplot(3,1,2),hold on,plot(time,Ca)
    subplot(3,1,3),hold on,plot(time,Na)
end

all_outputs
% columns: N outputs
% rows: N trials
toc

%% Saving
%save SA_output_matrix_control_1Hz_300s all_outputs output_names output_units
