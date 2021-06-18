% Generate ICs in the population

clear
close all
clc

%% Initial conditions
load yf_ham_ina_ran_1Hz % steady-state 1-Hz, w/out Acetylcholine

y0 = yfinal;
N_state_vars = length(y0);

%% Parameters
% load matrix all_parameters (columns: N parameters, rows: N trials)
load SA_par_matrix_1000_s0p1

[N_trials N_par] = size(all_parameters);

%% Input parameters
prot_index = 1; % 1) 'pace_cc';

prot_rate = 1; % (Hz) used only with 'pace_cc' and 'v_step'
prot_interval = 0; % (ms) used only for "recovery" protocols (and ERP) (and 17)
prot_vm = 0; % (mV)

% Drug parameters
drug_index = 0; drug_conc = 0;               % Drug Free
%drug_index = 1; drug_conc = 5 * (1E-6);    % Ranolazine (M)

% Other experimental conditions
exp_Temp = 310; % [K]
exp_Nao = 140; % [Na]o 130 mM in optimization, 140 otherwise

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
%p_SA = ones(1,19);

duration = 300e3;
tspan = [0 duration];
options = odeset('RelTol',1e-5,'MaxStep',1,'Stats','on'); 

%% Run cycle
all_ICs = zeros(N_trials,N_state_vars);

tic
parfor ii=1:N_trials
    X = sprintf('Run %d on %d',ii,N_trials); disp(X)
    p_SA = all_parameters(ii,:); % 19 parameters
    
    [t,y] = ode15s(@morotti_et_al_ham_ina_ran_model_SA,tspan,y0,options,p,p_SA);
    all_ICs(ii,:) = y(end,:);
end
toc

all_ICs;
% columns: N state variables
% rows: N trials

%% Saving
%save SA_ICs_matrix_1000_s0p1 all_ICs % Control
