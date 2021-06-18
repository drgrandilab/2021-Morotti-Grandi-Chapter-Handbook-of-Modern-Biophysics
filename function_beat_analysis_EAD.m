function outputs = function_beat_analysis_DAD(time,Vm,Ca,CaSR,Na,dVm,period,AP_index)
% Beat analysis first or last AP (index = 1 or 2)

% Outputs for logistic regression (1-5):
% 1) flag_abnormal_rep 2) flag_ead 3) flag_normal_AP 4) run_Vm_min 5) run_delta_Na
% Outputs for linear regression (6-22) - only for normal beat:
% 1) dVm_max 2) Vm_max 3) -Vm_min 4) AP_amp 5) APD90 6) APD70 7) APD50 8) APD30
% 9) Ca_max 10) Ca_min 11) CaT_amp 12) CaT_rise 13) CaT_decay_50 14) CaT_decay_63
% 15) Na_min 16) CaSR_max 17) CaSR_min
    
% Check Na stability
run_delta_Na = Na(end)-Na(1);

% Check for impaired repolarization
flag_abnormal_rep = 0;

run_Vm_min = min(Vm);
if run_Vm_min > -50
    flag_abnormal_rep = 1;
else
    for i = 1:time(end)/period
        ti_roi = find(time>period*(i-1)); ti_index = ti_roi(1);
        te_roi = find(time<period*i); te_index = te_roi(end);
       
        time_roi = time(ti_index:te_index);
        Vm_roi = Vm(ti_index:te_index);
        %Ca_roi = Ca(ti_index:te_index);
        %CaSR_roi = CaSR(ti_index:te_index);
        %Na_roi = Na(ti_index:te_index);
        %dVm_roi = dVm(ti_index:te_index);
        
        Vm_end_array(i) = Vm_roi(end);
    end
    if max(Vm_end_array) > -50
        flag_abnormal_rep = 1;
    end
end

% Check for EADs
flag_ead = 0;
if flag_abnormal_rep == 0
    for i = 1:time(end)/period
        ti_roi = find(time>period*(i-1)); ti_index = ti_roi(1);
        te_roi = find(time>period*(i-1)+500); te_index = te_roi(1);

        time_roi = time(ti_index:te_index);
        Vm_roi = Vm(ti_index:te_index);
        %Ca_roi = Ca(ti_index:te_index);
        %CaSR_roi = CaSR(ti_index:te_index);
        %Na_roi = Na(ti_index:te_index);
        dVm_roi = dVm(ti_index:te_index);
        dVm_dVm = dVm_roi(1:end-1).*dVm_roi(2:end);
        
        [max_dVm idx_dVm] = max(dVm_roi);
        dVm_0 = find(dVm_dVm(idx_dVm:end)<0);
        dVm_0_count = length(dVm_0);

        if dVm_0_count > 1
            ead_array(i) = 1;
        else
            ead_array(i) = 0;
        end
    end
    if sum(ead_array)>0
        flag_ead = 1;
    end
end

% Beat analysis
if flag_abnormal_rep == 1 || flag_ead == 1
    flag_normal_AP = 0;
    outputs = [flag_abnormal_rep flag_ead flag_normal_AP run_Vm_min run_delta_Na...
        zeros(1,17)];
else
    flag_normal_AP = 1;
    if AP_index == 1
        t_in = 0; t_fin = t_in+period-5;
    else
        t_in = time(end)-period; t_fin = t_in+period-5;
    end
    t_in_roi = find(time>t_in); t_in_index = t_in_roi(1);%-1;
    t_fin_roi = find(time>t_fin); t_fin_index = t_fin_roi(1);

    %% Em
    [dVm_max, index] = max(dVm(t_in_index:t_fin_index)); index = index-1; % max slope

    [Vm_max, index_max] = max(Vm(t_in_index:t_fin_index)); index_max = index_max-1; % peak
    Vm_min = min(Vm(t_in_index:t_fin_index)); % Em resting
    AP_amp = Vm_max-Vm_min;

    APfraction = 0.9;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    APD90 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);

    APfraction = 0.7;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    APD70 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);

    APfraction = 0.5;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    APD50 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);

    APfraction = 0.3;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    APD30 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);

    %% Ca
    Ca_min = min(Ca(t_in_index:t_fin_index)); % diast Ca

    [Ca_max, index_ca] = max(Ca(t_in_index:t_fin_index)); index_ca = index_ca-1; % peak CaT
    CaT_amp = Ca_max-Ca_min;

    CaT_rise = time(t_in_index+index_ca)-time(t_in_index);

    Cafraction = 0.5;
    Ca_roi = find(Ca(t_in_index+index_ca:t_fin_index)<Ca_max-Cafraction*CaT_amp)-1;
    CaT_decay_50 = time(t_in_index+index_ca+Ca_roi(1))-time(t_in_index+index_ca);

    Cafraction = 0.632;
    Ca_roi = find(Ca(t_in_index+index_ca:t_fin_index)<Ca_max-Cafraction*CaT_amp)-1;
    CaT_decay_63 = time(t_in_index+index_ca+Ca_roi(1))-time(t_in_index+index_ca);
    
    %% Na
    [Na_min, index_na] = min(Na(t_in_index:t_fin_index)); index_na = index_na-1; % diast Na
    
    %% CaSR
    CaSR_max = max(CaSR(t_in_index:t_fin_index)); % diast CaSR
    CaSR_min = min(CaSR(t_in_index:t_fin_index)); % syst CaSR
    
    %% Collect all outputs
    outputs_beat = [dVm_max Vm_max -Vm_min AP_amp APD90 APD70 APD50 APD30...
        Ca_max Ca_min CaT_amp CaT_rise CaT_decay_50 CaT_decay_63 Na_min...
        CaSR_max CaSR_min];
    
    outputs = [flag_abnormal_rep flag_ead flag_normal_AP run_Vm_min run_delta_Na...
        outputs_beat];
end
