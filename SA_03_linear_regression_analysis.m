% This file performs the linear regression analysis.

clear
close all
clc

color = [0 0 0]; % BLACK

%% Load parameters
% load matrix all_parameters (columns: N parameters, rows: N trials)
load SA_par_matrix_1000_s0p1 % sigma 0.1
% all_parameters
% parameter_names

% SA_par: 1) GNa 2) GNaB 3) IbarNaK 4) Gtof 5) GKr
% 6) GKs 7) GKur 8) GKp 9) GK1 10) GK,ACh
% 11) GClCa 12) GClB 13) GCa 14) GCaB 15) IbarPMCA
% 16) IbarNCX 17) VmaxSERCA 18) RyR 19) SR_leak

[N_trials N_pars] = size(all_parameters);

%% Load outputs
% load all_outputs output_names output_units

load SA_output_matrix_control_1Hz_300s

all_outputs_imported = all_outputs;
output_names_imported = output_names;
output_units_imported = output_units;

%% Selection of output of interest
% Outputs for logistic regression (1-5):
% 1) flag_abnormal_rep 2) flag_ead 3) flag_normal_AP 4) run_Vm_min 5) run_delta_Na
% Outputs for linear regression (6-22) - only for normal beat:
% 1) dVm_max 2) Vm_max 3) -Vm_min 4) AP_amp 5) APD90 6) APD70 7) APD50 8) APD30
% 9) Ca_max 10) Ca_min 11) CaT_amp 12) CaT_rise 13) CaT_decay_50 14) CaT_decay_63
% 15) Na_min 16) CaSR_max 17) CaSR_min

logistic_outputs = all_outputs_imported(:,1:3);

all_outputs = all_outputs_imported(:,6:end);
output_names = output_names_imported(:,6:end);
output_units = output_units_imported(:,6:end);

N_outputs = length(output_names); 

%% Linear Regression (normal beat only)
% Check basic properties - AP amplitude > 10 mV
crit_1 = find(all_outputs(:,4)>10);

good_count = length(crit_1);
all_good_parameters = all_parameters(crit_1,:);
all_good_outputs = all_outputs(crit_1,:);

%% Plot histogram APD90 & CaTamp
ind_hist = 5; % APD90
out_hist = all_good_outputs(:,ind_hist);

name_out_hist = output_names(ind_hist)
mean_out_hist = mean(out_hist)
std_out_hist = std(out_hist)

figure,set(gcf,'color','w')%,'Position',[50,100,1500,750])
subplot(2,1,1),histogram(out_hist,'BinWidth',5)
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel(name_out_hist)

ind_hist = 11; % CaTamp
out_hist = all_good_outputs(:,ind_hist);

name_out_hist = output_names(ind_hist)
mean_out_hist = mean(out_hist)
std_out_hist = std(out_hist)

subplot(2,1,2),histogram(out_hist)%histogram(out_hist,'BinWidth',5)
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel(name_out_hist)

%% Define X and Y
X = log(all_good_parameters);
Y = log(all_good_outputs);

%% Call the PLS routine
% PLS - nipals algorithm (2003)
[T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y]=...
    PLS_nipals(X,Y,rank(X));

% % PLS - svds algorithm (2010)
% [T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,...
%           Yhat,Yjack,R2x,R2y,RESSy,PRESSy,Q2,r2y_random,rv_random,...
%           Yhat4Press,Yhat4Ress]=PLS_jack_svds(X,Y,rank(X));

N_pars, N_outputs, N_trials, good_count

%% Goodness of fit
% R^2
% Fraction of the variance in the dependent variable which is explained by the model
% Calculate agreement of values predicted by regression (Yhat = Bpls*X) with original outputs (Y)
SSYT = sum((Y-ones(good_count,1)*mean(Y)).^2);
SSYR = sum((Yhat-ones(good_count,1)*mean(Y)).^2);
R2each = SSYR./SSYT;

% Residual Standard Deviation
% Standard deviation for a normal distribution, centered on the predicted regression line,
% representing the distribution of actually observed values
oSD = zeros(1,N_outputs);
rSD = zeros(1,N_outputs);
for dex = 1:N_outputs
    oSD(dex) = std(exp(Y(:,dex)));
    rSD(dex) = sqrt(sum((exp(Yhat(:,dex)) - exp(Y(:,dex)) ).^2) / (good_count-2));
end

%% Plot
N_figures = ceil(N_outputs/6);

dex1 = 1;
for figdex = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex = 1:6
        if dex1 <= N_outputs
            subplot(2,3,subdex)
            % Plot data points
            plot(exp(Y(:,dex1)),exp(Yhat(:,dex1)),'Marker','o','LineStyle','none','Color',color);
            xlabel(['Actual ', output_names{dex1}])
            ylabel(['Predicted ', output_names{dex1}])
            title(['R^2 = ',num2str(R2each(dex1),4)])
            set(gca,'box','off','tickdir','out','fontsize',10)
            % Plot identity line
            ylim_ind = get(gca,'ylim') ;
            xlim_ind = get(gca,'xlim') ;
            minpoint = min([ylim_ind(1),xlim_ind(1)]);
            maxpoint = max([ylim_ind(2),xlim_ind(2)]);
            hold on
            plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')      
            xlim([minpoint, maxpoint])
            ylim([minpoint, maxpoint])
            dex1 = dex1+1;
        end
    end
end

figure,set(gcf,'color','w')%,'Position',[50,100,1500,750])
bar(R2each,'FaceColor',color)
set(gca,'box','off','tickdir','out','fontsize',12)
set(gca,'XTick',1:N_outputs)
set(gca,'XTickLabel',output_names)
set(gca,'XLim',[0 N_outputs+1])
ylim([0 1])
rotateXLabels( gca(), 90)
title('R^2 values')

% dex1 = 1;
% for figdex = 1:N_figures
%     figure
%     set(gcf,'color','w','Position',[50,100,1500,750])
%     for subdex = 1:6
%         if dex1 <= N_outputs
%             subplot(2,3,subdex)
%             % Plot data points
%             plot(exp(Yhat(:,dex1)),(exp(Y(:,dex1))-exp(Yhat(:,dex1)))/rSD(dex1),'Marker','o','LineStyle','none','Color',color);
%             xlabel(['Predicted ', output_names{dex1}])
%             ylabel('Standardized residuals')
%             title(['rSD/oSD = ',num2str(rSD(dex1)/oSD(dex1),4)])
%             set(gca,'box','off','tickdir','out','fontsize',10)
%             % Plot identity line
%             xlim_ind = get(gca,'xlim') ;
%             hold on
%             plot([xlim_ind(1), xlim_ind(2)],[0,0],'--k') 
%             xlim([xlim_ind(1), xlim_ind(2)])
%             dex1 = dex1+1;
%         end
%     end
% end
% 
% figure,set(gcf,'color','w')%,'Position',[50,100,1500,750])
% subplot(2,2,1),bar(R2each,'FaceColor',color)
% set(gca,'box','off','tickdir','out','fontsize',12)
% set(gca,'XTick',1:N_outputs)
% set(gca,'XTickLabel',output_names)
% set(gca,'XLim',[0 N_outputs+1])
% ylim([0 1])
% rotateXLabels( gca(), 90)
% title('R^2 values')
% 
% subplot(2,2,2),bar(rSD,'FaceColor',color)
% set(gca,'box','off','tickdir','out','fontsize',12)
% set(gca,'XTick',1:N_outputs)
% set(gca,'XTickLabel',output_names)
% set(gca,'XLim',[0 N_outputs+1])
% %ylim([0 1])
% rotateXLabels( gca(), 90)
% title('rSD values')
% 
% subplot(2,2,4),bar(oSD,'FaceColor',color)
% set(gca,'box','off','tickdir','out','fontsize',12)
% set(gca,'XTick',1:N_outputs)
% set(gca,'XTickLabel',output_names)
% set(gca,'XLim',[0 N_outputs+1])
% %ylim([0 1])
% rotateXLabels( gca(), 90)
% title('oSD values')
% 
% subplot(2,2,3),set(gcf,'color','w')%,'Position',[50,100,1500,750])
% bar(rSD./oSD,'FaceColor',color)
% set(gca,'box','off','tickdir','out','fontsize',12)
% set(gca,'XTick',1:N_outputs)
% set(gca,'XTickLabel',output_names)
% set(gca,'XLim',[0 N_outputs+1])
% %ylim([0 1])
% rotateXLabels( gca(), 90)
% title('rSD/oSD values')

%brexit
%% Regression coeffcients
dex2 = 1;
for figdex2 = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex2 = 1:6
        if dex2 <= N_outputs
            subplot(2,3,subdex2)
            bar(Bpls(:,dex2),'FaceColor',color)
            title(output_names(dex2))
            set(gca,'box','off','tickdir','out','fontsize',10)
            set(gca,'XTick',1:N_pars)
            set(gca,'XTickLabel',parameter_names)
            set(gca,'XLim',[0 N_pars+1])
            rotateXLabels( gca(), 90)
            dex2 = dex2 + 1;
        end
    end
end

% N_figures_p = ceil(N_pars/6);
% dex3 = 1;
% for figdex3 = 1:N_figures_p
%     figure
%     set(gcf,'color','w','Position',[50,100,1500,750])
%     for subdex3 = 1:6
%         if dex3 <= N_pars
%             subplot(2,3,subdex3)
%             bar(Bpls(dex3,:),'FaceColor',color)
%             title(parameter_names(dex3))
%             set(gca,'box','off','tickdir','out','fontsize',10)
%             set(gca,'XTick',1:N_outputs)
%             set(gca,'XTickLabel',output_names)
%             set(gca,'XLim',[0 N_outputs+1])
%             rotateXLabels( gca(), 90)
%             dex3 = dex3 + 1;
%             %ylim([-1 1])
%         end
%     end
% end

%% Plot X, Y, B matrices
% % all_parameters
% figure; set(gcf,'color','w')
% imagesc(all_parameters);% colormap jet;
% set(gca,'box','off','tickdir','out','fontsize',10)
% set(gca,'YDir','normal')
% title('Parameters - all parameters');
% xlabel('Parameters');
% ylabel('Trials');
% %set(gca,'YTick',(1:N_trials))
% set(gca,'XTick',(1:N_pars))
% set(gca,'XTickLabel',parameter_names)
% rotateXLabels( gca(), 90)
% colorbar
    
% X = log(all_parameteres)
figure; set(gcf,'color','w')
imagesc(log(all_parameters));% colormap jet;
set(gca,'box','off','tickdir','out','fontsize',10)
set(gca,'YDir','normal')
title('Parameters - log(all parameters)');
xlabel('Parameters');
ylabel('Trials');
%set(gca,'YTick',(1:N_trials))
set(gca,'XTick',(1:N_pars))
set(gca,'XTickLabel',parameter_names)
rotateXLabels( gca(), 90)
colorbar

% % all_outputs
% figure; set(gcf,'color','w')
% imagesc(all_outputs);% colormap jet;
% %imagesc(norm_outputs);% colormap jet;
% set(gca,'box','off','tickdir','out','fontsize',10)
% set(gca,'YDir','normal')
% title('Outputs - all outputs');
% xlabel('Outputs');
% ylabel('Trials');
% %set(gca,'YTick',(1:N_trials))
% set(gca,'XTick',(1:N_outputs))
% set(gca,'XTickLabel',output_names)
% rotateXLabels( gca(), 90)
% colorbar

% Y = log(all_outputs)
figure; set(gcf,'color','w')
imagesc(log(all_outputs));% colormap jet;
%imagesc(norm_outputs);% colormap jet;
set(gca,'box','off','tickdir','out','fontsize',10)
set(gca,'YDir','normal')
title('Outputs - log(all outputs)');
xlabel('Outputs');
ylabel('Trials');
%set(gca,'YTick',(1:N_trials))
set(gca,'XTick',(1:N_outputs))
set(gca,'XTickLabel',output_names)
rotateXLabels( gca(), 90)
colorbar

% norm_outputs = all_outputs;
% mean_outputs = mean(all_outputs);
% for ii = 1:N_trials
%     norm_outputs(ii,:) = all_outputs(ii,:)./mean_outputs;
% end
% 
% figure; set(gcf,'color','w')
% %imagesc(all_outputs);% colormap jet;
% imagesc(norm_outputs);% colormap jet;
% set(gca,'box','off','tickdir','out','fontsize',10)
% set(gca,'YDir','normal')
% title('Outputs - norm(all outputs)');
% xlabel('Outputs');
% ylabel('Trials');
% %set(gca,'YTick',(1:N_trials))
% set(gca,'XTick',(1:N_outputs))
% set(gca,'XTickLabel',output_names)
% rotateXLabels( gca(), 90)
% colorbar

% figure; set(gcf,'color','w')
% %imagesc(all_outputs);% colormap jet;
% imagesc(log(norm_outputs));% colormap jet;
% set(gca,'box','off','tickdir','out','fontsize',10)
% set(gca,'YDir','normal')
% title('Outputs - log(norm(all outputs))');
% xlabel('Outputs');
% ylabel('Trials');
% %set(gca,'YTick',(1:N_trials))
% set(gca,'XTick',(1:N_outputs))
% set(gca,'XTickLabel',output_names)
% rotateXLabels( gca(), 90)
% colorbar

figure; set(gcf,'color','w')
imagesc(Bpls');% colormap jet;
set(gca,'box','off','tickdir','out','fontsize',10)
set(gca,'YDir','normal')
title('Regression coefficients (B)');
xlabel('Parameters');
ylabel('Outputs');
set(gca,'YTick',(1:N_outputs))
set(gca,'YTickLabel',output_names)
set(gca,'XTickLabel',parameter_names)
set(gca,'XTick',(1:N_pars))
rotateXLabels( gca(), 90)
colorbar


% %% Figure with 6 outputs
% output_6 = [5 7 3 11 14 15];
% 
% figure
% set(gcf,'color','w','Position',[50,100,1500,750])
% for dex4 = 1:length(output_6)
%     output_dex4 = output_6(dex4);
%     subplot(2,3,dex4)
%     bar(Bpls(:,output_dex4),'FaceColor',color)
%     title(output_names(output_dex4))
%     set(gca,'box','off','tickdir','out','fontsize',10)
%     set(gca,'XTick',1:N_pars)
%     set(gca,'XTickLabel',parameter_names)
%     set(gca,'XLim',[0 N_pars+1])
%     rotateXLabels( gca(), 90)
% end

%% APD90
output_index = 5;
B_bars = Bpls(:,output_index);
[descBars descI] = sort(abs(B_bars),'descend');
descNames = parameter_names(descI);

figure; set(gcf,'color','w')
bar(B_bars)
set(gca,'box','off','tickdir','out','fontsize',12)
title(output_names(output_index))
ylabel('Rregression coefficients')
set(gca,'XTick',1:N_pars)
set(gca,'XTickLabel',parameter_names)
set(gca,'XLim',[0 N_pars+1])
rotateXLabels( gca(), 90)

% % Figure with descending bars
% figure; set(gcf,'color','w')
% bar(descBars,'FaceColor',color)
% set(gca,'box','off','tickdir','out','fontsize',10)
% title(output_names(output_index))
% ylabel('Abs(regression coefficients)')
% set(gca,'XTick',1:N_pars)
% set(gca,'XTickLabel',descNames)
% set(gca,'XLim',[0 N_pars+1])
% rotateXLabels( gca(), 90)
% %xtickangle(90)
% 
% figure; set(gcf,'color','w')
% bar(B_bars(descI),'FaceColor',color)
% set(gca,'box','off','tickdir','out','fontsize',10)
% title(output_names(output_index))
% ylabel('Regression coefficients')
% set(gca,'XTick',1:N_pars)
% set(gca,'XTickLabel',descNames)
% set(gca,'XLim',[0 N_pars+1])
% rotateXLabels( gca(), 90)
% %xtickangle(45)

% Plot Perturbations
scale = (1-0.2:0.01:1+0.30);

% APD90
var_index = 5;
var_value = mean(all_good_outputs(:,var_index)); % ms

B1 = Bpls(13,var_index); % ICa
mod1 = var_value.*scale.^B1;
B2 = Bpls(9,var_index); % IK1
mod2 = var_value.*scale.^B2;
B3 = Bpls(5,var_index); % IKr
mod3 = var_value.*scale.^B3;

figure, set(gcf,'color','w'), hold on, grid on
hold on, plot(scale,mod1,scale,mod2,scale,mod3)
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('Scale factor (-)');
ylabel('APD90 (ms)');
legend('GCa','GK1','GKr')
xlim([0.8 1.4])

% %% Test with one model from the population
% 
% % Model index in the population
% model_index = 1;
% % Scaling factors
% example_parameters = all_good_parameters(model_index,:);
% % Output
% var_index = 5; % APD90
% % Output (mean in the population)
% mean_output_population = mean(all_good_outputs(:,var_index))
% % Output (simulated value)
% example_simulated_output = all_good_outputs(model_index,var_index)%;
% % Output (regression coefficients)
% output_bars = Bpls(:,var_index);
% 
% example_x = log(example_parameters);
% example_xz = (example_x-mean(X))./std(X);
% 
% example_yz = example_xz*output_bars;
% example_y = example_yz.*std(Y(:,var_index))+mean(Y(:,var_index));
% example_predicted_output = exp(example_y)%;
% 
% % Test with specified scaling factors
% 
% % Scaling factors
% example_parameters = ones(1,19);
%     example_parameters(9) = 0.12; % IK1
%     example_parameters(13) = 0.8; % ICaL
% % Output
% var_index = 5; % APD90
% % Output (mean in the population)
% mean_output_population = mean(all_good_outputs(:,var_index))
% % Output (regression coefficients)
% output_bars = Bpls(:,var_index);
% 
% example_x = log(example_parameters);
% example_xz = (example_x-mean(X))./std(X);
% 
% example_yz = example_xz*output_bars;
% example_y = example_yz.*std(Y(:,var_index))+mean(Y(:,var_index));
% example_predicted_output = exp(example_y)%;
