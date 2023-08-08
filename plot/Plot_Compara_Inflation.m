% Plotting the results of solving the modified Lagos-Wright model
% Comparing results for different levels of trend inflation
% Part of "The Long-Run Redistributive Effects of Monetary Policy"
% Christian Bustamante
% Last modified: 07/08/2023 at 22:46
clear; clc
rgbcol  = get(groot,'DefaultAxesColorOrder');
fntsize = 14;
xcolor;


%% Load data
mu_num = [1.01, 1.02, 1.10, 1.20];
mu_str = {'010','020','100','200'};
Nmu    = length(mu_num);
cd ..
cd Output
mkdir Figs
for i = 1:length(mu_num)
    eval(['Out_Data_',mu_str{i},' = Read_HDF5(''OutData_Steady_Chi050_Mu',mu_str{i},'.h5'');';]);   
end
cd ..


%% Some parameters
KnotsM = Out_Data_010.KnotsM;
Grid_M = Out_Data_010.Grid_M;
beta   = Out_Data_010.Par_Beta;
Nm     = length(KnotsM);
Ng     = length(Grid_M);
mlow   = min(KnotsM);
mhigh  = max(KnotsM);
inom   = Out_Data_010.Par_rReal*(1+inf) - 1;

% Allocating
Mu_Phi   = zeros(Nmu,1);
Mu_Velo  = zeros(Nmu,1);
Mu_Inf   = zeros(Nmu,1);
Mu_MBar  = zeros(Nmu,1);
Mu_P     = zeros(Nmu,1);
Mu_P_Std = zeros(Nmu,1);
Mu_MBar_Std = zeros(Nmu,1);
Mu_CAgg_DM  = zeros(Nmu,1);
Mu_CAgg_CM  = zeros(Nmu,1);


%% Arranging the data of differerent stationary equilibria
for m = 1:length(mu_str)
    
    % Temp
    eval(['Out_Data = Out_Data_',mu_str{m},';']);   
    Grid_M = Out_Data.Grid_M;
    Dist_F = Out_Data.Dist_F;
    Dist_P = Out_Data.Dist_P;
    Grid_P = Out_Data.Grid_P;    
    
    % Obtaining mean and variances of distributions of M and P
    Ng = length(Grid_M);
    MBar = 0.0;
    MBar_Var = 0.0;
    for i = 1:Ng
        MBar = MBar + Dist_F(i)*(Grid_M(i));
        MBar_Var  = MBar_Var  + Dist_F(i)*(Grid_M(i)^2);
    end
    MBar_Var = MBar_Var - MBar^2;
    Mean_P = 0.0;
    Var_P  = 0.0;
    for i = 1:Ng
        Mean_P = Mean_P + Dist_P(i)*(Grid_P(i));
        Var_P  = Var_P  + Dist_P(i)*(Grid_P(i)^2);
    end
    Var_P = Var_P - Mean_P^2;
    
    % Saving data in consolidated vector
    Mu_Phi(m)  = Out_Data.phim;
    Mu_Velo(m) = Out_Data.Velo;
    Mu_Inf(m)  = Out_Data.Inf_Rate;
    Mu_MBar(m) = MBar;
    Mu_P(m)    = Mean_P;
    Mu_MBar_Std(m) = sqrt(MBar_Var);
    Mu_P_Std(m)    = sqrt(Var_P);
    Mu_CAgg_DM(m)  = Out_Data.CAgg_DM_Real;
    Mu_CAgg_CM(m)  = Out_Data.CAgg_CM_Real;   
end


%% Table in levels
% Seven columns: 
% 1. Inflation, 2. Output in DM, 3. Price money, 4. Avg price in DM, 
% 5. Std Price in DM 
Tab_Comp = zeros(length(mu_num),5);
Tab_Comp(:,1) = 100*(mu_num-1)';
Tab_Comp(:,2) = Mu_CAgg_DM;
Tab_Comp(:,3) = Mu_Phi;
Tab_Comp(:,4) = Mu_P;
Tab_Comp(:,5) = Mu_P_Std;
Tab_Comp = Tab_Comp';

dfile  = 'Output/Table_D6.txt';
fileID = [1,fopen(dfile,'w')];

rLab = {'Inflation     ',...
        'Output_DM     ',...
        'Price_M       ',...
        'Mean_Price    ',...
        'SD_Price      '};
printf(fileID,'---------------------------------------------------------------------------------------------\n')
printf(fileID,'Table D9\n')
printf(fileID,'Stationary equilibrium for different levels of trend inflation\n')
printf(fileID,'One-asset model\n')
printf(fileID,'---------------------------------------------------------------------------------------------\n')
fmt = [ 'Inflation    ', ' %8.2f \t', ' %8.2f \t', ' %8.2f \t', ' %8.2f \t\n']; printf(fileID, fmt, Tab_Comp(1,:))

for i = 2:size(Tab_Comp,1)
    fmt = [rLab{i}, '%8.2f \t', ' %8.2f \t', ' %8.2f \t', ' %8.2f \t\n']; printf(fileID, fmt, Tab_Comp(i,:))
end
printf(fileID,'---------------------------------------------------------------------------------------------\n\n')
fclose(fileID(2));


%% Normalizing
Mu_Inf = Mu_Inf + 1;
inf_normal = 1.02;
index_normal  = find(Mu_Inf>inf_normal-1e-6 & Mu_Inf<inf_normal+1e-6);
Mu_N_Velo     = Mu_Velo    /Mu_Velo(index_normal);
Mu_N_MBar     = Mu_MBar    /Mu_MBar(index_normal);
Mu_N_MBar_Std = Mu_MBar_Std/Mu_MBar_Std(index_normal);
Mu_N_P        = Mu_P       /Mu_P(index_normal);
Mu_N_P_Std    = Mu_P_Std   /Mu_P_Std(index_normal);
Mu_N_CAgg_CM  = Mu_CAgg_CM /Mu_CAgg_CM(index_normal);
Mu_N_CAgg_DM  = Mu_CAgg_DM /Mu_CAgg_DM(index_normal);
Mu_N_Phi      = Mu_Phi     /Mu_Phi(index_normal);


%% Subplot
fig2 = figure;


% Subplot fot total output
subplot(2,2,1)
plot(100*(Mu_Inf-1),Mu_N_Phi,'Linewidth',2,'Marker','o','Markersize',7); hold on
plot(2,1,'MarkerFaceColor','r','Marker','o','Markersize',7,'MarkerEdgeColor','k')
xlim([min(100*(Mu_Inf-1)),max(100*(Mu_Inf-1))])
title('\textbf{(a)} Price of money','Interpreter','Latex')
xlabel('$\mu$, percent','Interpreter','Latex')
ylabel('$\phi$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')
set(gca,'Fontsize',fntsize)
grid on; grid minor; set(gca,'GridLineStyle',':')

% Subplot for output in DM
subplot(2,2,2)
plot(100*(Mu_Inf-1),Mu_N_CAgg_DM,'Linewidth',2,'Marker','o','Markersize',7); hold on
plot(2,1,'MarkerFaceColor','r','Marker','o','Markersize',7,'MarkerEdgeColor','k')
xlim([min(100*(Mu_Inf-1)),max(100*(Mu_Inf-1))])
title('\textbf{(b)} Output in DM','Interpreter','Latex')
xlabel('$\mu$, percent','Interpreter','Latex')
ylabel('$y_{DM}$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')
set(gca,'Fontsize',fntsize)
grid on; grid minor; set(gca,'GridLineStyle',':')

% Subplot for price in DM
subplot(2,2,3)
plot(100*(Mu_Inf-1),Mu_N_P,'Linewidth',2,'Marker','o','Markersize',7); hold on
plot(2,1,'MarkerFaceColor','r','Marker','o','Markersize',7,'MarkerEdgeColor','k')
xlim([min(100*(Mu_Inf-1)),max(100*(Mu_Inf-1))])
title('\textbf{(c)} Avg. price in DM','Interpreter','Latex')
xlabel('$\mu$, percent','Interpreter','Latex')
ylabel('$\bar{p_{DM}}$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')
set(gca,'Fontsize',fntsize)
grid on; grid minor; set(gca,'GridLineStyle',':')

% Subplot for standard deviation of price in DM
subplot(2,2,4)
plot(100*(Mu_Inf-1),Mu_N_P_Std,'Linewidth',2,'Marker','o','Markersize',7); hold on
plot(2,1,'MarkerFaceColor','r','Marker','o','Markersize',7,'MarkerEdgeColor','k')
xlim([min(100*(Mu_Inf-1)),max(100*(Mu_Inf-1))])
title('\textbf{(d)} Stdev. price in DM','Interpreter','Latex')
xlabel('$\mu$, percent','Interpreter','Latex')
ylabel('$\sigma_{pDM}$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')
set(gca,'Fontsize',fntsize)
grid on; grid minor; set(gca,'GridLineStyle',':')

% Exporting pdf
set(fig2,'PaperOrientation','landscape');
set(fig2,'PaperUnits','normalized');
set(fig2,'PaperPosition', [0 0 1 1]);
print(fig2,'Output/Figs/D8_Comp_Inflation','-dpdf')

