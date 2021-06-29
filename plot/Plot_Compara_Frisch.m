% Plotting the results of solving the modified Lagos-Wright model
% Comparing results for different Frisch elasticities
% Christian Bustamante
% Last modified: 24/03/2019 at 14:14
clear; clc
rgbcol  = get(groot,'DefaultAxesColorOrder');
fntsize = 14;
xcolor;


%% Load data
chivec_num = [0.0, 0.25, 0.50, 1.0, 2.0, 3.0];
chivec_str = {'000','025','050','100','200','300'};
cd ..
addpath('plot')
cd out
mkdir Figs
for i = 1:length(chivec_num)
    eval(['Out_Data_',chivec_str{i},' = Read_HDF5(''OutData_Steady_Chi',chivec_str{i},'_Mu020.h5'');';]);   
end
cd ..


%% Some parameters
KnotsM = Out_Data_000.KnotsM;
Grid_M = Out_Data_000.Grid_M;
beta   = Out_Data_000.Par_Beta;
Nm     = length(KnotsM);
Ng     = length(Grid_M);
mlow   = min(KnotsM);
mhigh  = max(KnotsM);
inf    = Out_Data_000.Inf_Rate;
inom   = Out_Data_000.Par_rReal*(1+inf) - 1;


%% Distributions of money F and G for values of chi
mhigh_plot = 3.0;
fig1 = figure;

% Subplot for chi = 0
subplot(2,2,1);
plot(Grid_M,Out_Data_000.Dist_F,'Linewidth',2); hold on
plot(Grid_M,Out_Data_000.Dist_G,'Linewidth',2);
xlim([0,mhigh_plot])
xlabel('$m$','Interpreter','Latex')
ylabel('Density','Interpreter','Latex')
title('\textbf{(a)} $\chi = 0$','Interpreter','Latex')
ybound = get(gca,'ylim');
plot([1,1],ybound,'--','Color',[.5,.5,.5]); hold on
ylim(ybound)
leg1 = legend('$F(m)$','$G(m)$');
set(leg1,'box','off');
set(leg1,'Interpreter','Latex','Fontsize',14);
set(leg1,'Location','Northeast');
set(gca,'TickLabelInterpreter','Latex')
set(gca,'Fontsize',fntsize)

% Subplot for chi = 0.25
subplot(2,2,2);
plot(Grid_M,Out_Data_025.Dist_F,'Linewidth',2); hold on
plot(Grid_M,Out_Data_025.Dist_G,'Linewidth',2);
xlim([0,mhigh_plot])
xlabel('$m$','Interpreter','Latex')
ylabel('Density','Interpreter','Latex')
title('\textbf{(b)} $\chi = 0.25$','Interpreter','Latex')
ybound = get(gca,'ylim');
plot([1,1],ybound,'--','Color',[.5,.5,.5]); hold on
ylim(ybound)
set(gca,'TickLabelInterpreter','Latex')
set(gca,'Fontsize',fntsize)

% Subplot for chi = 0.5
subplot(2,2,3);
plot(Grid_M,Out_Data_050.Dist_F,'Linewidth',2); hold on
plot(Grid_M,Out_Data_050.Dist_G,'Linewidth',2);
xlim([0,mhigh_plot])
xlabel('$m$','Interpreter','Latex')
ylabel('Density','Interpreter','Latex')
title('\textbf{(c)} $\chi = 0.5$','Interpreter','Latex')
ybound = get(gca,'ylim');
plot([1,1],ybound,'--','Color',[.5,.5,.5]); hold on
ylim(ybound)
set(gca,'TickLabelInterpreter','Latex')
set(gca,'Fontsize',fntsize)

% Subplot for chi = 3
subplot(2,2,4);
plot(Grid_M,Out_Data_300.Dist_F,'Linewidth',2); hold on
plot(Grid_M,Out_Data_300.Dist_G,'Linewidth',2);
xlim([0,mhigh_plot])
xlabel('$m$','Interpreter','Latex')
ylabel('Density','Interpreter','Latex')
title('\textbf{(d)} $\chi = 3$','Interpreter','Latex')
ybound = get(gca,'ylim');
plot([1,1],ybound,'--','Color',[.5,.5,.5]); hold on
ylim(ybound)
set(gca,'TickLabelInterpreter','Latex')
set(gca,'Fontsize',fntsize)

% Exporting pdf
set(fig1,'PaperOrientation','landscape');
set(fig1,'PaperUnits','normalized');
set(fig1,'PaperPosition', [0 0 1 1]);
% print(fig1,'out/Figs/Comp_Frisch_Dist','-dpdf')


%% Decision rules for m'
fig2 = figure;
plot(KnotsM,Out_Data_000.Gm,':s', 'Linewidth',2); hold on; 
plot(KnotsM,Out_Data_025.Gm,'-',  'Linewidth',2); hold on; 
plot(KnotsM,Out_Data_050.Gm,':x', 'Linewidth',2); hold on; 
plot(KnotsM,Out_Data_300.Gm,'-.', 'Linewidth',2,'Color',xcolo.green); hold on; 
plot([mlow,mhigh],[mlow,mhigh],'--','Color',xcolo.gray)
xlim([0,mhigh])
leg1 = legend('$\chi=0$','$\chi=0.25$','$\chi=0.5$','$\chi=3$','45 degree line');
set(leg1,'box','off');
set(leg1,'Interpreter','Latex','Fontsize',fntsize);
set(leg1,'Location','Northwest');
xlabel('$m$','Interpreter','Latex')
ylabel('Money holdings, $m^\prime(m)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',fntsize)
% print(fig2,'out/Figs/Comp_Frisch_DR_Money','-dpdf')


%% Decision rules for labor
fig3 = figure;
plot(KnotsM,Out_Data_000.Hstar,':s','Linewidth',2); hold on; 
plot(KnotsM,Out_Data_025.Hstar,'-', 'Linewidth',2); hold on; 
plot(KnotsM,Out_Data_050.Hstar,':x','Linewidth',2); hold on; 
plot(KnotsM,Out_Data_300.Hstar,'-.','Linewidth',2,'Color',xcolo.green); hold on; 
plot([mlow,mhigh],[0,0],'--','Color',xcolo.gray)
xlim([0,mhigh])
leg1 = legend('$\chi=0$','$\chi=0.25$','$\chi=0.5$','$\chi=3$');
set(leg1,'box','off');
set(leg1,'Interpreter','Latex','Fontsize',fntsize);
set(leg1,'Location','Southwest');
xlabel('$m$','Interpreter','Latex')
ylabel('Labor, $h(m)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',fntsize)
% print(fig3,'out/Figs/Comp_Frisch_DR_Labor','-dpdf')


%% Decision rules for consumption
fig4 = figure;
plot(KnotsM,Out_Data_000.Cstar,':s','Linewidth',2); hold on; 
plot(KnotsM,Out_Data_025.Cstar,'-', 'Linewidth',2); hold on; 
plot(KnotsM,Out_Data_050.Cstar,':x','Linewidth',2); hold on; 
plot(KnotsM,Out_Data_300.Cstar,'-.','Linewidth',2,'Color',xcolo.green); hold on; 
xlim([0,mhigh])
leg1 = legend('$\chi=0$','$\chi=0.25$','$\chi=0.5$','$\chi=3$');
set(leg1,'box','off');
set(leg1,'Interpreter','Latex','Fontsize',fntsize);
set(leg1,'Location','Southeast');
xlabel('$m$','Interpreter','Latex')
ylabel('Consumption, $c(m)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',fntsize)
% print(fig4,'out/Figs/Comp_Frisch_DR_Consumption','-dpdf')

