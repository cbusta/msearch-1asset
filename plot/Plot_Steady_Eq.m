% Plotting the results of solving the modified Lagos-Wright model
% Christian Bustamante
% Last modified: 24/03/2019 at 14:14
clear; clc
rgbcol  = get(groot,'DefaultAxesColorOrder');
fntsize = 14;


%% Load data
cd ..
addpath('plot')
cd out
mkdir Figs
Out = Read_HDF5('OutData_Steady_Chi050_Mu020.h5');
cd ..


%% Some parameters
beta   = Out.Par_Beta;
Nm     = length(Out.KnotsM);
Ng     = length(Out.Grid_M);
mlow   = min(Out.KnotsM);
mhigh  = max(Out.KnotsM);
inf    = Out.Inf_Rate;
inom   = Out.Par_rReal*(1+inf) - 1;


%% Mean and variance of money holdings before DM
M_Mean_F = 0.0;
M_Mean_G = 0.0;
M_Var_F  = 0.0;
M_Var_G  = 0.0;
for i = 1:Ng
    M_Mean_F = M_Mean_F + Out.Dist_F(i)*(Out.Grid_M(i));
    M_Mean_G = M_Mean_G + Out.Dist_G(i)*(Out.Grid_M(i));
    M_Var_F  = M_Var_F  + Out.Dist_F(i)*(Out.Grid_M(i)^2);
    M_Var_G  = M_Var_G  + Out.Dist_G(i)*(Out.Grid_M(i)^2);
end
M_Var_F = M_Var_F - M_Mean_F^2;
M_Var_G = M_Var_G - M_Mean_G^2;
MBar    = M_Mean_F;


%% Plot decision rules for M'
fig1 = figure;
plot(Out.KnotsM,Out.Gm,'Linewidth',2); hold on; 
plot([mlow,mhigh],[mlow,mhigh],'--')
xlim([0,mhigh])
leg1 = legend('Decision rule for $m^\prime, \ g(m)$','45 degree line');
set(leg1,'box','off');
set(leg1,'Interpreter','Latex','Fontsize',fntsize);
set(leg1,'Location','Best');
xlabel('$m$','Interpreter','Latex')
ylabel('Optimal choice of money, $m^\prime$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',14)
% print(fig1,'out/Figs/Decision_Rule_M','-dpdf')


%% Plot terms of trade
fig2 = figure;
mesh(Out.KnotsM,Out.KnotsM,Out.Dsc0)
xlim([0,mhigh]);ylim([0,mhigh])
xlabel('$m_s$','Interpreter','Latex')
ylabel('$m_b$','Interpreter','Latex')
zlabel('$d(m_b,m_s)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',14)
% print(fig2,'out/Figs/Terms_Trade_D','-dpdf')

fig3 = figure;
mesh(Out.KnotsM,Out.KnotsM,Out.Qsc0)
xlim([0,mhigh]);ylim([0,mhigh])
xlabel('$m_s$','Interpreter','Latex')
ylabel('$m_b$','Interpreter','Latex')
zlabel('$q(m_b,m_s)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',14)
% print(fig3,'out/Figs/Terms_Trade_Q','-dpdf')

fig4 = figure;
mesh(Out.KnotsM,Out.KnotsM,Out.Psc)
xlim([0,mhigh]);ylim([0,mhigh])
xlabel('$m_s$','Interpreter','Latex')
ylabel('$m_b$','Interpreter','Latex')
zlabel('$p_d(m_b,m_s)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',14)
% print(fig4,'out/Figs/Terms_Trade_P','-dpdf')


%% Plot distributions
fig5 = figure;
plot(Out.Grid_M,Out.Dist_G,'Linewidth',2); hold on
plot(Out.Grid_M,Out.Dist_F,'Linewidth',2); hold on
xlim([0,mhigh])
ybound = get(gca,'ylim');
plot([Out.MBar,Out.MBar],ybound,'--','Color',[.5,.5,.5]); hold on
ylim(ybound)
leg1 = legend('$G(m)$','$F(m)$');
set(leg1,'box','off');
set(leg1,'Interpreter','Latex','Fontsize',14);
set(leg1,'Location','Northeast');
xlabel('$m$','Interpreter','Latex')
ylabel('Density','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',fntsize)
% print(fig5,'out/Figs/Dist_M_Nolab','-dpdf')

% Create textbox
str_mbar  = strcat('$\bar{m} \ \  = ',num2str(Out.MBar,'%.4f'),'$');
str_mstdF = strcat('$\sigma_{F} \ = ',num2str(sqrt(M_Var_F),'%.4f'),'$');
str_mstdG = strcat('$\sigma_{G} \ = ',num2str(sqrt(M_Var_G),'%.4f'),'$');
str_infr  = strcat('$\pi \ \ \ = ',num2str(inf,'%.4f'),'$');
str_inom  = strcat('$i \ \ \ \: = ',num2str(inom,'%.4f'),'$');
annotation(fig5,'Textbox',[0.68 0.5 0.1 0.1],...
    'String',{str_mbar,str_mstdF,str_mstdG,str_infr,str_inom},...
    'Interpreter','latex','Fontsize',fntsize,'EdgeColor','none');
% print(fig5,'out/Figs/Dist_M','-dpdf')



%% Value functions
fig6 = figure;
plot(Out.KnotsM,Out.V,'r','Linewidth',2); hold on
plot(Out.KnotsM,Out.W,'Linewidth',2);
xlim([0,mhigh])
xlabel('$m$','Interpreter','Latex')
ylabel('Value','Interpreter','Latex')
leg1 = legend('Value at DM, $V(m)$','Value at CM, $W(m)$','Location','Northwest');
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',fntsize)
set(leg1,'box','off')
set(leg1,'Interpreter','Latex','Fontsize',fntsize);
set(leg1,'Location','Best');
% print(fig6,'out/Figs/Value_Fns','-dpdf')


%% Recision rules* for c2 and h2
fig7 = figure;
subplot(2,1,1)
plot(Out.KnotsM,Out.Cstar,'r','Linewidth',2);
xlim([0,mhigh])
xlabel('$m$','Interpreter','Latex')
ylabel('Consumption','Interpreter','Latex')
leg1 = legend('Consumption in CM, $c(m)$','Location','Southeast');
set(leg1,'Interpreter','Latex','Fontsize',fntsize);
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',fntsize)
set(leg1,'box','off')

subplot(2,1,2)
plot(Out.KnotsM,Out.Hstar,'Linewidth',2);
xlim([0,mhigh])
xlabel('$m$','Interpreter','Latex')
ylabel('Labor','Interpreter','Latex')
leg1 = legend('Labor in CM, $h(m)$','Location','Northeast');
set(leg1,'Interpreter','Latex','Fontsize',fntsize);
set(leg1,'box','off')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',fntsize)
% print(fig7,'out/Figs/Decision_Rule_C+','-dpdf')


%% Computing price distribution in the DM
Mean_P = 0.0;
Var_P  = 0.0;
for i = 1:Ng
    Mean_P = Mean_P + Out.Dist_P(i)*(Out.Grid_P(i));
    Var_P  = Var_P  + Out.Dist_P(i)*(Out.Grid_P(i)^2);
end
Var_P = Var_P - Mean_P^2;


%% Plotting price distribution in the DM
fig8 = figure;
plot(Out.Grid_P,Out.Dist_P,'Linewidth',2); hold on
xlabel('$p_d$','Interpreter','Latex')
ylabel('Density','Interpreter','Latex')
ybound = get(gca,'ylim');
plot([Mean_P,Mean_P],ybound,'--','Color',[.5,.5,.5]); hold on
ylim(ybound)
xlim([Out.Grid_P(1),Out.Grid_P(Ng)])
plow  = 0.95*min(Out.Grid_P(Out.Dist_P>0));
phigh = 1.05*max(Out.Grid_P(Out.Dist_P>0));
xlim([plow,phigh])
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',fntsize)
% print(fig8,'out/Figs/Dist_P_Nolab','-dpdf')

% Create textbox
str_pbar = strcat('$\bar{p}_d \; \, = ',num2str((Mean_P),'%.4f'),'$');
str_pstd = strcat('$\sigma_{pd} = ',num2str(sqrt(Var_P),'%.4f'),'$');
annotation(fig8,'Textbox',[0.68 0.75 0.1 0.1],...
    'String',{str_pbar,str_pstd},...
    'Interpreter','latex','Fontsize',fntsize,'EdgeColor','none');
% print(fig8,'out/Figs/Dist_P','-dpdf')



%% Markup
fig9 = figure;
mesh(Out.KnotsM,Out.KnotsM,Out.Mkup)
xlim([0,mhigh]);ylim([0,mhigh])
xlabel('$m_s$','Interpreter','Latex')
ylabel('$m_b$','Interpreter','Latex')
zlabel('$\theta_\mu(m_b,m_s;F)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',fntsize)
% print(fig9,'out/Figs/Terms_Trade_Mkup','-dpdf')


%% Computing markup distribution in the DM
Mean_Mkup = 0.0;
Var_Mkup  = 0.0;
for i = 1:Ng
    Mean_Mkup = Mean_Mkup + Out.Dist_Mkup(i)*(Out.Grid_Mkup(i));
    Var_Mkup  = Var_Mkup  + Out.Dist_Mkup(i)*(Out.Grid_Mkup(i)^2);
end
Var_Mkup = Var_Mkup - Mean_Mkup^2;


%% Plotting markup distribution in the DM
fig10 = figure;
plot(Out.Grid_Mkup,Out.Dist_Mkup,'Linewidth',2); hold on
xlabel('$\theta_\mu$','Interpreter','Latex')
ylabel('Density','Interpreter','Latex')
ybound = get(gca,'ylim');
plot([Mean_Mkup,Mean_Mkup],ybound,'--','Color',[.5,.5,.5]); hold on
ylim(ybound)
xlim([Out.Grid_Mkup(1),Out.Grid_Mkup(Ng)])
plow  = 0.95*min(Out.Grid_Mkup(Out.Dist_Mkup>0));
phigh = 1.05*max(Out.Grid_Mkup(Out.Dist_Mkup>0));
xlim([plow,phigh])

% Create textbox
str_pbar = strcat('$\bar{\theta}_\mu \; \, = ',num2str((Mean_Mkup),'%.4f'),'$');
str_pstd = strcat('$\sigma_{\theta\mu} = ',num2str(sqrt(Var_Mkup),'%.4f'),'$');
annotation(fig10,'Textbox',[0.68 0.75 0.1 0.1],...
    'String',{str_pbar,str_pstd},...
    'Interpreter','latex','Fontsize',fntsize,'EdgeColor','none');
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',fntsize)
% print(fig10,'out/Figs/Dist_Mkup','-dpdf')




