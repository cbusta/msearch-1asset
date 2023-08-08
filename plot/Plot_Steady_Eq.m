% Plotting the results of solving the modified Lagos-Wright model
% Part of "The Long-Run Redistributive Effects of Monetary Policy"
% Christian Bustamante
% Last modified: 07/08/2023 at 20:49
clear; clc
rgbcol  = get(groot,'DefaultAxesColorOrder');
fntsize = 14;


%% Load data
cd ..
cd Output
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


%% Plot terms of trade
fig2 = figure;
mesh(Out.KnotsM,Out.KnotsM,Out.Dsc0)
xlim([0,mhigh]);ylim([0,mhigh])
xlabel('$m_s$','Interpreter','Latex')
ylabel('$m_b$','Interpreter','Latex')
zlabel('$d(m_b,m_s)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',14)
fixpdf(fig2,'Output/Figs/D4a_Terms_Trade_D')

fig3 = figure;
mesh(Out.KnotsM,Out.KnotsM,Out.Qsc0)
xlim([0,mhigh]);ylim([0,mhigh])
xlabel('$m_s$','Interpreter','Latex')
ylabel('$m_b$','Interpreter','Latex')
zlabel('$q(m_b,m_s)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',14)
fixpdf(fig3,'Output/Figs/D4b_Terms_Trade_Q')

fig4 = figure;
mesh(Out.KnotsM,Out.KnotsM,Out.Psc)
xlim([0,mhigh]);ylim([0,mhigh])
xlabel('$m_s$','Interpreter','Latex')
ylabel('$m_b$','Interpreter','Latex')
zlabel('$p_d(m_b,m_s)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',14)
fixpdf(fig4,'Output/Figs/D4c_Terms_Trade_P')


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
% print(fig5,'Output/Figs/Dist_M_Nolab','-dpdf')
fixpdf(fig5,'Output/Figs/D5a_Dist_M_Nolab')


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
fixpdf(fig8,'Output/Figs/D5b_Dist_P_Nolab')





