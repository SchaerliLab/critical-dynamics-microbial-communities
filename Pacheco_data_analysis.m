%% Pacheco data analysis
% https://www.nature.com/articles/s41467-026-73686-w#Sec20

clc
clear
close all

%% Download data

filename = "41467_2026_73686_MOESM8_ESM.xlsx";

if ~isfile(filename)

    url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-73686-w/MediaObjects/41467_2026_73686_MOESM8_ESM.xlsx";

    try
        websave(filename,url);
        fprintf("Downloaded %s\n",filename);
    catch ME
        error("Could not download supplementary data:\n%s",ME.message);
    end

end

%% read data

% load the data
M = readmatrix("41467_2026_73686_MOESM8_ESM.xlsx", ...
            "Sheet","Figure3a");
% Leaf257-Leaf68, Figure3a (oscillation)
% oscillation data set:
data_osc = M(3:end, 2:19);
% time
t_osc = M(3:end, 1);
% start time
t0_osc = 10; % h
t0_pos_osc = find(t_osc>t0_osc,1,'first');

% Leaf257-Leaf265, Figure3a
data = M(3:end, 33:42);
% time
t = M(3:end, 32);
% start time
t0 = 21; % h
t0_pos = find(t>t0,1,'first');

% temporal derivative
dt = diff(t(t0_pos:end));
a = abs(diff(data(t0_pos:end,:)))./dt;

a = a(:);
a = a(~isnan(a) & a > 0); % Remove NaNs and zeros
a = sort(a);
N = length(a);
P = (1:N)' / N;


%% Visualization

log_a = log10(a);
log_P = log10(1-P);

figure
hold on
scatter(log_a,log_P, 15, 'filled', 'MarkerFaceAlpha', 0.5);

% Simple Linear Regression
lb = -1.5; %prctile(log_a, 10);
ub = -0.9; %prctile(log_a, 90);
fit_idx = (log_a > lb) & (log_a < ub);

p = polyfit(log_a(fit_idx), log_P(fit_idx), 1);

% sift the line parallel for better visualization
shift = .3;
p(2) = p(2);
y_fit = polyval(p, [lb,ub]);

plot([lb,ub] + shift, y_fit, 'r:', 'LineWidth', 2);


xlabel('\boldmath $\log_{10}(w)$', 'Interpreter', 'latex','FontSize',18);
ylabel('\boldmath $\log_{10}(P(W > w))$', 'Interpreter', 'latex','FontSize',18);
set(gca,'LineWidth',1.5,'FontSize',14)

grid on;

axes('Position',[.28 .3 .45 .45])
hold on
plot(t,data)
xline(t0,'k:','LineWidth',2)
xlabel('\boldmath \bf Time (h)','FontSize',18,'Interpreter','latex')
ylabel('\bf Relative area','FontSize',18,'Interpreter','latex')
set(gca,'LineWidth',1.5,'FontSize',14)

exportgraphics(gcf,'SOC_Pacheco.pdf','Resolution',300,'ContentType','vector')

%% visualize oscillation

figure
hold on
plot(t_osc,data_osc)
plot(t,mean(data_osc,2,'omitnan'),'k','LineWidth',3)
xline(t0_osc,'k:','LineWidth',2)
xlabel('\boldmath \bf Time (h)','FontSize',18,'Interpreter','latex')
ylabel('\bf Relative area','FontSize',18,'Interpreter','latex')
set(gca,'LineWidth',1.5,'FontSize',14)

exportgraphics(gcf,'oscillation.pdf','Resolution',300,'ContentType','vector')

%% ACF comparison

nLag = 250;

% Equidistant sampling
dt = t(2)-t(1);
t_lag = dt*(0:nLag);

% Compute ACFs
[acf_mean, acf_sem] = computeMeanACF(data,t0_pos,nLag);
[acf_mean_osc, acf_sem_osc] = computeMeanACF(data_osc,t0_pos_osc,nLag);

% Plot

figure
hold on

% Transition (red)
fill([t_lag fliplr(t_lag)], ...
     [acf_mean+acf_sem; flipud(acf_mean-acf_sem)]', ...
     [0.85 0.2 0.2], ...
     'FaceAlpha',0.25, ...
     'EdgeColor','none');

plot(t_lag,acf_mean,...
    'Color',[0.75 0 0],...
    'LineWidth',2.5)

% Oscillation (blue)
fill([t_lag fliplr(t_lag)], ...
     [acf_mean_osc+acf_sem_osc; flipud(acf_mean_osc-acf_sem_osc)]', ...
     [0.2 0.4 0.9], ...
     'FaceAlpha',0.25, ...
     'EdgeColor','none');

plot(t_lag,acf_mean_osc,...
    'Color',[0 0.2 0.8],...
    'LineWidth',2.5)

xlabel('\bf Lag time (h)','Interpreter','latex','FontSize',18)
ylabel('\bf Autocorrelation','Interpreter','latex','FontSize',18)

legend({'Criticality \pm SEM','Criticality',...
        'Oscillation \pm SEM','Oscillation'},...
        'Location','northeast')

grid on
box on
set(gca,'LineWidth',1.5,'FontSize',14)
xlim([0 max(t_lag)])

exportgraphics(gcf,'ACF.pdf','Resolution',300,'ContentType','vector')


% helper function for ACF
function [acf_mean,acf_sem] = computeMeanACF(data,t0_pos,nLag)

nRep = size(data,2);
acf = nan(nLag+1,nRep);

for i = 1:nRep

    x = data(t0_pos:end,i);
    x = x(~isnan(x));

    acf(:,i) = autocorr(x,'NumLags',nLag);

end

acf_mean = mean(acf,2,'omitnan');
acf_sem  = std(acf,0,2,'omitnan') ./ sqrt(sum(~isnan(acf),2));

end