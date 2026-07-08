%% ODE simulations
clc
clear
close all

%% Parameters

% 'oscillation', 'SOC', 'SOB'
behavior = 'oscillation'; 

switch behavior
    case 'oscillation'
        p.kAa = 0.001;
        p.kAb = 0.01;
        p.kAB = 1.5;
    case 'SOC'
        p.kAa = 0.001;
        p.kAb = 0.01;
        p.kAB = 0.8;
    case 'SOB'
        p.kAa = 0.05;
        p.kAb = 0.1;
        p.kAB = 1.2;
    otherwise
        error('Unknown parameter settings for "behavior".')
end

p.kBb = p.kAa;
p.kBa = p.kAb;
p.kBA = p.kAB;

p.kA = 1;
p.kB = 1;
p.n = 1;

p.source = 1e-20;
p.s = p.source;

% transient part of the stochastic and deterministic simulations
transient_det = 0.9;
transient_stoch = 0.1;

% show the full deterministic trajectory
show_trajectory = false;

%% fixed point for initial condition

syms A B a b
% we are looking for the non zero solution
assume([a>0,A>0,B>0,b>0])

% the model
eq = [p.kA*A*((1-A)*b^p.n/(1+b^p.n)-p.kAB*B*a^p.n/(1+a^p.n)),...
    p.kB*B*(a^p.n/(1+a^p.n)-B*a^p.n/(1+a^p.n)-p.kBA*A*b^p.n/(1+b^p.n)),...
    p.kAa*A-p.kBa*a*B,...
    p.kBb*B-p.kAb*b*A];

% find non trivial fixed point
s = solve(eq);

% c: A, B, a, b
c0_eq = [double(s.A),double(s.B),double(s.a),double(s.b)];

%% Simulation

% dA/dt = kA*A*((1-A)*b^n/(1+b^n)-kAB*B*a^n/(1+a^n));
% dB/dt = kB*B*((1-B)*a^n/(1+a^n)-kBA*A*b^n/(1+b^n));
% da/dt = kAa*A - kBa*a*B
% db/dt = kBb*B - kAb*b*A

% c: A, B, a, b
c0 = [.01,.1,.1,.1];

tspan = [0 1e5];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,c] = ode45(@(t,c) odefun_full(t,c,p),tspan,c0,options);

%% plot
% show the full deterministic trajectory

if show_trajectory == true
    figure
    set(gcf,'Position',[524   380   944   357])
    tiledlayout(1,2)
    nexttile
    plot(t,c(:,1:2))
    legend('A','B')
    xlabel('\boldmath $t$','FontSize',12,'Interpreter','latex')
    ylabel('\boldmath $c$','FontSize',12,'Interpreter','latex')
    nexttile
    % plot(t,c(:,3:4))
    hold on
    plot(t,c(:,3)./(1+c(:,3)))
    plot(t,c(:,4)./(1+c(:,4)))
    legend('a','b')
    xlabel('\boldmath $t$','FontSize',12,'Interpreter','latex')
    ylabel('\boldmath $c$','FontSize',12,'Interpreter','latex')
end

%% Explicit solver with Euler

% number of output time points
n_out = 1e6;

% number of time steps to solve the equation
nt = 1e6;

if rem(nt,n_out)~=0
    error('We ue round numbers for substeps.')
else
    % sub steps between output generation
    n_substeps = nt/n_out;
end
t_out = linspace(0,tspan(2),n_out+1);
c_out = zeros(n_out+1,4);

% time step
dt = tspan(2)/nt;

% output position
out_pos = 1;
% initial condition
c_out(out_pos,:) = c0_eq;
% concentration in the given time step
ct = c0_eq.';

% start for a given seed
rng(1)

% volume
omega = 1e4;
% scaling the stochasticity
scale = sqrt(dt/omega);

% iteration
for i = 1:nt

    % dA/dt = kA*A*((1-A)*b^n/(1+b^n)-kAB*B*a^n/(1+a^n));
    % dB/dt = kB*B*((1-B)*a^n/(1+a^n)-kBA*A*b^n/(1+b^n));
    % da/dt = kAa*A - kBa*a*B
    % db/dt = kBb*B - kAb*b*A

    % degradation terms will not contain the "*c" term
    Ha = ct(3)^p.n/(1+ct(3)^p.n);
    Hb = ct(4)^p.n/(1+ct(4)^p.n);
    r1 = p.kA*ct(1).*Hb+p.source;
    r1r1 = p.kA*ct(1).*Hb; % *ct(1)
    r1r2 = p.kA*p.kAB*ct(2)*Ha; % *ct(1)
    r2 = p.kB*ct(2)*Ha+p.source;
    r2r1 = p.kB*ct(2)*Ha; % *ct(2)
    r2r2 = p.kB*p.kBA*ct(1).*Hb; % *ct(2)
    r3 = p.kAa*ct(1);
    r3r = p.kBa*ct(2); % *ct(3)
    r4 = p.kBb*ct(2);
    r4r = p.kAb*ct(1); % *ct(4)

    f1 = r1*dt;
    d1 = (r1r1+r1r2)*dt;
    G1 = scale*(sqrt(r1)*randn - sqrt(r1r1*ct(1))*randn - sqrt(r1r2*ct(1))*randn);
    ct(1) = Patankar_Euler(ct(1), f1, d1, G1);
    f2 = r2*dt;
    d2 = (r2r1+r2r2)*dt;
    G2 = scale*(sqrt(r2)*randn - sqrt(r2r1*ct(2))*randn - sqrt(r2r2*ct(2))*randn);
    ct(2) = Patankar_Euler(ct(2), f2, d2, G2);
    f3 = r3*dt;
    d3 = r3r*dt;
    G3 = scale*(sqrt(r3)*randn - sqrt(r3r*ct(3))*randn);
    ct(3) = Patankar_Euler(ct(3), f3, d3, G3);
    f4 = r4*dt;
    d4 = r4r*dt;
    G4 = scale*(sqrt(r4)*randn - sqrt(r4r*ct(4))*randn);
    ct(4) = Patankar_Euler(ct(4), f4, d4, G4);

    % save the data
    if rem(i,n_substeps) == 0
        out_pos = out_pos+1;
        c_out(out_pos,:) = ct;
    end

end

%% cumulative distribution

if ~exist('output','dir')
    mkdir('output')
end

% get rid of the transient
start = find(t_out>tspan(2)*transient_stoch,1,'first');

% above threshold state
activated_state = c_out(start:end,2) > c0_eq(2);
% calculate the time above this state
d = diff([0; activated_state; 0]);
widths = (find(d==-1) - find(d==1))*dt;
clear d activated_state

% unique widths
[widths, ~, ic] = unique(widths);
% how many times each appears
counts = accumarray(ic, 1);   
% CCDF
ccdf = 1 - cumsum(counts)/sum(counts);


% log-log
figure
% plot(log10(widths),log10(ccdf),'Linewidth',2)
loglog(widths,ccdf,'o-','Linewidth',2)
% set(gca,'XScale','log','YScale','log')

% fit a line
start = 1;%find(widths>10^-0.2,1,'first');
switch behavior
    case 'SOC'
        stop = find(widths<51,1,'last');
    case 'SOB'
        stop = find(widths<5.3,1,'last');
    case 'oscillation'
        % no linear part
    otherwise 
        % determine linear range automatically
        % stop = min([find(widths<10^1.8,1,'last'),length(widths)-1]);
        stop = round(0.1*length(ccdf));

        % initial rough fit
        fp = polyfit(log10(widths(start:stop)),log10(ccdf(start:stop)),1);

        % residuals
        r = log10(ccdf(start:end-1)) - polyval(fp,log10(widths(start:end-1)));

        % detect variance change point
        stop = findchangepts(abs(r), ...
            'Statistic','std', ...
            'MaxNumChanges',1);
        stop = start-1+stop;
end

switch behavior
    case 'oscillation'
        % no linear part
    otherwise
        % fit only the initial linear part
        f = fit(log10(widths(start:stop)),log10(ccdf(start:stop)),'poly1');
        hold on
        loglog(widths(start:stop),10.^f(log10(widths(start:stop))),'Linewidth',2)
%         title(['slope = ' num2str(f.p1)])
end

xlabel('\boldmath \bf $w$','FontSize',16,'Interpreter','latex')
% WR: duration of intervals where R(t)>r
ylabel('\boldmath \bf $P(W_B > w)$','FontSize',16,'Interpreter','latex')
% title(['k_{AB} = ' num2str(p.kAB) ', k_{BA} = ' num2str(p.kBA)])
set(gca,'LineWidth',1.5,'FontSize',16)

% stochastic simulation in time on inset
max_pos = ceil(0.9*length(t_out));
axes('Position',[.28 .3 .32 .32])
plot(t_out(max_pos:end),c_out(max_pos:end,2))
xlabel('\boldmath $t$','FontSize',12,'Interpreter','latex')
ylabel('\boldmath $B$','FontSize',12,'Interpreter','latex')
set(gca,'LineWidth',1.5,'FontSize',12)

exportgraphics(gcf,['output/Crossfeeding_' behavior '_powerlaw.pdf'],'resolution',300)

%% phase space plot in 3D

figure
hold on

% plot the trajectory
% stochastic simulation
start = find(t_out>tspan(2)*transient_stoch,1,'first');
transparency = 0.5;
hp = plot3(c_out(start:end,3),c_out(start:end,4), ...
    c_out(start:end,2)./(c_out(start:end,2)+c_out(start:end,1)), ...
    'Color',[0.49, 0.18, 0.56, transparency]);
% plot the deterministic simulations
start = find(t>tspan(2)*transient_det,1,'first');
if std(c(start:end,3)) < 1e-2
    % fixed point
    hdet = scatter3(c(start:end,3),c(start:end,4),c(start:end,2)./(c(start:end,1)+c(start:end,2)), ...
        100,'filled','MarkerFaceColor',[0.47, 0.67, 0.19]);
else
    % limit cycle
    hdet = plot3(c(start:end,3),c(start:end,4),c(start:end,2)./(c(start:end,1)+c(start:end,2)), ...
        'LineWidth',6,'Color',[0.47, 0.67, 0.19]);
end

switch behavior
    case 'SOC'
        axlim = [0.05,0.15,0.05,0.15];
    case 'SOB'
        axlim = [0,1,0,1];
    case 'oscillation'
        axlim = axis;
end

% syms A B Ha Hb
% assume([A,b,Ha,Hb],'real')
syms A B a b
assume([A,b,a,b],'real')

% the model
eq_sub = [p.kA*A*((1-A)*b.^p.n/(1+b.^p.n)-p.kAB*B*a.^p.n/(1+a.^p.n)),...
    p.kB*B*((1-B)*a.^p.n/(1+a.^p.n)-p.kBA*A*b.^p.n/(1+b.^p.n))];

% find fixed points
sol = solve(eq_sub,[A,B]);
fsolA = matlabFunction(sol.A(4),'Vars',[a,b]);
fsolB = matlabFunction(sol.B(4),'Vars',[a,b]);

% Jacobian
J = jacobian(eq_sub,[A,B]);

% eigenvalues
eigs = eig(J);
feigs = matlabFunction(eigs,'Vars',{A,B,a,b});

% number of grid points in one direction
n_grid = 100;
[a_grid,b_grid] = meshgrid(linspace(axlim(1),axlim(2),n_grid),linspace(axlim(3),axlim(4),n_grid));

% check stability
fstab = @(A,B,a,b) all(real(feigs(A,B,a,b)) < 0);

% real part of leading eigenvalue
fmaxeig = @(A,B,a,b) max(real(feigs(A,B,a,b)));

% leading eigenvalue for three different fixed points of the subsystem
maxeig0 = zeros(size(a_grid)); % 1,0
maxeig1 = zeros(size(a_grid)); % 0,1
maxeig = zeros(size(a_grid)); % mixed

% value of the mixed fixed point
fix_mixed_A = fsolA(a_grid,b_grid);
fix_mixed_B = fsolB(a_grid,b_grid);


% calculate the eigenvals on the grid
for i = 1:n_grid
    for j = 1:n_grid
        maxeig0(i,j) = fmaxeig(fix_mixed_A(i,j),fix_mixed_B(i,j),a_grid(i,j),b_grid(i,j));
        maxeig1(i,j) = fmaxeig(0,1,a_grid(i,j),b_grid(i,j));
        maxeig(i,j) = fmaxeig(1,0,a_grid(i,j),b_grid(i,j));
    end
end

ratio = fix_mixed_B./(fix_mixed_B+fix_mixed_A);
ratio(ratio<-10) = nan;
ratio(ratio>10) = nan;
surf(a_grid,b_grid,ratio,maxeig0,...
    'EdgeColor','none','FaceAlpha',0.9)

surf(a_grid,b_grid,ones(size(a_grid)),maxeig1,...
    'EdgeColor','none','FaceAlpha',0.8)

surf(a_grid,b_grid,zeros(size(a_grid)),maxeig,...
    'EdgeColor','none','FaceAlpha',0.8)

% trajectories to the top
uistack(hp, 'top')
uistack(hdet, 'top')

% nullclines of the fast subsystem
fhill = @(x) x.^p.n./(1+x.^p.n);
fb1 = @(b) (p.kBA*fhill(b)./(1-p.kBA*fhill(b))).^(1/p.n);
fb2 = @(b) (p.kAB^-1*fhill(b)./(1-p.kAB^-1*fhill(b))).^(1/p.n);

fplot3(@(b)b,fb1,@(b)b./b,[axlim(1),axlim(2)],'LineWidth',4,'Color','k','LineStyle',':')
fplot(fb2,[axlim(3),axlim(4)],'LineWidth',4,'Color','k','LineStyle',':')

axis([axlim(1),axlim(2),axlim(3),axlim(4),0,1])

% symmetric color axis around zero
cmax = max(abs(maxeig0(:)));
clim([-cmax cmax])
% blue-white-red colormap
n = 256;
% softer blue tones
neg = [linspace(0.15,0.95,n/2)', ...
       linspace(0.35,0.98,n/2)', ...
       ones(n/2,1)];
% softer red tones
pos = [ones(n/2,1), ...
       linspace(0.98,0.35,n/2)', ...
       linspace(0.95,0.15,n/2)'];
cmap = [neg; pos];

colormap(cmap)

% colorbar
hc = colorbar;
ylabel(hc,'\boldmath $\max(\mathrm{Re}(\lambda))$', ...
    'FontSize',16, ...
    'Interpreter','latex')
xlabel('\boldmath $a$','FontSize',16,'Interpreter','latex')
ylabel('\boldmath $b$','FontSize',16,'Interpreter','latex')
zlabel('\boldmath $x_B$','FontSize',16,'Interpreter','latex')
set(gca,'LineWidth',1.5,'FontSize',16)
if p.kAB < 1
    view(-62,23)
else
    view(3)
end
box on
grid on

% colorbar ylabel is out, scale the axes
axpos = get(gca,'Position');
axpos(3) = axpos(3)*.9;
set(gca,'Position', axpos)

exportgraphics(gcf,['output/Crossfeeding_' behavior '_3D.pdf'],'resolution',300)


%% phase space 2D

figure
hold on

% stochastic simulation
transparency = 0.3;
start = find(t_out>tspan(2)*transient_stoch,1,'first');
plot(c_out(start:end,1),c_out(start:end,2),'Color',[0.49, 0.18, 0.56, transparency])

% nullclines
fB1 = @(A,Ha,Hb) 1-p.kBA*Hb/Ha*A;
fB2 = @(A,Ha,Hb) Hb/(p.kAB*Ha)*(1-A);
fplot(@(A) fB1(A,c0_eq(3)/(1+c0_eq(3)),c0_eq(4)/(1+c0_eq(4))),[0,1],'Linewidth',2,'Color',[0 0.4470 0.7410])
fplot(@(A) fB2(A,c0_eq(3)/(1+c0_eq(3)),c0_eq(4)/(1+c0_eq(4))),[0,1],'Linewidth',2,'Color',[0.8500 0.3250 0.0980])

% plot the deterministic simulations
start = find(t>tspan(2)*transient_det,1,'first');
if std(c(start:end,3)) < 1e-2
    % fixed point
    hdet = scatter(c(start:end,1),c(start:end,2), ...
        100,'filled','MarkerFaceColor',[0.47, 0.67, 0.19]);
else
    % limit cycle
    hdet = plot(c(start:end,1),c(start:end,2), ...
        'LineWidth',3,'Color',[0.47, 0.67, 0.19,1]);
end

axis([0,1,0,1])
xlabel('\boldmath $A$','FontSize',16,'Interpreter','latex')
ylabel('\boldmath $B$','FontSize',16,'Interpreter','latex')
set(gca,'LineWidth',1.5,'FontSize',16)

exportgraphics(gcf,['output/Crossfeeding_' behavior '_2D.pdf'],'resolution',300)
