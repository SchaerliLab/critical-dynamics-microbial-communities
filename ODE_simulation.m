%% ODE simulations
clc
clear
close all

%% Parameters

% 'oscillation', 'SOC', 'SOB'
behavior = 'oscillation';

switch behavior
    case 'oscillation'
        p.kAa = 0.001; % 0.1
        p.kBa = 0.01; % 0.01
        p.kAB = 1.5;
    case 'SOB'
        p.kAa = 0.05; % 0.1
        p.kBa = 0.1; % 0.01
        p.kAB = 1.4; %1.5;
    case 'SOC'
        p.kAa = 0.001; % 0.1
        p.kBa = 0.01; % 0.01
        p.kAB = 0.9;
end

p.kBA = 0.8;
p.kA = 1;
p.kB = 1;

p.n = 1;

% in (0,0) the system can die out in stochastic simulations
% noise induced extinction
% lets use a small source
p.source = 1e-20;
% just a different name for ODEs to change separately:
p.s = p.source;

% transient part of the stochastic and deterministic simulations
transient_det = 0.9;
transient_stoch = 0.1;

% show simplified, approximating dynamics
show_simplified = false;

%% fixed point for initial condition

syms A B a
% we are looking for the non zero solution
assume([a>0,A>0,B>0])

% the model
eq = [p.kA*A*(1-A-p.kAB*B*a^p.n/(1+a^p.n)),...
    p.kB*B*(a^p.n/(1+a^p.n)-B*a^p.n/(1+a^p.n)-p.kBA*A),...
    p.kAa*A-p.kBa*a*B];

s = solve(eq);

% c: A, B, a
c0_eq = [double(s.A),double(s.B),double(s.a)];
% c0 = c0_eq.*(1+randn(1,3)/10);
c0 = [0.1,0.1,0.1];

%% Simulation

% dA/dt = kA*A*((1-A)*b^n/(1+b^n)-kAB*B*a^n/(1+a^n)) + source;
% dB/dt = kB*B*((1-B)*a^n/(1+a^n)-kBA*A*b^n/(1+b^n)) + source;
% da/dt = kAa*A - kBa*a*B

% cp: p, a
c0p = [c0(2)/(c0(1)+c0(2)),c0(3)];

% cpN: p, N, a
c0pN = [c0(2)/(c0(1)+c0(2)),c0(1)+c0(2),c0(3)];

tspan = [0 1e5];

% original ode
options = odeset('RelTol',1e-12,'AbsTol',1e-16);
[t,c] = ode45(@(t,c) odefun(t,c,p),tspan,c0,options);

if show_simplified == true
    % transformed ode
    [tpN,cpN] = ode45(@(t,c) odefun_pN(t,c,p),tspan,c0pN,options);

    % simplified just with molar ratio
    [tp,cp] = ode45(@(t,c) odefun_p(t,c,p),tspan,c0p,options);
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
c_out = zeros(n_out+1,3);

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

    % degradation terms will not contain the "*c" term
    Ha = ct(3)^p.n/(1+ct(3)^p.n);
    r1 = p.kA*ct(1)+p.source;
    r1r1 = p.kA*ct(1); % *ct(1)
    r1r2 = p.kA*p.kAB*ct(2)*Ha; % *ct(1)
    r2 = p.kB*ct(2)*Ha+p.source;
    r2r1 = p.kB*ct(2)*Ha; % *ct(2)
    r2r2 = p.kB*p.kBA*ct(1); % *ct(2)
    r3 = p.kAa*ct(1);
    r3r = p.kBa*ct(2); % *ct(3)

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
loglog(widths,ccdf,'o-','Linewidth',2)

% fit a line
start = 1;%find(widths>10^-0.2,1,'first');
switch behavior
    case 'SOC'
        stop = find(widths<135,1,'last');
    case 'SOB'
        stop = find(widths<3.2,1,'last');
    case 'oscillation'
    otherwise 
        % determine linear range automatically
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
        f = fit(log10(widths(start:stop)),log10(ccdf(start:stop)),'poly1');
        hold on
        loglog(widths(start:stop),10.^f(log10(widths(start:stop))),'Linewidth',2)
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

exportgraphics(gcf,['output/unifeeding_' behavior '_powerlaw.pdf'],'resolution',300)

%% show the stochastic simulations on the reduced phase space

syms A B a
% we are looking for the non zero solution
assume(a>0)

% sub system, if we fix "a"
eq_sub = [p.kA*A*(1-A-p.kAB*B*a^p.n/(1+a^p.n)),...
    p.kB*B*(a^p.n/(1+a^p.n)-B*a^p.n/(1+a^p.n)-p.kBA*A)];

% find the fixed points
s = solve(eq_sub,[A,B]);

% calculate the Jacobian
J = jacobian(eq_sub,[A,B]);

% a values
n_avals = 1000;
a_max = 10;
a_vals = linspace(0,a_max,n_avals);

% stability of fixed points
% true: stable, false: instable
instab = false(numel(a_vals),numel(s.A));

% values of the fixed point for the B ratio
x_B = zeros(size(instab));

for n_fix = 1 : numel(s.A)
    % evaluate it at the fixed point
    JJ = subs(J,A,s.A(n_fix));
    JJ = subs(JJ,B,s.B(n_fix));

    % eigen values
    e = eig(JJ);
    % e = simplify(eig(JJ));
    % disp(e)

    % eigenvalues in the function of the wavenumber
    eigenvalues = matlabFunction(e,"Vars","a");

    %  evaluate in the given points
    instab(:,n_fix) = arrayfun(@(a) any(real(eigenvalues(a)) > 0), a_vals);

    % fix point value
    fix_A = matlabFunction(s.A(n_fix),"Vars","a");
    fix_B = matlabFunction(s.B(n_fix),"Vars","a");
    Aval = fix_A(a_vals);
    Bval = fix_B(a_vals);
    x_B(:,n_fix) = Bval./(Bval + Aval);

    % get rid of the negative values
    x_B(Aval<0 | Bval<0,n_fix) = nan;

end 

figure;
hold on

for i = 1:numel(s.A)
    
    % stable points
    idx_stable = ~instab(:,i);
    plot(a_vals(idx_stable), x_B(idx_stable,i), 'b', 'LineWidth', 2)
    
    % unstable points
    idx_unstable = instab(:,i);
    plot(a_vals(idx_unstable), x_B(idx_unstable,i), '--r', 'LineWidth', 2)
    
end
xlabel('\boldmath $a$','FontSize',16,'Interpreter','latex')
ylabel('\boldmath $x_B$','FontSize',16,'Interpreter','latex')
% title(['k_{AB} = ' num2str(p.kAB) ', k_{BA} = ' num2str(p.kBA)])
set(gca,'LineWidth',1.5,'FontSize',16)

% plot the stochastic simulations
start = find(t_out>tspan(2)*transient_stoch,1,'first');
hp = plot(c_out(start:end,3),c_out(start:end,2)./(c_out(start:end,1)+c_out(start:end,2)),'Color',[0.49, 0.18, 0.56]);
uistack(hp, 'bottom')
% plot the deterministic simulations
start = find(t>tspan(2)*transient_det,1,'first');
if std(c(start:end,3)) < 1e-2
    % fixed point
    scatter(c(start:end,3),c(start:end,2)./(c(start:end,1)+c(start:end,2)),100,'filled','MarkerFaceColor',[0.47, 0.67, 0.19])
else
    % limit cycle
    plot(c(start:end,3),c(start:end,2)./(c(start:end,1)+c(start:end,2)),'LineWidth',3,'Color',[0.47, 0.67, 0.19])
end

exportgraphics(gcf,['output/unifeeding_' behavior '_2D.pdf'],'resolution',300)

%% plot
% comparing the approximating descriptions

if show_simplified == true
    figure
    set(gcf,'Position',[524   380   944   357])
    tiledlayout(1,2)
    nexttile
    plot(t,c(:,1:2))
    xlabel('\boldmath $t$','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath $c$','FontSize',16,'Interpreter','latex')
    legend('A','B')
    nexttile
    % plot(t,c(:,3:4))
    hold on
    plot(t,c(:,3)./(1+c(:,3)))
    xlabel('\boldmath $t$','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath $H_a$','FontSize',16,'Interpreter','latex')

    figure
    set(gcf,'Position',[524   380   944   357])
    tiledlayout(1,2)
    nexttile
    plot(tpN,cpN(:,1:2))
    legend('p','N')
    xlabel('\boldmath $t$','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath $c$','FontSize',16,'Interpreter','latex')
    nexttile
    % plot(tpN,cpN(:,3:4))
    hold on
    plot(tpN,cpN(:,3)./(1+cpN(:,3)))
    xlabel('\boldmath $t$','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath $H_a$','FontSize',16,'Interpreter','latex')

    figure
    set(gcf,'Position',[524   380   944   357])
    tiledlayout(1,2)
    nexttile
    plot(tp,cp(:,1))
    xlabel('\boldmath $t$','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath $H_a$','FontSize',16,'Interpreter','latex')
    legend('p')
    nexttile
    % plot(tp,cp(:,2))
    hold on
    plot(tp,cp(:,2)./(1+cp(:,2)))
    xlabel('\boldmath $t$','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath $H_a$','FontSize',16,'Interpreter','latex')
end

%% show the corresponding potential
% assume N = 1
% n = 1
if p.n ~= 1
    return;
end

% the potential
V = @(Ha,x) ((p.kB*p.kBA)/4 - (Ha*p.kB)/4 - p.kA/4 + (Ha*p.kA*p.kAB)/4).*x.^4 + ...
(p.kA/3 + (2*Ha*p.kB)/3 - (2*p.kB*p.kBA)/3 - (Ha*p.kA*p.kAB)/3).*x.^3 + ...
(- (Ha*p.kB)/2 + (p.kB*p.kBA)/2 + p.s).*x.^2 - p.s*x;

% figure
% fsurf(@(Ha,p)V(Ha,p),[0,1,-0.2,1.2])
% view(2)

switch behavior
    case 'SOC'
        a = 10;
    case 'SOB'
        a = 3;
    case 'oscillation'
        a = 2.5;
end

figure
fplot(@(p)V(a/(1+a),p),[-0.2,1.2],'LineWidth',2)
xlabel('\boldmath $p$','FontSize',16,'Interpreter','latex')
ylabel('\boldmath $V$','FontSize',16,'Interpreter','latex')
% title(['\boldmath \bf $a$ = ' num2str(a)], 'Interpreter','latex','FontSize',18)

set(gca,'LineWidth',1.5,'FontSize',12)

exportgraphics(gcf,['output/unifeeding_' behavior '_potential.pdf'],'resolution',300)