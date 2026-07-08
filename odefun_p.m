function dcdt = odefun_p(~,c,p)

dcdt = zeros(2,1);

% dp/dt = -p*(p - 1)*(Ha*kB - kB*kBA - kA*p - Ha*kB*p + kB*kBA*p +
% Ha*kA*kAB*p) + s*(1-2*p)
% da/dt = kAa - kAa*p - a*kBa*p

Ha = c(2)^p.n/(1+c(2)^p.n);
dcdt(1) = -c(1)*(c(1) - 1)*(Ha*p.kB - p.kB*p.kBA - p.kA*c(1) - Ha*p.kB*c(1) + p.kB*p.kBA*c(1) + Ha*p.kA*p.kAB*c(1)) + ...
    + p.s*(1-2*c(1));
dcdt(2) = p.kAa - p.kAa*c(1) - c(2)*p.kBa*c(1);

end

