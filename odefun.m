function dcdt = odefun(~,c,p)
dcdt = zeros(3,1);

% dA/dt = kA*A*(1-A-kAB*B*a^n/(1+a^n)) + source;
% dB/dt = kB*B*((1-B)*a^n/(1+a^n)-kBA*A) + source;
% da/dt = kAa*A - kBa*a*B

Ha = c(3)^p.n/(1+c(3)^p.n);
dcdt(1) = p.kA*c(1)*(1-c(1)-p.kAB*c(2)*Ha) + p.source;
dcdt(2) = p.kB*c(2)*((1-c(2))*Ha-p.kBA*c(1)) + p.source;
dcdt(3) = p.kAa*c(1) - p.kBa*c(3)*c(2);

end

