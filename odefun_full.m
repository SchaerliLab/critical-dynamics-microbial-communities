function dcdt = odefun_full(~,c,p)
dcdt = zeros(4,1);

% dA/dt = kA*A*((1-A)*b^n/(1+b^n)-kAB*B*a^n/(1+a^n)) + s;
% dB/dt = kB*B*((1-B)*a^n/(1+a^n)-kBA*A*b^n/(1+b^n)) + s;
% da/dt = kAa*A - kBa*a*B
% db/dt = kBb*B - kAb*b*A

Ha = c(3)^p.n/(1+c(3)^p.n) + p.s;
Hb = c(4)^p.n/(1+c(4)^p.n) + p.s;
dcdt(1) = p.kA*c(1)*((1-c(1))*Hb-p.kAB*c(2)*Ha);
dcdt(2) = p.kB*c(2)*((1-c(2))*Ha-p.kBA*c(1)*Hb);
dcdt(3) = p.kAa*c(1) - p.kBa*c(3)*c(2);
dcdt(4) = p.kBb*c(2) - p.kAb*c(4)*c(1);

end

