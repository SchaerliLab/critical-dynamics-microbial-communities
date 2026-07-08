function dcdt = odefun_pN(~,c,p)

dcdt = zeros(3,1);

% dpdt = p*(p - 1)*(kA - Ha*kB - N*kA + N*kB*kBA + N*kA*p + Ha*N*kB*p -
% N*kB*kBA*p - Ha*N*kA*kAB*p) + (1 - 2*p)*s/N
% dNdt = -N*(N*kA - kA + kA*p - Ha*kB*p - 2*N*kA*p + N*kA*p^2 + Ha*N*kB*p^2
% - N*kB*kBA*p^2 + N*kB*kBA*p + Ha*N*kA*kAB*p - Ha*N*kA*kAB*p^2) + 2*s
% dadt = -N*(kAa*p - kAa + a*kBa*p)

Ha = c(3)^p.n/(1+c(3)^p.n);
dcdt(1) = c(1)*(c(1) - 1)*(p.kA - Ha*p.kB - c(2)*p.kA + c(2)*p.kB*p.kBA + ...
    c(2)*p.kA*c(1) + Ha*c(2)*p.kB*c(1) - c(2)*p.kB*p.kBA*c(1) - Ha*c(2)*p.kA*p.kAB*c(1)) + ...
    + (1 - 2*c(1))*p.s/c(2);
dcdt(2) = -c(2)*(c(2)*p.kA - p.kA + p.kA*c(1) - Ha*p.kB*c(1) - 2*c(2)*p.kA*c(1) + ...
    c(2)*p.kA*c(1)^2 + Ha*c(2)*p.kB*c(1)^2 - c(2)*p.kB*p.kBA*c(1)^2 + c(2)*p.kB*p.kBA*c(1) + ...
    Ha*c(2)*p.kA*p.kAB*c(1) - Ha*c(2)*p.kA*p.kAB*c(1)^2) + ...
    2*p.s;
dcdt(3) = p.kAa - p.kAa*c(1) - c(3)*p.kBa*c(1);

end

