function c_new = Patankar_Euler(c, f, d, G)
% Patankar-Euler method: https:%doi.org/10.1063/5.0146502
% f: sum of the deterministic production terms x dt
% d: sum of the deterministic degradation terms x dt
% the corresponding reaction rate is in the form: r_i =  d_i x c_i
% so d_i = r_i/c_i and we do not have to divide with small c_i
% G: sum of the diffusion (stochastic) terms x dW

% Euler–Maruyama
c_new = c + f - d*c + G;

if c_new < 0.0

    % Patankar-Euler
    c_new = (c + f + G)/(1. + d);

    % stochastic Patankar-Euler
    if c_new < 0.0
        c_new = (c + f)/(1. + d - G/c + (G/c).^2);
    end

end


end
