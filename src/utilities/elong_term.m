function elong_cor = elong_term(elong,alph,tau)
%Calculates the elongation correlation term
%   Detailed explanation goes here
    elong_cor = 0;
    if tau < elong
        for i=1:(elong-tau)
            elong_cor = elong_cor + alpha(elong-i,alph)*alpha(elong-i-tau,alph);
        end
    end
end

