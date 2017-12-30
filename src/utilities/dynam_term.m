function dynam_cor = dynam_term(elong,alph,tau,b)
%Calculates the dynamic correlation term
%   Detailed explanation goes here
    dynam_cor = 0;
    if tau < 0
        display(elong);
    end
    if tau < elong
        for i = 1:tau
            term1_add = 0;
            for j = 1:(elong-tau+i)
                term1_add = term1_add + alpha(elong - j, alph) * ...
                    alpha(elong - j - (tau - i),alph);
            end
            dynam_cor = dynam_cor + term1_add * exp(-b*i);
        end
        for i = 1:(elong-tau)
            term2_add = 0;
            for j = 1:(elong - tau - i)
                term2_add = term2_add + alpha(elong - j,alph) * ...
                    alpha(elong - j - tau - i, alph);
            end
            dynam_cor = dynam_cor + term2_add * exp(-b*i);
        end
    else
        for i =(tau - elong + 1):tau
            term1_add = 0;
            for j = 1:(elong - tau + i)
                term1_add = term1_add + alpha(elong - j, alph) * ...
                    alpha(elong - j - (tau - i),alph);
            end
            dynam_cor = dynam_cor + term1_add * exp(-b*i);
        end  
    end
    for i = (tau + 1):(elong + tau)
        term3_add = 0;
        for j = 1:(elong + tau - i)
            term3_add = term3_add + alpha(elong-j, alph) * ...
                alpha(elong - j - (i - tau),alph);
        end
        dynam_cor = dynam_cor + term3_add * exp(-b*i);
    end
end

