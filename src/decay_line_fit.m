function x = decay_line_fit(auto_cor,x0)
%Fits the elongation correlation term plus an exponential decay term for
%the dynamics correlation to the given auto correlation with a nonlinear
%least squares fitting
%  auto_cor: the auto correlation function that we are trying to fit
%  x0: intial values for [a,b,c,el,al]
%  el: estimated elongation time
%  al: estimated rise time
%  a,b,c: miscellaneous other parameters

    addpath('utilities/');
    f = @(vars) [0]; %vars = [a,b,c,el,al]
    for i = 1:length(auto_cor)
        f = @(vars) [f(vars) vars(1)*exp(-vars(2)*(i-1)) + ...
            vars(3)*elong_term(vars(4),vars(5),i-1) - auto_cor(i)];
    end
    x = lsqnonlin(f,x0,[0,0,0,0,0]);
    %plot result
    approx = zeros(1,length(auto_cor));
    for i = 1:length(approx)
        approx(i) = x(1) * exp(-x(2)*(i-1)) + x(3)*elong_term(x(4),x(5),i-1);
    end
    figure();
    grid on
    plot(0:length(approx)-1,approx);
end

