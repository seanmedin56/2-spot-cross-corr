function opt = decay_line_fit(auto_cor,upper_limits,lower_limits)
%Fits the elongation correlation term plus an exponential decay term for
%the dynamics correlation to the given auto correlation with a nonlinear
%least squares fitting
%  auto_cor: the auto correlation function that we are trying to fit
%  x0: intial values for [a,b,c,el,al]
%  el: estimated elongation time
%  al: estimated rise time
%  a,b,c: miscellaneous other parameters

    addpath('utilities/');
    
    %generate function which is a combination of an expoential decay term
    %and an elongation term
    f = @(vars) [0]; %vars = [a,b,c,el,al]
    for i = 1:length(auto_cor)
        f = @(vars) [f(vars) vars(1)*exp(-vars(2)*(i-1)) + ...
            vars(3)*elong_term(vars(4),vars(5),i-1) - auto_cor(i)];
    end
    
    % run non linear least squares on function with multiple random
    % starting points and choose the one with the lowest error
    low_err = 10000;
    opt = [0,0,0,0,0];
    for j = 1:100
        x0 = zeros(1,5);
        for i =1:5
            x0(i) = rand() * (upper_limits(i) - lower_limits(i)) + lower_limits(i);
        end
        [x,err] = lsqnonlin(f,x0,lower_limits,upper_limits);
        if err < low_err
            low_err = err;
            opt = x;
        end
    end
    %plot result
    approx = zeros(1,length(auto_cor));
    for i = 1:length(approx)
        approx(i) = opt(1) * exp(-opt(2)*(i-1)) + opt(3)*elong_term(opt(4),opt(5),i-1);
    end
    figure();
    grid on
    plot(0:length(approx)-1,approx);
end

