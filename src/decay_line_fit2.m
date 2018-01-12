function opt = decay_line_fit2(auto_cor,upper_limits,lower_limits)
%Fits the elongation correlation term plus an exponential decay term for
%the dynamics correlation to the given auto correlation with a nonlinear
%least squares fitting
%  auto_cor: the auto correlation function that we are trying to fit
%  x0: intial values for [b,el,al]
%  el: estimated elongation time
%  al: estimated rise time
%  b: parameter representing the dynamics

    addpath('utilities/');
    num_vars = 3;
    %generate function which is a combination of an expoential decay term
    %and an elongation term
    f = @(vars) [0]; %vars = [b,el,al]
    for i = 1:length(auto_cor)
        f = @(vars) [f(vars) full_func_cor(vars(2),vars(3),i-1,vars(1)) / ...
            full_func_cor(vars(2),vars(3),0,vars(1)) - auto_cor(i)];
    end
    
    % run non linear least squares on function with multiple random
    % starting points and choose the one with the lowest error
    low_err = 10000;
    opt = zeros(1,num_vars);
    options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8, 'OptimalityTolerance', 1e-8, 'display', 'off');
    for j = 1:10
        x0 = zeros(1,num_vars);
        for i =1:num_vars
            x0(i) = rand() * (upper_limits(i) - lower_limits(i)) + lower_limits(i);
        end
        [x,err] = lsqnonlin(f,x0,lower_limits,upper_limits,options);
        if err < low_err
            low_err = err;
            opt = x;
            disp(low_err);
            disp(opt);
        end
    end
    display(low_err);
    %plot result
    approx = zeros(1,length(auto_cor));
    for i = 1:length(approx)
        approx(i) = full_func_cor(opt(2),opt(3),i-1,opt(1)) / ...
            full_func_cor(opt(2),opt(3),0,opt(1));
    end
    figure();
    plot(0:length(approx)-1,approx);
    grid on
    deriv1 = approx(2:end) - approx(1:end-1);
    figure();
    plot(deriv1);
    grid on
    deriv2 = deriv1(2:end) - deriv1(1:end-1);
    figure();
    plot(deriv2);
    grid on
end

