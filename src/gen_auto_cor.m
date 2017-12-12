function gen_auto_cor(traces, auto, first, second, bootstraps, max_delay)
% Generates an autocorrelation and/or the derivatives of an
% autocorrelation for the traces
% traces: The traces we are taking the autocorrelation of
% auto: Boolean saying whether or not to include the autocorrelation
% first: Boolean determining if to include 1st derivative
% second: Boolean determining if to include 2nd derivative
% bootstraps: Boolean determining if to include bootstraps
% max_delay: How many points to take of the autocorrelation

corr = auto_corr_m_calc_norm(traces, max_delay);

corr_1st = corr(2:max_delay) - corr(1:max_delay-1);

corr_2nd = corr_1st(2:max_delay-1) - corr_1st(1:max_delay-2);

if auto
    h = figure;
    plot(0:max_delay-1, corr);
    title('Central Moment');
    grid on
end

if first
    h = figure;
    plot(corr_1st);
    title('Central Moment 1st Derivative');
    grid on
end

if second
    h = figure;
    plot(corr_2nd);
    title('Central Moment 2nd Derivative');
    grid on

end

