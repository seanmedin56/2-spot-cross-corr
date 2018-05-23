% analyzes the cross correlation and related statistics of two traces
%saves a plot of R and M for the autocorrelation (with bootstraps)
%also saves a plot for the first and second derivative (without bootstraps,
%but can easily be changed to include he boostraps)

addpath('utilities/');

if exist('project') ~= 1
    warning('Define "project" var (str): project identifier');
    return
end

% create directories for outputing plots
dirs = {['../out/' project '/']};

for i = 1:numel(dirs)
    if (exist(dirs{i}, 'dir') ~= 7)
        mkdir(dirs{i});
    end
end

%checks if traces exist in the workspace
if exist('trace1') ~= 1
    warning('Define trace1 (cell array of arrays)');
    return
end

%checks if traces exist in the workspace
if exist('trace2') ~= 1
    warning('Define trace2 (cell array of arrays)');
    return
end

%checks if max_delay exist in the workspace
if exist('max_delay') ~= 1
    warning('Define how many time delay points to be analyzed (int)');
    return
end

% -------------plots the autocovariance and the raw moment----------------

corr_r = smooth(cross_corr_r_calc(trace1, trace2, max_delay))';
if length(trace1) > 5
    [bootstrap_r, std1, std2] = corr_bootstraps(trace1, trace2, max_delay, 1000, 'r');
end

h = figure;
if length(trace1) > 5
    errorbar(0:max_delay-1, corr_r, bootstrap_r);
else
    plot(0:max_delay-1, corr_r);
end
title(['raw moment: ' project]);
xlabel('time delay');
savefig([dirs{1} 'raw_moment.fig']);
close(h);

corr_m = cross_corr_m_calc(trace1, trace2, max_delay);
if length(trace1) > 5
    [bootstrap_m, std1, std2] = corr_bootstraps(trace1, trace2, max_delay, 1000, 'm');
end

h = figure;
if length(trace1) > 5
    errorbar(0:max_delay-1, corr_m, bootstrap_m);
else
    plot(0:max_delay-1, corr_m);
end
title(['central moment: ' project]);
xlabel('time delay');
savefig([dirs{1} 'central_moment.fig']);
close(h);

% -----------plots the first and second derivatives ----------------------
% -----------of the autocovariance and the raw moment --------------------

r_1st_deriv = corr_r(2:max_delay) - corr_r(1:max_delay - 1);
m_1st_deriv = corr_m(2:max_delay) - corr_m(1:max_delay - 1);

r_2nd_deriv = r_1st_deriv(2:max_delay - 1) - r_1st_deriv(1:max_delay - 2);
m_2nd_deriv = m_1st_deriv(2:max_delay - 1) - m_1st_deriv(1:max_delay - 2);

h = figure;
plot(r_1st_deriv);
title(['first derivative of raw moment: ' project]);
xlabel('time delay');
savefig([dirs{1} 'raw_moment_1st_deriv.fig']);
close(h);

h = figure;
plot(m_1st_deriv);
title(['first derivative of central moment: ' project]);
xlabel('time delay');
savefig([dirs{1} 'central_moment_1st_deriv.fig']);
close(h);

h = figure;
plot(r_2nd_deriv);
title(['second derivative of raw moment: ' project]);
xlabel('time delay');
savefig([dirs{1} 'raw_moment_2nd_deriv.fig']);
close(h);

h = figure;
plot(m_2nd_deriv);
title(['second derivative of central moment: ' project]);
xlabel('time delay');
savefig([dirs{1} 'central_moment_2nd_deriv.fig']);
close(h);

