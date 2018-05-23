function stds = corr_2nd_deriv_bootstraps(trace1, trace2, max_delay, num_times, type)

% takes random selections of traces and calculates the standard deviation
% of the autocorrelations
% can be done for cross correlation (trace1 != trace2) or for
% autocorrelation (trace1 == trace2)
%   max_delay: number of time delay points to take in the auto correlation
%   num_times: number of times to sample the traces
%   type: 'r' = raw moment, anything else = central moment

    vals = cell([1 max_delay - 2]);

    for i = 1:num_times
        sample_idx = randi([1 length(trace1)], 1, length(trace1));
        samples1 = cell([1 length(trace1)]);
        samples2 = cell([1 length(trace1)]);
        for j = 1:length(samples1)
            samples1{j} = trace1{sample_idx(j)};
            samples2{j} = trace2{sample_idx(j)};
        end

        if type == 'r'
            corr = auto_corr_r_calc(samples1, samples2, max_delay);
        else
            corr = auto_corr_m_calc(samples1, samples2, max_delay);
        end
        c1st_deriv = corr(2:max_delay) - corr(1:max_delay - 1);
        c2nd_deriv = c1st_deriv(2:max_delay-1) - c1st_deriv(1:max_delay-2);
        for j = 1:length(c2nd_deriv)
            vals{j}(i) = c2nd_deriv(j);
        end
    end
    stds = zeros([1 length(vals)]);
    for i = 1:length(vals)
        stds(i) = std(vals{i});
    end
end