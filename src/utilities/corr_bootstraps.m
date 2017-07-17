function stds = corr_bootstraps(trace1, trace2, max_delay, num_times, type)

% takes random selections of traces and calculates the standard deviation
% of the autocorrelations

    vals = cell([1 max_delay]);

    for i = 1:num_times
        sample_idx = randi([1 length(trace1)], 1, length(trace1));
        sample1 = cell([1 length(trace1)]);
        sample2 = cell([1 length(trace2)]);
        for j = 1:length(sample1)
            sample1{j} = trace1{sample_idx(j)};
            sample2{j} = trace2{sample_idx(j)};
        end

        if type == 'r'
            corr = cross_corr_r_calc(sample1, sample2, max_delay);
        else
            corr = cross_corr_m_calc(sample1, sample2, max_delay);
        end
        for j = 1:length(corr)
            vals{j}(i) = corr(j);
        end
    end
    stds = zeros([1 length(vals)]);
    for i = 1:length(vals)
        stds(i) = std(vals{i});
    end
end