function [best_elong,best_rise,all] = grid_monte_carlo(time_res, points_per_trace, ...
                            num_traces, num_states, trans_mat, ...
                            rna_per_sec, fluo_per_rna, ...
                            init_dist,auto_cor,elong_min,elong_max, ...
                            elong_step, rise_min, rise_max, rise_step)
%Calculates the most likely elongation time and rise time for the
%autocorrelation with all the other parameters given
%   Detailed explanation goes here
    auto_cor1st = auto_cor(2:end) - auto_cor(1:end-1);
    auto_cor2nd = auto_cor1st(2:end) - auto_cor1st(1:end-1);
    elongs = elong_min:elong_step:elong_max;
    rises = rise_min:rise_step:rise_max;
    all = zeros(length(elongs),length(rises));
    best_err = 2^31;
    best_rise = rise_min;
    best_elong = elong_min;
    i = 1;
    for elong = elongs
        j = 1;
        for rise = rises
            traces = gen_data(elong,time_res,points_per_trace, ...
                            num_traces, num_states, trans_mat, ...
                            rna_per_sec, fluo_per_rna, rise, ...
                            init_dist, 0);
            corr = auto_corr_m_calc_norm(traces, length(auto_cor));
            corr1st = corr(2:end) - corr(1:end-1);
            corr2nd = corr1st(2:end) - corr1st(1:end-1);
            corr2nd = corr2nd * (auto_cor2nd(i) / corr2nd(i));
            err_arr = corr2nd(2:end) - auto_cor2nd(2:end);
            err = err_arr * err_arr';
            all(i,j) = err;
            if err < best_err
                best_err = err;
                best_rise = rise;
                best_elong = elong;
            end
            j = j + 1;
        end
        i = i + 1;
    end
end

