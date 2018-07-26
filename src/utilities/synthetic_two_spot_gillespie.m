function fluos = synthetic_two_spot_gillespie(num_points, alpha, ...
                                K, elong_times, R, deltaT, r_emission, noise, pi0, ...
                                fluo_per_rna, flor)

    % Generates a fluorescence sequence for the given model parameters
    % using the Gillespie algorithm.
    % 
    % INPUTS
    % num_points: length of the sequence in time steps
    % alpha: length of the MS2 loop in time steps
    % K: number of naive states
    % elong_times: cell array with elongation time for each construct
    % (must be in descending order)
    % R: transition rate matrix [sec^{-1}]
    % deltaT: time resolution [sec]
    % r_emission: emission rates [sec^{-1}]
    % noise: Gaussian noise [a.u.]
    % pi0: initial state probability mass function (pmf)
    % flor: mandatory waiting time between polymerases being loaded
    % 
    % OUTPUTS
    % fluos: the two traces--one for each color
    
    % uniformly distributed time points with resolution deltaT 
    % and length num_points
    times_unif = (1:num_points) * deltaT;
    
    % arrays for storing the naive states, times of transitions and 
    % cumulative fluorescence values during the Gillespie simulation
    naive_states = [];
    times = [];
    fluo_cum = [];
    arrival_times = [];
    
    % duration of the simulated process
    t_max = num_points * deltaT;
    
    % first state obtained by sampling the initial state pmf
    naive_states(1) = randsample(1:K, 1, true, pi0);
    
    % assign 0 values to the time and fluo_cum arrays at t = 0
    times(1) = 0;
    
    % variable to keep track of the current reaction time
    t = 0;
    
    % indexing variable to keep track of the current transition
    ind = 1;
    
    % keeps track of how much more wait time is required before loading
    dead = 0;
    
    % generate a sequence of naive states using the Gillespie algorithm
    while (t < t_max)
        
        lambda = -R(naive_states(ind),naive_states(ind));
        dt = exprnd(1/lambda);
        t = t + dt;
        
        rates = R(:,naive_states(ind));
        rates(naive_states(ind)) = 0;
        
        probs = rates / lambda;
        
        naive_states(ind + 1) = randsample(1:K, 1, true, probs);
        ind = ind + 1;
        times(ind) = t;
        
        avg_time = 1 / r_emission(naive_states(ind-1));
        ddt = exprnd(avg_time) + dead;
        while ddt < dt
            dead = flor;
            arrival_times = [arrival_times t - dt + ddt];
            next = exprnd(avg_time) + dead;
            ddt = ddt + next;
        end
        dead = max(0, flor - (ddt - dt));
    end
    
    % number of transitions in the simulated trace
    transition_count = ind-1;
    
    % loops through the different memories to generate the respective
    % traces
    fluos = cell([1 length(elong_times)]);
    for wi = 1:length(elong_times)

        % find the fluorescence at the uniformly distributed points taking the
        % MS2 loops into account
        arrival_times = sort(arrival_times);
        arrival_times_r = arrival_times + (elong_times{1} - elong_times{wi}) * deltaT;
        fluos{wi} = zeros(1, num_points);
        for k = 1:num_points
            t_end = times_unif(k);
            t_start = max([0, t_end - elong_times{wi}*deltaT]);

            ind_start_before = find(arrival_times_r >= t_start);
            if isempty(ind_start_before)
                continue;
            end
            i_start = ind_start_before(1);

            ind_end_before = find(arrival_times_r <= t_end);
            if isempty(ind_end_before)
                continue;
            end
            i_end = ind_end_before(end);

            times_window = arrival_times_r(i_start:i_end);

            for i = 1:(length(times_window))

                t2 = t_end - times_window(i);
                if t2 > alpha * deltaT
                    fluos{wi}(k) = fluos{wi}(k) + fluo_per_rna;
                else
                    fluos{wi}(k) = fluos{wi}(k) + fluo_per_rna * ...
                              t2 / (alpha * deltaT);
                end
            end
        end

        % Gaussian noise
        gauss_noise = normrnd(0,noise,1,num_points);

        % add a Gaussian noise
        fluos{wi} = fluos{wi} + gauss_noise;
    end
     
