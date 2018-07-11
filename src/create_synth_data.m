% creates synthetic two spot data (should be run from src folder)

addpath('utilities/');

% create directories for outputing data and the results
dirs = {'../out/synth_dat/', '../dat/synth_dat/'};
     
for i = 1:numel(dirs)
    if (exist(dirs{i}, 'dir') ~= 7)
        mkdir(dirs{i});
    end
end

% number of promoter states
K = 3;

% number of groups
nGR = 1;

% time resolution
deltaT = {10, 20, 10};

% system memory
w = {12, 3};

% time it takes to transcribe the MS2 loops [sec]
t_MS2 = 30;

% alpha: length of the MS2 loop in time steps
alpha = cell(1,nGR);
for n = 1:nGR
    alpha{n} = t_MS2 / deltaT{n};
end

% transition rates [sec^-1]
R = cell([1, nGR]);
R{1} = [-0.0115,  0.0095,  0.0000; ...
         0.0115, -0.0180,  0.0440; ...
         0.0000,  0.0085, -0.0440];
R{2} = R{1};
R{3} = R{1};

% transition probabilities
A = cell([1, nGR]);
for i = 1:nGR
    A{i} = rate_to_prob(R{i}, deltaT{i});
end

% initial state pmf
pi0 = cell([1, nGR]);
pi0{1} = [0.15, 0.25, 0.6];
pi0{2} = pi0{1};
pi0{3} = pi0{1};

% emission rates [a.u. / sec]
r_emission = cell([1, nGR]);
r_emission{1} = [2, 65, 130];
r_emission{2} = r_emission{1};
r_emission{3} = r_emission{2};

% background noise [a.u.]
noise = cell([1, nGR]);
for i = 1:nGR
    noise{i} = 0;
end

% state conversion (which states are active in the two loci)
conv = cell([1, nGR]);
conv{1} = [1, 0; 0, 2; 3, 0];
conv{2} = conv{1};
conv{3} = conv{1};

% fluorescence per rna [a.u. / rna]
fluo_per_rna = 350;

% time it takes for another polymerase to load [sec]
promoter_pause = 5;
loading_rate = 1;
wait_time = promoter_pause + loading_rate;

% mandatory wait time between arrivals
flor = 0;

% total number of time points in a pooled data set
n_points_total = [ 3000, 3000, 3000];

% average number of points per trace;
seq_length = [50, 50, 50];

% number of traces in each AP position
n_traces = round(n_points_total./seq_length);

% ---------------- Save parameters into a structure -------------------

% structure to store all synthetic parameter values
synthetic_parameters = struct;
synthetic_parameters.K = K;
synthetic_parameters.w = w;
synthetic_parameters.t_MS2 = t_MS2;
synthetic_parameters.deltaT = deltaT;
synthetic_parameters.alpha = alpha;

synthetic_parameters.R = R;
synthetic_parameters.A = A;
synthetic_parameters.pi0 = pi0;
synthetic_parameters.noise = noise;
synthetic_parameters.r_emission = r_emission;

synthetic_parameters.n_points_total = n_points_total;
synthetic_parameters.seq_length = seq_length;
synthetic_parameters.n_traces = n_traces;

% structure to store synthetic data sets
data = struct;

% index variable for the current generated trace
set_ind = 1;
n_sets = 1;

% maximum number of EM iterations
n_steps_max = 1000;

% tolerance parameter for inference convergence
eps = 10^(-4);

% matrix to store the time resolution and set number pairs
set_count_mat = zeros(n_sets, nGR);

for i = 1:nGR
    for set = 1:n_sets
        fluo_data = cell([1 length(w)]);
        old_fluo_data = cell([1 length(w)]);
        for wi = 1:length(w)
            fluo_data{wi} = cell([1, n_traces(i)]);
            old_fluo_data{wi} = cell([1, n_traces(i)]);
        end
        state_data = cell([1, n_traces(i)]);
        for tr = 1:n_traces(i)
            fluo_gill = synthetic_two_spot_gillespie(seq_length(i), alpha{i}, ...
                K, w, R{i}, deltaT{i}, r_emission{i}, noise{i}, pi0{i}, ...
                fluo_per_rna, flor, conv{i});
            
            for wi = 1:length(w)
                fluo_data{wi}{tr} = fluo_gill(wi).fluo_MS2;
                old_fluo_data{wi}{tr} = fluo_gill(wi).fluo_MS2_orig_noise;
            end
            state_data{tr} = fluo_gill.naive_states;
        end

        data(set_ind).fluo_data = fluo_data;
        data(set_ind).old_fluo_data = old_fluo_data;
        data(set_ind).deltaT = deltaT{i};
        data(set_ind).deltaT_ind = i;
        data(set_ind).n_traces = n_traces(i);
        data(set_ind).alpha = alpha{i};
        data(set_ind).naive_states = state_data;
        data(set_ind).set = set;
        
        set_count_mat(set, i) = set_ind;
        set_ind = set_ind + 1;
        
        
    end
end


% extract the current date in a string format
file_name = ['eve_like_dT' int2str(deltaT{1}) 'w' int2str(w{1}) 'w' int2str(w{2}) ...
    'n1k']; 

% save the generated data into a '.mat' file
save(['../dat/synth_dat/' file_name '_data.mat'], 'data');
save(['../dat/synth_dat/' file_name '_params.mat'], 'synthetic_parameters');
