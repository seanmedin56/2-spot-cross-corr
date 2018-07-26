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
K = 2;

% time resolution
deltaT = 20;

% elongation times of each construct (in descending order)
elong_times = {7, 3};

% time it takes to transcribe the MS2 loops [sec]
t_MS2 = 30;

% alpha: length of the MS2 loop in time steps
alpha = t_MS2 / deltaT;

% transition rates [sec^-1]
R = [-0.001,  0.001; ...
      0.001, -0.001];

% initial state pmf
pi0 = [0, 1];

% emission rates [mRNA / sec]
r_emission = [.0001, .2];


% background noise [a.u.]
noise = 0;

% fluorescence per rna [a.u. / rna]
fluo_per_rna = 50;

% mandatory wait time between arrivals
flor = 0;

% number of traces
num_traces = 100;

% how many seconds each trace is taken for
time_per_trace = 1200;

% how many time points per trace
num_points = floor(time_per_trace / deltaT);

% ---------------- Save parameters into a structure -------------------

% structure to store traces and synthetic parameters
data = struct;
data.K = K;
data.w = elong_times;
data.t_MS2 = t_MS2;
data.deltaT = deltaT;
data.alpha = alpha;

data.R = R;
data.pi0 = pi0;
data.noise = noise;
data.r_emission = r_emission;
red_traces = cell(1, num_traces);
green_traces  = cell(1, num_traces);

for tr = 1:num_traces
    fluos = synthetic_two_spot_gillespie(num_points, alpha, ...
        K, elong_times, R, deltaT, r_emission, noise, pi0, ...
        fluo_per_rna, flor); 
    red_traces{tr} = fluos{2};
    green_traces{tr} = fluos{1};
end

data.red_traces = red_traces;
data.green_traces = green_traces;

% extract the current date in a string format
file_name = ['w7w3dt20test1']; 

% save the generated data into a '.mat' file
save(['../dat/synth_dat/' file_name '_data.mat'], 'data');
