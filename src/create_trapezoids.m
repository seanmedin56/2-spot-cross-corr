% generates 2 trapezoids with various properties

%project name
project = 'slope diff';

% properties of PP7
intercept1 = 5;
slope1 = 50;
rise_time1 = 4;


% properties of MS2
intercept2 = 10;
slope2 = 300;
rise_time2 = 4;

time_steps = 60;
noise1 = normrnd(0, 0, [1 time_steps]);
noise2 = normrnd(0, 0, [1 time_steps]);

%creates trapezoid 1
trap1 = zeros([1 time_steps]);
for i = intercept1 + 1:rise_time1 + intercept1
    trap1(i) = slope1 + trap1(i-1);
end
trap1(rise_time1 + intercept1 + 1:end) = trap1(rise_time1 + intercept1);

% creates trapezoid 2
trap2 = zeros([1 time_steps]);
for i = intercept2 + 1:rise_time2 + intercept2
    trap2(i) = slope2 + trap2(i-1);
end
trap2(rise_time2 + intercept2 + 1:end) = trap2(rise_time2 + intercept2);
%trap1 = trap1 + noise1;
%trap2 = trap2 + noise2;
data = struct;
data.trap1 = trap1;
data.trap2 = trap2;
data.intercept1 = intercept1;
data.intercept2 = intercept2;
data.slope1 = slope1;
data.slope2 = slope2;
data.rise_time1 = rise_time1;
data.rise_time2 = rise_time2;

save(['../dat/' project], 'data');

figure;
plot(trap1, 'r');
hold on;
plot(trap2, 'g');
title('Trace');
xlabel('time steps');
ylabel('Spot Intensity (A.U)');

trace1 = cell([1 1]);
trace2 = cell([1 1]);
trace1{1} = trap1;
trace2{1} = trap2;
max_delay = 30;
cross_corr_analysis;
