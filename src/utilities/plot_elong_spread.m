function plot_elong_spread(elongs,intervals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
max_elong = max(elongs{1});
min_elong = min(elongs{1});
for i = 2:length(elongs)
    max_elong = max(max_elong,max(elongs{i}));
    min_elong = min(min_elong,min(elongs{i}));
end
interv = (max_elong - min_elong) / intervals;
edges = min_elong:interv:max_elong;
edges_for_plot = edges(1:end-1) + interv / 2;
disc_elongs = cell([1,length(elongs)]);
vals = cell([1,length(edges)-1]);
figure();
for i = 1:length(elongs)
    disc_elongs{i} = discretize(elongs{i},edges);
    actual_vals = zeros(1,length(edges) - 1);
    for val = 1:length(actual_vals)
        actual_vals(val) = length(disc_elongs{i}(disc_elongs{i} == val));
        vals{val} = [vals{val} actual_vals(val) / length(disc_elongs{i})];
    end
    plot(edges_for_plot,actual_vals);
    hold on
end
xlabel('elongation times');
ylabel('frequency');
title('distribution of elongation times');

figure();
means = zeros(1,length(vals));
stds = zeros(1,length(vals));
for i = 1:length(vals)
    means(i) = mean(vals{i});
    stds(i) = std(vals{i});
end
errorbar(edges_for_plot,means,stds);
xlabel('elongation times');
ylabel('frequency');
title('distribution of elongation times');
end
