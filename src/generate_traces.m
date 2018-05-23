% processes traces from 2 spot experiment

addpath('utilities/');

dataElong = open('C:\Users\seanm\Desktop\garcia_lab_stuff\elong_search\dat\real_data\20sSnippetStructures.mat');
s3 = dataElong.s3;
s5 = dataElong.s5;

% generates traces from data from 3' end of the gene that has more than
% minFrames
% frames of data collected

minFrames = 30;

preprocessed = struct;

% actually generates traces

idx = 1;
for k = 1:length(s3)
    if length(s3(k).newFrames) >= minFrames
        a = s3(k).SnippetGreen;
        b = s3(k).SnippetRed;
        t = s3(k).newFrames;
        
        greensum = [];
        rawredsum = [];
        greenmax = [];
        rawredmax = [];
        green98 = [];
        red98 = [];
        green95 = [];
        red95 = [];
        greengauss = [];
        redgauss = [];

        for i = 1:length(a)
            greensum(i) = sum(sum(a{i}));
            greenmax(i) = max(max(a{i}));
            rawredsum(i) = sum(sum(b{i}));
            rawredmax(i) = max(max(b{i}));
            a2 = reshape(a{i}, 1, numel(a{i}));
            b2 = reshape(b{i}, 1, numel(b{i}));
            green98(i) = prctile(a2, 98);
            red98(i) = prctile(b2, 98);
            green95(i) = prctile(a2, 95);
            red95(i) = prctile(b2, 95);
            
            %use gaussain technique taken from pipeline
            [greengauss(i), ~] = gauss_integral(double(a{i}));
            [redgauss(i), ~] = gauss_integral(double(b{i}));
            
        end
        
        preprocessed(idx).greensum = greensum;
        preprocessed(idx).greenmax = greenmax;
        preprocessed(idx).green98 = green98;
        preprocessed(idx).green95 = green95;
        preprocessed(idx).redsum = rawredsum;
        preprocessed(idx).redmax = rawredmax;
        preprocessed(idx).red98 = red98;
        preprocessed(idx).red95 = red95;
        preprocessed(idx).greengauss = greengauss;
        preprocessed(idx).redgauss = redgauss;
        preprocessed(idx).times = t * 20;
        idx = idx + 1;
    end
end

save('../dat/real_data/processed_data.mat', 'preprocessed');
