dataElong = open('C:\Users\seanm\Desktop\garcia_lab_stuff\elong_search\dat\real_data\20sSnippetStructures.mat');
s3 = dataElong.s3;
s5 = dataElong.s5;

% generates traces from data from 3' end of the gene that has more than
% minFrames
% frames of data collected

minFrames = 35;

preprocessed = struct;

% actually generates traces

idx = 1;
for k = 1:length(s3)
    if length(s3(k).newFrames) > minFrames
        a = s3(k).SnippetGreen;
        b = s3(k).SnippetRed;
        t = s3(k).newFrames;
        
        greensum = [];
        rawredsum = [];
        greenmax = [];
        rawredmax = [];

        for i = 1:length(a)
            greensum(i) = sum(sum(a{i}));
            greenmax(i) = max(max(a{i}));
            rawredsum(i) = sum(sum(b{i}));
            rawredmax(i) = max(max(b{i}));
        end
        
        preprocessed(idx).greensum = greensum;
        preprocessed(idx).greenmax = greenmax;
        preprocessed(idx).redsum = rawredsum;
        preprocessed(idx).redmax = rawredmax;
        idx = idx + 1;
    end
end

save('../dat/real_data/processed_data.mat', 'preprocessed');
