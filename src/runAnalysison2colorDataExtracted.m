%file from Aramando for getting traces from his 2 spot data

mcp5frameinfo = load('D:\Data\Armando\livemRNA\Data\Dropbox\DynamicsResults\2017-04-25-P2P_MS2-PP7-LacZ_20secc2\FrameInfo.mat');
mcp5frameinfo = mcp5frameinfo.FrameInfo;
dt = mcp5frameinfo(2).Time;

exppath = 'C:\Users\ArmandoReimer\Documents\Personal-Repository\Teaching\physio2017';
data = open([exppath, '\20sCalibrationSnippetStructures.mat']);
s3 = data.s3;
s5 = data.s5;

%%%%
%Make calibration curve.
greenvals = [];
redvals = [];
greenvalss = [];
redvalss = [];
for k = 1:length(s3)
    if length(s3(k).Frame) > 20
        
        a = s3(k).SnippetGreen;
        b = s3(k).SnippetRed;
        greensum = [];
        redsum = [];
        greenmax = [];
        redmax = [];

        for i = 1:length(a)
            greenvalss = [greenvalss, sum(sum(a{i}))];
            greenvals = [greenvals, max(max(a{i}))];
            redvalss = [redvalss,sum(sum(b{i}))];
            redvals = [redvals,max(max(b{i}))];
        end
    end
end


figure(1)
scatter(redvals, greenvals)
xlabel('red')
ylabel('green')
title('calibration curve (all points from long traces). max')
figure(2)
scatter(redvalss, greenvalss)
xlabel('red')
ylabel('green')
title('calibration curve (all points from long traces). sum')


%%%%

dataElon = open([exppath, '\20sSnippetStructures.mat']);
s3 = dataElon.s3;
s5 = dataElon.s5;

gaps = [];
for k = 1:length(s3)
    if length(s3(k).newFrames) > 35
        
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
                     
        %do the calibration
        
        redsum = (3.*rawredsum) - 900;
        redmax = 3.*rawredmax - 18;      
        
        figure(3)
        plot(t, greensum, 'g')
        hold on
        plot(t, redsum, 'r');
        title(['Sum trace:', num2str(k)]);
        legend('3''', '5''')
        hold off
        
        figure(4)
        plot(t, greenmax, 'g')
        hold on
        plot(t, redmax, 'r');
        title(['Max trace:', k]);
        legend('3''', '5''')
        hold off
        
        %%%%Do a cross-correlation?
        
        [ccor, lag] = xcorr(rawredsum, greenmax);
        %smooth it in preparation for two derivatives. probably pchip
        %but will compare spline too
        oversamplingrate = 1/10;
        xq1 = lag(1):oversamplingrate:lag(end);
        p = pchip(lag,ccor,xq1);
        p = smooth(p, 75);
        figure(5)
        plot(xq1, p);
        title('Cross-correlation of red (cal.) and green max values');
        xlabel('Lag time (frames)')
        ylabel('Cross-correlation')
        hold on
        plot(lag, ccor, 'rx');
        legend('pchip, 10x sample', 'raw data')
        hold off
        
        figure(6)
        dc = diff(p);
        dxq1 = (xq1(1)+(oversamplingrate/2)):oversamplingrate:(xq1(end)-(oversamplingrate/2)); %we lost 1 sample with the derivative. 
        %this now samples in between the previous timepoints.
        plot(dxq1, dc, 'rx')
        title('First derivative of smoothed correlation')
        hold on
        plot(dxq1, smooth(dc, 75))
        hold off
        
        figure(7)
        ddc = diff(dc);
        ddxq1 = dt/60*(xq1(1)+oversamplingrate:oversamplingrate:xq1(end)-oversamplingrate); %we lost 1 sample with the derivative. 
        plot(ddxq1, ddc, 'rx')
        title('Second derivative of smoothed correlation')
        hold on
        ddc = smooth(ddc, 50);
        plot(ddxq1, ddc)
        xlabel('Time lag (min)')
        [m, i] = max(ddc);
        scatter(ddxq1(i), m, 'k', 'linewidth', 1);
        text(ddxq1(i), m, [num2str(ddxq1(i)), ' mins'],'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
        hold off
        
        elt = ddxq1(i);
        genelength = 1.4 + 3; 
        elr = elt / genelength; %kbs/min

        %%%%
        newgap = 0;
        gaps = [gaps, newgap];
    end
end


