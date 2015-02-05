close all
% Load the file
dirname = '\\smdnas\Paige-Lab\Paige Lab\Labs\A-Lab\Experiments\GSM\Cycles Extinction\';
[filename, pathname] = uigetfile([dirname , '*_trace.mat']);
load([pathname, filename]);
trialList = fieldnames(cyclestraces);
if not(length(trialList) == 24)
    error('Wrong number of traces in file')
end

conditionsList = {
    1,'Visual','Yes',20,20,-20,40,29500;
    2,'Visual','Yes',20,-40,40,80,69880;
    3,'Visual','Yes',20,10,-10,20,21200;
    4,'Visual','Yes',20,-20,-20,40,34300;
    5,'Visual','Yes',20,40,-40,80,52280;
    6,'Visual','Yes',20,-10,10,20,24000;
    7,'Visual','Yes',20,20,-20,40,39100;
    8,'Visual','Yes',20,-10,10,20,18400;
    9,'Visual','Yes',20,40,-40,80,61080;
    10,'Auditory','Yes',20,20,-20,40,29500;
    11,'Auditory','Yes',20,-40,40,80,69880;
    12,'Auditory','Yes',20,10,-10,20,21200;
    13,'Auditory','Yes',20,-20,20,40,34300;
    14,'Auditory','Yes',20,40,-40,80,52280;
    15,'Auditory','Yes',20,-10,10,20,24000;
    16,'Auditory','Yes',20,20,-20,40,39100;
    17,'Auditory','Yes',20,-10,10,20,18400;
    18,'Auditory','Yes',20,40,-40,80,61080;
    19,'Imaginary','Yes',20,20,-20,40,3000;
    20,'Imaginary','Yes',20,-40,40,80,3000;
    21,'Imaginary','Yes',20,10,-10,20,3000;
    22,'Imaginary','Yes',20,-20,20,40,3000;
    23,'Imaginary','Yes',20,40,-40,80,3000;
    24,'Imaginary','Yes',20,-10,10,20,3000
    };

figure('name', filename(1:end-15))

subplot(3,3,1)
title('Visual with Extinction - 0.5 Hz')
subplot(3,3,2)
title('Auditory with Extinction - 0.5 Hz')
subplot(3,3,3)
title('Imagined - 0.5 Hz')

subplot(3,3,4)
title('Visual with Extinction - 0.25 Hz')
subplot(3,3,5)
title('Auditory with Extinction - 0.25 Hz')
subplot(3,3,6)
title('Imagined - 0.25 Hz')

subplot(3,3,7)
title('Visual with Extinction - 0.125 Hz')
subplot(3,3,8)
title('Auditory with Extinction - 0.125 Hz')
subplot(3,3,9)
title('Imagined - 0.125 Hz')

Ysum = zeros(9,3001);
Yeach = zeros(24,3001);
xrange = [-1500:1500];

for i = 1:length(trialList)
    
    thistrace = cyclestraces.(trialList{i});
    
    if strcmp(conditionsList{i,2}, 'Auditory')
        col = 2;
    elseif strcmp(conditionsList{i,2}, 'Visual')
        col = 1;
    else
        col = 3; % imaginary
    end
    
    if conditionsList{i,7} ==20
        row = 1;
    elseif conditionsList{i,7} ==40
        row = 2;
    else
        row = 3; % imaginary
    end
    
    
    position = col+3*(row-1);
    
    xmin = conditionsList{i,8}-1500;
    xmax = conditionsList{i,8}+1500;
    if xmax > length(thistrace(:,4))
        xmax = length(thistrace(:,4));
        X = [xmin:xmax]';
        Y = thistrace(X,4);
        Y(length(X)+1:3001) = 0;
    else
        X = [xmin:xmax]';
        Y = thistrace(X,4);
    end
    
    Y = Y - Y(1); % adjust for the Y location at time = zero/extinction
    % Y = Y-outFile(i).Y(1); % adjust for starting position
    Y = Y *sign(conditionsList{i,5}); % adjust for direction
    % Y = Y/outFile(i).velocity; % adjust for velocity
    
    subplot(3,3, position)
    hold on
    plot(xrange,Y)
    
    %     offset = floor(X(1)); % say it starts at 5000
    %     dur = length(X);
    %     before = zeros(1,offset+5000);
    %     after = zeros(1,-offset+10000-dur);
    %     snip = [before, Y, after];
    Ysum(position,:) = Ysum(position,:) + Y';
    Yeach(i,:) = Y';
end
for position = 1:9
    subplot(3,3, position)
    hold on
    plot(xrange, Ysum(position,:)/3, 'r')
end
excelfilename = [pathname,filename(1:end-15), '_ext.xls'];
matfilename = [pathname,filename(1:end-15), '_ext.mat'];
xlswrite(excelfilename, Ysum',1)
save(matfilename, 'Yeach')
