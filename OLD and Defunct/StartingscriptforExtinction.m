close all
dirname = '\\delphi\Paige-Lab\Paige Lab\Users\Christy\GSM_extinction\Analyzed Data';
filelist = dir([dirname, '\*.mat']);

for j = 1:length(filelist)
    load([dirname , '/', filelist(j).name])
    figure('name', filelist(j).name)
    
    subplot(2,4,1)
    title('Visual with Extinction - 10 deg/sec')
    subplot(2,4,2)
    title('Auditory with Extinction - 10 deg/sec')
    subplot(2,4,3)
    title('Visual No Extinction - 10 deg/sec')
    subplot(2,4,4)
    title('Auditory No Extinction - 10 deg/sec')
        subplot(2,4,5)
    title('Visual with Extinction - 20 deg/sec')
    subplot(2,4,6)
    title('Auditory with Extinction - 20 deg/sec')
    subplot(2,4,7)
    title('Visual No Extinction - 20 deg/sec')
    subplot(2,4,8)
    title('Auditory No Extinction - 20 deg/sec')

    Ysum = zeros(8,15000);

    for i = 1:length(outFile)
              
        if strcmp(outFile(i).stimtype, 'Warm-up')
            continue
        end
        
        if isempty(outFile(i).Y)
            continue
        else
            Y = outFile(i).Y;
        end
        
        if strcmp(outFile(i).stimtype, 'Auditory')
            col = 2;
        else
            col = 1;
        end
        if outFile(i).extinctiontime == -1
            row = 2;
            X = outFile(i).X;
        else
            row = 0;
            X = outFile(i).X-outFile(i).extinctiontime;
        end
        if outFile(i).velocity == 10
        position = col+row;
        else
            position = col+row+4;
        end
        
        I = find(abs(X) < 1,1);
        if isempty(I)
            I = 1;
        end
        
        Y = Y - Y(I); % adjust for the Y location at time = zero/extinction
        % Y = Y-outFile(i).Y(1); % adjust for starting position
        Y = Y *sign(outFile(i).startpos); % adjust for direction
       % Y = Y/outFile(i).velocity; % adjust for velocity
             
        subplot(2,4, position)
        hold on
        plot(X,Y)
        
        offset = floor(X(1)); % say it starts at 5000
        dur = length(X);
        before = zeros(1,offset+5000);
        after = zeros(1,-offset+10000-dur);
        snip = [before, Y, after];
        Ysum(position,:) = Ysum(position,:) + snip;
    end
    for position = 1:8
        subplot(2,4, position)
        hold on
        plot((-5000:9999), Ysum(position,:)/9, 'r')
    end
    excelfilename = [filelist(j).name(1:end-3), 'xls'];
    xlswrite(excelfilename, Ysum',1)
end
