function EOG_Cycles_Analysis()
% Function to extract EGO smooth pursuit data from ramps
% Release version 1/14/2014
% Paul Allen email: pdallen73@gmail.com
% Update 1/27/2014 - reads in ILD-ITD files - assumes always the same
% 'angles' for a given trial number
% Appears also to work for visual pursuit
% Lines 589 and 590 assignin('base','TargetMotion',positions.speakerAZ); assignin('base','EyeMotion',positions.eyesAZ)
% This code sends the data for the current plot of Target and Eye motion to the workspace (panel 2).
% where it can be extracted for plotting. It doesn't have the pursuit segments, but those can be added if needed for plots
% Update 1/28/2014 Can read Imaginary Pursuit
% Update 2/11/2014 Fixed bugs where crashed if trial duration too short or if ITD/ILD/imagined longer than 50
% Update 2/11/2014.2 Imagined and Headphone trials use %
%  [ data_struc.trial(#).state_var.TrialIndex -100 ] to identify trial number in lookup table
% Update 2/11/2014 Cludge fix on look up trial number, using modulo 50
% Update 5/20/2014 Open S-Lab data 
% Update 5/21/2014 Added function xml2struct in body of function
% Update 5/22/2014 Modify Rescale button to use ginput to get a mouseclick
%    to control the rescale
% Update 6/6/2014 Return the 'previous' button
% Branch 6/10/2014 for the Cycles program
% Updated smoothing window for velocities to be 40
% Update 6/17/2014 to port the new load XML code that checks the XML for number of channels and variable names
%                  Add last segment to smooth pursuit - and plot all if no saccades
% Update 6/24/2014 This version pretty much nails desaccading and allows
%   PursuitMotion to be accessed from the workspace for manipulation outside
%   the program
% Update 7/8/2014 Modify desaccading to remove 50 ms before and after
% saccade. Send snippet "start_time" to work space
% Update 12/16/2014 Add the 'saccade' at start of movement


mydata.MainFigure = figure('name', 'Data Viewer', 'numbertitle', 'off', 'menubar', 'none');

allData.num_trials=0;
set(mydata.MainFigure,'userdata',allData);

mydata.selecting=uicontrol(mydata.MainFigure,'style','text',...
    'string','Measuring and Selecting Trials...',...
    'position',[210 180 100 50], 'visible','off');

mydata.reading=uicontrol(mydata.MainFigure,'style','text',...
    'string','Reading Data from Labview file...',...
    'position',[210 180 100 50], 'visible','off');

mydata.filemenu = uimenu(mydata.MainFigure, 'Label', 'File');
mydata.openh5 = uimenu(mydata.filemenu, 'Label', 'Open A-Lab file', 'callback', {@LoadData mydata});

LoadData([],[],mydata)

end

function ListCallback(cbo,~,mydata)
plotThis=get(cbo,'Value');
allData=get(mydata.MainFigure,'userdata');
allData.currentPlot=plotThis;
set(mydata.MainFigure,'userdata',allData);
PlotTrial([],[],mydata,plotThis)
end

function Rescale(cbo,~,mydata)
scale = get(cbo,'String');
allData=get(mydata.MainFigure,'userdata');
plotThis=allData.currentPlot;
f=allData.trialList;
trialname=f{plotThis}{1};
allData.data.(trialname).scaleF = str2double(scale);
set(mydata.MainFigure,'userdata',allData);
%PlotTrial([],[],mydata,plotThis)
end

function CopyGainOffset(~,~,mydata)
allData=get(mydata.MainFigure,'userdata');
plotThis=allData.currentPlot;
f=allData.trialList;
trialname=f{plotThis}{1};
prevtrialname=f{plotThis-1}{1};
allData.data.(trialname).scaleF = allData.data.(prevtrialname).scaleF;
allData.data.(trialname).offsetF = allData.data.(prevtrialname).offsetF ;
set(mydata.MainFigure,'userdata',allData);
set(mydata.offsetFactdisp,'string',num2str(allData.data.(trialname).offsetF))
set(mydata.scaleFactdisp,'string',num2str(allData.data.(trialname).scaleF))
end

function ClickRescale(~,~,mydata) % redo this as a rescale function
allData=get(mydata.MainFigure,'userdata');
plotThis=allData.currentPlot;
f=allData.trialList;
trialname=f{plotThis}{1};

%msgbox('Click on the top figure near the saccade just prior to the movement','How to rescale')
[x,~] = ginput(1);

x = round(x);
spkstart = 0;
spkend = mean(allData.data.(trialname).speaker_pos(x-250:x+250));

eyesAZ = allData.data.(trialname).eye_pos;
eyestart = mean(eyesAZ(x-500:x-250));
eyeend = mean(eyesAZ(x+250:x+500));

scaleF = (spkend-spkstart)/(eyeend-eyestart);

offsetF = eyeend*scaleF-spkend;

allData.data.(trialname).scaleF = scaleF;
allData.data.(trialname).offsetF = offsetF ;


set(mydata.MainFigure,'userdata',allData);
set(mydata.offsetFactdisp,'string',num2str(allData.data.(trialname).offsetF))
set(mydata.scaleFactdisp,'string',num2str(allData.data.(trialname).scaleF))
end

function Reoffset(cbo,~,mydata)
offset = get(cbo,'String');
allData=get(mydata.MainFigure,'userdata');
plotThis=allData.currentPlot;
f=allData.trialList;
trialname=f{plotThis}{1};
allData.data.(trialname).offsetF = str2double(offset);
set(mydata.MainFigure,'userdata',allData);
%PlotTrial([],[],mydata,plotThis)
end

function LoadData(~,~,mydata)

choice = questdlg('What type of session will you be loading?', ...
    'Session Type', ...
    'Imaginary','Headphone (ITD or ILD)','Free Field','Free Field');


set(mydata.reading,'visible','on');drawnow

cd '\\smdnas\Paige-Lab\Paige Lab\Labs\A-Lab\Experiments\GSM';
%cd 'C:\Users\pallen\Dropbox\1312 Cloninger EOG reader' cd 'C:\Users\paul
%allen\Dropbox\1312 Cloninger EOG reader'
[xmlfile,PathName] = uigetfile('*hedr.xml','Select the S-Lab output file');

set(mydata.reading,'visible','on');drawnow
disp('Loading file')

out = xml2struct([PathName, xmlfile]);
assignin('base', 'var',out)
if ~strcmp(out.children(4).children(2).name, 'Sample_Rate')
    error('Check XML file - Sample_Rate not at Correct line')
end

scalefilename = [PathName, xmlfile(1:end-4),'_scale.mat'];
scalefile = dir(scalefilename);

if length(scalefile) == 1
    load(scalefilename) % variable scale is loaded
    disp('Loading previously saved scaling factors')
else
    disp('No previous scaling factors found')
    scale = zeros(500,2);
    scale(:,1) = 1;
end

sampleRate = str2double(out.children(4).children(2).children.data);

numChannels = length([out.children(4).children(8).children.children])+length([out.children(4).children(4).children.children]);

binFile = [PathName, xmlfile(1:end-8), 'smpl'];
fid = fopen(binFile,'r','b');
analogdigitaldata = fread(fid, [numChannels,inf], 'int16')';
fclose(fid); 

%[FileName,PathName] = uigetfile('C:\Users\pallen\Desktop\Cloninger - Read
%data into Matlab\_sample ramp data\*.0*','Select the A-Lab output file');
%data_struc = read_raw(PathName,FileName);




[B,A] = butter(12,0.1,'low');
trialstartsample = 1;

varName = cell((length(out.children(8).children)-1)/2-1,2);
j = 1;
for i = 2:2:length(out.children(8).children)-1
    varName{j,2} = i;
    varName{j,1} = out.children(8).children(i).name;
    j = j+1;
end

samplecountIndex = varName{find(strcmp('Sample_Count',varName)),2};
trialnumberindex = varName{find(strcmp('Trial_Index',varName)),2};
speaker_preAZIndex = varName{find(strcmp('Pre_Inr_P2',varName)),2};
speaker_stimAZIndex = varName{find(strcmp('Inr_Mov_P',varName)),2};
LVtrialnumberindex = varName{find(strcmp('LV_Trial',varName)),2};
numeventsindex = varName{find(strcmp('Event_Count',varName)),2};

for trialNum = 8:2:size(out.children,2)-1 % number of trials
    try
    %% try scale and offset so that mean of first and last 500ms matches eyes to speaker location
    trialendsample = str2double(out.children(trialNum).children(samplecountIndex).children.data);
    lookuptrial = str2double(out.children(trialNum).children(trialnumberindex).children.data);
    speaker_preAZ = str2double(out.children(trialNum).children(speaker_preAZIndex).children.data);
    speaker_destinationAZ = str2double(out.children(trialNum).children(speaker_stimAZIndex).children.data);
    
    Excel_Cluster(trialNum/2-3,:) = [...
        lookuptrial,... % Trial Number
        str2double(out.children(trialNum).children(LVtrialnumberindex).children.data),... %LV trial number
        trialendsample,... % Number of samples
        str2double(out.children(trialNum).children(numeventsindex).children.data),... % Number of events
        speaker_preAZ,... % where the speaker was at the start of the trial for calibration
        speaker_destinationAZ... % Where the inner speaker heads to
        ];
    if trialendsample-trialstartsample < 1000
       warning('Something Wierd with a short duration trial');
       disp(trialNum);
        continue
    end
      eyesAZ = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 6))*100/32768;
    %eyesEL = filtfilt(B,A,data_struc.trial(trialNum).data.LEEL);
    eyesEL = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 7))*100/32768;
    %speaker = filtfilt(B,A,data_struc.trial(trialNum).data.OSAZ*1.429-4.6);
    speaker = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 3))*100/32768/1.0282;
    
    trialstartsample = trialendsample+1;
  
 
    
    switch choice
        case 'Headphone (ITD or ILD)'
            % disp('Inferred angle from ITD/ILD will be used for Target Motion')
            [startITD, endITD, ITDrate] = lookupheadphoneITD(lookuptrial-100);
            duroftraceSamp = length(speaker);
            durofstartSec = 1; % assume a second of dwell time at the start
            durofmoveSec = abs(endITD-startITD)/ITDrate;
            durofmvtSamp = floor(durofmoveSec*sampleRate);
            startofmvtSamp = floor(durofstartSec*sampleRate);
            endofmvtSamp = startofmvtSamp+durofmvtSamp;
            speaker(1:startofmvtSamp) = startITD;
            speaker(startofmvtSamp:endofmvtSamp) = startITD:(endITD-startITD)/durofmvtSamp:endITD;
            speaker(endofmvtSamp:duroftraceSamp) = endITD;
        case 'Imaginary'
            % disp('Inferred angle from Imagined pursuit')
            [startpos, endpos, rate] = lookupimaginaryITD(mod(lookuptrial,50));
            duroftraceSamp = length(speaker);
            durofstartSec = 1; % assume a second of dwell time at the start
            durofmoveSec = abs(endpos-startpos)/rate;
            durofmvtSamp = floor(durofmoveSec*sampleRate);
            startofmvtSamp = floor(durofstartSec*sampleRate);
            endofmvtSamp = startofmvtSamp+durofmvtSamp;
            speaker(1:startofmvtSamp) = startpos;
            speaker(startofmvtSamp:endofmvtSamp) = startpos:(endpos-startpos)/durofmvtSamp:endpos;
            speaker(endofmvtSamp:duroftraceSamp) = endpos;
    end
    
        spkstart = mean(speaker(500:1000));
        spkend = 0;%mean(speaker(end-550:end-50));
        
        eyestart = mean(eyesAZ(500:1000));
        eyeend = mean(eyesAZ(2000:2500));%mean(eyesAZ(end-550:end-50));
        
        scaleF = (spkend-spkstart)/(eyeend-eyestart);
        if scaleF <0
            scaleF = 0.6;
        end
        
        offsetF = eyeend*scaleF-spkend;
    catch
        disp('problem with auto fitting')
        continue
    end
    
    data.(['trial_', num2str(trialNum/2-3)]).eye_pos= resample(eyesAZ,1000,sampleRate)  ;
    
    data.(['trial_', num2str(trialNum/2-3)]).eye_posEL= resample(eyesEL,1000,sampleRate)  ;
    %  data.(['trial_',
    %  num2str(trialNum)]).eye_pos=data_struc.trial(trialNum).data.LEAZ;
    data.(['trial_', num2str(trialNum/2-3)]).speaker_pos=resample(speaker, 1000, sampleRate);
    % data.(['trial_',
    % num2str(trialNum)]).speaker_pos=data_struc.trial(trialNum).data.OSAZ*1.429;
    data.(['trial_', num2str(trialNum/2-3)]).scaleF = scale(trialNum/2-3,1);
    data.(['trial_', num2str(trialNum/2-3)]).offsetF = scale(trialNum/2-3,2);
end

set(mydata.reading,'visible','off');drawnow

allData=get(mydata.MainFigure,'userdata');

allData.data=data;
allData.Subject = xmlfile(1:6);
allData.filename = xmlfile;
allData.currentPlot=1;
allData.num_trials=num2cell(1:trialNum/2-3);
set(mydata.MainFigure,'userdata',allData);
SelectTrials([],[],mydata)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
set(mydata.MainFigure,'menubar','figure')
end

function RemoveTrial(~,~,mydata)
allData=get(mydata.MainFigure,'userdata');
trialList=allData.trialList;
removethis=allData.currentPlot;
trialList(removethis)=[];
allData.trialList=trialList;
set(mydata.MainFigure,'userdata',allData);
% if length(trialList)>=removethis
%     set(mydata.listbox,'string',[trialList{1:end}],'value',removethis);
%     PlotTrial([],[],mydata,removethis)
% else
%     set(mydata.listbox,'string',[trialList{1:end}],'value',1);
%     PlotTrial([],[],mydata,1)
% end
end

function SelectTrials(~,~,mydata)
% set(mydata.selecting,'visible','on');drawnow
allData=get(mydata.MainFigure,'userdata');
f=fields(allData.data);
trialList= [];

for i = 1:length(f)
    meas=MeasureTrial(mydata,f{i},1,0);
    allData.data.(f{i}).meas=meas;
    trialList=[trialList {f(i)}];
end
if isempty(trialList)
    errordlg('No trials meet criteria','No Trials');
end
uicontrol('style','text','position',[10,180,120,20],'string',['Trials: ',num2str(length(trialList))])
allData.trialList = trialList;
allData.num_trials=length(trialList);

mydata.flagthistrace=uicontrol(mydata.MainFigure,'style','checkbox',...
    'position',[10,620,50,20],...
    'callback',{@flagthistrace mydata});

% mydata.flaglabel = uicontrol(mydata.MainFigure,'Style','text', ...
%     'Position', [70,620,80,20], 'String', 'Flagged Trial', 'visible', 'off');

mydata.removethis=uicontrol(mydata.MainFigure,'style','push','string',...
    'Remove Trial','position',[10,150,120,20],...
    'TooltipString','Removes the current trial from ths list',...
    'callback',{@RemoveTrial mydata});

mydata.Export=uicontrol('style','push','string','Export Data',...
    'TooltipString','Export the smooth and saccade segments',...
    'position',[15 50 100 20],'callback', {@Export mydata});

set(mydata.MainFigure,'userdata',allData);

mydata.scaleFactdisp= uicontrol(mydata.MainFigure,'Style','edit',...
    'Position', [10,590,50,20], 'Callback', {@Rescale mydata});

mydata.scalelabel = uicontrol(mydata.MainFigure,'Style','text', ...
    'Position', [70,590,80,20], 'String', 'Scale Factor');

mydata.offsetFactdisp= uicontrol(mydata.MainFigure,'Style','edit',...
    'Position', [10,560,50,20], 'Callback', {@Reoffset mydata});

mydata.offsetlabel = uicontrol(mydata.MainFigure,'Style','text', ...
    'Position', [70,560,80,20], 'String', 'Offset Factor');

mydata.copyGainOffset= uicontrol(mydata.MainFigure,'Style','pushbutton',...
    'Position', [100,530,80,20], 'String', 'Click Rescale', 'Callback', {@ClickRescale mydata});

mydata.copyGainOffset= uicontrol(mydata.MainFigure,'Style','pushbutton',...
    'Position', [10,530,80,20], 'String', 'Previous Scale', 'Callback', {@CopyGainOffset mydata});

mydata.listbox = uicontrol('Style', 'listbox',...
    'Position', [10,200,120,320], 'Callback', {@ListCallback mydata});


set(mydata.listbox,'string',[trialList{1:end}],...
    'Callback', {@ListCallback mydata},'visible','on');
PlotTrial([],[],mydata,1)

end

function Export(~,~,mydata)
directoryname = uigetdir('', 'Please select the folder to save the Excel analysis output into');
cd(directoryname)
allData=get(mydata.MainFigure,'userdata');
f=allData.trialList;
filename = [allData.filename(1:end-4),'_Slopes.xls'];
scalefilename = [allData.filename(1:end-4),'_scale.mat'];
%sheet = 1;% allData.filename(end-2:end);
% Export to Excel
xlswrite(filename, { 'Subject'	'Trial#'	'Start angle'	'End angle' 	'Velocity'	'Peak Pursuit'	'Mean Pursuit'	'Weighted Pursuit'	'# Pursuit segments'	'Fraction pursuit' 'Scale Factor' 'Offset'}, 'Sheet1', 'A1')

scale = zeros(length(f),2);
scale(:,1) = 1;

for i = 1: length(f)
    xlcell = ['A', num2str(i+1), ':L', num2str(i+1)];
    trialnumber=f{i}{1};
    meas=MeasureTrial(mydata,trialnumber,allData.data.(trialnumber).scaleF,allData.data.(trialnumber).offsetF );
    disp(['Analyzing trial ' num2str(i)])
    
    try %Summary Page
        validSegments = meas.pursuitDurations > 250;
        pursuitVelocitiesList = abs(meas.pursuitMeanVelocities(validSegments));
        
        pursuitfraction = sum(meas.pursuitAmplitudes(validSegments))/(meas.end_position-meas.start_position);
        numSegments = sum(validSegments);
        
        weightedMean = abs(sum(dot(meas.pursuitDurations(validSegments),pursuitVelocitiesList))/sum(meas.pursuitDurations(validSegments)));
        
                scale(i, 1) = allData.data.(trialnumber).scaleF;
        scale(i, 2) = allData.data.(trialnumber).offsetF;
        
        % disp({ allData.Subject,	trialnumber,	meas.start_position, 	meas.end_position,	...
        %     meas.ramp_speed,	max(pursuitVelocitiesList), mean(pursuitVelocitiesList)	,weightedMean,...
        %     numSegments,	pursuitfraction, allData.data.(trialnumber).scaleF, allData.data.(trialnumber).offsetF})
        
        xlswrite(filename, { allData.Subject,	trialnumber,	meas.start_position, ...
            meas.end_position,	meas.ramp_speed,	max(pursuitVelocitiesList), mean(pursuitVelocitiesList)	...
            ,weightedMean,	numSegments,	pursuitfraction, scale(i, 1),...
            scale(i, 2)}, 'Sheet1', xlcell)
    catch
        disp(['No valid data for ' trialnumber])
        xlswrite(filename, { allData.Subject,	trialnumber, 'No valid segments'},  'Sheet1', ['A', num2str(i+1), ':C', num2str(i+1)])
    end
    
    try % Smooth Pursuit Page
        outputPursuitsegments = {};
        outputPursuitsegments{1} = allData.Subject;
        outputPursuitsegments{2} = trialnumber;
        
        startList = meas.pursuitMovements_start(validSegments);
        endList = meas.pursuitMovements_end(validSegments);
        
        if isempty(endList)
            error('No smooth pursuit')
        end
        
        for j = 1:length(pursuitVelocitiesList)
            outputPursuitsegments{(j-1)*3+3} = startList(j);
            outputPursuitsegments{(j-1)*3+4} = endList(j);
            outputPursuitsegments{(j-1)*3+5} = pursuitVelocitiesList(j);
        end
        
        %disp(outputPursuitsegments)
        if j > 8
            xlcell = ['A', num2str(i),':A',char(64+(j-8)*3),num2str(i)];
        else
            xlcell = ['A', num2str(i),':',char(64+j*3+2),num2str(i)];
        end
        xlswrite(filename, outputPursuitsegments, 'Sheet2', xlcell)
    catch
        disp(['No valid pursuit data for ' trialnumber])
        xlswrite(filename, { allData.Subject,	trialnumber, 'No valid segments'},  'Sheet2', ['A', num2str(i+1), ':C', num2str(i+1)])
    end
    
    try % Smooth pursuit shorter than 250 ms
        outputPursuitsegments = {};
        outputPursuitsegments{1} = allData.Subject;
        outputPursuitsegments{2} = trialnumber;
        
        startList = meas.pursuitMovements_start(~validSegments);
        endList = meas.pursuitMovements_end(~validSegments);
        tooshortpursuits = abs(meas.pursuitMeanVelocities(~validSegments));
        
        if isempty(endList)
            error('No short smooth pursuit')
        end
        
        for j = 1:length(tooshortpursuits)
            outputPursuitsegments{(j-1)*3+3} = startList(j);
            outputPursuitsegments{(j-1)*3+4} = endList(j);
            outputPursuitsegments{(j-1)*3+5} = tooshortpursuits(j);
        end
        
        % disp(outputPursuitsegments)
        if j > 8
            xlcell = ['A', num2str(i),':A',char(64+(j-8)*3),num2str(i)];
        else
            xlcell = ['A', num2str(i),':',char(64+j*3+2),num2str(i)];
        end
        xlswrite(filename, outputPursuitsegments, 'Sheet3', xlcell)
    catch
        disp(['No valid short pursuit data for ' trialnumber])
        xlswrite(filename, { allData.Subject,	trialnumber, 'No valid segments'},  'Sheet3', ['A', num2str(i+1), ':C', num2str(i+1)])
    end
    
    try % Saccades
        outputPursuitsegments = {};
        outputPursuitsegments{1} = allData.Subject;
        outputPursuitsegments{2} = trialnumber;
        
        startList = meas.Gshifts_start;
        endList = meas.Gshifts_end;
        saccadeSpeed = meas.GSpeakVelocity;
        
        if isempty(endList)
            error('No saccades')
        end
        
        for j = 1:length(saccadeSpeed)
            outputPursuitsegments{(j-1)*3+3} = startList(j);
            outputPursuitsegments{(j-1)*3+4} = endList(j);
            outputPursuitsegments{(j-1)*3+5} = saccadeSpeed(j);
        end
        
        % disp(outputPursuitsegments)
        if j > 8
            xlcell = ['A', num2str(i),':A',char(64+(j-8)*3),num2str(i)];
        else
            xlcell = ['A', num2str(i),':',char(64+j*3+2),num2str(i)];
        end
        xlswrite(filename, outputPursuitsegments, 'Sheet4', xlcell)
    catch
        disp(['No valid saccade data for ' trialnumber])
        xlswrite(filename, { allData.Subject,	trialnumber, 'No valid segments'},  'Sheet4', ['A', num2str(i+1), ':C', num2str(i+1)])
    end
    
    
end

xlswrite(filename, { 'Export Complete'},  'Sheet1', ['A', num2str(i+2), ':A', num2str(i+2)])
disp([ 'Export finished. Check Excel file ' filename ])
save(scalefilename, 'scale')

end

function [m]= MeasureTrial(mydata,trialname, scaleF, offsetF)
allData=get(mydata.MainFigure,'userdata');

%horizontal
positions.eyesAZ=allData.data.(trialname).eye_pos*scaleF-offsetF;
positions.eyesEL=allData.data.(trialname).eye_posEL;
positions.speakerAZ=allData.data.(trialname).speaker_pos;

[velocities, ~]= calcva(positions,1);

try % calculate the trail start positions and velocity
%     m.start_position = round(mean(positions.speakerAZ(200:300)));
%     m.end_position = round(mean(positions.speakerAZ(end-300:end)));
    
    start_time = find(abs(velocities.vspeakerAZ(3500:end))>5,1)+3500-500;
    
    end_time = length(positions.speakerAZ); %find(abs(positions.speakerAZ-m.end_position)<2,1);
    
    rampVels = abs(velocities.vspeakerAZ(start_time:end_time));
    
    m.ramp_speed = round(mean(rampVels(rampVels>0.8*max(rampVels)))/5)*5;
    m.start_time = start_time;
   % m.end_time = end_time;
    
    % disp(m.ramp_speed)
    
    if start_time < 501
        startplot = 1;
    else
        startplot=start_time-500; % starts the plor 20 samples prior to the arm movement
    end
    
    if end_time+1500 > length(positions.speakerAZ)
        endplot=length(positions.speakerAZ);
    else
        endplot=end_time+1500; % ends the plot a short time after the end of the arm movement
    end
    
    if startplot>endplot
        disp('bad movement')
        startplot = 1;
        endplot = length(positions.speakerAZ);
    end
    
 
    
    %horizontal
    positions.eyesAZ=positions.eyesAZ(startplot:endplot); %eye posiitons
    positions.eyesEL=positions.eyesEL(startplot:endplot);
    positions.speakerAZ=positions.speakerAZ(startplot:endplot); %speaker positions
    
    [velocities, accelerations]= calcva(positions,1);
    
    eyesAZvelon = 50; % saccades if faster then 50
    eyesAZaccon = 1500;
    
   % PURvelon = 0; % Adam had 20 here - smooth pursuit goes down to slower here
   % PURaccon = 0;
    
    Gshifts_ind= union(find(abs(velocities.veyesAZ) > eyesAZvelon),...
        find(abs(accelerations.aeyesAZ) > eyesAZaccon)); %indicies of saccades
    %assignin('base', 'Gshifts_ind', Gshifts_ind)
    
    extrawidthofsaccade = 10; % assume 50 ms
    saccademask = 1:extrawidthofsaccade;
    for i = 1:length(Gshifts_ind)
        saccademask = [saccademask, (Gshifts_ind(i)-extrawidthofsaccade):(Gshifts_ind(i)+extrawidthofsaccade)];
    end
   % assignin('base', 'saccademask', saccademask)
    Gshifts_ind = unique(saccademask);
    
    Gshifts_ind = sort([475:525,Gshifts_ind]);
    
    pursuitMovements_ind = 1:length(velocities.veyesAZ);%union(find(abs(velocities.veyesAZ) > PURvelon),...
        %find(abs(accelerations.aeyesAZ) > PURaccon)); % indicies of pursuit
    
    pursuitMovements_ind = setxor(pursuitMovements_ind,Gshifts_ind); %cleans up
    
    %eliminate movements before 100ms
   % pursuitMovements_ind=pursuitMovements_ind(pursuitMovements_ind>100);
    
    Gshifts=find(diff(Gshifts_ind)>2);% reduced this number from 20, allows shorter saccades to be counted
    
    pursuitMovements= find(diff(pursuitMovements_ind)>2); % reduced this number from 20, allows shorter smooth to be counted
    
    numGshifts=length(Gshifts)+1;
    
    numPursuitMovements = length(pursuitMovements)+1;
    m.positions=positions;
    
    if ~isempty(Gshifts_ind)
        Gshifts_start=zeros(1,numGshifts);
        Gshifts_end=zeros(1,numGshifts);
        Gshifts_start(1)=Gshifts_ind(1);
        if numGshifts == 1
            Gshifts_end(1)= Gshifts_ind(end);
        else
            Gshifts_end(1)= Gshifts_ind(Gshifts(1));
        end
        if numGshifts > 1
            Gshifts_start(numGshifts)=Gshifts_ind(Gshifts(numGshifts-1)+1);
            Gshifts_end(numGshifts)=Gshifts_ind(end);
            if numGshifts > 2
                for i = 2:numGshifts-1
                    Gshifts_start(i)=Gshifts_ind(Gshifts(i-1)+1);
                    Gshifts_end(i)=Gshifts_ind(Gshifts(i));
                    
                end
            end
            if Gshifts_start(1) < 75
                Gshifts_start=Gshifts_start(2:end);
                Gshifts_end=Gshifts_end(2:end);
                numGshifts=numGshifts-1;
            end
        end
        
        m.numGshifts = numGshifts;
        m.Gshifts_start=Gshifts_start;
        m.Gshifts_end=Gshifts_end;
        m.Gdurations= Gshifts_end-Gshifts_start;
        m.Gamplitudes= positions.eyesAZ(Gshifts_end)-positions.eyesAZ(Gshifts_start);
        
        for i = 1:numGshifts
            if m.Gamplitudes(i) > 0
                m.GSpeakVelocity(i)=max(velocities.veyesAZ(Gshifts_start:Gshifts_end));
            else
                m.GSpeakVelocity(i)=min(velocities.veyesAZ(Gshifts_start:Gshifts_end));
            end
        end
        
    else
        m.numGshifts = 0;
    end
    
    if ~isempty(pursuitMovements)
        m.pursuitMovements_ind=pursuitMovements_ind;
        pursuitMovements_start=zeros(1,numPursuitMovements);
        pursuitMovements_end=zeros(1,numPursuitMovements);
        pursuitMovements_start(1)=pursuitMovements_ind(1);
        if ~isempty(pursuitMovements)
            pursuitMovements_end(1)= pursuitMovements_ind(pursuitMovements(1));
        else
            pursuitMovements_end(1)=pursuitMovements_ind(end);
        end
        if numPursuitMovements > 1
            pursuitMovements_start(numPursuitMovements)=pursuitMovements_ind(pursuitMovements(numPursuitMovements-1)+1);
            pursuitMovements_end(numPursuitMovements)=pursuitMovements_ind(end);
            if numPursuitMovements > 2
                for i = 2:numPursuitMovements-1
                    pursuitMovements_start(i)=pursuitMovements_ind(pursuitMovements(i-1)+1);
                    pursuitMovements_end(i)=pursuitMovements_ind(pursuitMovements(i));
                end
            end
        end
        m.numPursuitMovements=numPursuitMovements;
        m.pursuitMovements_start=pursuitMovements_start;
        m.pursuitMovements_end=pursuitMovements_end;
        for i =1:numPursuitMovements
            m.pursuitMeanVelocities(i)=mean(velocities.veyesAZ(pursuitMovements_start(i):pursuitMovements_end(i)));
            x = pursuitMovements_start(i):pursuitMovements_end(i);
            y = positions.eyesAZ(pursuitMovements_start(i):pursuitMovements_end(i))';
            p = polyfit(x,y,1);
            m.pursuitSlope(i) = p(1)*1000;
            %             m.pursuitHMeanVelocities(i)=round(mean(velocities.vhh(pursuitMovements_start(i):pursuitMovements_end(i))));
            %             m.pursuitEMeanVelocities(i)=round(mean(velocities.veh(pursuitMovements_start(i):pursuitMovements_end(i))));
        end
        m.pursuitDurations=pursuitMovements_end-pursuitMovements_start;
        m.pursuitAmplitudes= round(positions.eyesAZ(pursuitMovements_end)-positions.eyesAZ(pursuitMovements_start));
        
        for i = 1:length(m.pursuitDurations)
            if m.pursuitDurations(1) > 50
                m.PursuitStart = pursuitMovements_start(1)-20;
                break;
            end
        end
        
        
    end
    
catch
    disp('Trace too disordered to calculate movement parameters')
end

end

function PlotTrial(~,~,mydata,plotThis,~)
allData=get(mydata.MainFigure,'userdata');

if nargin<4
    plotThis=allData.currentPlot;
end

f=allData.trialList;
trialname=f{plotThis}{1};

%if length(trialList)>=removethis
set(mydata.listbox,'string',[allData.trialList{1:end}],'value',plotThis);
%    PlotTrial([],[],mydata,removethis)
%else
%     set(mydata.listbox,'string',[trialList{1:end}],'value',1);
%     PlotTrial([],[],mydata,1)
% end

scaleF = allData.data.(trialname).scaleF;
offsetF = allData.data.(trialname).offsetF;
set(mydata.offsetFactdisp,'string',num2str(offsetF))
set(mydata.scaleFactdisp,'string',num2str(scaleF))

meas=MeasureTrial(mydata,trialname, scaleF, offsetF);
assignin('base','MeasuresTrial',meas)

positions=meas.positions;
%[velocities, ~]= calcva(meas.positions,1);

%Make sure all figs have black background
whitebg('k')
% top plot- full record with EyeAZ in green, Speaker in blue, and zero in
% red
subplot(3,1,1) % top plot, the raw record

hold off
plot(allData.data.(trialname).eye_pos*scaleF-offsetF,'g')
hold on
plot(allData.data.(trialname).speaker_pos,'m')
plot([0,length(allData.data.(trialname).speaker_pos)], [0 0] , 'r')
plot([0,length(allData.data.(trialname).speaker_pos)], [4 4] , 'r:')
plot([0,length(allData.data.(trialname).speaker_pos)], [-4 -4] , 'r:')
title('Arm motion and EOG for full record')

% Second plot with the eye EL channel for blinks and sleep
subplot(3,1,2)
hold off
plot(positions.eyesEL*scaleF,'g')
title('EOG elevation channel during movement period')
hold on

% Main movement plot with the smooth pursuit and the saccades
subplot(3,1,3)

hold off
plot(positions.eyesAZ,'g')
title('Smooth pursuit and saccades extracted')
hold on
plot(positions.speakerAZ,'m')

assignin('base','TargetMotion',positions.speakerAZ)
assignin('base','EyeMotion',positions.eyesAZ)

netSaccade = 0;
shortSeg = 0;
PursuitMotion = [];
try % mark saccades in red
    i = -1;
    if meas.numGshifts > 0
        for i = 1:meas.numGshifts
            plot(meas.Gshifts_start(i):meas.Gshifts_end(i),positions.eyesAZ(meas.Gshifts_start(i):meas.Gshifts_end(i)),'linewidth',2.5,'color','r')
            plot(meas.pursuitMovements_start(i-shortSeg):meas.pursuitMovements_end(i-shortSeg),positions.eyesAZ(meas.pursuitMovements_start(i-shortSeg):meas.pursuitMovements_end(i-shortSeg))-netSaccade,'linewidth',2.5,'color','b')
            xes =  meas.pursuitMovements_start(i-shortSeg):meas.pursuitMovements_end(i-shortSeg);
            yes = positions.eyesAZ(meas.pursuitMovements_start(i-shortSeg):meas.pursuitMovements_end(i-shortSeg))-netSaccade;
            segment = [xes',yes];
            PursuitMotion = [PursuitMotion, segment'];
            % desaccade
            % get smooth segments before and after the saccade
            prevSmooth = positions.eyesAZ(meas.pursuitMovements_end(i-shortSeg)-100:meas.pursuitMovements_end(i-shortSeg));
            nextSmooth = positions.eyesAZ(meas.pursuitMovements_start(i+1-shortSeg):(meas.pursuitMovements_start(i+1-shortSeg)+100));
            
            % calculate their slopes and take average
            fit1 = polyfit(1:length(prevSmooth), prevSmooth',1);
            fit2 = polyfit(1:length(nextSmooth), nextSmooth',1);
            
            % multiply that slope by duration of saccade to get AZ drift during saccade
            drift = (meas.Gshifts_end(i)-meas.Gshifts_start(i))*(fit1(1)+fit2(1))/2;
            netSaccade = netSaccade + positions.eyesAZ(meas.Gshifts_end(i))- positions.eyesAZ(meas.Gshifts_start(i)) - drift ;

            if i == meas.numGshifts
                plot(meas.pursuitMovements_start(end):meas.pursuitMovements_end(end),positions.eyesAZ(meas.pursuitMovements_start(end):meas.pursuitMovements_end(end))-netSaccade,'linewidth',2.5,'color','b')
            end
        end
    else
        plot(positions.eyesAZ,'linewidth',2.5,'color','b')%% Plot whole trace because no saccades
    end
    assignin('base','PursuitMotion',PursuitMotion)
    assignin('base','start_time', meas.start_time)
catch
    disp(i)
    disp(meas.numGshifts)
    assignin('base','PursuitMotion',PursuitMotion)
    assignin('base','start_time', meas.start_time)
    disp('error Tried and failed to plot smooth pursuit')
end


end

function vels = ParabolicDiff(pos,n)
%define fcn Velocity - calculates velocity from position array (pos)and
%uses 2n +1 points in the parabola n must be an interger >= 1.
if nargin < 2;
    n = 1;
end

if ~isempty(find(size(pos)==0,1))
    vels=zeros(length(pos),1);
    return
end

%Calculate length of rows and columns in array pos
rows = length(pos(:,1));
cols = length(pos(1,:));
if cols > 1
    pos=pos(1,:)';
end

%discard first and last n points and for all columns calculate velocity
%using the k point parabola fcn and put the result into the array
%vels(i,j).

vels=zeros(length(pos),1);

for i = (n+1):(rows-n);
    vels(i) = 0;
    q = 0;
    for m = 1:n;
        qtemp = 2*(m^2);
        q = q+qtemp;
        velstemp =  (-m*pos(i-m)+m*pos(i+m));
        vels(i) = vels(i) + velstemp;
    end
    vels(i) = (vels(i)/q)*1000;
end

vels(1:n)=ones(n,1)*vels(n+1);
vels=vels';
end

function [velocities, accelerations]= calcva(positions,~)
n=40;
velocities.veyesAZ=ParabolicDiff(positions.eyesAZ,n);

if nargin> 1
    velocities.vspeakerAZ=ParabolicDiff(positions.speakerAZ,n);
end

accelerations.aeyesAZ=ParabolicDiff(velocities.veyesAZ,n);
end

function y=read_raw(pathname, filename)
% [filename,pathname] = uigetfile('*.*'); Matlab program to work with
%Host App: Client for PXI 11_08_02 PXI App: Global 11_06_02; Record module
%11_06_02 Adapted from Scott Seidman's version readLV and Anand Joshi's
%version LVdatareader_v1 Modified by JR: 09.12.05

% INPUT:       y=LVdatareader_JR1('H:\VN-Lab\ct0830L.000')
%
% OUTPUT:      label: 'Labview data acquisition program'
%              version: 3 labname: ''
%             userinfo: [1x66 char] username: 't c'
%                 date: '8/30/2005' time: '11:56 AM'
%             filename: 'ct0830L.000'
%          configfiles: 'Glab_reynolds_TRQ,Exptcfg LV7Glab1_llb v2'
%             expttype: ''
%        stickycomment: ''
%              comment: ''
%             scanrate: 1000
%         interchdelay: 3.0000e-006
%            numofchan: 9
%             channame: {1x9 cell}
%              channum: [0 1 2 5 6 8 9 14 15]
%       chansamplerate: [500 500 500 500 500 500 500 500 500]
%                 gain: [1x9 double]
%      scalemultiplier: [1x9 double]
%               offset: [0 0 0 0 0 0 0 0 0]
%          scalefactor: [18 36 40 4.8000 4.8000 7.7220 7.7220 1 2.8000]
%         typeofsignal: [0 1 0 0 0 0 0 0 0]
%      State_variables: {11x1 cell}
%           eventnames: {48x1 cell}
%     eventdescription: {48x1 cell}
%        num_of_trials: 4 triallocation: [4x1 double]
%                trial: [1x4 struct]
%to access events: y.trial(m).event.'eventname' to access state variables:
%y.trial(m).state_var.'state variable name' Each Event will have {lowfreq,
%high freq, lowres_occ, highres_occ} Each State variable will further have
%"value" as its field

if nargin == 0
    [filename, pathname] = uigetfile('*.0*');
end



% function y=LVdatareader_v1(filename)
tic
%Program to read the binary header and the data stored.
% filename=input('Input filename with path:   ','s');
fid=fopen([pathname,filename],'r','b');					% big-endian format read.
% Read the header
headerlength = fread(fid,1,'int32');            %read header length (4 byte= int32)
y.label = char((fread(fid,32,'char'))');		%Label info "Labview data acquisition program"
y.version = fread(fid,1,'float32');				%Version info to keep track of datafiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Head Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Lab name and version information
stringlgth1 = fread(fid,1,'int32');
y.labname = char((fread(fid,stringlgth1,'char'))');

status = fseek(fid,94,-1);
% read user supplied header, (subj name,date, time, filename)
userinfolength = fread(fid,1,'int32');
y.userinfo=char((fread(fid,userinfolength,'char'))');

rem = y.userinfo;
str = char(13); % ascill equivalent of carriage return '\r'
while true
    [token rem] = strtok(rem,str);
    
    if(strcmp(strtok(token),'Name'))
        y.username = token(8:end);
    elseif(strcmp(strtok(token),'Date'))
        y.date = token(8:end);
    elseif(strcmp(strtok(token),'Time'))
        y.time = token(8:end);
    elseif(strcmp(strtok(token),'Filename'))
        y.filename = token(12:end);
    end
    
    if isempty(token),break; end
end

%Configuration files information
status = fseek(fid,298,-1);
configlgth = fread(fid,1,'int32');
y.configfiles = char((fread(fid,configlgth,'char'))');

% Experiment type information
status = fseek(fid,452,-1);
stringlgth2 = fread(fid,1,'int32');
y.expttype = char((fread(fid,stringlgth2,'char'))');

%Sticky comments for the entire experiment
status = fseek(fid,506,-1);
stringlgth3 = fread(fid,1,'int32');
y.stickycomment = char((fread(fid,stringlgth3,'char'))');

% read user supplied comments from header block
status = fseek(fid,610,-1);
commentlength = fread(fid,1,'int32');
y.comment = char((fread(fid,commentlength,'char'))');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hardware Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
status = fseek(fid,1114,-1);
y.scanrate = fread(fid,1,'float32');			%Scan rate stored from the clk config
%assignin('base','scanrate',y.scanrate) ASSIGNIN sends scanrate to current
%workspace

y.interchdelay = fread(fid,1,'float32');		%Interchannel delay acqrd. from clk config.
y.numofchan = fread(fid,1,'int32');
chanlistlength = fread(fid,1,'int32');
chanacqd = char((fread(fid,chanlistlength,'char'))');

%parse channel names into cell string
startchanname=[ 1 findstr(',',chanacqd)+1];	%finds starting point of each chan name in the list
endchanname=[ findstr(',',chanacqd)-1 length(chanacqd)];	%finds end point of each chan name in the list

for i=1:y.numofchan
    channame(i)={chanacqd(startchanname(i):endchanname(i))};	%split chan names into array
    y.channame(i)=channame(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Channel Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
status = fseek(fid,1450,-1);

for i=1:y.numofchan
    y.channum(i)=fread(fid,1,'int32');
    y.chansamplerate(i) = fread(fid,1,'float32');
    y.gain(i)=fread(fid,1,'float32');
    %note, gain seems wrong!!  Says 0.5, Should be 1.0?5
    y.scalemultiplier(i)=fread(fid,1,'float32');
    %note, scalemultiplier seems wrong!!  Says 0.3052, Should be 1.0?
    y.offset(i)=fread(fid,1,'float32');
    y.scalefactor(i)=fread(fid,1,'float32');
    y.typeofsignal(i)=fread(fid,1,'int32');		%0 = position , 1 = Velocity, 2 = acceleration
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State variable declaration and description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_of_state_variables = fread(fid,1,'int32');
lgth_of_state_variable_string = fread(fid,1,'int32');
State_variables = char((fread(fid,lgth_of_state_variable_string,'char'))');
y.State_variables = parse_string(State_variables,',');  % parsing out the comma separated names into string array

% Max number of State variables possi
status = fseek(fid,2000-lgth_of_state_variable_string,0);%offset to read the event names

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Event/Epoch Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_of_eventstype = fread(fid,1,'int32');
lgth_of_event_string = fread(fid,1,'int32');
eventnames = char((fread(fid,lgth_of_event_string,'char'))');
y.eventnames = parse_string(eventnames,',');    %converting eventnames string into array of event names strings

status = fseek(fid,500-lgth_of_event_string,0);         %Offset to read the descriptions

%Events/ Epoch description
lgth_of_eventdescription = fread(fid,1,'int32');
eventdescription = char((fread(fid,lgth_of_eventdescription,'char'))'); % event description is newline separated
% converting the carriage return string description into array of string
y.eventdescription = parse_string(eventdescription,sprintf('\r'));

status = fseek(fid,2500-lgth_of_eventdescription,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trial information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y.num_of_trials = fread(fid,1,'int32');
% get starting location of the trial data
y.triallocation = fread(fid,y.num_of_trials,'int32');
offset = 4 + 4*y.num_of_trials;

for m = 1:y.num_of_trials
    
    %Trial location and reading
    status = fseek(fid,(y.triallocation(m)+offset),-1);
    %reading the trial
    total_trial_hdrlength = fread(fid,1,'int32');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %State Variable Structure: y.trial(m).state_var.'state var name' Assign
    %State variable names as fields of state_var structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_of_state_variables = fread(fid,1,'int32');
    y.trial(m).state_var = [];
    for j = 1:num_of_state_variables
        
        state_var_name_lgth(j) = fread(fid,1,'int32');
        tmpstatevarname = char((fread(fid,state_var_name_lgth(j),'char'))');
        if (isempty(tmpstatevarname))
            tmpstatevarname = ('empty string');
        end
        tmpstatevarname = strrep(tmpstatevarname,' ','_');
        tmpstatevarname = strrep(tmpstatevarname,'#','_');
        status = fseek(fid,(20 - state_var_name_lgth(j)),0); % Max character for State variable name  = 20.
        value = fread(fid,1,'float32');
        y.trial(m).state_var = setfield(y.trial(m).state_var,tmpstatevarname,value);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Events Structure: y.trial(m).event.'event name' Struct Event ={ Low
    %freq, High freq, lowres_occ, highres_occ} Each individual event
    %structure identified by its own name will be fields of the
    %trial.event. e.g. y.trial.event.STROT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    total_num_of_event_types = fread(fid,1,'int32');
    y.trial(m).event = [];
    for k = 1:total_num_of_event_types
        trialevent_name_lgth(k) = fread(fid,1,'int32');
        
        tmpeventname = char((fread(fid,trialevent_name_lgth(k),'char'))');
        if (isempty(tmpeventname))
            tmpeventname = ('empty string');
        end
        tmpeventname = strrep(tmpeventname,' ','_');
        tmpeventname = strrep(tmpeventname,'#','_');
        status = fseek(fid,(20 - trialevent_name_lgth(k)),0); % Max character for event name  = 20.
        eventdata.lowfreq = fread(fid,1,'float32');
        eventdata.highfreq = fread(fid,1,'float32');
        num_of_events_of_one_type = fread(fid,1,'int32');
        %Read the Lowres occurences and Highres occurences in array format
        for n=(1:num_of_events_of_one_type)
            eventdata.lowres_occ(n,1) = fread(fid,1,'int32');
            eventdata.highres_occ(n,1) = fread(fid,1,'int32');
        end
        y.trial(m).event = setfield(y.trial(m).event,tmpeventname,eventdata);
        clear eventdata; %reset eventdata within for loop for next event
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DATA STORAGE Data is stored in 2D array: Rows = # channels, Cols =
    %data/channel Datasize is a 2D array equal to {# channels, size of each
    %channel} For the current version, each channel will have the same
    %size. read data data=fread(fid,'int16');
    % to get into array of proper name in base workspace
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y.trial(m).datasize = fread(fid,y.numofchan,'int32');
    
    
    y.trial(m).data=[];
    for i=1:y.numofchan
        % parse data into channels
        datasplit =(fread(fid,y.trial(m).datasize(i),'int16')).*y.scalemultiplier(i).*y.scalefactor(i);
        % strip off leading blanks
        tmpname=char(channame(i));
        tmpname=tmpname(abs(tmpname)~=32);
        % create channel variable in base workspace
        %keyboard;
        y.trial(m).data = setfield(y.trial(m).data,tmpname,datasplit);
        
        %y.trial(m) = struct(tmpname,datasplit); assignin('caller',[tmpname
        %num2str(m)], datasplit);
        
    end
    %      y.trial(m).time =
    %      (1:y.trial(m).datasize(1))/y.chansamplerate(1);
    %y.trial(m).data = struct({[tmpname]},{[datasplit]});
    % pause;
end
t = toc;
%%%%%%%%%%%%%%%

fclose(fid);
end

function parsed_str=parse_string(str,delimiter)
% creates a 2D cell matrix of a string parsed by spaces and new lines
parsed_str={};
rownum=1;

while(not(isempty(str)))
    [line,str]=strtok(str,delimiter);
    colnum=1;
    parsed_str(rownum,colnum)={line};
    rownum=rownum+1;
end
end

function [startITD, endITD, ITDrate]= lookupheadphoneITD(trialNum)
conditionsList = [0,0,1;0,0,1;0,0,1;0,0,1;0,300,100;300,100,300;100,-300,300;...
    -300,300,100;300,-150,300;-150,50,100;50,300,300;300,-300,100;-300,300,100;...
    300,-150,300;-150,50,100;50,-350,100;-350,250,300;250,-150,100;-150,-300,300;...
    -300,300,300;300,-100,300;-100,300,100;300,-300,300;-300,0,300;0,0,1;0,0,1;...
    0,0,1;0,-300,100;-300,-100,300;-100,300,300;300,-300,100;-300,150,300;150,-50,100;...
    -50,-300,300;-300,300,100;300,-300,100;-300,150,300;150,-50,100;-50,350,100;...
    350,-250,300;-250,150,100;150,300,300;300,-300,300;-300,100,300;100,-300,100;...
    -300,300,300;300,0,300;0,0,1;0,0,1;0,0,1];
if trialNum > length(conditionsList)
    disp(trialNum)
    startITD = 0;
    endITD = 0;
    ITDrate = 1;
else
    startITD = conditionsList(trialNum, 1)/10;
    endITD = conditionsList(trialNum, 2)/10;
    ITDrate = conditionsList(trialNum, 3)/10;
end
end

function [startpos, endpos, rate]= lookupimaginaryITD(trialNum)
conditionsList = [0,0,1;20,20,1;-20,-20,1;0,0,1;30,-30,10;-30,30,10;30,-30,10;-30,30,10;...
    30,-30,10;-30,30,10;30,-30,10;-30,30,10;30,-30,10;-30,30,10;30,-30,30;-30,30,30;30,-30,30;...
    -30,30,30;30,-30,30;-30,30,30;30,-30,30;-30,30,30;30,-30,30;-30,30,30;30,-30,10;-30,30,10;...
    30,-30,10;-30,30,10;30,-30,10;-30,30,30;30,-30,30;-30,30,30;30,-30,30;-30,30,30;...
    30,-30,30;0,0,1;20,20,1;-20,-20,1];
if trialNum > length(conditionsList)
    disp(trialNum)
    startpos = 0;
    endpos = 0;
    rate = 1;
else
    startpos = conditionsList(trialNum, 1);
    endpos = conditionsList(trialNum, 2);
    rate = conditionsList(trialNum, 3);
end
end

function out = xml2struct(xmlfile)
% XML2STRUCT Read an XML file into a MATLAB structure.

% written By Douglas M. Schwarz, douglas.schwarz@kodak.com
xml = xmlread(xmlfile);

children = xml.getChildNodes;
for i = 1:children.getLength
   out(i) = node2struct(children.item(i-1));
end
end

function s = node2struct(node)

s.name = char(node.getNodeName);

if node.hasAttributes
   attributes = node.getAttributes;
   nattr = attributes.getLength;
   s.attributes = struct('name',cell(1,nattr),'value',cell(1,nattr));
   for i = 1:nattr
      attr = attributes.item(i-1);
      s.attributes(i).name = char(attr.getName);
      s.attributes(i).value = char(attr.getValue);
   end
else
   s.attributes = [];
end

try
   s.data = char(node.getData);
catch
   s.data = '';
end

if node.hasChildNodes
   children = node.getChildNodes;
   nchildren = children.getLength;
   c = cell(1,nchildren);
   s.children = struct('name',c,'attributes',c,'data',c,'children',c);
   for i = 1:nchildren
      child = children.item(i-1);
      s.children(i) = node2struct(child);
   end
else
   s.children = [];
end
end