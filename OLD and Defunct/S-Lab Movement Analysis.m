function HeadFree_Ramp_Analysis()
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
% Update 6/12/2014 Real variables in XML to robustly encode spreadsheet columns, not hardcoded
% Update 6/17/2014 Read the config portion of the XML to determine number of Analog channels
% Update 10/1/2014 Branch to work with head-free gaze pursuit ramps
%   Labview channel 3 is inner Az; 7 is Elmar right horizontal;
%                   8 is Elmar Right vertical; 9 is Elmar Left horizontal;
%                   10 is Elmar Left vertical;13 is Head Yaw
%                   14 is Head Pitch; 15 is Head Roll;
%   Goal is export analysed data
%   Removed headphone option
% Update 11/11/2014 Calibrate head motion from first several blocks
% 11-19-2014 this line added to test git
% Update 12/2/2014 added calculation of head start move
% Update 12/3/2014 add in gaze start and export the times

%%\\smdnas\paige-lab\Paige Lab\Labs\A-Lab\RunTables\Gaze Speaker Movment

% Update 12/9/14 Rewrote the saccade detection section and improved
% stability of the ploting
% Update 12/11/2014 Major update to Excel Export and graphical interface
% Update 12/16/2014 Add .scale file to save the scale parameters; add head
% start and gaze start times to the summary page. Removed old comments

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
PlotTrial([],[],mydata,plotThis);
end

function Rescale(cbo,~,mydata)
scale = get(cbo,'String');
allData=get(mydata.MainFigure,'userdata');
plotThis=allData.currentPlot;
f=allData.trialList;
trialname=f{plotThis}{1};
allData.data.(trialname).scaleF = str2double(scale);
set(mydata.MainFigure,'userdata',allData);
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
set(mydata.offsetFactdisp,'string',num2str(allData.data.(trialname).offsetF));
set(mydata.scaleFactdisp,'string',num2str(allData.data.(trialname).scaleF));
end

function ClickRescale(~,~,mydata) % redo this as a rescale function
allData=get(mydata.MainFigure,'userdata');
plotThis=allData.currentPlot;
f=allData.trialList;
trialname=f{plotThis}{1};

[x,~] = ginput(1);

x = round(x);
spkstart = 0;
spkend = mean(allData.data.(trialname).speaker(x-250:x+250));

scaleAZ = allData.data.scaleAZ;
offsetAZ = allData.data.offsetAZ;
headpos = mean(allData.data.(trialname).Head_Yaw(x-500:x+250))*scaleAZ+offsetAZ;

eyesAZ = (allData.data.(trialname).eyesAZ_R+allData.data.(trialname).eyesAZ_L)/2;
eyestart = mean(eyesAZ(x-500:x-250));
eyeend = mean(eyesAZ(x+250:x+500));

scaleF = (spkend-spkstart)/(eyeend-eyestart);

offsetF = -headpos-eyestart*scaleF;

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
end

function LoadData(~,~,mydata)
set(mydata.reading,'visible','on');drawnow
cd '\\smdnas\Paige-Lab\Paige Lab\Labs\A-Lab\Experiments\GSM';
[xmlfile,PathName] = uigetfile('*Head_Free.hedr.xml','Select the S-Lab output file');

set(mydata.reading,'visible','on');drawnow
disp('Loading file')

out = xml2struct([PathName, xmlfile]);
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

numChannels = length([out.children(4).children(8).children.children])+length([out.children(4).children(4).children.children]);
binFile = [PathName, xmlfile(1:end-8), 'smpl'];

fid = fopen(binFile,'r','b');
analogdigitaldata = fread(fid, [numChannels,inf], 'int16')';
fclose(fid);

[B,A] = butter(12,0.1,'low');
trialstartsample = 1;

varName = cell((length(out.children(8).children)-1)/2-1,2);
j = 1;
for i = 2:2:length(out.children(8).children)-1
    varName{j,2} = i;
    varName{j,1} = out.children(8).children(i).name;
    j = j+1;
end

samplecountIndex = varName{strcmp('Sample_Count',varName),2};
trialnumberindex = varName{strcmp('Trial_Index',varName),2};
speaker_preAZIndex = varName{strcmp('Pre_Inr_P2',varName),2};
speaker_stimAZIndex = varName{strcmp('Inr_Mov_P',varName),2};
speaker_velAZIndex = varName{strcmp('Inr_Mov_V',varName),2};
LVtrialnumberindex = varName{strcmp('LV_Trial',varName),2};
numeventsindex = varName{strcmp('Event_Count',varName),2};
inneroffset = -0.7109;% varName{find(strcmp('Inr_Offset',varName)),2}/100;

Excelinfo = zeros((size(out.children,2)-1)/2-4, 7);

for trialNum = 8:2:size(out.children,2)-1 % number of trials
    trialendsample = str2double(out.children(trialNum).children(samplecountIndex).children.data);
    lookuptrial = str2double(out.children(trialNum).children(trialnumberindex).children.data);
    speaker_preAZ = str2double(out.children(trialNum).children(speaker_preAZIndex).children.data);
    speaker_destinationAZ = str2double(out.children(trialNum).children(speaker_stimAZIndex).children.data);
    speaker_velocityAZ = str2double(out.children(trialNum).children(speaker_velAZIndex).children.data);
    
    Excelinfo(trialNum/2-3,:) = [...
        lookuptrial,... % Trial Number
        str2double(out.children(trialNum).children(LVtrialnumberindex).children.data),... %LV trial number
        trialendsample,... % Number of samples
        str2double(out.children(trialNum).children(numeventsindex).children.data),... % Number of events
        speaker_preAZ,... % where the speaker was at the start of the trial for calibration
        speaker_destinationAZ,... % Where the inner speaker heads to
        speaker_velocityAZ... % how fast it went
        ];
    if trialendsample-trialstartsample < 1000
        warning('Something wierd with a short duration trial');
        disp(trialNum);
        continue
    end
    speaker = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 3))*100/32768/1.02286 +inneroffset;
    eyesAZ_R = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 8))*100/32768/3;
    eyesAZ_L = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 10))*100/32768/3;
    eyesEL_R = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 9))*100/32768/3;
    eyesEL_L = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 11))*100/32768/3;
    
    Head_Yaw = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 14))*0.01;
    Head_Pitch = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 15))*0.01;
    Head_Roll = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 16))*0.01;
    
    trialstartsample = trialendsample+1;
    
    data.(['trial_', num2str(trialNum/2-3)]).speaker = speaker;
    
    data.(['trial_', num2str(trialNum/2-3)]).eyesAZ_R = eyesAZ_R;
    data.(['trial_', num2str(trialNum/2-3)]).eyesAZ_L = eyesAZ_L;
    data.(['trial_', num2str(trialNum/2-3)]).eyesEL_R = eyesEL_R;
    data.(['trial_', num2str(trialNum/2-3)]).eyesEL_L = eyesEL_L;
    
    data.(['trial_', num2str(trialNum/2-3)]).Head_Yaw = Head_Yaw;
    data.(['trial_', num2str(trialNum/2-3)]).Head_Pitch = Head_Pitch;
    data.(['trial_', num2str(trialNum/2-3)]).Head_Roll = Head_Roll;
    
    data.(['trial_', num2str(trialNum/2-3)]).scaleF = scale(trialNum/2-3,1);
    data.(['trial_', num2str(trialNum/2-3)]).offsetF = scale(trialNum/2-3,2);
end

iOne = ['trial_',num2str(find(Excelinfo(:,1) == 1, 1, 'last' ))];
iTwo = ['trial_',num2str(find(Excelinfo(:,1) == 2, 1, 'last' ))];
iThree = ['trial_',num2str(find(Excelinfo(:,1) == 3, 1, 'last' ))];
iFour = ['trial_',num2str(find(Excelinfo(:,1) == 4, 1, 'last' ))];
iFive = ['trial_',num2str(find(Excelinfo(:,1) == 5, 1, 'last' ))];

data.scaleAZ = 40/(mean(data.(iThree).Head_Yaw(end-1000:end))-mean(data.(iOne).Head_Yaw(end-1000:end)));
data.offsetAZ = -mean(data.(iTwo).Head_Yaw(end-1000:end))*data.scaleAZ;
data.scaleEL = 20/(mean(data.(iFour).Head_Pitch(end-1000:end))-mean(data.(iFive).Head_Pitch(end-1000:end)));
data.offsetEL = -mean(data.(iTwo).Head_Pitch(end-1000:end))*data.scaleEL;

set(mydata.reading,'visible','off');drawnow
assignin('base', 'Excelinfo',Excelinfo)
allData=get(mydata.MainFigure,'userdata');
allData.Excelinfo = Excelinfo;
allData.data=data;
allData.Subject = xmlfile(1:6);
allData.filename = xmlfile;
allData.currentPlot=6;
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
end

function SelectTrials(~,~,mydata)
allData=get(mydata.MainFigure,'userdata');
f=fields(allData.data);
trialList= cell(length(f)-4,1);
trialkeyList = cell(length(f)-4,1);

for i = 1:length(f)-4
    trialList(i) ={f(i)};
    trialkeyList{i} = ['Trial_', num2str(allData.Excelinfo(i,1)), '_LVtrial_', num2str(allData.Excelinfo(i,2))];
end

if isempty(trialList)
    errordlg('No trials meet criteria','No Trials');
end

uicontrol('style','text','position',[10,180,120,20],'string',['Trials: ',num2str(length(trialList))])
allData.trialList = trialList;
allData.trialkeyList = trialkeyList;
allData.num_trials=length(trialList);

mydata.flagthistrace=uicontrol(mydata.MainFigure,'style','checkbox',...
    'position',[10,620,50,20],...
    'callback',{@flagthistrace mydata});

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


set(mydata.listbox,'string',[trialkeyList{1:end}],...
    'Callback', {@ListCallback mydata},'visible','on');
PlotTrial([],[],mydata)

end

function Export(~,~,mydata)
directoryname = uigetdir('', 'Please select the folder to save the Excel analysis output into');
disp('Starting Excel Export')
cd(directoryname)
allData=get(mydata.MainFigure,'userdata');
f=allData.trialList;
filename = [allData.filename(1:end-4),'_Slopes.xls'];
scalefilename = [allData.filename(1:end-4),'_scale.mat'];

% Export to Excel
xlswrite(filename, { 'Subject'	'LVTrial#'  'XLTrial'	'Start angle'	'End angle' 	'Velocity'	'Peak Gaze Pursuit'	'Mean Gaze Pursuit'	'Weighted Gaze Pursuit'	'# Pursuit segments'	'Fraction pursuit > 250ms'  'Fraction all pursuit' 'Head Start' 'Gaze Start' 'Scale Factor' 'Offset'}, 'Sheet1', 'A1')
xlswrite(filename, { 'Subject'	'LVTrial#'  'XLTrial'	'Segment Start'  'Segment End'  'Segment Velocity' 	'Pursuit Segments > 250 ms'	}, 'Sheet2', 'A1')
xlswrite(filename, { 'Subject'	'LVTrial#'  'XLTrial'	'Segment Start'  'Segment End'  'Segment Velocity' 	'Pursuit Segments < 250 ms'	}, 'Sheet3', 'A1')
xlswrite(filename, { 'Subject'	'LVTrial#'  'XLTrial'	'Segment Start'  'Segment End'  'Segment Velocity' 	'Saccade Segments'	}, 'Sheet4', 'A1')

scale = zeros(length(f),2);
scale(:,1) = 1;

for i = 1:length(f)
    
    xlcell = ['A', num2str(i+1), ':P', num2str(i+1)];
    trialnumber=f{i}{1};
    xltrialnum = allData.Excelinfo(i,1);
    if allData.Excelinfo(i,7) == 0 % Skip if zero velocity
        disp(['No speaker movement Trial ' num2str(i)])
        xlswrite(filename, { allData.Subject,	trialnumber, xltrialnum, 'No speaker movement this trial'},  'Sheet1', ['A', num2str(i+1), ':D', num2str(i+1)])
        xlswrite(filename, { allData.Subject,	trialnumber, xltrialnum, 'No speaker movement this trial'},  'Sheet2', ['A', num2str(i+1), ':D', num2str(i+1)])
        xlswrite(filename, { allData.Subject,	trialnumber, xltrialnum, 'No speaker movement this trial'},  'Sheet3', ['A', num2str(i+1), ':D', num2str(i+1)])
        xlswrite(filename, { allData.Subject,	trialnumber, xltrialnum, 'No speaker movement this trial'},  'Sheet4', ['A', num2str(i+1), ':D', num2str(i+1)])
        continue
    end
    
    meas=MeasureTrial(mydata,trialnumber,allData.data.(trialnumber).scaleF,allData.data.(trialnumber).offsetF );
    disp(['Analyzing trial ' num2str(i)])
    startAngle = allData.Excelinfo(i,5);
    endAngle = allData.Excelinfo(i,6);
    velocity = allData.Excelinfo(i,7);
    
    try %Summary Page - Sheet 1
        validSegments = meas.pursuitDurations > 250;
        pursuitVelocitiesList = abs(meas.pursuitMeanVelocities(validSegments));
        pursuitfraction = sum(meas.pursuitAmplitudes(validSegments))/(endAngle-startAngle);
        allpursuitfraction = sum(meas.pursuitAmplitudes)/(endAngle-startAngle);
        numSegments = sum(validSegments);
        headstart = meas.HeadStart-meas.TargetStart;
        gazestart = meas.GazeStart-meas.TargetStart;
        scale(i, 1) = allData.data.(trialnumber).scaleF;
        scale(i, 2) = allData.data.(trialnumber).offsetF;
        
        weightedMean = abs(sum(dot(meas.pursuitDurations(validSegments),pursuitVelocitiesList))/sum(meas.pursuitDurations(validSegments)));
        
        xlswrite(filename, { allData.Subject,	trialnumber,	xltrialnum, startAngle, ...
            endAngle,	velocity,	max(pursuitVelocitiesList), mean(pursuitVelocitiesList(2:end))	...
            ,weightedMean,	numSegments,	pursuitfraction, allpursuitfraction, headstart, gazestart, scale(i, 1),...
            scale(i, 2)}, 'Sheet1', xlcell)
    catch
        disp(['Problem writing sheet 1 (summary) data for ' trialnumber])
        xlswrite(filename, { allData.Subject,	trialnumber,xltrialnum, 'No valid segments'},  'Sheet1', ['A', num2str(i+1), ':D', num2str(i+1)])
    end
    
    try % Smooth Pursuit Page - Sheet 2
        outputPursuitsegments = {length(pursuitVelocitiesList),1};
        outputPursuitsegments{1} = allData.Subject;
        outputPursuitsegments{2} = trialnumber;
        outputPursuitsegments{3} = xltrialnum;
        
        startList = meas.pursuitMovements_start(validSegments);
        endList = meas.pursuitMovements_end(validSegments);
        
        if isempty(endList)
            error('No smooth pursuit greater than 250ms duration')
            xlswrite(filename, { allData.Subject,	trialnumber, xltrialnum, 'No smooth pursuit > 250 ms duration'},  'Sheet2', ['A', num2str(i+1), ':D', num2str(i+1)])
        end
        
        for j = 1:length(pursuitVelocitiesList)
            outputPursuitsegments{(j-1)*3+4} = startList(j);
            outputPursuitsegments{(j-1)*3+5} = endList(j);
            outputPursuitsegments{(j-1)*3+6} = pursuitVelocitiesList(j);
        end
        
        %disp(outputPursuitsegments)
        rowLen = length(outputPursuitsegments);
        if rowLen > 26
            xlcell = ['A', num2str(i+1),':A',char(64+rowLen-26),num2str(i+1)];
        else
            xlcell = ['A', num2str(i+1),':',char(64+rowLen),num2str(i+1)];
        end
        xlswrite(filename, outputPursuitsegments, 'Sheet2', xlcell)
    catch
        disp(['No valid pursuit data for ' trialnumber])
        xlswrite(filename, { allData.Subject,	trialnumber, xltrialnum, 'No smooth pursuit < 250ms'},  'Sheet2', ['A', num2str(i+1), ':D', num2str(i+1)])
    end
    
    try % Smooth pursuit shorter than 250 ms - Sheet 3
        tooshortpursuits = abs(meas.pursuitMeanVelocities(~validSegments));
        
        outputPursuitsegments = {length(tooshortpursuits),1};
        outputPursuitsegments{1} = allData.Subject;
        outputPursuitsegments{2} = trialnumber;
        outputPursuitsegments{3} = xltrialnum;
        
        startList = meas.pursuitMovements_start(~validSegments);
        endList = meas.pursuitMovements_end(~validSegments);
        
        
        if isempty(endList)
            error('No smooth pursuit shorter than 250ms')
            xlswrite(filename, { allData.Subject,	trialnumber, xltrialnum, 'No smooth pursuit < 250 ms duration'},  'Sheet3', ['A', num2str(i+1), ':D', num2str(i+1)])
            
        end
        
        for j = 1:length(tooshortpursuits)
            outputPursuitsegments{(j-1)*3+4} = startList(j);
            outputPursuitsegments{(j-1)*3+5} = endList(j);
            outputPursuitsegments{(j-1)*3+6} = tooshortpursuits(j);
        end
        
        rowLen = length(outputPursuitsegments);
        if rowLen > 26
            xlcell = ['A', num2str(i+1),':A',char(64+rowLen-26),num2str(i+1)];
        else
            xlcell = ['A', num2str(i+1),':',char(64+rowLen),num2str(i+1)];
        end
        xlswrite(filename, outputPursuitsegments, 'Sheet3', xlcell)
    catch
        disp(['No valid short pursuit data for ' trialnumber])
        xlswrite(filename, { allData.Subject,	trialnumber, xltrialnum, 'No valid segments'},  'Sheet3', ['A', num2str(i+1), ':D', num2str(i+1)])
    end
    
    try % Saccades - Sheet 4
        outputPursuitsegments = {meas.numGshifts+3,1};
        outputPursuitsegments{1} = allData.Subject;
        outputPursuitsegments{2} = trialnumber;
        outputPursuitsegments{3} = xltrialnum;
        
        startList = meas.Gshifts_start;
        endList = meas.Gshifts_end;
        saccadeSpeed = meas.Gamplitudes./meas.Gdurations'*1000;
        
        if isempty(endList)
            error('No saccades')
        end
        
        for j = 1:meas.numGshifts
            outputPursuitsegments{(j-1)*3+4} = startList(j);
            outputPursuitsegments{(j-1)*3+5} = endList(j);
            outputPursuitsegments{(j-1)*3+6} = saccadeSpeed(j);
        end
        
        rowLen = length(outputPursuitsegments);
        if rowLen > 26
            xlcell = ['A', num2str(i+1),':A',char(64+rowLen-26),num2str(i+1)];
        else
            xlcell = ['A', num2str(i+1),':',char(64+rowLen),num2str(i+1)];
        end
        xlswrite(filename, outputPursuitsegments, 'Sheet4', xlcell)
    catch
        disp(['No valid saccade data for ' trialnumber])
        xlswrite(filename, { allData.Subject,	trialnumber,xltrialnum, 'No valid saccades'},  'Sheet4', ['A', num2str(i+1), ':D', num2str(i+1)])
    end
    
end
xlswrite(filename, { 'Export Complete'},  'Sheet1', ['A', num2str(i+2), ':A', num2str(i+2)])
save(scalefilename, 'scale')
disp([ 'Export finished. Check Excel file ' filename ])

end

function [m]= MeasureTrial(mydata,trialname, scaleF, offsetF)
allData=get(mydata.MainFigure,'userdata');
scaleAZ = allData.data.scaleAZ;
offsetAZ = allData.data.offsetAZ;
%horizontal
positions.eyesAZ=(allData.data.(trialname).eyesAZ_R+allData.data.(trialname).eyesAZ_L)/2*scaleF+offsetF;
positions.eyesEL=allData.data.(trialname).eyesEL_R;
positions.speakerAZ=allData.data.(trialname).speaker;
positions.headYAW=allData.data.(trialname).Head_Yaw*scaleAZ+offsetAZ;
positions.gaze = positions.eyesAZ + positions.headYAW;

[velocities, ~]= calcva(positions,20);

try % calculate the trail start and end positions
    m.start_position = round(mean(positions.speakerAZ(200:300)));
    m.end_position = round(mean(positions.speakerAZ(end-300:end)));
    
    % Determine when the speaker starts moving by a linear fit to 80% of
    % the motion
    ramptrail = (positions.speakerAZ'-m.start_position)*sign(m.end_position-m.start_position);
    x1 = find(ramptrail > 0.1*max(ramptrail),1);
    x2 = find(ramptrail > 0.9*max(ramptrail),1);
    X = ramptrail(x1:x2);
    Y = x1:x2;
    start_time =floor(roots(polyfit(X,Y,1))+x1);
    end_time = length(positions.speakerAZ);
    
    %speakerpath = (positions.speakerAZ(start_time-500:end)'-m.start_position)*sign(m.end_position-m.start_position);
    %  assignin('base','speakerpath', speakerpath)
    % Determine when the head starts moving by a linear fit to 80% of
    % the motion
    
    headstartpos = round(mean(positions.headYAW(start_time-500:start_time)));
    headpath = (positions.headYAW(start_time-500:end)'-headstartpos)*sign(m.end_position-m.start_position);
    %assignin('base','headpath', headpath)
    %     x1 = find(headpath > 0.1*max(headpath),1);
    %     x2 = find(headpath > 0.9*max(headpath),1);
    %     X = headpath(x1:x2);
    %     Y = x1:x2;
    %     start_head =floor(roots(polyfit(X,Y,1))+x1);
    
    start_head =find(headpath > 0.1*max(headpath),1);
    
    % Find the start of the head movement by interpolating back from the
    % midpoint of the start of the movement
    gazestartpos = round(mean(positions.gaze(start_time-500:start_time)));
    gazepath = (positions.gaze(start_time-500:end_time)'-gazestartpos)*sign(m.end_position-m.start_position);
    %assignin('base','gazepath', gazepath)
    %     x1 = find(gazepath > 0.25*max(gazepath),1);
    %     x2 = find(gazepath > 0.75*max(gazepath),1);
    %     X = gazepath(x1:x2);
    %     Y = x1:x2;
    %     start_gaze =floor(roots(polyfit(X,Y,1))+x1);
    start_gaze =find(gazepath > 0.1*max(gazepath),1);
    
    
    
    rampVels = abs(velocities.vspeakerAZ(start_time:end_time));
    m.ramp_speed = round(mean(rampVels(rampVels>0.8*max(rampVels)))/5)*5;
    
    if start_time < 501
        startplot = 1;
    else
        startplot=start_time-500; % starts the plot 500ms prior to the arm movement
    end
    
    if end_time+1500 > length(positions.speakerAZ)
        endplot=length(positions.speakerAZ);
    else
        endplot=end_time+1500; % ends the plot a short time after the end of the arm movement
    end
    
    if startplot>endplot
        disp('Start of plot after end of plot')
        startplot = 1;
        endplot = length(positions.speakerAZ);
    end
    
    %horizontal
    positions.eyesAZ=positions.eyesAZ(startplot:endplot); %eye posiitons
    positions.eyesEL=positions.eyesEL(startplot:endplot);
    positions.speakerAZ=positions.speakerAZ(startplot:endplot); %speaker positions
    positions.headYAW=allData.data.(trialname).Head_Yaw(startplot:endplot)*scaleAZ+offsetAZ;
    positions.gaze = positions.eyesAZ + positions.headYAW;
    numSamples = length(positions.eyesAZ);
    
    [velocities, accelerations]= calcva(positions,20);
    
    eyesAZvelon = 40; % saccades if faster then 50
    eyesAZaccon = 1500;
    
    Gshifts_ind= union(find(abs(velocities.gazeAZ) > eyesAZvelon),...
        find(abs(accelerations.gazeAZ) > eyesAZaccon)); %indicies of saccades
    
    % adds a 'saccade' where the start of the speaker movement occurs
    Gshifts_ind = sort([475:525,Gshifts_ind]);
    
    % anything not saccade is pursuit
    pursuitMovements_ind = setxor(1:length(positions.eyesAZ),Gshifts_ind);
    
    %eliminate movements before 100ms
    pursuitMovements_ind=pursuitMovements_ind(pursuitMovements_ind>100);
    
    Gshifts=find(diff(Gshifts_ind)>20);
    numGshifts=length(Gshifts)+1;
    numPursuitMovements = numGshifts+1 ;
    
    if ~isempty(Gshifts_ind)
        Gshifts_start=zeros(1,numGshifts);
        Gshifts_end=zeros(1,numGshifts);
        pursuitMovements_start=zeros(1,numGshifts+1);
        pursuitMovements_end=zeros(1,numGshifts+1);
        
        Gshifts_start(1)=Gshifts_ind(1);
        pursuitMovements_start(1) = 100;
        pursuitMovements_end(1) = Gshifts_start(1)-1;
        
        if numGshifts == 1
            Gshifts_end(1)= Gshifts_ind(end);
            pursuitMovements_start(2) = Gshifts_end(1)+1;
            pursuitMovements_end(2) = numSamples;
        end
        if numGshifts > 1
            Gshifts_end(1)= Gshifts_ind(Gshifts(1));
            pursuitMovements_start(2) = Gshifts_ind(Gshifts(1))+1;
            pursuitMovements_end(2) = Gshifts_ind(Gshifts(1)+1)-1;
            % this is because there Gshifts(1) is the end of the first
            % saccade and Gshifts(2) is the beginning of the second saccade
            
            %the last shift
            Gshifts_start(numGshifts)=Gshifts_ind(Gshifts(numGshifts-1)+1);
            Gshifts_end(numGshifts)=Gshifts_ind(end);
            pursuitMovements_start(numPursuitMovements) = Gshifts_ind(end)+1;
            pursuitMovements_end(numPursuitMovements) = numSamples;
            
            if numGshifts > 2
                for i = 2:numGshifts-1
                    Gshifts_start(i)=Gshifts_ind(Gshifts(i-1)+1);
                    Gshifts_end(i)=Gshifts_ind(Gshifts(i));
                    pursuitMovements_start(i+1) = Gshifts_ind(Gshifts(i))+1;
                    pursuitMovements_end(i+1) = Gshifts_ind(Gshifts(i)+1)-1;
                end
            end
            
        end
        
    else
        m.numGshifts = 0;
        pursuitMovements_start(1) = 100;
        pursuitMovements_end(1) = numSamples;
    end
    
    for i =1:numPursuitMovements
        m.pursuitMeanVelocities(i)=mean(velocities.gazeAZ(pursuitMovements_start(i):pursuitMovements_end(i)));
        x = pursuitMovements_start(i):pursuitMovements_end(i);
        y = positions.gaze(pursuitMovements_start(i):pursuitMovements_end(i))';
        p = polyfit(x,y,1);
        m.pursuitSlope(i) = p(1)*1000;
    end
    
    m.positions=positions;
    
    m.TargetStart = start_time-startplot; % ie always at 500ms
    m.HeadStart = start_head; % 500 is the buffer for plotting before the start of the motion
    m.GazeStart = start_gaze;
    
    % the saccades - measure eyes
    m.numGshifts = numGshifts;
    m.Gshifts_start=Gshifts_start;
    m.Gshifts_end=Gshifts_end;
    m.Gdurations= Gshifts_end-Gshifts_start;
    m.Gamplitudes= positions.eyesAZ(Gshifts_end)-positions.eyesAZ(Gshifts_start);
    
    % Pursuit - measure gaze
    m.pursuitMovements_ind=pursuitMovements_ind;
    m.numPursuitMovements=numPursuitMovements;
    m.pursuitMovements_start=pursuitMovements_start;
    m.pursuitMovements_end=pursuitMovements_end;
    m.pursuitDurations=pursuitMovements_end-pursuitMovements_start;
    m.pursuitAmplitudes= round(positions.gaze(pursuitMovements_end)-positions.gaze(pursuitMovements_start));
    m.headpursuitAmplitude = round(positions.headYAW(pursuitMovements_end)-positions.headYAW(pursuitMovements_start));
    
catch
    disp('Problem extracting saccades and smooth pursuit')
    m.positions=positions;
    m.numGshifts = 0;
end

end

function PlotTrial(~,~,mydata,plotThis,~)
allData=get(mydata.MainFigure,'userdata');
assignin('base','allData',allData);
if nargin<4
    plotThis=allData.currentPlot;
end

f=allData.trialList;
trialname=f{plotThis}{1};

scaleAZ = allData.data.scaleAZ;
offsetAZ = allData.data.offsetAZ;
scaleEL = allData.data.scaleEL;
offsetEL = allData.data.offsetEL;

set(mydata.listbox,'string',[allData.trialList{1:end}],'value',plotThis);

scaleF = allData.data.(trialname).scaleF;
offsetF = allData.data.(trialname).offsetF;
set(mydata.offsetFactdisp,'string',num2str(offsetF));
set(mydata.scaleFactdisp,'string',num2str(scaleF));

meas=MeasureTrial(mydata,trialname, scaleF, offsetF);
assignin('base','MeasuresTrial',meas);

positions=meas.positions;

%Make sure all figs have black background
whitebg('k')
% top plot- full record with EyeAZ in green, Speaker in blue, and zero in red
subplot(3,1,1)

hold off
plot(scaleF*allData.data.(trialname).eyesAZ_R+offsetF,'g')
hold on
plot(scaleF*allData.data.(trialname).eyesAZ_L+offsetF,'g:')
plot(allData.data.(trialname).speaker,'m')
plot(allData.data.(trialname).Head_Yaw*scaleAZ+offsetAZ,'y')

gaze = allData.data.(trialname).Head_Yaw*scaleAZ+offsetAZ + scaleF*(allData.data.(trialname).eyesAZ_R+allData.data.(trialname).eyesAZ_L)/2+offsetF;
plot(gaze,'b')

plot([0,length(allData.data.(trialname).speaker)], [0 0] , 'r')
plot([0,length(allData.data.(trialname).speaker)], [4 4] , 'r:')
plot([0,length(allData.data.(trialname).speaker)], [-4 -4] , 'r:')
title({'Arm motion, Eyes and Head Yaw for full record',allData.filename(1:end-9)})

if allData.Excelinfo(plotThis, 7) == 0
    subplot(3,1,2)
    hold off
    plot([0,length(allData.data.(trialname).speaker)], [0 0] , 'r')
    title('No Speaker Movement on this Trial')
    hold on
    subplot(3,1,3)
    hold off
    plot([0,length(allData.data.(trialname).speaker)], [0 0] , 'r')
    title('No Speaker Movement on this Trial')
    hold on
    return
end
% Plot of head motions
subplot(3,1,2)
hold off
plot(allData.data.(trialname).Head_Yaw*scaleAZ+offsetAZ,'y')
title('Head Motion: Yaw is Yellow; Pitch is Red; Roll is Blue (not calibrated)')
hold on
plot(allData.data.(trialname).Head_Pitch*scaleEL+offsetEL,'r')
plot(allData.data.(trialname).Head_Roll,'b')

% Main movement plot with the smooth pursuit and the saccades
subplot(3,1,3)

hold off
plot(positions.eyesAZ,'g')
title({'Smooth pursuit and saccades extracted',['Excel trial ',num2str(allData.Excelinfo(plotThis,1)) ,' : ',sprintf('Az1 = %d Az2 = %d Vel = %d',allData.Excelinfo(plotThis,5:7))]})
hold on
plot(positions.speakerAZ,'m')
plot(positions.headYAW,'y')
plot(positions.gaze, 'b')

try
    plot([meas.TargetStart,meas.TargetStart],[-10,10], 'LineStyle', ':','LineWidth', 1, 'Color', 'm');
    plot([meas.HeadStart,meas.HeadStart],[-10,10], 'LineStyle', ':','LineWidth', 1, 'Color', 'y');
    plot([meas.GazeStart,meas.GazeStart],[-10,10], 'LineStyle', ':','LineWidth', 1, 'Color', 'b');
catch
    disp('Tried and failed - meas.TargetStart')
end

netSaccade = 0;
shortSeg = 0;
PursuitMotion = [];
GazeMotion = [];
try % mark saccades in red
    i = -1;
    if meas.numGshifts > 0
        for i = 1:meas.numGshifts
            plot(meas.Gshifts_start(i):meas.Gshifts_end(i),positions.eyesAZ(meas.Gshifts_start(i):meas.Gshifts_end(i)),'linewidth',2.5,'color','r')
            
            xes =  meas.pursuitMovements_start(i-shortSeg):meas.pursuitMovements_end(i-shortSeg);
            yes = positions.eyesAZ(meas.pursuitMovements_start(i-shortSeg):meas.pursuitMovements_end(i-shortSeg))-netSaccade;
            yesGAZE = positions.gaze(meas.pursuitMovements_start(i-shortSeg):meas.pursuitMovements_end(i-shortSeg))-netSaccade;
            
            plot(xes,yes,'linewidth',1,'color','w')
            plot(xes,yesGAZE,'linewidth',2.5,'color','c')
            
            segment = [xes',yes];
            gazesegment = [xes',yesGAZE];
            PursuitMotion = [PursuitMotion, segment'];
            GazeMotion = [GazeMotion, gazesegment'];
            
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
                plot(meas.pursuitMovements_start(end):meas.pursuitMovements_end(end),positions.eyesAZ(meas.pursuitMovements_start(end):meas.pursuitMovements_end(end))-netSaccade,'linewidth',1,'color','w')
                plot(meas.pursuitMovements_start(end):meas.pursuitMovements_end(end),positions.eyesAZ(meas.pursuitMovements_start(end):meas.pursuitMovements_end(end))+positions.headYAW(meas.pursuitMovements_start(end):meas.pursuitMovements_end(end))-netSaccade,'linewidth',2.5,'color','c')
            end
        end
    else
        plot(positions.eyesAZ,'linewidth',2.5,'color','b')%% Plot whole trace because no saccades
    end
    
catch
    disp('error Tried and failed to plot smooth pursuit')
    disp(i)
    disp(meas.numGshifts)
    assignin('base','PursuitMotion',PursuitMotion)
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

function [velocities, accelerations]= calcva(positions,n)

if nargin == 1
    n = 20;
end

velocities.veyesAZ=ParabolicDiff(positions.eyesAZ,n);
velocities.vspeakerAZ=ParabolicDiff(positions.speakerAZ,n);
velocities.gazeAZ=ParabolicDiff(positions.gaze,n);
velocities.headYAW=ParabolicDiff(positions.headYAW,2*n);

accelerations.aeyesAZ=ParabolicDiff(velocities.veyesAZ,n);
accelerations.gazeAZ=ParabolicDiff(velocities.gazeAZ,n);
accelerations.headYAW=ParabolicDiff(velocities.headYAW,2*n);
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