function SLab_Movement()
%% Function to extract smooth pursuit movement data from S-Lab data files
%
% Based on previously developed code written by Paul Allen in 2014 and 2015
% This program will stand alone and be able to read all of the S-Lab Movement
% Experiment Data types and provide a uniform analysis of all the data
%
%% Launch the function
mydata.MainFigure = figure('name', 'Data Viewer', 'numbertitle', 'off', 'menubar', 'none');

allData.num_trials=0;
set(mydata.MainFigure,'userdata',allData);

mydata.selecting=uicontrol(mydata.MainFigure,'style','text',...
    'string','Measuring and Selecting Trials...',...
    'position',[210 180 100 50], 'visible','off');

mydata.reading=uicontrol(mydata.MainFigure,'style','text',...
    'string','Reading Data from Labview XML file...',...
    'position',[210 180 100 50], 'visible','off');

mydata.filemenu = uimenu(mydata.MainFigure, 'Label', 'File');
mydata.openh5 = uimenu(mydata.filemenu, 'Label', 'Open S-Lab file', 'callback', {@LoadData mydata});
LoadData([],[],mydata)
end

function ListCallback(cbo,~,mydata)
plotThis=get(cbo,'Value');
allData=get(mydata.MainFigure,'userdata');
allData.currentPlot=plotThis;
set(mydata.MainFigure,'userdata',allData);
PlotTrial([],[],mydata,plotThis);
end

function Rethresh(cbo,~,mydata)
thresh = get(cbo,'String');
allData=get(mydata.MainFigure,'userdata');
plotThis=allData.currentPlot;
f=allData.trialList;
trialname=f{plotThis}{1};
allData.data.(trialname).saccThresh = str2double(thresh);
set(mydata.MainFigure,'userdata',allData);
end

function startWin(cbo,~,mydata)
startWin = get(cbo,'String');
allData=get(mydata.MainFigure,'userdata');
plotThis=allData.currentPlot;
f=allData.trialList;
trialname=f{plotThis}{1};
allData.data.(trialname).startWin = str2double(startWin);
set(mydata.MainFigure,'userdata',allData);
end

function endWin(cbo,~,mydata)
endWin = get(cbo,'String');
allData=get(mydata.MainFigure,'userdata');
plotThis=allData.currentPlot;
f=allData.trialList;
trialname=f{plotThis}{1};
allData.data.(trialname).endWin = str2double(endWin);
set(mydata.MainFigure,'userdata',allData);
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
% This function copies the previous trial scaling parameters into the
% current trial
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
whitebg([0.3,0.3,0.3])
[x,~] = ginput(1);
whitebg('k')
x = round(x);
spkstart = 0;
spkend = mean(allData.data.(trialname).speaker(x-250:x+250));

scaleAZ = allData.data.scaleAZ;
offsetAZ = allData.data.offsetAZ;
eyes = allData.data.(trialname).eyes;

headpos = mean(allData.data.(trialname).Head_Yaw(x-500:x+250))*scaleAZ+offsetAZ;

if eyes == 1
    if allData.data.(trialname).speaker(1) >= 0
        eyesAZ = allData.data.(trialname).eyesAZ_R;
    else
        eyesAZ = allData.data.(trialname).eyesAZ_L;
    end
elseif eyes == 2
    eyesAZ = allData.data.(trialname).eyesAZ_R;
else
    eyesAZ = allData.data.(trialname).eyesAZ_L;
end
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
% This function loads S-Lab XML data for movement experiements and converts
% to a data structure which is saved in mydata.

set(mydata.reading,'visible','on');drawnow
cd '\\smdnas\Paige-Lab\Paige Lab\Labs\A-Lab\Experiments\GSM';
[xmlfile,PathName] = uigetfile('*.hedr.xml','Select the S-Lab output file');
assignin('base', 'xmlfile',xmlfile)
set(mydata.reading,'visible','on');drawnow
disp('Loading file')

datafilename = [PathName, xmlfile(1:end-4),'_allData.mat'];
datafile = dir(datafilename);

if length(datafile) == 1
    load(datafilename) % variable scale is loaded
    disp('Found previously saved session')
    disp(['Loading file ' xmlfile])
    try
        disp(length(allData.data.trial_1.eyesAZ))
    catch
        for i=1:allData.num_trials
            allData.data.(['trial_', num2str(i)]).eyesAZ = (allData.data.(['trial_', num2str(i)]).eyesAZ_R+allData.data.(['trial_', num2str(i)]).eyesAZ_L)/2;
        end
    end
else
    disp('No previous session found')
    disp(['Loading file ' xmlfile])
    scale = zeros(500,2);
    scale(:,1) = 1;
    scale(:,3) = 0;
    scale(:,4) = 0;
    scale(:,5) = 0;
    
    out = xml2struct([PathName, xmlfile]);
    % NB use the the xml2struct function of Falkena, Wanner, Smirnov & Mo
    assignin('base', 'out',out)
    
    numAnalogChans = length(out.GXML_Root.Device_Config.Analog_Names.String);
    numDigitalChans = length(out.GXML_Root.Device_Config.Digital_Names.String);
    numChannels = numAnalogChans+numDigitalChans;
    
    binFile = [PathName, xmlfile(1:end-8), 'smpl'];
    
    fid = fopen(binFile,'r','b');
    analogdigitaldata = fread(fid, [numChannels,inf], 'int16')';
    fclose(fid);
    
    [B,A] = butter(12,0.1,'low');
    
    inneroffset = str2double(out.GXML_Root.Excel_Cluster{1, 1}.Inr_Offset.Text);
    
    totalNumTrials = length(out.GXML_Root.Excel_Cluster);
    Excelinfo = zeros(totalNumTrials, 7);
    trialstartsample = 1;
    for trialNum = 1:totalNumTrials % number of trials
        
        Excelinfo(trialNum,:) = [...
            str2double(out.GXML_Root.Excel_Cluster{1, trialNum}.LV_Trial.Text),... %LV trial number
            str2double(out.GXML_Root.Excel_Cluster{1, trialNum}.Trial_Index.Text),... %LV trial number
            str2double(out.GXML_Root.Excel_Cluster{1, trialNum}.Sample_Count.Text),... % Number of samples
            str2double(out.GXML_Root.Excel_Cluster{1, trialNum}.Event_Count.Text),... % Number of events
            str2double(out.GXML_Root.Excel_Cluster{1, trialNum}.Pre_Inr_P2.Text),... % where the speaker was at the start of the trial for calibration
            str2double(out.GXML_Root.Excel_Cluster{1, trialNum}.Inr_Mov_P.Text),... % Where the inner speaker heads to
            str2double(out.GXML_Root.Excel_Cluster{1, trialNum}.Inr_Mov_V.Text)... % how fast it went
            ];
        
        trialendsample = str2double(out.GXML_Root.Excel_Cluster{1, trialNum}.Sample_Count.Text);
        Inner_Az_index = 3;
        
        if numAnalogChans == 7
            eyesAZ_R_index = 6;% this is spoofed due to EOG
            eyesAZ_L_index = 6;
            eyesEL_R = 5;
            eyesEL_L = 5; % this is spoofed due to EOG
            
            Head_Yaw = zeros(trialendsample-trialstartsample+1,1);
            Head_Pitch = zeros(trialendsample-trialstartsample+1,1);
            Head_Roll = zeros(trialendsample-trialstartsample+1,1);
            
        else
            eyesAZ_R_index = 8;
            eyesAZ_L_index = 10;
            eyesEL_R = 9;
            eyesEL_L = 11;
            
            Head_Yaw = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 14))*0.01;
            Head_Pitch = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 15))*0.01;
            Head_Roll = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, 16))*0.01;
            
        end
        
        speaker = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, Inner_Az_index))*100/32768/1.02286 +inneroffset;
        eyesAZ_R = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, eyesAZ_R_index))*100/32768/3;
        eyesAZ_L = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, eyesAZ_L_index))*100/32768/3;
        eyesEL_R = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, eyesEL_R))*100/32768/3;
        eyesEL_L = filtfilt(B,A,analogdigitaldata(trialstartsample:trialendsample, eyesEL_L))*100/32768/3;
        
        if (max(Head_Yaw) - min(Head_Yaw) > 100)
            Head_Yaw = zeros(length(Head_Yaw),1);
        end
        
        trialstartsample = trialendsample+1;
        
        data.(['trial_', num2str(trialNum)]).speaker = speaker;
        
        data.(['trial_', num2str(trialNum)]).eyesAZ_R = eyesAZ_R;
        data.(['trial_', num2str(trialNum)]).eyesAZ_L = eyesAZ_L;
        data.(['trial_', num2str(trialNum)]).eyesEL_R = eyesEL_R;
        data.(['trial_', num2str(trialNum)]).eyesEL_L = eyesEL_L;
        data.(['trial_', num2str(trialNum)]).eyesAZ = (eyesAZ_R+eyesAZ_L)/2;
        
        data.(['trial_', num2str(trialNum)]).Head_Yaw = Head_Yaw;
        data.(['trial_', num2str(trialNum)]).Head_Pitch = Head_Pitch;
        data.(['trial_', num2str(trialNum)]).Head_Roll = Head_Roll;
        
        data.(['trial_', num2str(trialNum)]).scaleF = scale(trialNum,1);
        data.(['trial_', num2str(trialNum)]).offsetF = scale(trialNum,2);
        data.(['trial_', num2str(trialNum)]).saccThresh = scale(trialNum,3);
        data.(['trial_', num2str(trialNum)]).startWin = scale(trialNum,4);
        data.(['trial_', num2str(trialNum)]).endWin = scale(trialNum,5);
    end
    
    iOne = ['trial_',num2str(find(Excelinfo(:,2) == 1, 1, 'last' ))];
    iTwo = ['trial_',num2str(find(Excelinfo(:,2) == 2, 1, 'last' ))];
    iThree = ['trial_',num2str(find(Excelinfo(:,2) == 3, 1, 'last' ))];
    iFour = ['trial_',num2str(find(Excelinfo(:,2) == 4, 1, 'last' ))];
    iFive = ['trial_',num2str(find(Excelinfo(:,2) == 5, 1, 'last' ))];
    
    if isempty(strfind(xmlfile , 'Free'))
        data.scaleAZ = 0;
        data.offsetAZ = 0;
        data.scaleEL = 0;
        data.offsetEL = 0;
    else
        try
            data.scaleAZ = 40/(mean(data.(iThree).Head_Yaw(end-1000:end))-mean(data.(iOne).Head_Yaw(end-1000:end)));
            data.offsetAZ = -mean(data.(iTwo).Head_Yaw(end-1000:end))*data.scaleAZ;
            data.scaleEL = 20/(mean(data.(iFour).Head_Pitch(end-1000:end))-mean(data.(iFive).Head_Pitch(end-1000:end)));
            data.offsetEL = -mean(data.(iTwo).Head_Pitch(end-1000:end))*data.scaleEL;
        catch
            disp('Improper presentation design - unscaled head movements')
            data.scaleAZ = 1;
            data.offsetAZ = 0;
            data.scaleEL = 1;
            data.offsetEL = 0;
        end
    end
    
    allData=get(mydata.MainFigure,'userdata');
    allData.Excelinfo = Excelinfo;
    allData.data=data;
    allData.Subject = xmlfile(1:6);
    allData.filename = xmlfile;
    allData.PathName = PathName;
    allData.currentPlot=6;
    allData.num_trials=num2cell(1:trialNum/2-3);
    
end

if strfind(xmlfile, 'Cycles') > 0
    allData.cycles = 1;
else
    allData.cycles = 0;
end

set(mydata.reading,'visible','off');drawnow
assignin('base', 'allData',allData)
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

mydata.removethis=uicontrol(mydata.MainFigure,'style','push','string',...
    'Remove Trial','position',[10,150,100,20],...
    'TooltipString','Removes the current trial from ths list',...
    'callback',{@RemoveTrial mydata});

mydata.SaveSession=uicontrol('style','push','string','Save Session',...
    'TooltipString','Save the current session to MAT file',...
    'position',[15 80 100 20],'callback', {@SaveSession mydata});

mydata.Export=uicontrol('style','push','string','Export Data',...
    'TooltipString','Export the smooth and saccade segments',...
    'position',[15 50 100 20],'callback', {@Export mydata});

set(mydata.MainFigure,'userdata',allData);

mydata.startWin= uicontrol(mydata.MainFigure,'Style','edit',...
    'Position', [10,680,50,20], 'Callback', {@startWin mydata});

mydata.startWinlabel = uicontrol(mydata.MainFigure,'Style','text', ...
    'Position', [70,680,80,20], 'String', 'Start of Window');

mydata.endWin= uicontrol(mydata.MainFigure,'Style','edit',...
    'Position', [10,650,50,20], 'Callback', {@endWin mydata});

mydata.endWinlabel = uicontrol(mydata.MainFigure,'Style','text', ...
    'Position', [70,650,80,20], 'String', 'End of Window');

mydata.saccThresh= uicontrol(mydata.MainFigure,'Style','edit',...
    'Position', [10,620,50,20], 'Callback', {@Rethresh mydata});

mydata.saccThreshlabel = uicontrol(mydata.MainFigure,'Style','text', ...
    'Position', [70,620,80,20], 'String', 'Sacc Thesh %');

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

mydata.eyes= uicontrol(mydata.MainFigure,'Style','popup',...
    'Position', [100,730,80,20], 'String', 'Alternating|Right|Left', 'Callback', {@SelectEyes mydata});

mydata.listbox = uicontrol('Style', 'listbox',...
    'Position', [10,200,120,320], 'Callback', {@ListCallback mydata});

set(mydata.listbox,'string',[trialkeyList{1:end}],...
    'Callback', {@ListCallback mydata},'visible','on');

PlotTrial([],[],mydata)

end

function SaveSession(~,~,mydata)
allData=get(mydata.MainFigure,'userdata');
datafilename = [allData.PathName, allData.filename(1:end-4),'_allData.mat'];
save(datafilename, 'allData')
end

function SelectEyes(cbo,~,mydata)
allData=get(mydata.MainFigure,'userdata');
plotThis=allData.currentPlot;
f=allData.trialList;
trialname=f{plotThis}{1};
allData.data.(trialname).eyes = get(cbo,'Value');
allData.data.(trialname).eyes
set(mydata.MainFigure,'userdata',allData)
end


function Export_Cycles(mydata)

allData=get(mydata.MainFigure,'userdata');
f=allData.trialList;

tracefilename = [allData.filename(1:end-4),'_trace.mat'];

for i = 1: length(f)
    %try % mark saccades in red
    trialnumber=f{i}{1};
    disp(['Analyzing trial ' num2str(i)])
    meas=MeasureTrial(mydata,trialnumber,allData.data.(trialnumber).scaleF,allData.data.(trialnumber).offsetF );
    
    outputtraces = double.empty(length(meas.positions.speakerAZ),0);
    outputtraces(:,1) = meas.positions.speakerAZ;
    outputtraces(:,2) = meas.positions.eyesAZ;
    outputtraces(:,3) = NaN(1,length(meas.positions.eyesAZ));
    outputtraces(:,4) = outputtraces(:,3);
    outputtraces(:,5) = meas.positions.headYAW;
    
    netSaccade = 0;
    % desacIndex = 1;
    if meas.numGshifts > 0
        for j = 1:meas.numGshifts
            
            xes =  meas.pursuitMovements_start(j):meas.pursuitMovements_end(j);
            outputtraces(xes,3) = meas.positions.eyesAZ(xes)-netSaccade;
            outputtraces(xes,4) = outputtraces(xes,3);
            
            sacdur = meas.pursuitMovements_start(j+1)-meas.pursuitMovements_end(j);
            lastpursuitpoint = outputtraces(xes(end),4);
            
            if (length(xes) > 100 && (meas.pursuitMovements_start(j+1)+100) < length(meas.positions.eyesAZ) )
                % desaccade
                % get smooth segments before and after the saccade
                prevSmooth = meas.positions.eyesAZ(meas.pursuitMovements_end(j)-100:meas.pursuitMovements_end(j));
                nextSmooth = meas.positions.eyesAZ(meas.pursuitMovements_start(j+1):(meas.pursuitMovements_start(j+1)+100));
                
                % calculate their slopes and take average
                fit1 = polyfit(1:length(prevSmooth), prevSmooth',1);
                fit2 = polyfit(1:length(nextSmooth), nextSmooth',1);
                
                % multiply that slope by duration of saccade to get AZ drift during saccade
                drift = sacdur*(fit1(1)+fit2(1))/2;
                driftsegment = (0:1/sacdur:1)*drift + lastpursuitpoint;
            else
                drift = 0;
                driftsegment = lastpursuitpoint * ones(1, sacdur+1);
            end
            netSaccade = netSaccade + meas.positions.eyesAZ(meas.Gshifts_end(j))- meas.positions.eyesAZ(meas.Gshifts_start(j)) - drift ;
            
            outputtraces(meas.pursuitMovements_end(j):meas.pursuitMovements_start(j+1),4) = driftsegment';
            
            if j == meas.numGshifts
                xes =  meas.pursuitMovements_start(end):meas.pursuitMovements_end(end);
                outputtraces(xes,3) = meas.positions.eyesAZ(xes)-netSaccade;
                outputtraces(xes,4) = outputtraces(xes,3);
            end
            
        end
        
    else
        outputtraces(:,3) = meas.positions.eyesAZ;
        outputtraces(:,4) = outputtraces(:,3);
    end
    cyclestraces.(trialnumber) = outputtraces;
end
save(tracefilename, 'cyclestraces')
disp([ 'Export finished. Check Mat file ' tracefilename ])
assignin('base','cyclestraces',cyclestraces)
end


function Export(~,~,mydata)
%Look here for the one place that you need to change to modify the length
%restriction on valid segments
durationofvalidsegments = 250;

allData=get(mydata.MainFigure,'userdata');

if allData.cycles == 1
    Export_Cycles(mydata)
    return
end

directoryname = allData.PathName;
if isempty(dir(directoryname))
    directoryname = uigetdir;
end
disp('Starting Excel Export')
f=allData.trialList;
filename = [directoryname, allData.filename(1:end-4),'_Slopes.xls'];

% Export to Excel
xlswrite(filename, { 'Subject'	'LVTrial#'  'XLTrial'	'Start angle'	'End angle' 	'Velocity'...
    'Central Slope' 'Central Slope Head YAW' 'Peak Gaze Pursuit'	'Time to Peak' 'Mean Gaze Pursuit'  '# Pursuit segments' 'Fraction all pursuit' ...
    '# Valid pursuit segments'	'Fraction valid pursuit'   'Head Start' 'Gaze Start'}, 'Sheet1', 'A1')
xlswrite(filename, { 'Subject'	'LVTrial#'  'XLTrial'	'Segment Start'  'Segment End'  'Segment Velocity' 	'Valid Pursuit Segments'	}, 'Sheet2', 'A1')
xlswrite(filename, { 'Subject'	'LVTrial#'  'XLTrial'	'Segment Start'  'Segment End'  'Segment Velocity' 	'Invalid Pursuit Segments'	}, 'Sheet3', 'A1')
xlswrite(filename, { 'Subject'	'LVTrial#'  'XLTrial'	'Segment Start'  'Segment End'  'Segment Velocity' 	'Saccade Segments'	}, 'Sheet4', 'A1')

for i = 1:length(f)
    
    xlcell = ['A', num2str(i+1), ':Q', num2str(i+1)];
    trialnumber=f{i}{1};
    xltrialnum = str2double(trialnumber(7:end)); %allData.Excelinfo(i,1);
    if allData.Excelinfo(xltrialnum,7) == 0 % Skip if zero velocity
        disp(['No speaker movement Trial ' num2str(i)])
        xlswrite(filename, { allData.Subject,	trialnumber, xltrialnum, 'No speaker movement this trial'},  'Sheet1', ['A', num2str(i+1), ':D', num2str(i+1)])
        xlswrite(filename, { allData.Subject,	trialnumber, xltrialnum, 'No speaker movement this trial'},  'Sheet2', ['A', num2str(i+1), ':D', num2str(i+1)])
        xlswrite(filename, { allData.Subject,	trialnumber, xltrialnum, 'No speaker movement this trial'},  'Sheet3', ['A', num2str(i+1), ':D', num2str(i+1)])
        xlswrite(filename, { allData.Subject,	trialnumber, xltrialnum, 'No speaker movement this trial'},  'Sheet4', ['A', num2str(i+1), ':D', num2str(i+1)])
        continue
    end
    
    disp(['Analyzing trial ' num2str(xltrialnum)])
    meas=MeasureTrial(mydata,trialnumber,allData.data.(trialnumber).scaleF,allData.data.(trialnumber).offsetF );
    
    startAngle = allData.Excelinfo(xltrialnum,5);
    endAngle = allData.Excelinfo(xltrialnum,6);
    velocity = allData.Excelinfo(xltrialnum,7);
    
    try %Summary Page - Sheet 1
        validSegments = meas.pursuitDurations > durationofvalidsegments; %<--- this is the number to change in order to include shorter saccades
        validSegments(1) = 0; % always exclude the first, flat segment
        pursuitVelocitiesList = abs(meas.pursuitMeanVelocities(validSegments));
        centralslope = abs(allData.data.(trialnumber).centralslope);
        centralslopeYAW = allData.data.(trialnumber).centralslopeYAW;
        [peakpursuit, peakindex] = max(pursuitVelocitiesList);
        SPstartlist = meas.pursuitMovements_start(validSegments);
        SPendlist = meas.pursuitMovements_end(validSegments);
        time2peak = (SPendlist(peakindex)-SPstartlist(peakindex))/2+SPstartlist(peakindex)-500;
        meanpursuit =  mean(pursuitVelocitiesList);
        numSegments = meas.numPursuitMovements;
        allpursuitfraction = sum(meas.pursuitAmplitudes)/(endAngle-startAngle);
        numvalidSegments = sum(validSegments);
        pursuitfraction = sum(meas.pursuitAmplitudes(validSegments))/(endAngle-startAngle);
        
        headstart = meas.HeadStart-meas.TargetStart;
        gazestart = meas.GazeStart-meas.TargetStart;
        
        xlswrite(filename, { allData.Subject,	trialnumber,	xltrialnum, startAngle, ...
            endAngle,	velocity,	centralslope, centralslopeYAW, peakpursuit, time2peak ,	...
            meanpursuit, numSegments ,allpursuitfraction,	numvalidSegments, pursuitfraction, headstart, gazestart}, 'Sheet1', xlcell)
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
disp([ 'Export finished. Check Excel file ' filename ])

end

function [m]= MeasureTrial(mydata,trialname, scaleF, offsetF)
allData=get(mydata.MainFigure,'userdata');
scaleAZ = allData.data.scaleAZ;
offsetAZ = allData.data.offsetAZ;
%horizontal
positions.eyesAZ=allData.data.(trialname).eyesAZ*scaleF+offsetF;
positions.eyesEL=allData.data.(trialname).eyesEL_R;
positions.speakerAZ=allData.data.(trialname).speaker;
positions.headYAW=allData.data.(trialname).Head_Yaw*scaleAZ+offsetAZ;

positions.gaze = positions.eyesAZ + positions.headYAW;
positions.headYAW = positions.headYAW;

velocities= calcva(positions,20);

try % calculate the trail start and end positions
    m.start_position = round(mean(positions.speakerAZ(200:300)));
    m.end_position = round(mean(positions.speakerAZ(end-300:end)));
    %
    % Determine when the speaker starts moving by abs velocity change
    ramptrail = abs(velocities.speakerAZ);
    start_time = find(ramptrail> 0.2*max(ramptrail),1);%floor(roots(polyfit(X,Y,1))+x1);
    end_time = length(positions.speakerAZ);
    
    
    rampVels = abs(velocities.speakerAZ(start_time:end_time));
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
    
    startWin = allData.data.(trialname).startWin;
    if startWin == 0
        set(mydata.startWin,'string',num2str(startplot));
        allData.data.(trialname).startWin = startplot;
        set(mydata.MainFigure,'userdata',allData);
    else
        startplot = startWin;
    end
    
    endWin = allData.data.(trialname).endWin;
    if endWin == 0
        set(mydata.endWin,'string',num2str(endplot));
        allData.data.(trialname).endWin = endplot;
        set(mydata.MainFigure,'userdata',allData);
    else
        endplot = endWin;
    end
    
    %horizontal
    positions.eyesAZ=positions.eyesAZ(startplot:endplot); %eye posiitons
    positions.eyesEL=positions.eyesEL(startplot:endplot);
    positions.speakerAZ=positions.speakerAZ(startplot:endplot); %speaker positions
    positions.headYAW=allData.data.(trialname).Head_Yaw(startplot:endplot)*scaleAZ+offsetAZ;
    positions.gaze = positions.eyesAZ + positions.headYAW;
    numSamples = length(positions.eyesAZ);
    
    velocities = calcva(positions,20);
    
    headtrail = abs(velocities.headYAW);
    start_head = find(headtrail > 0.2*max(headtrail),1);
    
    gazetrail = abs(velocities.gaze);
    start_gaze =find(gazetrail > 0.2*max(gazetrail),1);
    
    smoothvels = calcva(positions,500);
    
    thisThresh = allData.data.(trialname).saccThresh/100;
    if thisThresh == 0
        thisThresh =  3*std(abs(velocities.gaze-smoothvels.gaze))/max(abs(velocities.gaze-smoothvels.gaze));
        set(mydata.saccThresh,'string',num2str(thisThresh*100));
        
        allData.data.(trialname).saccThresh = thisThresh*100;
        set(mydata.MainFigure,'userdata',allData);
    end
    
    eyesAZvelon = thisThresh*max(abs(velocities.gaze-smoothvels.gaze));
    
    Gshifts_ind = find(abs(velocities.gaze-smoothvels.gaze) > eyesAZvelon);
    skirtwidth = 20; % The saccades are surrounded by 20ms of excluded trace
    Gshifts_withSkirt = [];
    for index = Gshifts_ind
        Gshifts_withSkirt = [ Gshifts_withSkirt, [index - skirtwidth: index + skirtwidth]];
    end
    Gshifts_withSkirt(Gshifts_withSkirt < 1) = 500;
    Gshifts_withSkirt(Gshifts_withSkirt > numSamples) = 500;
    
    % adds a 'saccade' where the start of the speaker movement occurs
    Gshifts_ind = unique([475:525,Gshifts_withSkirt]);
    
    % anything not saccade is pursuit
    pursuitMovements_ind = setxor(1:length(positions.eyesAZ),Gshifts_ind);
    
    %eliminate movements before 100ms
    
    Gshifts=find(diff(Gshifts_ind)>20);
    numGshifts=length(Gshifts)+1;
    
    if pursuitMovements_ind(end) == numSamples
        numPursuitMovements = numGshifts+1;
    else
        numPursuitMovements = numGshifts;
    end
    
    if ~isempty(Gshifts_ind)
        Gshifts_start=zeros(1,numGshifts);
        Gshifts_end=zeros(1,numGshifts);
        pursuitMovements_start=zeros(1,numGshifts+1);
        pursuitMovements_end=zeros(1,numGshifts+1);
        
        Gshifts_start(1)=Gshifts_ind(1);
        pursuitMovements_start(1) = 1;
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
        m.pursuitMeanVelocities(i)=mean(velocities.gaze(pursuitMovements_start(i):pursuitMovements_end(i)));
        x = pursuitMovements_start(i):pursuitMovements_end(i);
        y = positions.gaze(pursuitMovements_start(i):pursuitMovements_end(i))';
        p = polyfit(x,y,1);
        m.pursuitSlope(i) = p(1)*1000;
    end
    
    m.positions=positions;
    m.velocities=velocities;
    m.smoothvels=smoothvels;
    
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
set(mydata.listbox,'string',[allData.trialList{1:end}],'value',plotThis);

scaleAZ = allData.data.scaleAZ;
offsetAZ = allData.data.offsetAZ;

scaleF = allData.data.(trialname).scaleF;
offsetF = allData.data.(trialname).offsetF;

meas=MeasureTrial(mydata,trialname, scaleF, offsetF);
assignin('base','meas',meas);
allData=get(mydata.MainFigure,'userdata');

startWin = allData.data.(trialname).startWin;
endWin = allData.data.(trialname).endWin;

if ~isfield(allData.data.(trialname), 'eyes')
    allData.data.(trialname).eyes = 1;
end
eyes = allData.data.(trialname).eyes;

set(mydata.endWin,'string',num2str(endWin))
set(mydata.startWin,'string',num2str(startWin))
set(mydata.offsetFactdisp,'string',num2str(offsetF))
set(mydata.scaleFactdisp,'string',num2str(scaleF))
set(mydata.eyes,'Value',eyes)

if eyes == 1
    if allData.data.(trialname).speaker(1) >= 0
        allData.data.(trialname).eyesAZ = allData.data.(trialname).eyesAZ_R;
    else
        allData.data.(trialname).eyesAZ = allData.data.(trialname).eyesAZ_L;
    end
elseif eyes == 2
    allData.data.(trialname).eyesAZ = allData.data.(trialname).eyesAZ_R;
else
    allData.data.(trialname).eyesAZ = allData.data.(trialname).eyesAZ_L;
end

positions=meas.positions;

%Make sure all figs have black background
whitebg('k')
% top plot- full record with EyeAZ in green, Speaker in blue, and zero in red
subplot(4,1,1)

hold off
plot(scaleF*allData.data.(trialname).eyesAZ_R+offsetF,'b')
hold on
plot(scaleF*allData.data.(trialname).eyesAZ_L+offsetF,'b:')
plot(allData.data.(trialname).speaker,'m')
plot(allData.data.(trialname).Head_Yaw*scaleAZ+offsetAZ,'y')

gaze = allData.data.(trialname).Head_Yaw*scaleAZ+offsetAZ + scaleF*allData.data.(trialname).eyesAZ+offsetF;
plot(gaze,'g')

plot([0,length(allData.data.(trialname).speaker)], [0 0] , 'r')
plot([0,length(allData.data.(trialname).speaker)], [4 4] , 'r:')
plot([0,length(allData.data.(trialname).speaker)], [-4 -4] , 'r:')

plot([startWin, startWin], [-20, 20], 'r--')
plot([endWin, endWin], [-20, 20], 'r--')
title({'Arm motion(M), Eyes(B), Gaze(G) and Head Yaw(Y) for full record',allData.filename(1:end-9)})

% What to do if the speaker doesn't move ???
if (meas.start_position - meas.end_position) == 0
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
if allData.cycles == 0
    % Plot of velocities
    subplot(4,1,2)
    hold off
    plot(meas.velocities.speakerAZ,'m')
    title('Velocity Plot: Arm motion(M), Eyes(B), Gaze(G) and Head Yaw(Y)')
    hold on
    plot(meas.velocities.eyesAZ,'b')
    plot(meas.velocities.gaze,'g')
    plot(meas.velocities.headYAW,'y')
    
    
    % Plot of saccade candidates
    subplot(4,1,3)
    hold off
    thisThresh = allData.data.(trialname).saccThresh*max(abs(meas.velocities.gaze-meas.smoothvels.gaze))/100;
    plot([0, length(meas.velocities.gaze)],[thisThresh,thisThresh],'m-')
    title('Candidate Saccades and threshold')
    hold on
    sacctrace = abs(meas.velocities.gaze-meas.smoothvels.gaze);
    plot(sacctrace,'w')
    sacctrace(sacctrace < thisThresh) = 0;
    plot(sacctrace,'r')
    
    % Main movement plot with the smooth pursuit and the saccades
    
    subplot(4,1,4)
else
    subplot(4,1,2)
    hold off
    thisThresh = allData.data.(trialname).saccThresh*max(abs(meas.velocities.gaze-meas.smoothvels.gaze))/100;
    plot([0, length(meas.velocities.gaze)],[thisThresh,thisThresh],'m-')
    title('Candidate Saccades and threshold')
    hold on
    sacctrace = abs(meas.velocities.gaze-meas.smoothvels.gaze);
    plot(sacctrace,'w')
    sacctrace(sacctrace < thisThresh) = 0;
    plot(sacctrace,'r')
    
    subplot(2,1,2)
end

hold off
plot(positions.eyesAZ,'g')
title({'Smooth pursuit and saccades extracted',['Excel trial ',num2str(allData.Excelinfo(plotThis,1)) ,' : ',sprintf('Az1 = %d Az2 = %d Vel = %d',allData.Excelinfo(plotThis,5:7))]})

hold on

plot(positions.speakerAZ,'m')
plot(positions.headYAW,'y')
plot(positions.gaze, 'b')

try
    plot([meas.TargetStart,meas.TargetStart],[-10,10], 'LineStyle', ':','LineWidth', 1, 'Color', 'm');
catch
    disp('Could not plot start of speaker motion')
end
try
    plot([meas.HeadStart,meas.HeadStart],[-10,10], 'LineStyle', ':','LineWidth', 1, 'Color', 'y');
catch
    disp('Could not plot start of head motion')
end
try
    plot([meas.GazeStart,meas.GazeStart],[-10,10], 'LineStyle', ':','LineWidth', 1, 'Color', 'b');
catch
    disp('Could not plot start of gaze motion')
end

netSaccade = 0;
shortSeg = 0;

try % mark saccades in red
    i = -1;
    desacpursuit = zeros(1,length(meas.pursuitMovements_ind));
    desacIndex = 1;
    if meas.numGshifts > 0
        for i = 1:meas.numGshifts
            plot(meas.Gshifts_start(i):meas.Gshifts_end(i),positions.eyesAZ(meas.Gshifts_start(i):meas.Gshifts_end(i)),'linewidth',2.5,'color','r')
            
            xes =  meas.pursuitMovements_start(i-shortSeg):meas.pursuitMovements_end(i-shortSeg);
            yes = positions.eyesAZ(meas.pursuitMovements_start(i-shortSeg):meas.pursuitMovements_end(i-shortSeg))-netSaccade;
            yesGAZE = positions.gaze(meas.pursuitMovements_start(i-shortSeg):meas.pursuitMovements_end(i-shortSeg))-netSaccade;
            
            plot(xes,yes,'linewidth',1,'color','w')
            plot(xes,yesGAZE,'linewidth',2.5,'color','c')
            
            desacpursuit(desacIndex:desacIndex+length(xes)-1) = yesGAZE;
            desacIndex = desacIndex+length(xes);
            
            if (length(xes) > 100 && (meas.pursuitMovements_start(i+1-shortSeg)+100) < length(positions.eyesAZ) )
                % desaccade
                % get smooth segments before and after the saccade
                prevSmooth = positions.eyesAZ(meas.pursuitMovements_end(i-shortSeg)-100:meas.pursuitMovements_end(i-shortSeg));
                nextSmooth = positions.eyesAZ(meas.pursuitMovements_start(i+1-shortSeg):(meas.pursuitMovements_start(i+1-shortSeg)+100));
                
                % calculate their slopes and take average
                fit1 = polyfit(1:length(prevSmooth), prevSmooth',1);
                fit2 = polyfit(1:length(nextSmooth), nextSmooth',1);
                
                % multiply that slope by duration of saccade to get AZ drift during saccade
                drift = (meas.Gshifts_end(i)-meas.Gshifts_start(i))*(fit1(1)+fit2(1))/2;
                
            else
                drift = 0;
            end
            netSaccade = netSaccade + positions.eyesAZ(meas.Gshifts_end(i))- positions.eyesAZ(meas.Gshifts_start(i)) - drift ;
            if i == meas.numGshifts
                
                xes =  meas.pursuitMovements_start(end):meas.pursuitMovements_end(end);
                yes = positions.eyesAZ(meas.pursuitMovements_start(end):meas.pursuitMovements_end(end))-netSaccade;
                yesGAZE = positions.gaze(meas.pursuitMovements_start(end):meas.pursuitMovements_end(end))-netSaccade;
                
                plot(xes,yes,'linewidth',1,'color','w')
                plot(xes,yesGAZE,'linewidth',2.5,'color','c')
                
                desacpursuit(desacIndex:desacIndex+length(xes)-1) = yesGAZE;
                
            end
        end
    else
        
        yes = positions.eyesAZ;
        xes =  1:length(yes);
        yesGAZE = positions.gaze;
        
        plot(xes,yes,'linewidth',1,'color','w')
        plot(xes,yesGAZE,'linewidth',2.5,'color','c')
        
        desacpursuit = yesGAZE;
    end
    allData.data.(trialname).desacindex = meas.pursuitMovements_ind;
    allData.data.(trialname).desacpursuit = desacpursuit;
    if allData.cycles == 0
        % Calculate the slope of the smooth pursuit during the middle half of the speaker movement
        spkStart = meas.start_position;
        spkEnd = meas.end_position;
        distance = spkEnd-spkStart;
        if distance >0
            ispk25 = find((meas.positions.speakerAZ > spkStart+0.25*distance),1);
            ispk75 = find((meas.positions.speakerAZ > spkStart+0.75*distance),1);
        else
            ispk25 = find((meas.positions.speakerAZ < spkStart+0.25*distance),1);
            ispk75 = find((meas.positions.speakerAZ < spkStart+0.75*distance),1);
        end
        
        hxes = ispk25:ispk75;
        hyes = positions.headYAW(hxes)';
        % try and fine the head movement in this range
        fitparh = polyfit(hxes, hyes, 1);
        allData.data.(trialname).centralslopeYAW = fitparh(1)*1000;
        
        ipursuitStart = find(meas.pursuitMovements_ind > ispk25,1);
        ipursuitEnd = find(meas.pursuitMovements_ind > ispk75,1);
        xes = meas.pursuitMovements_ind(ipursuitStart:ipursuitEnd);
        yes = desacpursuit(ipursuitStart:ipursuitEnd);
        
        plot(xes,yes, 'go')
        fitpar = polyfit(xes, yes, 1);
        text(xes(1),yes(1)-10,{['Central Slope = ' num2str(fitpar(1)*1000, '%0.1f') ' deg/s'];['Head C. Slope = ' num2str(fitparh(1)*1000, '%0.1f') ' deg/s']})
        
        allData.data.(trialname).centralslope = fitpar(1)*1000;
    end
    set(mydata.MainFigure,'userdata',allData);
catch
    disp('Error - Tried and failed to plot smooth pursuit')
    disp(['It broke on iteration ' num2str(i) ' of ' num2str(meas.numGshifts)])
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
    pos=pos';
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

function velocities = calcva(positions,n)

if nargin == 1
    n = 20;
end

velocities.eyesAZ=ParabolicDiff(positions.eyesAZ,n);
velocities.speakerAZ=ParabolicDiff(positions.speakerAZ,n);
velocities.gaze=ParabolicDiff(positions.gaze,n);
velocities.headYAW=ParabolicDiff(positions.headYAW,n);
end

function [ s ] = xml2struct( file )
%Convert xml file into a MATLAB structure
% [ s ] = xml2struct( file )
%
% A file containing:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</Element>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Will produce:
% s.XMLname.Attributes.attrib1 = "Some value";
% s.XMLname.Element.Text = "Some text";
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
% s.XMLname.DifferentElement{1}.Text = "Some more text";
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
% s.XMLname.DifferentElement{2}.Text = "Even more text";
%
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012

if (nargin < 1)
    clc;
    help xml2struct
    return
end

if isa(file, 'org.apache.xerces.dom.DeferredDocumentImpl') || isa(file, 'org.apache.xerces.dom.DeferredElementImpl')
    % input is a java xml object
    xDoc = file;
else
    %check for existance
    if (exist(file,'file') == 0)
        %Perhaps the xml extension was omitted from the file name. Add the
        %extension and try again.
        if (isempty(strfind(file,'.xml')))
            file = [file '.xml'];
        end
        
        if (exist(file,'file') == 0)
            error(['The file ' file ' could not be found']);
        end
    end
    %read the xml file
    xDoc = xmlread(file);
end

%parse xDoc into a MATLAB structure
s = parseChildNodes(xDoc);

end

% ----- Subfunction parseChildNodes -----
function [children,ptext,textflag] = parseChildNodes(theNode)
% Recurse over node children.
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012
children = struct;
ptext = struct; textflag = 'Text';
if hasChildNodes(theNode)
    childNodes = getChildNodes(theNode);
    numChildNodes = getLength(childNodes);
    
    for count = 1:numChildNodes
        theChild = item(childNodes,count-1);
        [text,name,attr,childs,textflag] = getNodeData(theChild);
        
        if (~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section'))
            %XML allows the same elements to be defined multiple times,
            %put each in a different cell
            if (isfield(children,name))
                if (~iscell(children.(name)))
                    %put existsing element into cell format
                    children.(name) = {children.(name)};
                end
                index = length(children.(name))+1;
                %add new element
                children.(name){index} = childs;
                if(~isempty(fieldnames(text)))
                    children.(name){index} = text;
                end
                if(~isempty(attr))
                    children.(name){index}.('Attributes') = attr;
                end
            else
                %add previously unknown (new) element to the structure
                children.(name) = childs;
                if(~isempty(text) && ~isempty(fieldnames(text)))
                    children.(name) = text;
                end
                if(~isempty(attr))
                    children.(name).('Attributes') = attr;
                end
            end
        else
            ptextflag = 'Text';
            if (strcmp(name, '#cdata_dash_section'))
                ptextflag = 'CDATA';
            elseif (strcmp(name, '#comment'))
                ptextflag = 'Comment';
            end
            
            %this is the text in an element (i.e., the parentNode)
            if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                    ptext.(ptextflag) = text.(textflag);
                else
                    %what to do when element data is as follows:
                    %<element>Text <!--Comment--> More text</element>
                    
                    %put the text in different cells:
                    % if (~iscell(ptext)) ptext = {ptext}; end
                    % ptext{length(ptext)+1} = text;
                    
                    %just append the text
                    ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                end
            end
        end
        
    end
end
end

% ----- Subfunction getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
% Create structure of node info.

%make sure name is allowed as structure name
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012
name = toCharArray(getNodeName(theNode))';
name = strrep(name, '-', '_dash_');
name = strrep(name, ':', '_colon_');
name = strrep(name, '.', '_dot_');

attr = parseAttributes(theNode);
if (isempty(fieldnames(attr)))
    attr = [];
end

%parse child nodes
[childs,text,textflag] = parseChildNodes(theNode);

if (isempty(fieldnames(childs)) && isempty(fieldnames(text)))
    %get the data of any childless nodes
    % faster than if any(strcmp(methods(theNode), 'getData'))
    % no need to try-catch (?)
    % faster than text = char(getData(theNode));
    text.(textflag) = toCharArray(getTextContent(theNode))';
end

end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
% Create attributes structure.
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012

attributes = struct;
if hasAttributes(theNode)
    theAttributes = getAttributes(theNode);
    numAttributes = getLength(theAttributes);
    
    for count = 1:numAttributes
        %attrib = item(theAttributes,count-1);
        %attr_name = regexprep(char(getName(attrib)),'[-:.]','_');
        %attributes.(attr_name) = char(getValue(attrib));
        
        %Suggestion of Adrian Wanner
        str = toCharArray(toString(item(theAttributes,count-1)))';
        k = strfind(str,'=');
        attr_name = str(1:(k(1)-1));
        attr_name = strrep(attr_name, '-', '_dash_');
        attr_name = strrep(attr_name, ':', '_colon_');
        attr_name = strrep(attr_name, '.', '_dot_');
        attributes.(attr_name) = str((k(1)+2):(end-1));
    end
end
end
