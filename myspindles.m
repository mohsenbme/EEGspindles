function myspindles()
% Spindle analysis 
% Programmed by Mohsen Naji; Aug 26 2017, Sleep and Cognition Lab, UC Riverside
% This program assumes the input edf file only contains Stage 2 or Stage 3
% sleep and the data were re-referenced to proper electrodes
% Spindle detection algorithm uses morlet wavelet and adopted from E. Wamsley et al. "Reduced Sleep Spindles and Spindle Coherence in 
% Schizophrenia: Mechanisms of Impaired Memory Consolidation", Biol Psychiatry 71, 2012, pp 154-161
% Frequency estimation is based on the average time differences between
% either of local maxima and minima of the detected spindles

[FileName,PathName] = uigetfile('*.edf','Select the Stage 2 edf file');
[X,Channels,fs]=edf2matsamefs([PathName FileName]);
stage_min=size(X,2)/(fs*60);
spindle_frq=cell(1,length(Channels));
spindle_activity=spindle_frq;
spindle_amplitude=spindle_frq;
spindle_count=zeros(length(Channels),1);
% spindle_density=zeros(1,length(Channels));
mean_frq=zeros(length(Channels),1);
mean_amp=mean_frq;
mean_dur=mean_frq;
mean_act=mean_frq;
for j=1:length(Channels)
    cx{1,1}=X(j,:);
    disp(['spindle analysis for ' Channels{1,j}]);
    [detection,bmn_s,bnd_s] = Mywamsley_spindledet(cx,X(j,:),256);
    spindle_count(j,1)=size(bnd_s,1);
    for i=1:size(bnd_s,1)
        x=X(j,bnd_s(i,1):bnd_s(i,2));
        [~,lcs]=findpeaks(x);
        ffmx=fs./(lcs(2:end)-lcs(1:end-1));
        [~,lcs]=findpeaks(-x);
        ffmn=fs./(lcs(2:end)-lcs(1:end-1));
        spindle_frq{1,j}(i,1)=(mean(ffmx(find(ffmx>11 & ffmx<17)))+mean...
            (ffmn(find(ffmn>11 & ffmn<17))))/2;
        spindle_amplitude{1,j}(i,1)=max(abs(hilbert(detrend(x))));
        spindle_activity{1,j}(i,1)=max(abs(hilbert(detrend(x))))*bmn_s(i);
    end
    sp_frq=spindle_frq{1,j};
    sp_amp=spindle_amplitude{1,j};
    sp_dur=bmn_s;
    sp_act=spindle_activity{1,j};
    ch_output=table(sp_frq,sp_amp,sp_dur,sp_act);
    writetable(ch_output,[PathName FileName(1:end-4) '_' Channels{1,j}(1:2) '.csv'],'Delimiter',',','QuoteStrings',true);
    mean_frq(j,1)=mean(spindle_frq{1,j});
    mean_amp(j,1)=mean(spindle_amplitude{1,j});
    mean_dur(j,1)=mean(bmn_s);
    mean_act(j,1)=mean(spindle_activity{1,j});
end
spindle_density=spindle_count./stage_min;
Channels=Channels';
disp('writing outputs');
output=table(Channels,spindle_count,spindle_density,mean_frq,mean_amp,mean_dur,mean_act);
writetable(output,[PathName FileName(1:end-4) '_averages.csv'],'Delimiter',',','QuoteStrings',true);
end
%%
function [X,Channels,fs]=edf2matsamefs(fileloc)
% [FileName,PathName] = uigetfile('/bazhlab/naji/home/EDFs_ACH_500Hz/*.edf','Select the edf data file');
% load([PathName FileName]);
fid=fopen(fileloc);% in format of [PathName FileName]
a=fread(fid,236,'*char');
ndr=fread(fid,8,'*char');
ndr=str2double(ndr'); %number of data records in sec
a=fread(fid,8,'*char');
drdur=str2double(a'); %duration of each data record in sec
ns=fread(fid,4,'*char'); ns=ns'; ns=str2double(ns);% number of signal channels
Channels=cell(1,ns);
for i=1:ns
    C=fread(fid,16,'*char');C=C';
    Channels{i}=C;
end
fread(fid,ns*80,'*char'); % channel transducer type can be extracted
fread(fid,ns*8,'*char'); %channel physical dimension can be extracted
phmn=zeros(1,ns); phmx=phmn;dmn=phmn;dmx=dmn;
for i=1:ns
    pm=fread(fid,8,'*char');pm=pm';
    phmn(i)=str2double(pm);
end                         %phys min
for i=1:ns
    pm=fread(fid,8,'*char'); pm=pm';
    phmx(i)=str2double(pm);%phys max
end
for i=1:ns
    dm=fread(fid,8,'*char');dm=dm'; 
    dmn(i)=str2double(dm);
end                         %dig min
for i=1:ns
    dx=fread(fid,8,'*char'); dx=dx';
    dmx(i)=str2double(dx);
end                         %dig max
scalefac=(phmx-phmn)./(dmx-dmn);
dc=phmx-scalefac.*dmx;

fread(fid,ns*80,'*char'); % prefilters
nr=zeros(1,ns);
for i=1:ns
    nrc=fread(fid,8,'*char'); nrc=nrc';
    nr(i)=str2double(nrc); %number of samples in each data record
end
if sum(ismember(nr,nr(1)))==length(nr)
    fs=nr(1);
else
    sprintf('Data cant be stored in a single matrix')
end
fread(fid,ns*32,'*char');
ch_fs=nr/drdur;
if mean(nr)==nr(1) && mean(ch_fs)==ch_fs(1)
X=zeros(ns,nr(1)*ndr);
% for i=1:ns
%     X{i,1}=zeros(1,nr(i)*ndr);
% end
fs=ch_fs(1);
end
disp('Reading EDF file...');
for i=1:ndr
    for j=1:ns
        s=fread(fid,nr(j),'int16').*scalefac(j)+dc(j);s=s';
        X(j,(i-1)*nr(j)+1:i*nr(j))=s;
    end


end

fclose(fid);
end
%%%
function [detection,bmn_s,bnd_s] = Mywamsley_spindledet(C3_N2,C3,fs)
% WAMSLEY Detect sleep spindles using the Wamsley algorithm.
% E. Wamsley et al. "Reduced Sleep Spindles and Spindle Coherence in 
% Schizophrenia: Mechanisms of Impaired Memory Consolidation", 
% Biol Psychiatry 71, 2012, pp 154-161
%


signalmean = threshold_wamsley(C3_N2,fs);
[detection,bmn_s,bnd_s] = wamsley(C3,fs,signalmean);

    function signalmean = threshold_wamsley(C3nrem2,fs)
        % THRESHOLD_WAMSLEY Calculates the amplitude criteria for spindle detection
        % using the wamsley method.
        % Input is a structure where each entry contains a continuous segment of
        % EEG data recorded at C3-M2 during S2. The sampling frequency is the
        % final input.
        
        %% Define parameters for the wavelet analysis
        fb = 13.5;
        fc = 0.5;
        scale = fs*fc/fb; %9.48 for 256Hz; %3.7 for 100 Hz
        
        Ltotal = 0;
        for k = 1:length(C3nrem2)
            signal = C3nrem2{k}; L = length(signal);
            %% Perform wavelet transformation
            EEGWave = cwt(signal,scale,['cmor' num2str(fb) '-' num2str(fc)]);
            EEGData = real(EEGWave.^2);
            
            %% Take Moving Average
            EEGData = EEGData.^2;
            window = ones(ceil(fs/10),1)/ceil(fs/10); % create 100ms window to convolve with
            EEGData2 = filter(window,1,EEGData); % take the moving average using the above window
            MA(Ltotal+1:Ltotal+L) = EEGData2;
            Ltotal = Ltotal+L;
        end
        
        %% Determine amplitude threshold
        signalmean = mean(MA);
    end

    function [detection,bmn_s,bnd_s] = wamsley(C3,fs,signalmean)
        % WAMSLEY Detect sleep spindles in EEG given the amplitude criteria.
        % Input is the EEG signal we wish to detect spindles in, the sampling
        % frequency and the amplitude criteria.
        % Output is a vector containing the detection of spindles.
        
        %% Define parameters for the wavelet analysis
        fb = 13.5;
        fc = 0.5;
        scale = fs*fc/fb; %9.48 for 256Hz; %3.7 for 100 Hz
        
        EEGWave = cwt(C3,scale,['cmor' num2str(fb) '-' num2str(fc)]);
        EEGData = real(EEGWave.^2);
        
        %% Take Moving Average
        EEGData = EEGData.^2;
        window = ones(ceil(fs/10),1)/ceil(fs/10); % create 100ms window to convolve with
        EEGData2 = filter(window,1,EEGData); % take the moving average using the above window
        
        %% Determine amplitude threshold
        threshold = signalmean.*4.5; % defines the threshold
        
        %% Find Peaks in the MS Signal
        current_data=EEGData2;
        
        over=current_data>threshold; % Mark all points over threshold as '1'
        detection = zeros(length(current_data),1);
        detection(over) = 1;
        [begins,ends] = find_spindles(detection);
        [detection,begins,ends] = maximum_duration(detection,begins,ends,3,fs);
        [detection,begins_03,ends_03] = minimum_duration(detection,begins,ends,0.4,fs);
        
        locs_03=(zeros(1,length(current_data)))';  % Create a vector of zeros the length of the MS signal
        for i=1:((length(current_data))-ceil(fs*0.4));  % for the length of the signal, if the sum of 30 concurrent points = Fs*0.3, mark a spindle
            if sum(over(i:(i+(ceil(fs*0.4)-1))))==ceil(fs*0.4);
                locs_03(i,1)=1;
            end
        end
        
        spin_03=zeros((length(locs_03)),1);  % only mark a spindle in vector 'spin' at the end of a 300ms duration peak
        for i=1:length(locs_03);
            if locs_03(i,1)==1 && locs_03(i+1,1)==0;
                spin_03(i,1)=1;
            end
        end
        
        for i=513:length(spin_03);%201-->513  % for every spindle marked in 'spin', delete the spindle if there is also a spindle within the second preceeding it
            if spin_03(i,1)==1 && sum(spin_03((i-fs):(i-1)))>0;
                spin_03(i,1)=0;
                idx = find(i>=begins_03 & i<=ends_03);
                if isempty(idx) == 0
                    detection(begins_03(idx):ends_03(idx)) = 0;
                else
                    error('Did not find spindle beginning and ending around a spin point')
                end
            end
        end
        t=find((detection(2:end)-detection(1:end-1))~=0);
        t=[0 t' length(detection)];
        clear bnd_s
        for i=1:length(t)-1
          bnd_s(i,:)=[t(i)+1 t(i+1)];
        end
        bmn_s=(bnd_s(:,2)-bnd_s(:,1)+1)./fs;
        bmn_s(find(detection(bnd_s(:,1))==0),:)=[];
        bnd_s(find(detection(bnd_s(:,1))==0),:)=[];
    end

    function [begins, ends] = find_spindles(bv)
        % FIND_SPINDLES - find start and end index' of spindles.
        % Input is a binary vector bv containing ones where spindles are detected.
        % Output is vectors containing the index' of spindle beginnings and ends
        % (first sample of spindle and last sample of spindle, respectively).
        
        sise = size(bv);
        E = bv(2:end)-bv(1:end-1); % Find start and end of intervals with spindles
        
        begins = find(E==1)+1;
        if bv(1) == 1
            if sise(1) > 1
                begins = [1; begins];
            elseif sise(2) > 1
                begins = [1 begins];
            else
                error('The input signal is not one dimensional')
            end
        elseif numel(begins) == 0 && bv(1) == 0
            begins = NaN;
        end
        
        ends = find(E==-1);
        if bv(end) == 1
            if sise(1) > 1
                ends = [ends; length(bv)];
            elseif sise(2) > 1
                ends = [ends length(bv)];
            else
                error('The input signal is not one dimensional')
            end
        elseif numel(ends) == 0 && bv(end) == 0
            ends = NaN;
        end
    end

    function [bv,begins,ends] = maximum_duration(bv,begins,ends,max_dur,fs)
        % MAXIMUM_DURATION - checks the sample duration of the spindles.
        % Input is a vector containing ones in the interval where the spindle is
        % and indexs describing the start and end of the spindle. The last two
        % inputs are the maximum duration given in seconds and the sampling
        % frequency given in Hz.
        % Output is a vector containing ones in the interval where the spindle with
        % duration shorter than or equal to the maximum duration is and indexs
        % describing the start and end of the spindle.
        
        duration_samples = ends-begins+1;
        for k = 1:length(begins)
            if duration_samples(k) > max_dur*fs
                bv(begins(k):ends(k)) = 0;
                begins(k) = 0;
                ends(k) = 0;
            end
        end
        begins = begins(begins~=0);
        ends = ends(ends~=0);
    end

    function [bv,begins,ends] = minimum_duration(bv,begins,ends,min_dur,fs)
        % MINIMUM_DURATION - checks the sample duration of the spindles.
        % Input is a vector containing ones in the interval where the spindle is
        % and indexs describing the start and end of the spindle. The last two
        % inputs are the minimum duration given in seconds and the sampling
        % frequency given in Hz.
        % Output is a vector containing ones in the interval where the spindle with
        % duration longer than or equal to the minimum duration is and indexs
        % describing the start and end of the spindle.
        
        duration_samples = ends-begins+1;
        for k = 1:length(begins)
            if duration_samples(k) < ceil(min_dur*fs)
                bv(begins(k):ends(k)) = 0;
                begins(k) = 0;
                ends(k) = 0;
            end
        end
        begins = begins(begins~=0);
        ends = ends(ends~=0);
    end

end
        
        