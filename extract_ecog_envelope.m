function [signal_filtered,signal_filtered_power_envelope,StimulusCode] = extract_ecog_envelope(param,datafile)

%% user changeable section

% clear variables;
% close all;
% opengl neverselect
% clc

%% enable distributed computing if toolbox is enabled 
%if license('test','distrib_computing_toolbox')
%    if matlabpool('size') == 0 
%        matlabpool 
%    end
%end
    
% % define output directory for plots
% param.outdir = './plots-AMC26';
% 
% % define subject name
% param.SubjectName = 'AMC026';
% 
% % define the input data
% %datafile{1} = '/Volumes/BCI_STORAGE/BCI_SYSTEM/data/SUBJECT033/DAY3/screening/ECOG001/ECOGS001R03.dat';
% 
% for idx = 1:1,%5,
% 
%     datafile{idx}                  = sprintf('/Volumes/ecog/AMC/SUBJECT026/BCI2000/data/DAY3/AuditoryAttentionTask/ECOG001/ECOGS001R%02d.dat',idx);
%     filename.audio_overt{idx}           = sprintf('/Volumes/BCI_MAC/AMC/SUBJECT026/data/DAY3/AuditoryAttentionTask/ECOG001/tasks/audio_overt/ECOGS001R%02d.wav',idx);
% %    filename.audio_overt_original{idx}  = sprintf('/Volumes/ecog/AMC/SUBJECT022/tasks/audio_overt_original/ECOGS001R%02d.wav',idx);
% 
% end
% 
% 
% %param.topo.ElectrodeLocationList = csvread('/Users/pbrunner/Documents/bird-server/ecog/MCH/SWEENEY/coordinates-sweeney.csv');
% 
% %% do not change code from her on
% % spectral parameters
% param.filter_car  = 1;     % 0 == off 
%                            % 1 == on
% 
% param.filter_type = 'IIR'; % IIR | FIR
%                                
% 
% % --- all bands ---
% % param.f_low_stop  = [ 1, 9,19,29,38,58, 78, 98,116,146,176,206];
% % param.f_low_pass  = [ 2,10,20,30,40,60, 80,100,120,150,180,210];
% % param.f_high_pass = [10,20,30,40,60,80,100,120,150,180,210,240];
% % param.f_high_stop = [11,21,31,41,62,82,102,124,154,184,214,244];
% % param.topo_size   = [4,3];
% 
% % --- gamma ---
% % param.f_low_stop  = [ 80];
% % param.f_low_pass  = [ 82];
% % param.f_high_pass = [ 98];
% % param.f_high_stop = [100];
% % param.topo_size   = [1,1];
% 
% % --- broadband gamma ---
% param.f_low_stop  = [ 50];
% param.f_low_pass  = [ 60];
% param.f_high_pass = [250];
% param.f_high_stop = [300];
% param.topo_size   = [1,1];
% 
% % --- highpass filter --- 
% param.highpass.Wp = 0.50; % Hz
% param.highpass.Ws = 0.05; % Hz
% param.highpass.Rp = 3;    % dB
% param.highpass.Rs = 30;   % dB
% 
% % --- highpass filter --- 
% param.lowpass.Wp = 3; % Hz
% param.lowpass.Ws = 4; % Hz
% param.lowpass.Rp = 3;    % dB
% param.lowpass.Rs = 30;   % dB
% 
% 
% 
% define the invalid channels
ch_overlap = [];
ch_noise   = [];
ch_other   = [];

% param.samplingrate = 20;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK MATLAB DEPENDENCIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Checking MATLAB dependencies \n');

% test if the following mex files and MATLAB licenses are available
dependency.load_bcidat           = ~(exist('load_bcidat','file') ~= 3) ; 
dependency.filter_design_toolbox = license('test','filter_design_toolbox');
dependency.signal_toolbox        = license('test','signal_toolbox');
dependency.signal_blocks         = license('test','signal_blocks');
dependency.image_toolbox         = license('test','image_toolbox');

% present an error message if any of the tested mex files and MATLAB licenses are not available
if ~dependency.load_bcidat,
    warning('error: The mex file load_bcidat is missing for this platform!'); %#ok<WNTAG>
end

if ~dependency.filter_design_toolbox, 
    warning('error: This matlab script needs the filter design toolbox installed!'); %#ok<WNTAG>
end

if ~dependency.signal_toolbox, 
    warning('error: This matlab script needs the signal processing toolbox installed!'); %#ok<WNTAG>
end

if ~dependency.signal_blocks, 
    warning('error: This matlab script needs the DSP system toolbox installed!'); %#ok<WNTAG>
end

if ~dependency.image_toolbox, 
    warning('error: This matlab script needs the image processing toolbox installed!'); %#ok<WNTAG>
end

% stop analysis if any of the tested mex files and MATLAB licenses are not available
if ~dependency.load_bcidat   || ~dependency.filter_design_toolbox || ~dependency.signal_toolbox || ...
   ~dependency.signal_blocks || ~dependency.image_toolbox,

%   error('error: Please install the missing toolbox'); 

end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CHECK DATA FILE DEPENDENCIES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Checking param consistentcy \n');

% present an error message if any of the tested parameters is not consistent
if length(param.f_low_stop) ~= length(param.f_low_pass) || length(param.f_high_pass) ~= length(param.f_high_stop) || length(param.f_low_pass) ~= length(param.f_high_pass)
    error('error: f_low_stop, f_low_pass, f_high_pass, f_high_stop must contain the same number of entries'); 
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CHECK DATA FILE DEPENDENCIES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Checking data file dependencies \n');

% go through all data files
for idx=1:length(datafile),
    
    % test if this file exists
    if ~(exist(datafile{idx},'file'))
        error('error: Data file %d/%d (%s) does not exist!',idx,length(datafile),datafile{idx});
    end
    
    % load just the header of this data file
    [ ~, ~, parameters{idx} ] = load_bcidat(datafile{idx},[0 0]); %#ok<*SAGROW>
end

% initialze dependency test variables
dependency.SamplingRate       = 1;
dependency.SourceChList       = 1;

% go through all combinations of data files
for idxA=1:length(parameters),
    for idxB=1:length(parameters),
        
        % compare these parameters between all files
        dependency.SamplingRate        =  dependency.SamplingRate       &      (parameters     {idxA}.SamplingRate.NumericValue   == parameters     {idxB}.SamplingRate.NumericValue);
        dependency.SourceChList        =  dependency.SourceChList       & (mean(parameters     {idxA}.SourceChList.NumericValue   == parameters     {idxB}.SourceChList.NumericValue)   == 1);
       
    end
end

% present an error message if there is a parameter mismatch across data files
if ~dependency.SamplingRate, 
    warning('error: Data files dont have the same SamplingRate!'); %#ok<WNTAG>
end

if ~dependency.SourceChList, 
    warning('error: Data files dont have the same SourceChList!'); %#ok<WNTAG>
end


% stopy analysis if there is a parameter mismatch across data files
if ~dependency.SamplingRate   %|| ~dependency.SourceChList,
   error('error: Data files are too different to be analyzed together!'); 
end

% check if common sampling rate is a fraction of the sampleblock size
decimation_factor = parameters{1}.SamplingRate.NumericValue / param.samplingrate;
if round(parameters{1}.SampleBlockSize.NumericValue / decimation_factor) ~= (parameters{1}.SampleBlockSize.NumericValue / decimation_factor),
 %   error('error: decimation factor (%d) has to be a fraction of the SampleBlockSize (%d)',decimation_factor,parameters{1}.SampleBlockSize.NumericValue);
end


if isfield(param,'TransmitChList'),
    parameters{1}.TransmitChList.NumericValue = param.TransmitChList;
end

param.channels = setdiff(parameters{1}.TransmitChList.NumericValue',unique([ch_overlap,ch_noise,ch_other]));
parameters = parameters{1};


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE HIGH-PASS FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define passband, stopband and attenuation
highpass{1}.Wp = param.highpass.Wp/(parameters.SamplingRate.NumericValue/2); 
highpass{1}.Ws = param.highpass.Ws/(parameters.SamplingRate.NumericValue/2);
highpass{1}.Rp = param.highpass.Rp; 
highpass{1}.Rs = param.highpass.Rs;

% calculate the minimum filter order for butterworth filter
[highpass{1}.n,highpass{1}.Wn] = buttord(highpass{1}.Wp,highpass{1}.Ws,highpass{1}.Rp,highpass{1}.Rs);
highpass{1}.n = highpass{1}.n + rem(highpass{1}.n,2);

% caclulate the filter coefficients in Zero-Pole-Gain design
[highpass{1}.z,highpass{1}.p,highpass{1}.k] = butter(highpass{1}.n,highpass{1}.Wn,'high');
[highpass{1}.sos,highpass{1}.g]=zp2sos(highpass{1}.z,highpass{1}.p,highpass{1}.k);
highpass{1}.h=dfilt.df2sos(highpass{1}.sos,highpass{1}.g);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE BAND-PASS FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(param.filter_type,'FIR'), 
    
    % FIR bandpass definition using Kaiser filter
    for idx=1:length(param.f_low_stop),

        % define passband, stopband and attenuation
        bandpass{idx}.fcuts = [param.f_low_stop(idx) param.f_low_pass(idx) param.f_high_pass(idx) param.f_high_stop(idx)];
        bandpass{idx}.mags  = [0 1 0];
        bandpass{idx}.devs  = [0.01 0.05 0.01];
        
        % calculate the minimum filter order        
        [bandpass{idx}.n,bandpass{idx}.Wn,bandpass{idx}.beta,bandpass{idx}.ftype] = kaiserord(bandpass{idx}.fcuts,bandpass{idx}.mags,bandpass{idx}.devs,parameters.SamplingRate.NumericValue);
        bandpass{idx}.n = bandpass{idx}.n + rem(bandpass{idx}.n,2);
        
        % caclulate the filter coefficients in a,b design        
        bandpass{idx}.b = fir1(bandpass{idx}.n,bandpass{idx}.Wn,bandpass{idx}.ftype,kaiser(bandpass{idx}.n+1,bandpass{idx}.beta),'noscale');
        [bandpass{idx}.H,bandpass{idx}.f] = freqz(bandpass{idx}.b,1,1024,parameters.SamplingRate.NumericValue);

    end
    
elseif strcmp(param.filter_type,'IIR'),

    % IIR bandpass definition using Butterworth filter
    for idx=1:length(param.f_low_stop),

        % define passband, stopband and ripples
        bandpass{idx}.Wp = [param.f_low_pass(idx) param.f_high_pass(idx)]/(parameters.SamplingRate.NumericValue/2); 
        bandpass{idx}.Ws = [param.f_low_stop(idx) param.f_high_stop(idx)]/(parameters.SamplingRate.NumericValue/2);
        bandpass{idx}.Rp = 3; 
        bandpass{idx}.Rs = 30;

        % calculate the minimum filter order
        [bandpass{idx}.n,bandpass{idx}.Wn] = buttord(bandpass{idx}.Wp,bandpass{idx}.Ws,bandpass{idx}.Rp,bandpass{idx}.Rs);
        bandpass{idx}.n = bandpass{idx}.n + rem(bandpass{idx}.n,2);

        % caclulate the filter coefficients in Zero-Pole-Gain design
        [bandpass{idx}.z, bandpass{idx}.p, bandpass{idx}.k] = butter(bandpass{idx}.n,bandpass{idx}.Wn,'bandpass');
        [bandpass{idx}.sos,bandpass{idx}.g]=zp2sos(bandpass{idx}.z,bandpass{idx}.p,bandpass{idx}.k);
        bandpass{idx}.h=dfilt.df2sos(bandpass{idx}.sos,bandpass{idx}.g);

    end    
    
else 
    error('error: invalid filter type in param.filter_type'); 
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE LOW-PASS FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define passband, stopband and ripples 
lowpass{1}.Wp = param.lowpass.Wp/(parameters.SamplingRate.NumericValue/2); 
lowpass{1}.Ws = param.lowpass.Ws/(parameters.SamplingRate.NumericValue/2);
lowpass{1}.Rp = param.lowpass.Rp; 
lowpass{1}.Rs = param.lowpass.Rs;

% calculate the minimum filter order
[lowpass{1}.n,lowpass{1}.Wn] = buttord(lowpass{1}.Wp,lowpass{1}.Ws,lowpass{1}.Rp,lowpass{1}.Rs);
lowpass{1}.n = lowpass{1}.n + rem(lowpass{1}.n,2);

% caclulate the filter coefficients in Zero-Pole-Gain design
[lowpass{1}.z,lowpass{1}.p,lowpass{1}.k] = butter(lowpass{1}.n,lowpass{1}.Wn,'low');
[lowpass{1}.sos,lowpass{1}.g]=zp2sos(lowpass{1}.z,lowpass{1}.p,lowpass{1}.k);
lowpass{1}.h=dfilt.df2sos(lowpass{1}.sos,lowpass{1}.g);


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE LINE-NOISE FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the line noise frequency and bandwidth
peak.fcenter = 60;
peak.bw      = 0.001;

% calculate the IIR-peak filter coefficients in a,b format 
peak.wo = peak.fcenter/(parameters.SamplingRate.NumericValue/2);  
peak.bw = peak.bw;
[peak.b,peak.a] = iirpeak(peak.wo,peak.bw);  

% define the harmonics of line noise frequency
param.filter.notch.fcenter = [60,120,180,240];
param.filter.notch.bw      = ones(1,length(param.filter.notch.fcenter)).*0.001;

% calculate the IIR-peak filter coefficients in a,b format 
for idx = 1:length(param.filter.notch.fcenter),
    notch{idx}.wo = param.filter.notch.fcenter(idx)/(parameters.SamplingRate.NumericValue/2);  
    notch{idx}.bw = param.filter.notch.bw(idx);
    [notch{idx}.b,notch{idx}.a] = iirnotch(notch{idx}.wo,notch{idx}.bw);  
end

%%

visualize_bandpass_filters(bandpass,parameters.SamplingRate.NumericValue,param)


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND CONCATENATE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize variables that will hold signals and states from all data files
signal = single([]);
states.StimulusCode = [];

% go through all data files
for idx_file = 1:length(datafile),

    fprintf(1, '> Processing data file %d/%d \n', idx_file,length(datafile));
    
    %% load data file
    fprintf(1, '>> Loading data file %s \n',datafile{idx_file});

    [ signal_loop, states_loop, void ] = load_bcidat(datafile{idx_file});
    signal_loop = signal_loop(:,param.channels);    
    
    signal_offset = signal_loop(1,:);
    
    %% highpass filter
    fprintf(1, '>> Highpass filtering signal \n');

    fprintf(1,'[');
    parfor idx_channel=1:size(signal_loop,2),
        warning('off', 'signal:filtfilt:ParseSOS');
        signal_loop(:,idx_channel) = single(filtfilt(highpass{1}.sos,highpass{1}.g,double(signal_loop(:,idx_channel)-signal_offset(idx_channel)))); %#ok<PFBNS>
        fprintf(1,'.');    
    end
    fprintf(1,'] done\n');

    % concatenate states and signals from data files
    signal              = [signal;signal_loop]; %#ok<AGROW>
    states.StimulusCode = [states.StimulusCode;states_loop.StimulusCode];
    clear signal_loop
    clear states_loop
    
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEASSURE LINE-NOISE POWER BEFORE SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Meassuring 60 Hz noise power before signal processing \n');

fprintf(1,'[');
parfor idx_channel=1:size(signal,2),
    % calculate average root-mean-square of the line-noise
    signal_noise_before(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
fprintf(1,'] done\n');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND DATA CHANNELS WITH SIGNIFICANT LINE-NOISE POWER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Searching for channels with significant 60 Hz noise \n');

fprintf(1,'[');
if isfield(parameters,'DeviceIDMaster'),
        

        % get the channels that actually had signals 
        list_channels = param.channels;

        % check if any channels are left and 
        if ~isempty(list_channels),
            % calculate the common average reference signal 
            signal_mean = mean(signal(:,list_channels),2);
        else
            % if no channels is left then the common average reference signal is zero
            signal_mean = zeros(size(signal,1),1);
        end
            
        % for each channel on this amp
        parfor idx_ch=list_channels';
            
            % subtract the common average signal from each channel of this amp          
            signal_preliminary = double(signal(:,idx_ch));% - double(signal_mean);
            
            % calculate the residual line-noise root-mean-square power
            signal_noise(idx_ch) = mean(sqrt(filtfilt(peak.b,peak.a,signal_preliminary).^2),1); %#ok<PFBNS>
            fprintf(1,'.');
            
        end
else
    warning('error: common average reference filtering is only supported for g.USBamp data');  %#ok<WNTAG>
end

clear signal_mean

fprintf(1,'] done\n');

% find those channels for which the line-noise power is 1.5 standard deviations higher than the average 
%param.channels_noise = find(signal_noise > (mean(signal_noise(param.channels))+1.5*std(signal_noise(param.channels))));
param.channels_noise = find(signal_noise > 14); %12

param.channels_noise = intersect(param.channels_noise,param.channels);
param.channels_selected = setdiff(param.channels,param.channels_noise);

fprintf(1, '> Found %d channels with significant 60 Hz noise: ',length(param.channels_noise));
fprintf(1, '%d ',param.channels_noise);
fprintf(1, '\n');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMOVE COMMON NOISE USING A COMMON-AVERAGE REFERENCE FILTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf(1, '> Common average filtering signal \n');

% now calculate the common-average reference signal again, without those channels that have significant line-noise
fprintf(1,'[');
if isfield(parameters,'DeviceIDMaster') && param.filter_car > 0,
    

        % exclude the channels that had signifiant line-noise
        list_channels = param.channels_selected;
        
        % check if any channels are left and 
        if ~isempty(list_channels),   
            
            % calculate the common average reference signal 
            signal_mean = mean(signal(:,param.channels_selected),2);
        
            % subtract the common average signal from each channel of this amp
%            for idx_ch=list_channels,
            for idx_ch=1:size(signal,2),
                signal(:,idx_ch) = signal(:,idx_ch) - signal_mean;
                fprintf(1,'.');
            end
        else
            % if no channel is left then the signal is not filtered
            for idx_ch=param.channels_selected,
                signal(:,idx_ch) = signal(:,idx_ch);
                fprintf(1,'.');
            end
        end
      
end 
fprintf(1,'] done\n');
 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEASSURE LINE-NOISE POWER AFTER SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Meassuring 60 Hz noise power after signal processing \n');

fprintf(1,'[');
% for each channels calculate the root-mean-square line-noise power
parfor idx_channel=1:size(signal,2),
    signal_noise_after(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
fprintf(1,'] done\n');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPORT LINE-NOISE POWER AFTER SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Reduced 60 Hz noise from %.2f to %.2f uV',mean(signal_noise_before(param.channels_selected)),mean(signal_noise_after(param.channels_selected)));
fprintf(1, '\n');



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTER OUT REMAINING LINE-NOISE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Notch filtering signal \n');
fprintf(1,'[');
% for each channel
parfor idx_channel=1:size(signal,2),
    
    % get the signal for this channel
    signal_preliminary = double(signal(:,idx_channel));
    
    % remove all harmonics of line-noise
    for idx = 1:length(param.filter.notch.fcenter), %#ok<PFBNS>
        signal_preliminary = filtfilt(notch{idx}.b,notch{idx}.a,signal_preliminary); %#ok<PFBNS>
    end 
    
    % return the signal
    signal(:,idx_channel) = single(signal_preliminary);
    
    fprintf(1,'.');
end
fprintf(1,'] done\n');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT AND DECIMATE BANDPOWER FOR EACH FREQUENCY BAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% for each band-pass filter process the signal 
for idx_bandpass = 1:length(bandpass),
    
fprintf(1, '> Extracting Frequency band %d/%d (%d-%d Hz)\n',idx_bandpass,length(bandpass),param.f_low_pass(idx_bandpass),param.f_high_pass(idx_bandpass));   
    
    %% band-pass filtering signal
    fprintf(1, '>> Band-pass filtering signal %d-%d Hz\n',param.f_low_pass(idx_bandpass),param.f_high_pass(idx_bandpass));

    warning('off', 'signal:filtfilt:ParseSOS');

    fprintf(1,'[');
    parfor idx_channel=1:size(signal,2),

        % check if a,b or zero-pole-gain filter design
        if isfield(bandpass{idx_bandpass},'b'),  %#ok<PFBNS>
            % filter signal using a,b coeffiienct 
            signal_loop(:,idx_channel) = single(filtfilt(bandpass{idx_bandpass}.b,1,double(signal(:,idx_channel))));
        else
            % filter signal using zero-pole-gain coefficients
            signal_loop(:,idx_channel) = single(filtfilt(bandpass{idx_bandpass}.sos,bandpass{idx_bandpass}.g,double(signal(:,idx_channel))));
        end
        fprintf(1,'.');
    end

    %% remove IIR artifact
    if ~isfield(bandpass{idx_bandpass},'b'),
        % set first and last few samples of the band-pass filtered signal to zero
        signal_loop(1:bandpass{idx_bandpass}.n*30      ,:) = 0;
        signal_loop(end-bandpass{idx_bandpass}.n*30:end,:) = 0;
    end
   
    signal_raw = signal_loop;
    
    fprintf(1,'] done\n');

    %% hilbert filter signal
    fprintf(1, '>> Extracting signal envelope \n');
    
    fprintf(1,'['); 
    parfor idx_channel=1:size(signal_loop,2),
        signal_loop(:,idx_channel) = abs(hilbert(double(signal_loop(:,idx_channel))));
        fprintf(1,'.');        
    end
    fprintf(1,'] done\n');
    
        
    %% low-pass filter signal 
    fprintf(1, '>> Low-pass filtering signal \n');

    warning('off', 'signal:filtfilt:ParseSOS');

    fprintf(1,'[');
    parfor idx_channel=1:size(signal_loop,2),
        % extract the envelope of the band-power by low-pass filtering the root-mean square signal 
        signal_loop(:,idx_channel) = single(filtfilt(lowpass{1}.sos,lowpass{1}.g,double(signal_loop(:,idx_channel)))); %#ok<PFBNS>
        fprintf(1,'.');
    end
    fprintf(1,'] done\n');

    %% down-sample signal
    fprintf(1, '>> Down-sampling signal \n');
       
    % calculate decimation factor for 12 Hz sampling rate
    decimation_factor = parameters.SamplingRate.NumericValue / param.samplingrate;

    % decimate StimulusCode to 12 Hz
    StimulusCode = downsample(double(states.StimulusCode),decimation_factor);

    % initialize variable for decimating band-power envelope to 12 Hz
    signal_decimated_loop = zeros(size(decimate(double(signal_loop(:,1)),decimation_factor),1),size(signal_loop,2),'single');
    
    % initialize variable for decimating band-power envelope to 12 Hz
    signal_raw_decimated_loop = zeros(size(decimate(double(signal_loop(:,1)),decimation_factor),1),size(signal_loop,2),'single');
    
    fprintf(1,'[');
    parfor idx_channel=1:size(signal,2),

        % decimate band-power envelope
        signal_decimated_loop(:,idx_channel) = single(decimate(double(signal_loop(:,idx_channel)),decimation_factor));
        
        % downsample filtered signal
        signal_raw_decimated_loop(:,idx_channel) = single(downsample(double(signal_raw(:,idx_channel)),decimation_factor));        
        
        fprintf(1,'.');
    end
    
    clear signal_loop;
    clear signal_raw;
    
    % return variable for each band-pass
    signal_filtered_power_envelope{idx_bandpass} = signal_decimated_loop;
    signal_filtered{idx_bandpass} = signal_raw_decimated_loop;
    clear signal_decimated_loop;
    clear signal_raw_decimated_loop;
    
    fprintf(1,'] done\n');

end 



clear signal; 
   





