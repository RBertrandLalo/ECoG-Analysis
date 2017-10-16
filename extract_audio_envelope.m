function envelope = extract_audio_envelope(param,filename)


% clear all
% close all
% 
% param.filter_type = 'IIR'; % IIR | FIR
% 
% % --- speech ---
% param.f_low_stop  = [ 100];
% param.f_low_pass  = [ 200];
% param.f_high_pass = [2000];
% param.f_high_stop = [2200];
% 
% param.samplingrate = 10;
% 
% for idx = 1:1, %5,
%     filename{idx}            = sprintf('/Volumes/BCI_MAC/AMC/SUBJECT026/data/DAY3/AuditoryAttentionTask/ECOG001/tasks/audio_overt/ECOGS001R%02d.wav',idx);
% end

%% enable distributed computing if toolbox is enabled 
if license('test','distrib_computing_toolbox')
    if matlabpool('size') == 0 
        matlabpool 
    end
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

signal = [];

for idx_file=1:length(filename),
    [signal_loop, samplingrate, nbits] = wavread(filename{idx_file});
    
    signal = [signal;signal_loop];
end

clear signal_loop

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
        [bandpass{idx}.n,bandpass{idx}.Wn,bandpass{idx}.beta,bandpass{idx}.ftype] = kaiserord(bandpass{idx}.fcuts,bandpass{idx}.mags,bandpass{idx}.devs,samplingrate);
        bandpass{idx}.n = bandpass{idx}.n + rem(bandpass{idx}.n,2);
        
        % caclulate the filter coefficients in a,b design        
        bandpass{idx}.b = fir1(bandpass{idx}.n,bandpass{idx}.Wn,bandpass{idx}.ftype,kaiser(bandpass{idx}.n+1,bandpass{idx}.beta),'noscale');
        [bandpass{idx}.H,bandpass{idx}.f] = freqz(bandpass{idx}.b,1,1024,samplingrate);

    end
    
elseif strcmp(param.filter_type,'IIR'),

    % IIR bandpass definition using Butterworth filter
    for idx=1:length(param.f_low_stop),

        % define passband, stopband and ripples
        bandpass{idx}.Wp = [param.f_low_pass(idx) param.f_high_pass(idx)]/(samplingrate/2); 
        bandpass{idx}.Ws = [param.f_low_stop(idx) param.f_high_stop(idx)]/(samplingrate/2);
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE LOW-PASS FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define passband, stopband and ripples 
lowpass{1}.Wp = 9/(samplingrate/2); 
lowpass{1}.Ws = 10/(samplingrate/2);
lowpass{1}.Rp = 3; 
lowpass{1}.Rs = 30;

% calculate the minimum filter order
[lowpass{1}.n,lowpass{1}.Wn] = buttord(lowpass{1}.Wp,lowpass{1}.Ws,lowpass{1}.Rp,lowpass{1}.Rs);
lowpass{1}.n = lowpass{1}.n + rem(lowpass{1}.n,2);

% caclulate the filter coefficients in Zero-Pole-Gain design
[lowpass{1}.z,lowpass{1}.p,lowpass{1}.k] = butter(lowpass{1}.n,lowpass{1}.Wn,'low');
[lowpass{1}.sos,lowpass{1}.g]=zp2sos(lowpass{1}.z,lowpass{1}.p,lowpass{1}.k);
lowpass{1}.h=dfilt.df2sos(lowpass{1}.sos,lowpass{1}.g);


%%

%visualize_bandpass_filters(bandpass,samplingrate,param);

%%

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

    fprintf(1,'] done\n');

    %% root-mean-square filter signal
%     fprintf(1, '>> Root-mean-square filtering signal \n');
% 
%     fprintf(1,'[');
%     parfor idx_channel=1:size(signal_loop,2),
%         % calculate band-power using bay applying the root-mean-square to the band-pass filtered signal
%         signal_loop(:,idx_channel) = sqrt(signal_loop(:,idx_channel).^2);
%         fprintf(1,'.');
%     end
%     fprintf(1,'] done\n');

    %%
    fprintf(1, '>> Extracting signal envelope \n');
    parfor idx_channel=1:size(signal_loop,2),
        signal_loop(:,idx_channel) = abs(hilbert(double(signal_loop(:,idx_channel))));
    end
        

    %% low-pass filter signal 
%     fprintf(1, '>> Low-pass filtering signal \n');
% 
%     warning('off', 'signal:filtfilt:ParseSOS');
% 
%     fprintf(1,'[');
%     parfor idx_channel=1:size(signal_loop,2),
%         % extract the envelope of the band-power by low-pass filtering the root-mean square signal 
%         signal_loop(:,idx_channel) = single(filtfilt(lowpass{1}.sos,lowpass{1}.g,double(signal_loop(:,idx_channel)))); %#ok<PFBNS>
%         fprintf(1,'.');
%     end
%     fprintf(1,'] done\n');

    %% down-sample signal
    fprintf(1, '>> Down-sampling signal \n');
       
    % calculate decimation factor for 12 Hz sampling rate
    decimation_factor = samplingrate / param.samplingrate;

    % initialize variable for decimating band-power envelope to 12 Hz
    signal_decimated_loop = zeros(size(decimate(double(signal_loop(:,1)),decimation_factor),1),size(signal_loop,2),'single');

    fprintf(1,'[');
    parfor idx_channel=1:size(signal,2),

        % decimate band-power envelope to 12 Hz
        signal_decimated_loop(:,idx_channel) = single(decimate(double(signal_loop(:,idx_channel)),decimation_factor));
        fprintf(1,'.');
    end
    
    clear signal_loop;
    
    % return variable for each band-pass
    envelope{idx_bandpass} = signal_decimated_loop;
    clear signal_decimated_loop;
    
    fprintf(1,'] done\n');

end    
