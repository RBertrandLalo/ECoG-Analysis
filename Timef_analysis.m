addpath ([cd,'/eeglab/functions/timefreqfunc/'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPOCHING THE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list_stimuluscodes = [1,2 , 3, 4, 7,8,9,10];
sampling_rate = void.SamplingRate.NumericValue ; 


list_baseline      = -0.2*sampling_rate:-1;
list_offset        = 0.05*sampling_rate:0.25*sampling_rate; 
list_signal        = -0.2*sampling_rate:0.6*sampling_rate;  
num_channels       = size(signal,2);

% Take the mean of task and baseline of gamma power for each trial 
for idx_stimcode=1:length(list_stimuluscodes),
	index_onset{idx_stimcode}     = find(diff(states.StimulusCode==list_stimuluscodes(idx_stimcode)) == 1);
    signal_matrix{idx_stimcode} = zeros(length(index_onset{idx_stimcode}),num_channels,length(list_signal)); 

    for idx_signal = 1:length(list_signal),
        signal_matrix{idx_stimcode}(:,:,idx_signal) = signal(index_onset{idx_stimcode}+list_signal(idx_signal),:);
    end
end





channel2plot0 = [3 4 52];
channel2plot1 = [12:14];
channel2plot2 = [20:24];
channel2plot3 = [28:30];
channel2plot4 = [34:38];
channel2plot5 = [70:74];
channel2plot6 = [86:88];

selected_ch =[channel2plot0  channel2plot1 channel2plot2 channel2plot3 channel2plot4  channel2plot5 channel2plot6]; 
for idx_ch=selected_ch
    % sig_ch_PS = [squeeze(signal_matrix{2}(:,idx_ch,:)); squeeze(signal_matrix{4}(:,idx_ch,:))];
    sig_ch = [squeeze(signal_matrix{2}(:,idx_ch,:)); squeeze(signal_matrix{4}(:,idx_ch,:)); squeeze(signal_matrix{6}(:,idx_ch,:)); squeeze(signal_matrix{8}(:,idx_ch,:))];
    sig_ch = [squeeze(signal_matrix{1}(:,idx_ch,:)); squeeze(signal_matrix{3}(:,idx_ch,:)); squeeze(signal_matrix{5}(:,idx_ch,:)); squeeze(signal_matrix{7}(:,idx_ch,:))];
    
    sig_ch_PS = [squeeze(signal_matrix{1}(:,idx_ch,:)); squeeze(signal_matrix{3}(:,idx_ch,:))];
    sig_ch_PD = [squeeze(signal_matrix{2}(:,idx_ch,:)); squeeze(signal_matrix{4}(:,idx_ch,:))];
    sig_ch_US = [squeeze(signal_matrix{5}(:,idx_ch,:)); squeeze(signal_matrix{7}(:,idx_ch,:))];
    sig_ch_UD = [squeeze(signal_matrix{6}(:,idx_ch,:)); squeeze(signal_matrix{8}(:,idx_ch,:))];
    
    sig_ch = sig_ch_PD';
    
    num_frame = size(sig_ch, 1);
    sig_ch_vec = sig_ch(:);
    
    [ersp,itc,powbase,times,freqs,erspboot,itcboot,itcphase] = timef(sig_ch_vec',num_frame,[-200, 600], sampling_rate, 0, 'itcmax',.4,'maxfreq', 170,'plotphase','off', 'erspmax', 2.3  );
    
    saveas(gcf,strcat('ch',num2str(idx_ch), '_timef_PD.jpg'))
    close all
end

%% PLOT ON BRAIN 


% plot for each electrode
ch2plot = 1:92
inset_width = .035;
inset_height = .035;
xl = xlim;
yl = ylim;
zl = zlim;
clear 'ersp'
idx_ch = 35

for idx_ch = 1:92,
    sig_ch_PS = [squeeze(signal_matrix{1}(:,idx_ch,:)); squeeze(signal_matrix{3}(:,idx_ch,:))]; 
   
    sig_ch = sig_ch_PS';
    num_frame = size(sig_ch, 1);
    sig_ch_vec = sig_ch(:);
    [ersp_PS{idx_ch},itc,powbase,times_PS,freqs_PS,erspboot,itcboot,itcphase] = timef(sig_ch_vec',num_frame,[-200, 600], sampling_rate, 0, 'maxfreq', 170,'plotphase','off', 'plotitc', 'off','plotersp','on');
    clear 'sig_ch'    
end

for idx_ch = 1:92,

    sig_ch_US = [squeeze(signal_matrix{5}(:,idx_ch,:)); squeeze(signal_matrix{7}(:,idx_ch,:))];

    sig_ch = sig_ch_US';
    num_frame = size(sig_ch, 1);
    sig_ch_vec = sig_ch(:);
    [ersp_US{idx_ch},itc,powbase,times_US,freqs_US,erspboot,itcboot,itcphase] = timef(sig_ch_vec',num_frame,[-200, 600], sampling_rate, 0, 'maxfreq', 170,'plotphase','off', 'plotitc', 'off','plotersp','on');
    clear 'sig_ch'    
end


for idx_ch = 1:92,

    sig_ch_PD = [squeeze(signal_matrix{2}(:,idx_ch,:)); squeeze(signal_matrix{4}(:,idx_ch,:))];

    sig_ch = sig_ch_PD';
    num_frame = size(sig_ch, 1);
    sig_ch_vec = sig_ch(:);
    [ersp_PD{idx_ch},itc,powbase,times_PD,freqs_PD,erspboot,itcboot,itcphase] = timef(sig_ch_vec',num_frame,[-200, 600], sampling_rate, 0, 'maxfreq', 170,'plotphase','off', 'plotitc', 'off','plotersp','on');
    clear 'sig_ch'    
end

for idx_ch = 1:92,

    sig_ch_UD = [squeeze(signal_matrix{6}(:,idx_ch,:)); squeeze(signal_matrix{8}(:,idx_ch,:))];

    sig_ch = sig_ch_PD';
    num_frame = size(sig_ch, 1);
    sig_ch_vec = sig_ch(:);
    [ersp_UD{idx_ch},itc,powbase,times_UD,freqs_UD,erspboot,itcboot,itcphase] = timef(sig_ch_vec',num_frame,[-200, 600], sampling_rate, 0, 'maxfreq', 170,'plotphase','off', 'plotitc', 'off','plotersp','on');
    clear 'sig_ch'    
end





