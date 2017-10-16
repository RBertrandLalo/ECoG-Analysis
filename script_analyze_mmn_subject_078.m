clear variables;
close all;
opengl neverselect
clc

addpath ([cd,'/eeglab/functions/sigprocfunc'])
addpath ([cd,'/eeglab/functions/adminfunc'])
addpath ([cd,'/eeglab/functions/guifunc'])
addpath ([cd,'/export_fig/']);


% spectral parameters
ecog.param.filter_car  = 1;     % 0 == off 
                           % 1 == on

ecog.param.filter_type = 'IIR'; % IIR | FIR

% --- all bands ---
% ecog.param.f_low_stop  = [0.1, 4,15,30, 64,220];
% ecog.param.f_low_pass  = [  1, 8,18,35, 70,250];
% ecog.param.f_high_pass = [ 30,12,25,50,170,350];
% ecog.param.f_high_stop = [ 50,16,30,55,200,380];
% ecog.param.topo_size   = [3,3];

ecog.param.f_low_stop  = [ 64];
ecog.param.f_low_pass  = [ 70];
ecog.param.f_high_pass = [ 170];
ecog.param.f_high_stop = [ 200];
ecog.param.topo_size   = [3,3];


% --- broadband ---
% ecog.param.f_low_stop  = [ 0.1 ];
% ecog.param.f_low_pass  = [ 1 ];
% ecog.param.f_high_pass = [30 ];
% ecog.param.f_high_stop = [50 ];
% ecog.param.topo_size   = [1,1];



% --- gerv's bands ---
% ecog.param.f_low_stop  = [ 4,15,30, 64,220];
% ecog.param.f_low_pass  = [ 8,18,35, 70,250];
% ecog.param.f_high_pass = [12,25,50,170,350];
% ecog.param.f_high_stop = [16,30,55,200,380];
% ecog.param.topo_size   = [2,3];


% --- highpass filter --- 
ecog.param.highpass.Wp = 0.50; % Hz
ecog.param.highpass.Ws = 0.05; % Hz
ecog.param.highpass.Rp = 3;    % dB
ecog.param.highpass.Rs = 30;   % dB

% --- lowpass filter --- 
ecog.param.lowpass.Wp = 30;    % Hz
ecog.param.lowpass.Ws = 40;    % Hz
ecog.param.lowpass.Rp = 3;     % dB
ecog.param.lowpass.Rs = 30;    % dB

ecog.param.samplingrate = 120;

ecog.param.TransmitChList = 1:92;
% ecog.param.TransmitChList = 1:64; ecog.param.bad_ch =   []; 
ecog.param.bad_ch = [15,27,43,77,93:134];
% ecog.param.bad_ch = [15,43,80,93:134];
ecog.param.remove_ch = intersect(ecog.param.TransmitChList,ecog.param.bad_ch);
ecog.param.good_ch = setdiff(ecog.param.TransmitChList,intersect(ecog.param.TransmitChList,ecog.param.bad_ch));


%ecog.param.good_ch = ecog.param.TransmitChList

sz_label{1} = 'match';
sz_label{2} = 'mismatch';

%%
datafile{1}  = '/Users/rzb06/Codes_MMN_ECoG/data/Su78/raw/ECOGSu78R01.dat';
datafile{2}  = '/Users/rzb06/Codes_MMN_ECoG/data/Su78/raw/ECOGSu78R02.dat';
datafile{3}  = '/Users/rzb06/Codes_MMN_ECoG/data/Su78/raw/ECOGSu78R03.dat';
datafile{4}  = '/Users/rzb06/Codes_MMN_ECoG/data/Su78/raw/ECOGSu78R04.dat';

idx = 1;

[signal{idx},envelope_ecog{idx},StimulusCode{idx}] = extract_ecog_envelope(ecog.param,datafile);

StimulusCode_corr{idx} = circshift(StimulusCode{idx} + circshift(StimulusCode{idx},[12,0]),[6,0]);

%%


i = 1;

condition{i}.txt{1}  = 'match';
condition{i}.txt{2}  = 'mismatch';
condition{i}.code{1} = 'StimulusCode_corr{idx} == 1';
condition{i}.code{2} = 'StimulusCode_corr{idx} == 2';
i = i + 1;


for i=1:length(condition),  

   corr_var{1} = zeros(size(StimulusCode{idx}));
   sz_eval_1 = ['corr_var{1}(',condition{i}.code{1},') =  1;'];
   sz_eval_2 = ['corr_var{1}(',condition{i}.code{2},') = -1;'];
   eval(sz_eval_1);
   eval(sz_eval_2);
   index = find(corr_var{1} ~= 0); 
   correlation{i}      = correlate_envelopes(envelope_ecog{idx},corr_var,index);
end


%% evoked-potential analysis

% list_stimuluscodes = [1,2,3,4,7,8,9,10];
list_stimuluscodes = [1,2,3,4,5,6,11,12];
list_offset        = -30:60; 
list_baseline      = -24:-13;
num_channels       = size(signal{1}{1},2);


for idx_band = 1:length(signal{1}),

    for idx_stimcode=1:length(list_stimuluscodes),

        index_onset{idx_stimcode}     = find(diff(StimulusCode{1}==list_stimuluscodes(idx_stimcode)) == 1);

        response_matrix{idx_band}{idx_stimcode} = zeros(length(index_onset{idx_stimcode}),num_channels,length(list_offset));       
        

        for idx_offset = 1:length(list_offset),
            response_matrix{idx_band}{idx_stimcode}(:,:,idx_offset) = signal{1}{idx_band}(index_onset{idx_stimcode}+list_offset(idx_offset),:);
        end
        
        
        % reject top and bottom 5th percentile 
        for idx_ch = 1:num_channels,
            for idx_offset = 1:length(list_offset),
                valid_range = prctile(response_matrix{idx_band}{idx_stimcode}(:,idx_ch,idx_offset),[5 95]);
                index_valid = find(response_matrix{idx_band}{idx_stimcode}(:,idx_ch,idx_offset) >= valid_range(1) & response_matrix{idx_band}{idx_stimcode}(:,idx_ch,idx_offset) <= valid_range(2));
                response_matrix_avg{idx_band}{idx_stimcode}(idx_ch,idx_offset) = mean(response_matrix{idx_band}{idx_stimcode}(index_valid,idx_ch,idx_offset));
            end
        end        

        
%        response_matrix_avg{idx_band}{idx_stimcode} = squeeze(mean(response_matrix{idx_band}{idx_stimcode},1));
        response_matrix_avg{idx_band}{idx_stimcode}(ecog.param.remove_ch,:) = NaN;
        
    end

end

%% evoked-power analysis
% list_stimuluscodes = [1,2];
list_stimuluscodes = [1,2,3,4,5,6,11,12];
list_offset        = -30:60; 
list_baseline      = -24:-13;
num_channels       = size(signal{1}{1},2);


for idx_band = 1:length(signal{1})

    for idx=1:length(list_stimuluscodes)
        
        index_onset{idx}     = find(diff(StimulusCode{1}==list_stimuluscodes(idx)) == 1);

        power_matrix{idx_band}{idx}    = zeros(length(index_onset{idx}),num_channels,length(list_offset));      
        baseline_matrix{idx_band}{idx} = zeros(length(index_onset{idx}),num_channels,length(list_baseline));       
        
        
    end  



    for idx_stimcode = 1:length(list_stimuluscodes)


        for idx_offset = 1:length(list_offset)
            power_matrix{idx_band}{idx_stimcode}(:,:,idx_offset) = envelope_ecog{1}{idx_band}(index_onset{idx_stimcode}+list_offset(idx_offset),:);
        end
        
        for idx_baseline = 1:length(list_baseline)
            baseline_matrix{idx_band}{idx_stimcode}(:,:,idx_baseline) = envelope_ecog{1}{idx_band}(index_onset{idx_stimcode}+list_baseline(idx_baseline),:);
        end

        baseline_mean{idx_band}{idx_stimcode} = mean(baseline_matrix{idx_band}{idx_stimcode},3); 
        baseline_std{idx_band}{idx_stimcode}  = zeros(size(baseline_mean{idx_band}));
        baseline_median{idx_band}{idx_stimcode} = median(baseline_matrix{idx_band}{idx_stimcode},3); 
        task_mean{idx_band}{idx_stimcode} = mean(power_matrix{idx_band}{idx_stimcode},3); 
        task_median{idx_band}{idx_stimcode} = median(power_matrix{idx_band}{idx_stimcode},3); 
        
        
        for idx_onset = 1:length(index_onset{idx_stimcode})
            
            for idx_ch = 1:num_channels
               
                baseline_std{idx_band}{idx_stimcode}(idx_onset,idx_ch) = std(baseline_matrix{idx_band}{idx_stimcode}(idx_onset,idx_ch,:));
                
                power_matrix{idx_band}{idx_stimcode}(idx_onset,idx_ch,:) = (squeeze(power_matrix{idx_band}{idx_stimcode}(idx_onset,idx_ch,:)) - baseline_mean{idx_band}{idx_stimcode}(idx_onset,idx_ch)) ./ baseline_std{idx_band}{idx_stimcode}(idx_onset,idx_ch);
                
            end
            
        end
        
        % reject top and bottom 5th percentile 
        for idx_ch = 1:num_channels,
            for idx_offset = 1:length(list_offset),
                valid_range = prctile(power_matrix{idx_band}{idx_stimcode}(:,idx_ch,idx_offset),[5 95]);
                index_valid = find(power_matrix{idx_band}{idx_stimcode}(:,idx_ch,idx_offset) >= valid_range(1) & power_matrix{idx_band}{idx_stimcode}(:,idx_ch,idx_offset) <= valid_range(2));
                power_matrix_avg{idx_band}{idx_stimcode}(idx_ch,idx_offset) = mean(power_matrix{idx_band}{idx_stimcode}(index_valid,idx_ch,idx_offset));
            end
        end
        
%        power_matrix_avg{idx_band}{idx_stimcode} = squeeze(mean(power_matrix{idx_band}{idx_stimcode},1));
        power_matrix_avg{idx_band}{idx_stimcode}(ecog.param.remove_ch,:) = NaN;
        
    end

end

%%

close all

%% load template

[ ~, ~, parameters ] = load_bcidat('/Users/rzb06/Google Drive/peter-raphaelle/data/PXB_session_2/NantesPredictiveTask/ECOG001/ECOGS001R04.dat',[0 0]);


%% render evoked power
topos.dimmax            = 512;
topos.outline           =   2;
topos.min_radius        =   1;
topos.max_radius        =  20;
topos.min_value         =   0;

topos.color_map         = jet(256);
topos.color_bar_ticks   = 4;

topos.signed            = 1;

topos.num_topos         = 1;
topos.num_topos_rows    = 1;
topos.num_topos_cols    = 1;
topos.size_width        = 660+1;
topos.size_height       = 700+1;

[topos.eloc void.labels void.Th void.Rd void.indices] = readlocs( 'eloc64.loc', 'filetype', 'loc');
% topos.ElectrodeLocation = parameters.ElectrodeLocation; 

power_matrix_avg_avg{idx_band}{1} = (power_matrix_avg{idx_band}{1}+power_matrix_avg{idx_band}{3} + power_matrix_avg{idx_band}{5} + power_matrix_avg{idx_band}{7})/4; 
power_matrix_avg_avg{idx_band}{2} = (power_matrix_avg{idx_band}{2}+power_matrix_avg{idx_band}{4} + power_matrix_avg{idx_band}{6} + power_matrix_avg{idx_band}{8})/4; 

for idx_band = 1:size(power_matrix_avg_avg,2),
    for idx_stimcode = 1:size(power_matrix_avg_avg{idx_band},2),

        topos.max_value = max(max(abs(power_matrix_avg_avg{idx_band}{idx_stimcode})));
        topos.min_value = 0;        
        
%         sz_filename = sprintf('./results/result_evoked_power_%d_%d_Hz_%s.mp4',ecog.param.f_low_pass(idx_band),ecog.param.f_high_pass(idx_band),sz_label{idx_stimcode});
        sz_filename = sprintf('./results_EEG/result_evoked_power_%d_%d_Hz_%s.mp4',ecog.param.f_low_pass(idx_band),ecog.param.f_high_pass(idx_band),sz_label{idx_stimcode});

        obj_video = VideoWriter(sz_filename);
        obj_video.FrameRate = 10; 
        
        open(obj_video);
        
        fprintf(1, '>> Rendering %s \n',sz_filename);   
        fprintf(1,'[');         
        
        for idx_offset = 1:size(power_matrix_avg_avg{idx_band}{idx_stimcode},2)

            fprintf(1,'.');                   
            
            topos.values = power_matrix_avg_avg{idx_band}{idx_stimcode}(:,idx_offset);
    
            results = framework_render_electrodes_eeg(topos);

            figure(1),
            imshow(results.image)       
            offset_time = list_offset(idx_offset) * 1/ecog.param.samplingrate * 1000;
            sz_title_top = sprintf('evoked power %d-%d Hz, %s',ecog.param.f_low_pass(idx_band),ecog.param.f_high_pass(idx_band),sz_label{idx_stimcode});
            sz_title_bottom = sprintf('%3.0fms',offset_time);
            annotation1 = annotation(gcf,'textbox','Position',[0 0.95 1 0],'LineStyle','none','FitHeightToText','off','FontSize',25,'FontWeight','bold','HorizontalAlignment','center','String',sz_title_top);   
            annotation2 = annotation(gcf,'textbox','Position',[0 0.17 1 0],'LineStyle','none','FitHeightToText','off','FontSize',25,'FontWeight','bold','HorizontalAlignment','center','String',sz_title_bottom);   
            pause(0.1)
            
            frame = getframe;
            writeVideo(obj_video,frame);
        
        end
        
        close(obj_video);
        
        fprintf(1,'] done\n');                    
        
    end
end

%% render evoked potential
topos.dimmax            = 512;
topos.outline           =   2;
topos.min_radius        =   1;
topos.max_radius        =  20;
topos.min_value         =   0;

topos.color_map         = jet(256);
topos.color_bar_ticks   = 4;

topos.signed            = 1;

topos.num_topos         = 1;
topos.num_topos_rows    = 1;
topos.num_topos_cols    = 1;
topos.size_width        = 660+1;
topos.size_height       = 700+1;

%[topos.eloc void.labels void.Th void.Rd void.indices] = readlocs( 'eloc64.loc', 'filetype', 'loc');
topos.ElectrodeLocation = parameters.ElectrodeLocation; 

for idx_band = 1:size(response_matrix_avg,2),
    for idx_stimcode = 1:size(response_matrix_avg{idx_band},2),
     
        topos.max_value = max(max(abs(response_matrix_avg{idx_band}{idx_stimcode})));
        topos.min_value = 0;        
        
        sz_filename = sprintf('./results/result_evoked_response_%d_%d_Hz_%s.mp4',ecog.param.f_low_pass(idx_band),ecog.param.f_high_pass(idx_band),sz_label{idx_stimcode});
        
        obj_video = VideoWriter(sz_filename);
        obj_video.FrameRate = 10; 
        
        open(obj_video);

        fprintf(1, '>> Rendering %s \n',sz_filename);   
        fprintf(1,'[');         
        
        for idx_offset = 1:size(response_matrix_avg{idx_band}{idx_stimcode},2)

            fprintf(1,'.');        
            
            topos.values = response_matrix_avg{idx_band}{idx_stimcode}(:,idx_offset);

            results = framework_render_electrodes_eeg(topos);

            figure(1),
            imshow(results.image)       
            offset_time = list_offset(idx_offset) * 1/ecog.param.samplingrate * 1000;
            sz_title_top = sprintf('evoked potential %d-%d Hz, %s',ecog.param.f_low_pass(idx_band),ecog.param.f_high_pass(idx_band),sz_label{idx_stimcode});
            sz_title_bottom = sprintf('%3.0fms',offset_time);
            annotation1 = annotation(gcf,'textbox','Position',[0 0.95 1 0],'LineStyle','none','FitHeightToText','off','FontSize',25,'FontWeight','bold','HorizontalAlignment','center','String',sz_title_top);   
            annotation2 = annotation(gcf,'textbox','Position',[0 0.17 1 0],'LineStyle','none','FitHeightToText','off','FontSize',25,'FontWeight','bold','HorizontalAlignment','center','String',sz_title_bottom);   
            pause(0.1)

            frame = getframe;
            writeVideo(obj_video,frame);        
        
        end
        
        close(obj_video);
        
        fprintf(1,'] done\n');                    
        
    end
end


%% render evoked power (brain)

addpath([cd,'/activeBrain']);
load([cd,'/brain_models/','AMC078.mat']);

viewstruct.viewvect     = [90, 0]; % view
viewstruct.lightpos     = [150, 0, 0]; % light

viewstruct.what2view = {'brain', 'activations'}; % remove 'activations' to only visualize the brain

% sets min and max of colorbar
cmapstruct.enablecolorbar = 1;
cmapstruct.enablecolormap = 1;
cmapstruct.fading         = 0;
c_steps                   = 256;
cmapstruct.cmap           = jet(c_steps);

for idx_band = 1:size(power_matrix_avg,2),
    for idx_stimcode = 1:size(power_matrix_avg{idx_band},2),

        cmapstruct.cmax = max(max(abs(power_matrix_avg{idx_band}{idx_stimcode})));
        cmapstruct.cmin = -cmapstruct.cmax;        
        
        sz_filename = sprintf('./results/result_evoked_power_brain_%d_%d_Hz_%s.mp4',ecog.param.f_low_pass(idx_band),ecog.param.f_high_pass(idx_band),sz_label{idx_stimcode});
        
        obj_video = VideoWriter(sz_filename);
        obj_video.FrameRate = 10; 
        
        open(obj_video);
        
        fprintf(1, '>> Rendering %s \n',sz_filename);   
        fprintf(1,'['); 
        
        for idx_offset = 1:size(power_matrix_avg{idx_band}{idx_stimcode},2)

            fprintf(1,'.');        
            
            tala.activations = power_matrix_avg{idx_band}{idx_stimcode}(:,idx_offset);
            tala.activations(isnan(tala.activations)) = 0;
            
            % plots viewstruct.what2view (i.e. brain w/wo activations)
            close(gcf)
            figure('Visible','Off','MenuBar','none','ToolBar','none','Position',[100,100,1024,768])

            cbar = activateBrain(cortex, vcontribs, tala, ix, cmapstruct, viewstruct); 
 %           set(cbar,'Location','south');
            set(cbar,'Ticks',[1,c_steps/4,c_steps/2,c_steps/4*3,c_steps+1]);
            set(cbar,'TickLabels',{sprintf('%2.2f',cmapstruct.cmin),sprintf('%2.2f',cmapstruct.cmin/2),'0',sprintf('%2.2f',cmapstruct.cmax/2),sprintf('%2.2f',cmapstruct.cmax)}) 
            set(cbar,'FontSize',12,'FontWeight','bold');
            
            offset_time = list_offset(idx_offset) * 1/ecog.param.samplingrate * 1000;
            sz_title_top = sprintf('evoked power %d-%d Hz, %s',ecog.param.f_low_pass(idx_band),ecog.param.f_high_pass(idx_band),sz_label{idx_stimcode});
            sz_title_bottom = sprintf('%3.0fms',offset_time);
            annotation1 = annotation(gcf,'textbox','Position',[0 0.95 1 0],'LineStyle','none','FitHeightToText','off','FontSize',25,'FontWeight','bold','HorizontalAlignment','center','String',sz_title_top);   
            annotation2 = annotation(gcf,'textbox','Position',[0 0.17 1 0],'LineStyle','none','FitHeightToText','off','FontSize',25,'FontWeight','bold','HorizontalAlignment','center','String',sz_title_bottom);   
           
            frame = getframe(gcf);
            
            writeVideo(obj_video,frame);
        
        end
        
        close(obj_video);

        fprintf(1,'] done\n');            

    end
end


%% render evoked potential (brain)

addpath([cd,'/activeBrain']);
load([cd,'/brain_models/','AMC078.mat']);

viewstruct.viewvect     = [90, 0]; % view
viewstruct.lightpos     = [150, 0, 0]; % light

viewstruct.what2view = {'brain', 'activations'}; % remove 'activations' to only visualize the brain

% sets min and max of colorbar
cmapstruct.enablecolorbar = 1;
cmapstruct.enablecolormap = 1;
cmapstruct.fading         = 0;
c_steps                   = 256;
cmapstruct.cmap           = jet(c_steps);

for idx_band = 1:size(power_matrix_avg,2),
    for idx_stimcode = 1:size(power_matrix_avg{idx_band},2),

        cmapstruct.cmax = max(max(abs(response_matrix_avg{idx_band}{idx_stimcode})));
        cmapstruct.cmin = -cmapstruct.cmax;        
        
        sz_filename = sprintf('./results/result_evoked_potential_brain_%d_%d_Hz_%s.mp4',ecog.param.f_low_pass(idx_band),ecog.param.f_high_pass(idx_band),sz_label{idx_stimcode});
        
        obj_video = VideoWriter(sz_filename);
        obj_video.FrameRate = 10; 
        
        open(obj_video);
        
        fprintf(1, '>> Rendering %s \n',sz_filename);   
        fprintf(1,'['); 
        
        for idx_offset = 1:size(power_matrix_avg{idx_band}{idx_stimcode},2)

            fprintf(1,'.');        
            
            tala.activations = response_matrix_avg{idx_band}{idx_stimcode}(:,idx_offset);
            tala.activations(isnan(tala.activations)) = 0;
            
            % plots viewstruct.what2view (i.e. brain w/wo activations)
            close(gcf)
            figure('Visible','Off','MenuBar','none','ToolBar','none','Position',[100,100,1024,768])

            cbar = activateBrain(cortex, vcontribs, tala, ix, cmapstruct, viewstruct); 
 %           set(cbar,'Location','south');
            set(cbar,'Ticks',[1,c_steps/4,c_steps/2,c_steps/4*3,c_steps+1]);
            set(cbar,'TickLabels',{sprintf('%2.2fuV',cmapstruct.cmin),sprintf('%2.2fuV',cmapstruct.cmin/2),'0uV',sprintf('%2.2fuV',cmapstruct.cmax/2),sprintf('%2.2fuV',cmapstruct.cmax)}) 
            set(cbar,'FontSize',12,'FontWeight','bold');
            
            offset_time = list_offset(idx_offset) * 1/ecog.param.samplingrate * 1000;
            sz_title_top = sprintf('evoked potential %d-%d Hz, %s',ecog.param.f_low_pass(idx_band),ecog.param.f_high_pass(idx_band),sz_label{idx_stimcode});
            sz_title_bottom = sprintf('%3.0fms',offset_time);
            annotation1 = annotation(gcf,'textbox','Position',[0 0.95 1 0],'LineStyle','none','FitHeightToText','off','FontSize',25,'FontWeight','bold','HorizontalAlignment','center','String',sz_title_top);   
            annotation2 = annotation(gcf,'textbox','Position',[0 0.17 1 0],'LineStyle','none','FitHeightToText','off','FontSize',25,'FontWeight','bold','HorizontalAlignment','center','String',sz_title_bottom);   
           
            frame = getframe(gcf);
            
            writeVideo(obj_video,frame);
        
        end
        
        close(obj_video);

        fprintf(1,'] done\n');            

    end
end

