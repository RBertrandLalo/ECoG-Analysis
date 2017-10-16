function visualize_bandpass_filters(bandpass,SamplingRate,param)

warning('off', 'signal:sigtools:fvtool:fvtool:useLegendMethod');

if isfield(bandpass{1},'b')
    
    vis.bandpass.sz_filter_string = 'fvtool(';
    
    for idx=1:length(bandpass),
        vis.bandpass.sz_filter_string = [vis.bandpass.sz_filter_string,sprintf('bandpass{%d}.b,[1],',idx)];
        vis.bandpass.legend{idx} = sprintf('%d-%d Hz',param.f_low_pass(idx),param.f_high_pass(idx));        
    end
    
    vis.bandpass.sz_filter_string = [vis.bandpass.sz_filter_string,'''FrequencyScale'',''linear'',''Fs'',SamplingRate);'];

    eval(vis.bandpass.sz_filter_string);
%    xlim([0,max(param.f_high_pass)+50]);
%    ylim([-80,10]);
    legend(vis.bandpass.legend);
    title('Band-pass filters','FontSize',18,'FontWeight','bold');    
    
elseif isfield(bandpass{1},'h')
    
    for idx=1:length(bandpass),
        vis.bandpass.h(idx) = bandpass{idx}.h;
        vis.bandpass.legend{idx} = sprintf('%d-%d Hz',param.f_low_pass(idx),param.f_high_pass(idx));
    end
    
    fvtool(vis.bandpass.h,'FrequencyScale','linear','Fs',SamplingRate);
%    xlim([0,max(param.f_high_pass)+50]);
%    ylim([-400,10]);
    legend(vis.bandpass.legend);
    title('Band-pass filters','FontSize',18,'FontWeight','bold');
end
