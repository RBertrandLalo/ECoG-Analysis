function visualize_notch_filters(bandpass,SamplingRate,param)

warning('off', 'signal:sigtools:fvtool:fvtool:useLegendMethod');

isfield(bandpass{1},'b')

vis.bandpass.sz_filter_string = 'fvtool(';

for idx=1:length(bandpass),
    vis.bandpass.sz_filter_string = [vis.bandpass.sz_filter_string,sprintf('bandpass{%d}.b,bandpass{%d}.a,',idx,idx)];
    vis.bandpass.legend{idx} = sprintf('%d Hz',param.filter.notch.fcenter(idx));        
end

vis.bandpass.sz_filter_string = [vis.bandpass.sz_filter_string,'''FrequencyScale'',''linear'',''Fs'',SamplingRate);'];

eval(vis.bandpass.sz_filter_string);
%    xlim([0,max(param.f_high_pass)+50]);
%    ylim([-80,10]);
legend(vis.bandpass.legend);
title('Band-pass filters','FontSize',18,'FontWeight','bold');    
    
