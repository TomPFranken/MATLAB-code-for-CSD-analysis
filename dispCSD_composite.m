function [f,smoothax,lfpax]=...
    dispCSD_composite(dataFile,doplot,ignorechans,interpolationmethod,yinterpolfactor,xlims,ylims_ch,CSD_max,dorectify,FsLFP)

if doplot
screenSize = get(0,'Screensize');

f=figure

numSubFigs = 5;
subPlotIndx = 0;

subPlotIndx = subPlotIndx+1;
h_sp = mysubplot(1,numSubFigs,subPlotIndx,0.1);
nchannels=dispCSD_LFP_traces(dataFile,1,h_sp,ignorechans);
lfpax=gca;
dispCSD_LFP_traces(dataFile,1,h_sp,ignorechans);

subPlotIndx = subPlotIndx+1;
h_sp = mysubplot(1,numSubFigs,subPlotIndx,0.1);
CSDtraceax=dispCSD_CSD_traces(dataFile,0,h_sp,ignorechans);
set(get(CSDtraceax,'Title'),'String','CSD traces');

l=load([dataFile '_LFP'],'goodchansCSD');
goodchansCSD=l.goodchansCSD;

goodchansCSD=sort(setdiff(goodchansCSD,ignorechans),'descend');

yticks=get(CSDtraceax,'YTick');
dist=mean(diff(yticks));
set(CSDtraceax,'YLim',[yticks(1)-dist*0.5 yticks(end)+dist*3]);
templabels=sort(goodchansCSD);
yticklabels=set(CSDtraceax,'YTickLabel',templabels(2:end-1))

subPlotIndx = subPlotIndx+1;
h_sp = mysubplot(1,numSubFigs,subPlotIndx,0.1);
dispCSD_CSD_cmap(dataFile,0,h_sp,ignorechans,CSD_max);
cmapax=gca;
nchansplotted=(numel(goodchansCSD)-2);

set(gca,'YTick',[min(goodchansCSD)+1:max(goodchansCSD)-1]-1,'YTickLabel',max(goodchansCSD)-1:-1:min(goodchansCSD)+1,'TickDir','out')

subPlotIndx = subPlotIndx+1;
h_sp = mysubplot(1,numSubFigs,subPlotIndx,0.1);
dispCSD_CSD_cmap_smooth(dataFile,0,h_sp,ignorechans,interpolationmethod,yinterpolfactor,xlims,ylims_ch,CSD_max,dorectify,FsLFP);
ylims=get(gca,'YLim');

smoothax=gca;

end


set(f,'Position',[1 60 screenSize(3)*0.9 screenSize(4)*0.3]);
end
