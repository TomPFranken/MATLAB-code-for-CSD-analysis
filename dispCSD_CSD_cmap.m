function dispCSD_CSD_cmap(dataFile,dispTitle,H,ignorechans,CSD_max)
%
%
fn = [dataFile '_LFP.mat'];
if exist(fn,'file')
    load(fn)
else
    fprintf(1,'File %s does not exist !!!\n',fn);
    return;
end

if ~exist('dispTitle','var') || isempty(dispTitle)
    dispTitle = 0;
end

if exist('H','var')
    try
        figure(H);
    catch
        subplot(H);
    end
else
    H = figure;
end
hold on

LFP = LFP_SEG_MEAN;
CSD = CSD_SEG_MEAN; %microV

nchans = size(LFP_SEG_MEAN,1);

ignorechans_corr=sort(nchans-ignorechans+1);
ignorechans_corr_forCSD=sort((nchans-1)-setdiff([min(ignorechans)-1 ignorechans max(ignorechans)+1],[0 1 nchans nchans+1])+1);

LFP=LFP(setdiff(1:size(LFP,1),ignorechans_corr),:);
CSD=CSD(setdiff(1:size(CSD,1),ignorechans_corr_forCSD),:);


x = -Twin:Twin;
cmap = colormap_redblackblue();

if isnan(CSD_max)
CSD_max = max(abs(CSD(:)));
end

imagesc(x,1:size(CSD,1),CSD);
set(gca,'Clim',[-CSD_max CSD_max]);
colormap(cmap);
axis tight
set(gca,'YDir','reverse');
set(gca,'Ytick',[]);
plot([0 STIM_DUR],[1 1],'w','LineWidth',2);
set(gca,'Xlim',[-30 Twin]);

if dispTitle
    titleStr = dataFile;
    title(titleStr,'interpreter','none');
end


