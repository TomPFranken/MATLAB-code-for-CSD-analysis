function numChannels=dispCSD_LFP_traces(dataFile,dispTitle,H,ignorechans)
%
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

LFP = LFP_SEG_MEAN; %uV
CSD = CSD_SEG_MEAN;

nchans = size(LFP_SEG_MEAN,1);

ignorechans_corr=sort(nchans-ignorechans+1);
ignorechans_corr_forCSD=sort((nchans-1)-setdiff([min(ignorechans)-1 ignorechans max(ignorechans)+1],[0 1 nchans nchans+1])+1);

LFP=LFP(setdiff(1:size(LFP,1),ignorechans_corr),:);
CSD=CSD(setdiff(1:size(CSD,1),ignorechans_corr_forCSD),:);

numChannels = size(LFP,1);

ytickpos=[];yticklabels=[];

x = -Twin:Twin;
SCALE = max(abs([min(LFP(:)) max(LFP(:))]));
SHIFT = 1.2;
SHIFT_OFFSET = 0;
for i=numChannels:-1:1
    y = LFP(i,:)/SCALE;
    
    SHIFT_OFFSET = SHIFT_OFFSET+SHIFT;
    y = y+SHIFT_OFFSET;
    plot(x,y,'k');
    
    ytickpos=[ytickpos y(x==0)];
    yticklabels=[yticklabels i];
end

%plot  a scale
scalevals=[-1 0 1].*200; %uV
% scalevals=[1 0 1].*200; %uV
y = scalevals./SCALE;
SHIFT_OFFSET = SHIFT_OFFSET+SHIFT;
y = y+SHIFT_OFFSET;
xoffset=50;
plot([0 0]+xoffset,y([1 3]),'k-')
plot(([-1 1].*5)+xoffset,y([1 1]),'k-')
plot(([-1 1].*5)+xoffset,y([3 3]),'k-')
plot(([-1 1].*5)+xoffset,y([2 2]),'k-')

set(gca,'Ytick',ytickpos,'YTickLabel',fliplr(yticklabels));
% set(gca,'YDir','reverse');
% plot([0 STIM_DUR],[1 1]*SHIFT/2,'k','LineWidth',2);

if dispTitle
    titleStr = dataFile;
    title(titleStr,'interpreter','none');
end

