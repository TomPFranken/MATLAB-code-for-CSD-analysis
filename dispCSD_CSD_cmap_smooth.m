function dispCSD_CSD_cmap_smooth...
    (dataFile,dispTitle,H,ignorechans,interpolationmethod,yinterpolfactor,xlims,ylims_ch,CSD_max,dorectify,FsLFP)
%
%
%
% yinterpolfactor=2;


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
CSD = CSD_SEG_MEAN;

nchans = size(LFP_SEG_MEAN,1);

ignorechans_corr=sort(nchans-ignorechans+1);
ignorechans_corr_forCSD=sort((nchans-1)-setdiff([min(ignorechans)-1 ignorechans max(ignorechans)+1],[0 1 nchans nchans+1])+1);

LFP=LFP(setdiff(1:size(LFP,1),ignorechans_corr),:);
CSD=CSD(setdiff(1:size(CSD,1),ignorechans_corr_forCSD),:);


x = -Twin:Twin; %samples
xtimes= (x./FsLFP).*1e3; %ms

xlims_samples=round(xlims.*1e-3.*FsLFP);

cmap = colormap_redblackblue();
I_FACTOR = yinterpolfactor;
M = calcInterpMap(CSD,size(CSD,2),size(CSD,1)*I_FACTOR,interpolationmethod);


goodxi=find(x>=min(xlims_samples) & x<=max(xlims_samples));
M.ZI=M.ZI(:,goodxi);
xtimes=xtimes(goodxi);

ylims_ch_inv=nchans-ylims_ch+1;
startpoint=(I_FACTOR/2)+0.5-I_FACTOR;
yticks=startpoint+((ylims_ch_inv-1)*I_FACTOR);
tempy=1:size(M.ZI,1);
goodyi=find(tempy>=(min(yticks)-(0.5*I_FACTOR)) & tempy<=(max(yticks)+(0.5*I_FACTOR)));
M.ZI=M.ZI(goodyi,:);
yticks_lim=yticks-(min(goodyi)-1);  %-min(goodyi...) because some rows have been removed

yticks_ch=[(nchans-1):-1:2];
yticks_ch_inv=[nchans-yticks_ch+1]; 
yticks=startpoint+((yticks_ch_inv-1)*I_FACTOR)-(min(goodyi)-1); %-min(goodyi...) because some rows have been removed

yticks_lim=[yticks_lim(1)-(yticks(2)-yticks(1)).*0.5 yticks_lim(2)+(yticks(2)-yticks(1)).*0.5];

if isnan(CSD_max)
CSD_max = max(abs(M.ZI(:)));
disp(['CSD_max: ' num2str(CSD_max)])
end

clims=[-CSD_max CSD_max];
if dorectify
clims=[-CSD_max 0];
cmap=cmap(1:ceil(size(cmap,1)/2),:);
end

imagesc(xtimes,1:size(M.ZI,1),M.ZI);
set(gca,'Clim',clims);
colormap(cmap);
axis tight
set(gca,'YDir','reverse');
set(gca,'Ytick',yticks,'YTickLabel',yticks_ch,'TickDir','out');
plot([0 STIM_DUR],[1 1]*I_FACTOR,'w','LineWidth',2);
set(gca,'Xlim',xlims);

if dispTitle
    titleStr = dataFile;
    title(titleStr,'interpreter','none');
end

cb=colorbar;
set(cb,'Units','centimeters','Position',[18 11.5 0.3 3],'TickDir','out')

