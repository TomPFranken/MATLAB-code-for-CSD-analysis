clear all;

datapath=['CSD_forVahed\CSD_20240222\'];

currdate='20240222';

filenameML=[currdate(3:end) '_Buzz_CSDtask.bhv2'];basefilenameIntan=[currdate '_' currdate(3:end)]; filename_np = ['Buzz' currdate '-CSD_g1_t0.imec0.lf.bin'];
dorectify=0; dozscore=0;

% load MonkeyLogic data and extract codes
[allcodes_ML,rising,code_start,code_end] = getMonkeyLogicCodes(datapath(1:end-1),filenameML);
% load Intan data, extract barcodes and eventcodes
[barcodes_intan, eventcodes_intan] = readIntanCodes(datapath,basefilenameIntan,rising);
if sum(mod(eventcodes_intan(:,1),128)==mod(allcodes_ML(:,1),128)) ~= size(allcodes_ML,1), error('event codes do not match'); end
eventcodes_intan(:,1)=allcodes_ML(:,1);
% load neuropixels data
fid = fopen([datapath,filename_np], 'r');
dat = double(fread(fid, [385 Inf], '*int16'));

lateralpositionwanted=0;

fclose(fid);
metadata=SGLX_readMeta.ReadMeta(filename_np,datapath);
%extract barcodes from neuropixels data
barcodedata_np=dat(end,:);
allbarcodes_np = getBarCodes(barcodedata_np,str2num(metadata.imSampRate));
%align eventcodes from MonkeyLogic/Intan to timing of Neuropixels data, using the barcodes
timestamps_np=alignTimestampsFromBarcodes(eventcodes_intan(:,2),allbarcodes_np,barcodes_intan);
eventcodedata=[eventcodes_intan(:,1) timestamps_np]; %first column is event codes, second column is time stamps in neuropixels time values (samples)
chRow=int32(0:383)'; %included channels (?)
probPos=[repmat([0;32],round(384/2),1) sort(repmat([0:20:20*(round(384/2)-1)]',2,1))]; % % column 1 = lateral position of electrode; column 2: z-position of electrode (0=tip)
posschannels=find(probPos(:,1)==lateralpositionwanted); %make sure all channels are from the same column on the probe
channels_selected=posschannels(1:3:numel(posschannels));
% channels_selected=posschannels(10:3:numel(posschannels));
unpitch=unique(diff(probPos(channels_selected,2)));
if numel(unpitch)~=1
    error('pitch not unique');
end
interpolfactor=unpitch/10; %make sure to interpolate every 10 um
nchannels_selected=numel(channels_selected);
FsLFP=str2num(metadata.imSampRate);
Fsorig=FsLFP;
probemap=[flipud([1:nchannels_selected]') [1:nchannels_selected]'];

%transform data to voltage
%(https://billkarsh.github.io/SpikeGLX/Sgl_help/Metadata_30.html)
Imax = str2num(metadata.imMaxInt);
Vmax =  str2num(metadata.imAiRangeMax); 
gain = str2num(metadata.imChan0lfGain); 
dat = dat.* (Vmax/Imax/gain) .*1e6; %microV?

%load task data to get Task Object per trial
b=mlread([datapath filenameML]);
alltaskobjects=[];b
for i=1:numel(b)
currtaskobject=b(i).TaskObject.Attribute{2}{2};
currtaskobject=currtaskobject(max(strfind(currtaskobject,'\'))+1:strfind(currtaskobject,'.')-1);
alltaskobjects=[alltaskobjects {currtaskobject}];
end
clear currtaskobject;

whichtaskobjects={'all'}; %if 'all', any object included 

clear data;
data.codeData.codes=(eventcodedata(:,1)');
data.codeData.times=(eventcodedata(:,2)');
tempdat=dat(channels_selected,:);
tempLFP=[];
for i=1:nchannels_selected
    disp(i)
tempLFP=[tempLFP {double(tempdat(i,:))}];
end
data.LFP=tempLFP;

%fixed code values
code_start=[9];code_end=[18];%start and end of trial
mincode=100;
rewardcode=117;

starttrials=strfind(eventcodedata(:,1)',code_start);
endtrials=strfind(eventcodedata(:,1)',code_end);

strobecodes=data.codeData.codes;
strobetimestamps=data.codeData.times./FsLFP; %sec

allcsdmaps=[];
for wi=1:numel(whichtaskobjects)

    lfpdata=[];
    timeinfo=[];
    condnos=[];
    totalflashes=0;
    for i=1:numel(starttrials)

        currtaskobject=alltaskobjects{i};

        starti=starttrials(i);
        endi=(endtrials(find((endtrials-starttrials(i)>0)==1, 1 )));
        curreventcodes=strobecodes(starti:endi);
        curreventtimes=strobetimestamps(starti:endi)-strobetimestamps(starti); %sec

        f=find(curreventcodes>mincode & curreventcodes<rewardcode);

        if sum(strcmpi(whichtaskobjects(wi),currtaskobject))==1 || sum(strcmpi(whichtaskobjects(wi),'all'))==1

            if numel(f)>=4

                flashon_times=curreventtimes(f(1:2:end-3)); %ignore last flash
                flashoff_times=curreventtimes(f(2:2:end-2));
                trialstarttime_absol=strobetimestamps(starti);
                trialend_time=curreventtimes(end);


                if numel(intersect(curreventcodes,6))>0

                    totalflashes=totalflashes+numel(flashon_times);

                    condnos=[condnos curreventcodes(f(1))-mincode];

                    clear temp;
                    temp.flashon_times=round(flashon_times.*FsLFP); %samples
                    temp.flashoff_times=round(flashoff_times.*FsLFP);
                    temp.trialend_time=round(trialend_time.*FsLFP);
                    timeinfo=[timeinfo {temp}];

                    sample_start=round(trialstarttime_absol.*FsLFP);
                    trial_dur_samples=round(trialend_time.*FsLFP);

                    temp=[];
                    for j=1:size(probemap,1)
                        %               data=nexFile.contvars{omneticstoatlas(j,1)}.data(sample_start:sample_start+trial_dur_samples-1);

                        currdata=double(data.LFP{probemap(j,1)}(sample_start:sample_start+trial_dur_samples-1));

                        temp=[temp;{currdata}];
                    end
                    lfpdata=[lfpdata temp];

                end
            end
        end
    end

    clear miga;
    miga.lfp_mat=lfpdata; %units? Need to be scaled to micro-V

    miga.event_mat=timeinfo;
    miga.channelnum=flipud(probemap(:,end));
    miga.condno=condnos;

    save([datapath 'miga'],'miga','-v7.3');
    delete([datapath 'miga_LFP.mat']);
    delete([datapath 'miga_LFP_lock.mat']);
    extractLFP_csdfltrg([datapath 'miga'],1,0,[],FsLFP,dorectify,dozscore,unpitch)

    CSD_max = NaN;

    [f,smoothax,lfpax]=dispCSD_composite([datapath 'miga'],1,...
        [],'cubic',interpolfactor,[-100 200],[nchannels_selected 1],CSD_max,dorectify,FsLFP);

    allcsdmaps=[allcsdmaps smoothax];
end

f_csd=figure;
C = copyobj(allcsdmaps,f_csd);

for i=1:numel(C)
    axes(C(i));
set(C(i),'Units','centimeters','Position',[2+(12*(i-1)) 2 10 24]);
cb=colorbar;
set(cb,'Units','centimeters','Position',[2+(12*(i-1))+10.1 22 0.3 3],'TickDir','out');
set(get(cb,'YLabel'),'String','CSD (uA.mm-3)')
end
colormap(get(f,'ColorMap'))
set(f_csd,'Position',[-2.5250    0.0417    2.4547    1.3193].*1e3);

close(f)

%esthetic changes
xticks_wanted=[0:50:150];
xticks_real=xticks_wanted;
set(gca,'XTick',xticks_real,'XTickLabel',xticks_wanted);
xlim([min(xticks_wanted) max(xticks_wanted)]);
set(gca,'TickLength',[0.005 0.005]);
curryticklabels=str2num(get(gca,'YTickLabel'));

%plot actual electrode row numbers: row 1: deepest pair of electrode
%channels. Top row is row 192. Use rows rather than electrode numbering
%because one can do interpolation for depth of rows, but not for individual
%electrodes
electroderows=[];
for ci=1:numel(channels_selected)
    electroderows=[electroderows;find(posschannels==(channels_selected(ci)))];
end
newyticklabels_electroderows=flipud(electroderows(2:end-1));
corryticks_electroderows=get(gca,'YTick');
set(gca,'YTickLabel',newyticklabels_electroderows);

%add in between electrode numbers
electroderowswanted=1:1:192;
rico=(corryticks_electroderows(2)-corryticks_electroderows(1))/(newyticklabels_electroderows(2)-newyticklabels_electroderows(1));
offset=corryticks_electroderows(1)-newyticklabels_electroderows(1)*rico;
ytickswanted=rico*electroderowswanted+offset;
[ytickswanted,sorti]=sort(ytickswanted);
electroderowswanted=electroderowswanted(sorti);
set(gca,'YTickLabel',electroderowswanted,'YTick',ytickswanted);

ylabel('Electrode row (1 = deep, near electrode tip)');
xlabel('Time re. stimulus onset (ms)')

% use barcodes to match event codes from CSD task to time stamps that correspond to the spike times derived from the ap.bin file
load([datapath 'barcodes_np.mat'])
timestamps_np=alignTimestampsFromBarcodes(eventcodes_intan(:,2),barcodes_np,barcodes_intan);
eventcodes_np=[eventcodes_intan(:,1) timestamps_np]; %first column is event codes, second column is time stamps in neuropixels time values (samples)

%get spike times
datapath_ksdata = [datapath 'kilosort4'];
spike_times_phy=readNPY([datapath_ksdata '\spike_times.npy']);
spike_clusters_phy=readNPY([datapath_ksdata '\spike_clusters.npy']);
spike_templates_phy=readNPY([datapath_ksdata '\spike_templates.npy']); %this will be the same as spike_clusters_phy if not manually curated in phy
channelpos=readNPY([datapath_ksdata '\channel_positions.npy']); %coordinates of electrodes on probe
t = readtable([datapath_ksdata '\cluster_KSLabel.tsv'], "FileType","text",'Delimiter', '\t'); %read KS cluster labels
spike_positions=readNPY([datapath_ksdata '\spike_positions.npy']);
selectedUnits_phyinds=t{strcmpi(t{:,2},'good'),1}; %only include units labeled as 'good' by KS
selectedUnits_depths=nan(size(selectedUnits_phyinds));
for tempi=1:numel(selectedUnits_phyinds)
    allspikes=(find(spike_clusters_phy==selectedUnits_phyinds(tempi)));
    currdepth=double(mean(spike_positions(allspikes,2))); %compute depth as average depth (?) across all spike occurrences. Probably these are not identical numbers because of drift correction?
    selectedUnits_depths(tempi)=currdepth;
end