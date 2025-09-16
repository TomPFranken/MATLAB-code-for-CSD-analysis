%edits
%20231209: TF add fs as input parameter; removed factors Fs/1e3 (not
%correct for FsLFP different from 1000)
%20240118: add rectify, add z-score (Bijanzadeh et al). Change conv2 to
%imgaussfilt
%20240419: change filters to same settings as eLife paper. Need to check if
%this works

function extractLFP_csdfltrg(dataFile,saveFlag,isMUA,ignorechans,Fs,dorectify,dozscore,electrodepitch)
%
%
%

if isMUA
outfn = [dataFile '_MUA.mat'];  
lockfn = [dataFile '_MUA_lock.mat'];
else
outfn = [dataFile '_LFP.mat'];
lockfn = [dataFile '_LFP_lock.mat'];
end

if (saveFlag==1) && exist(outfn,'file')
    fprintf(1,'File %s already exists !!!\n',outfn);
    return
end

% attempt to create a lock on the file
if ~exist(lockfn,'file')
    save(lockfn,'lockfn','-v7.3');
else
    return; % some other process has the file locked
end


% load the miga file
fn = [dataFile '.mat'];

if exist(fn,'file')
    load(fn)
else
    fprintf(1,'File %s does not exist !!!\n',fn);
    return;
end


numTrials = length(miga.condno);

numChannels = length(miga.channelnum);

LFP_SEG_MEAN = [];
LFP_SEG_STD  = [];
CSD_SEG_MEAN = [];
STIM_DUR     = 0;

Twin_sec = 0.2; %in s

N_SAMPLES = 0;

Twin = round(Twin_sec*Fs); %samples

% Create filters for filtering and smoothing the LFPs
% For now these parameters are copied from the Plexcsd code
[NL DL] = butter(4,88/(Fs/2),'low');
[NH DH] = butter(2,3.3/(Fs/2),'high');

for trial=1:numTrials

    E = miga.event_mat{trial};
    if isempty(E)
        continue;
    end

    numFlash = numel(E.flashon_times);
    if numFlash == 0
        continue;
    end

    % extarct LFP segment for each flash
    for i=1:numFlash
        tWin_start = round(E.flashon_times(i))-Twin; %samples
        tWin_stop  = round(E.flashon_times(i))+Twin;

        if (tWin_start < 1) || (tWin_stop > round(E.trialend_time))
            continue;
        end

        N_SAMPLES = N_SAMPLES+1;

        LFP_TRACES = zeros(numChannels,tWin_stop-tWin_start+1);
        for c=1:numChannels
            % extract the LFP
            if isMUA
            LFP = miga.mua_mat{c,trial};    
            else
            LFP = miga.lfp_mat{c,trial};
            end
            
            try
            LFPtrace = LFP(tWin_start:tWin_stop);
            catch ME
disp('');
            end
            % subtract mean and scale
            LFPtrace = LFPtrace-mean(LFPtrace); %micro-V

            % filter %probably better to do this before windowing. Then
            % subtracting mean may not be necessary?
            LFPtrace = filter(NH,DH,LFPtrace); 
            LFPtrace = filter(NL,DL,LFPtrace);

            LFP_TRACES(c,:) = LFPtrace;
        end
        % smooth
        LFP_AUG = [LFP_TRACES(1,:); LFP_TRACES; LFP_TRACES(end,:)];

        LFP_AUG = imgaussfilt(LFP_AUG,((3*20)/electrodepitch)); %seems like 2 is a good number

        LFP_TRACES = LFP_AUG(2:end-1,:); %microV

        LFP_SEG(:,:,N_SAMPLES) = LFP_TRACES; %micro-V

        STIM_DUR = STIM_DUR+ round(E.flashoff_times(i))- round(E.flashon_times(i));
    end

end % trials

LFP_SEG_MEAN = mean(LFP_SEG,3); %micro-V
LFP_SEG_STD  = std(LFP_SEG,0,3);



STIM_DUR = STIM_DUR/N_SAMPLES;

if isMUA
    MUA_SEG=LFP_SEG;
    MUA_SEG_MEAN=LFP_SEG_MEAN;
    MUA_SEG_STD=LFP_SEG_STD;
save(outfn,'MUA_SEG','MUA_SEG_MEAN','MUA_SEG_STD', ...
           'STIM_DUR', ...
           'N_SAMPLES','Twin','-v7.3');
else
    
    currnchan=size(LFP_SEG,1);

    ignorechan_inds=sort(currnchan-ignorechans+1);
    goodchaninds=setdiff(1:currnchan,ignorechan_inds);
    goodchansCSD=sort(setdiff([currnchan:-1:1],ignorechans),'descend');
    
    CSD_SEG = -diff(LFP_SEG(goodchaninds,:,:),2); %micro-V
    CSD_SEG_MEAN = mean(CSD_SEG,3); %micro-V

    % transform to uA/mm-3 (same procedure as in Poort et al. 2016)
    CSD_SEG_MEAN= CSD_SEG_MEAN.*1e-6; %V
    CSD_SEG_MEAN = (CSD_SEG_MEAN./(100*1e-6))./(100*1e-6); %V/(m^2)
    CSD_SEG_MEAN = CSD_SEG_MEAN.*0.4; %A/m3
    CSD_SEG_MEAN = CSD_SEG_MEAN.*1e6; %uA/m3
    CSD_SEG_MEAN = CSD_SEG_MEAN./1e9; %uA/mm3

    if dozscore
        % %across channels z score
        % CSDb_mean=mean(CSD_SEG_MEAN(:,Twin-(0.1*Fs)));
        % CSDb_std=std(CSD_SEG_MEAN(:,Twin-(0.1*Fs)));
        % CSD_SEG_MEAN=(CSD_SEG_MEAN-CSDb_mean)./CSDb_std;

        % %per channel Z score: result seems to be less noisy than across
        % %channels z score
        for currchani=1:size(CSD_SEG_MEAN,1)
        CSDb_mean=mean(CSD_SEG_MEAN(currchani,(Twin-(0.1*Fs)):Twin));
        CSDb_std=std(CSD_SEG_MEAN(currchani,(Twin-(0.1*Fs)):Twin));
        CSD_SEG_MEAN(currchani,:)=(CSD_SEG_MEAN(currchani,:)-CSDb_mean)./CSDb_std;
        end
    end

    if dorectify
    CSD_SEG_MEAN(CSD_SEG_MEAN>0)=0; %rectify
    end

save(outfn,'LFP_SEG','LFP_SEG_MEAN','LFP_SEG_STD', ...
           'CSD_SEG_MEAN', ...
           'STIM_DUR', ...
           'N_SAMPLES','Twin','goodchansCSD','-v7.3');
end
% delete the lock file
delete(lockfn);

