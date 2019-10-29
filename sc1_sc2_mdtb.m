function [ varargout ] = sc1_sc2_mdtb( what, varargin )
%[VARARGOUT] = SC1_SC2(WHAT, VARARGIN) This functions does every analysis
%that can be done on the mdtb dataset.
%   Detailed explanation goes here

numDummys = 3;   % number of dummy scans per run
numTRs    = 601; % number of scans per run

%%% setting path for the working directories
baseDir = '/Users/ladan/Documents/Project-Cerebellum/Cerebellum_Data';
% baseDir = '/home/ladan/Documents/Data/Cerebellum-MDTB';

%%% setting directory names
behavDir     ='/data';                  %% behavioral data directory.
suitDir      = 'suit';                  %% directory where the anatomicals used in creating flatmaps are stored.
regDir       = 'RegionOfInterest';      %% The ROI directory 
wbDir        = 'surfaceWb';

suitToolDir  = '/Users/ladan/Documents/MATLAB/suit';

runLst  = 1:16;    %% run numbers to use in GLM
ntp     = 598;     %% number of time points after eliminating the dummy scans

% cd(baseDir)

% Hard-coding some variables
%%% subjects
subj_name = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11',...
    's12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24',...
    's25','s26','s27','s28','s29','s30','s31'};
returnSubjs = [2, 3, 4, 6, 8, 9, 10, 12, 14, 15, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31]; % "good" subjects
%========================================================================================================================
% (3) GLM.

% GLM Directories - change glmNum when appropriate

funcRunNum = [51,66];  % first and last behavioural run numbers (16 runs per subject)
run{1}  = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16'};
run{2}  = {'17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32'};
runB = [51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66];  % Behavioural labelling of runs
sess = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2];                  % session number

runs{1}{1}=1:8;   % study 1, session 1
runs{1}{2}=9:16;  % study 1, session 2
runs{2}{1}=17:24; % study 2, session 1
runs{2}{2}=25:32; % study 2, session 2

hemI    = {'L', 'R'}; % left and right hemisphere
hemName = {'CortexLeft', 'CortexRight'};

warning('off')
switch what
    case 'PHYS:mdtb:get_log' % creates log files for RESP and PULS
        % use the function extractCMRRPhysio.m to get the log files from
        % the dcm file.
        % Example: sc1_sc2_mdtb('PHYS:mdtb:get_log', 'sn', 26, 'sess', 1, 'scan', 3)
        
        sn   = returnSubjs;
        sess = 1:2;    %% set it to either 1 or 2
        scan = 1:8;    %% check out the dicom folder for subject and find the scanning session with all the files
        
        vararginoptions(varargin, {'sn', 'sess', 'scan'});
        
        dicomDir = fullfile(baseDir, 'sc1', 'imaging_data_dicom');
        
        for s = sn
            for ss = sess
                for sca = scan
                    fprintf('\n********** getting the log files for %s session %d scan %d **********\n', subj_name{s}, ss, sca)
                    physSubjDir = fullfile(dicomDir, sprintf('%s_%d', subj_name{s}, ss), 'physio', num2str(sca));
                    cd(physSubjDir)
                    
                    % convert the dicom file into log files
                    %%% get the dicom file name
                    dcm_struct = dir('*.dcm');
                    dcmName    = dcm_struct.name;
                    
                    extractCMRRPhysio(dcmName);
                end % sca (scan)
            end % ss (sess)
        end % s (sn)
    case 'PHYS:mdtb:get_reg' % creates a text file for PULS and RESP regressors
        % using tapas latest version to get the regressors for
        % physiological signals using the example code that Suzanne
        % provided. Have to run 'PHYS:mdtb:get_log' if you don't already
        % have the log files!
        % Example: sc1_sc2_mdtb('PHYS:mdtb:get_reg', 'sn', 26)
        
        sn   = [24, 25, 26, 27, 28, 30, 31];
        sess = 1:2;
        scan = 1:8;
        
        vararginoptions(varargin, {'sn', 'sess', 'scan'})
        
        % all the log files are stored in the dicom directory
        dicomDir = fullfile(baseDir, 'sc1', 'imaging_data_dicom');
        
        % create directory for 
        PhysioDir = fullfile(baseDir, 'Physio');
        dircheck(PhysioDir);
        
        % get the regressors
        for s = sn
            outDir = fullfile(PhysioDir, subj_name{s});
            dircheck(outDir);
            for ss = sess
                for sca = scan
                    % Create default parameter structure with all fields
                    physioName = sprintf('%s_sess%d_scan%d_physio.mat', subj_name{s}, ss, sca);
                    regTxtName = sprintf('%s_sess%d_scan%d_reg.txt', subj_name{s}, ss, sca);
                    physio = tapas_physio_new();
                    
                    % Individual Parameter settings. Modify to your need and remove default settings
                    physio.save_dir = {outDir};
                    physio.log_files.vendor = 'Siemens_Tics';
                    
                    %%% check if PULS and RESP exist
                    cd(fullfile(dicomDir, sprintf('%s_%d', subj_name{s}, ss), 'physio', num2str(sca)))
                    PULS_struct = dir('*PULS.log');
                    RSP_struct  = dir('*RESP.log');
                    log_struct  = dir('*Info.log');
                    if ~isempty(PULS_struct)
                        PULS_name   = PULS_struct.name;
                    else
                        fprintf('WARNING: There is no PULS log file for sess %d scan %d of subject %s\n', ss, sca, subj_name{s});
                        PULS_name = '';
                    end % if there is PULS log file
                    if ~isempty(RSP_struct)
                        RSP_name    = RSP_struct.name;
                    else
                        fprintf('WARNING: There is no RESP log file for sess %d scan %d of subject %s\n', ss, sca, subj_name{s});
                        RSP_name = '';
                    end % if there is PULS log file
                    log_name = log_struct.name;
                    
                    % setting up the parameters for tapas!
                    %%% I will be using the default parameters as in
                    %%% Suzannes's example code, except for number of
                    %%% slices and number scans, etc. which are set
                    %%% according to the data acquisition parameters
                    physio.log_files.cardiac     = {PULS_name};
                    physio.log_files.respiration = {RSP_name};
                    physio.log_files.scan_timing = {log_name};
                    
                    physio.log_files.relative_start_acquisition = 0;
                    physio.log_files.align_scan                 = 'last';
                    physio.scan_timing.sqpar.Nslices            = 48;      % feel free to change for future analysis!
                    physio.scan_timing.sqpar.TR                 = 1;       % feel free to change for future analysis!
                    physio.scan_timing.sqpar.Ndummies           = 0;       % feel free to change for future analysis!
                    physio.scan_timing.sqpar.Nscans             = 601;     % feel free to change for future analysis!
                    physio.scan_timing.sqpar.onset_slice        = 24;      % is set to the middle slice, feel free to change for future analysis!
                    physio.scan_timing.sync.method              = 'scan_timing_log';
                    physio.preproc.cardiac.modality             = 'PPU';
                    physio.preproc.cardiac.filter.include       = false;   % should I set to true?
                    physio.preproc.cardiac.filter.type          = 'butter';
                    physio.preproc.cardiac.filter.passband      = [0.3 9];
                    
                    physio.preproc.cardiac.initial_cpulse_select.method             = 'auto_matched';
                    physio.preproc.cardiac.initial_cpulse_select.max_heart_rate_bpm = 90;
                    physio.preproc.cardiac.initial_cpulse_select.file               = 'initial_cpulse_kRpeakfile.mat';
                    physio.preproc.cardiac.initial_cpulse_select.min                = 0.4;
                    physio.preproc.cardiac.posthoc_cpulse_select.method             = 'off';
                    physio.preproc.cardiac.posthoc_cpulse_select.percentile         = 80;
                    physio.preproc.cardiac.posthoc_cpulse_select.upper_thresh       = 60;
                    physio.preproc.cardiac.posthoc_cpulse_select.lower_thresh       = 60;
                    
                    physio.model.orthogonalise                         = 'none';
                    physio.model.censor_unreliable_recording_intervals = false; % should I set it to true?
                    physio.model.output_multiple_regressors            = regTxtName;
                    physio.model.output_physio                         = physioName;
                    physio.model.retroicor.include                     = false;
                    physio.model.retroicor.order.c                     = 3;
                    physio.model.retroicor.order.r                     = 4;
                    physio.model.retroicor.order.cr                    = 1;
                    physio.model.rvt.include                           = true;
                    physio.model.rvt.delays                            = 0;
                    physio.model.hrv.include                           = true;
                    physio.model.hrv.delays                            = 0;
                    physio.model.noise_rois.include                    = true;
                    physio.model.noise_rois.thresholds                 = 0.9;
                    physio.model.noise_rois.n_voxel_crop               = 0;
                    physio.model.noise_rois.n_components               = 1;
                    physio.model.noise_rois.force_coregister           = 1;
                    physio.model.movement.include                      = false;
                    physio.model.movement.order                        = 6;
                    physio.model.movement.censoring_threshold          = 0.5;
                    physio.model.movement.censoring_method             = 'FD';
                    physio.model.other.include                         = false;
                    physio.verbose.level                               = 2;
                    physio.verbose.process_log                         = cell(0, 1);
                    physio.verbose.fig_handles                         = zeros(0, 1);
                    physio.verbose.use_tabs                            = false;
                    physio.ons_secs.c_scaling                          = 1;
                    physio.ons_secs.r_scaling                          = 1;
                    
                    % Run physiological recording preprocessing and noise modeling
                    [~, R, ~] = tapas_physio_main_create_regressors(physio);
                    save(fullfile(outDir, sprintf('%s_sess%d_scan%d_R.mat', subj_name{s}, ss, sca)), 'R', '-v7.3')
                    close all
                    keyboard;
                    fprintf('\n********** Physio regressors extracted for %s sess %d scan %d **********\n', subj_name{s}, ss, sca);
                end % sca (scan)
            end % ss (sess)
        end % s (sn)
    case 'PHYS:mdtb:lin_reg' % regresses HRV and RVT on the regressor(s) for instructions
        %%% uses one subject
        % Example: sc1_sc2_mdtb('PHYS:mdtb:lin_reg')
        sn         = 24;
        experiment_num = 1; 
        ppmethod   = '';
        glm        = 7;
        sess       = 1:2;
        scan       = 1:8;
        
        vararginoptions(varargin, {'sn', 'experiment_num', 'ppmethod', 'glm', 'sess', 'scan'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
        end
        PhysioDir = fullfile(baseDir, 'Physio');
        
        for s = sn
            
            % load in the SPM file (This could take a while)
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'));
            
            % load in SPM_info file
            T      = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            
            phys_cell = cell(length(sess), length(scan)); % preallocating the cell mat that will have all the linear reg results for a subject
            for ss = sess
                for sca = scan
                    
                    X     = SPM.xX.X(1:598, 1:end - 16);                                        % discarding the intercepts
                    ind   = ((T.sess == ss) & (T.run == sca) & (T.inst == 1) & (T.deriv == 0)); % get the indices for run1, instructions, non-derivatives
                    X_ind = X(:, ind);                                                          
                    X_reg = [zeros(3, 16); X_ind];                                              % My design matrix, adding zeros for dummies!
                    % load in HRV and RVT regressors
                    load(fullfile(PhysioDir, subj_name{s}, sprintf('%s_sess%d_scan%d_physio.mat', subj_name{s}, ss, sca)));
                    
                    HRV = physio.model.R(:, 1);
                    RVT = physio.model.R(:, 2);
                    
                    % Do OLS regression for HRV
                    [b_hrv,bint_hrv,r_hrv,rint_hrv,stats_hrv] = regress(HRV,X_reg);
                    % Do OLS regression for RVT
                    [b_rvt,bint_rvt,r_rvt,rint_rvt,stats_rvt] = regress(RVT,X_reg);
                    
                    physio_reg.b_hrv     = b_hrv;
                    physio_reg.b_rvt     = b_rvt;
                    physio_reg.bint_hrv  = bint_hrv;
                    physio_reg.bint_rvt  = bint_rvt;
                    physio_reg.r_hrv     = r_hrv;
                    physio_reg.r_rvt     = r_rvt;
                    physio_reg.rint_hrv  = rint_hrv;
                    physio_reg.rint_rvt  = rint_rvt;
                    physio_reg.stats_hrv = stats_hrv;
                    physio_reg.stats_rvt = stats_rvt;
                    
                    % represent instructions with a single regressor:
                    Xreg_uni = sum(X_reg, 2);
                    [b1_hrv,bint1_hrv,r1_hrv,rint1_hrv,stats1_hrv] = regress(HRV,Xreg_uni);
                    [b1_rvt,bint1_rvt,r1_rvt,rint1_rvt,stats1_rvt] = regress(RVT,Xreg_uni);
                    
                    physio_reg.b1_hrv     = b1_hrv;
                    physio_reg.bint1_hrv  = bint1_hrv;
                    physio_reg.r1_hrv     = r1_hrv;
                    physio_reg.rint1_hrv  = rint1_hrv;
                    physio_reg.stats1_hrv = stats1_hrv;
                    physio_reg.b1_rvt     = b1_rvt;
                    physio_reg.bint1_rvt  = bint1_rvt;
                    physio_reg.r1_rvt     = r1_rvt;
                    physio_reg.rint1_rvt  = rint1_rvt;
                    physio_reg.stats1_rvt = stats1_rvt;
                    
                    save(fullfile(PhysioDir, subj_name{s}, sprintf('%s_sess%d_scan%d_lin_reg.mat', subj_name{s}, ss, sca)), 'physio_reg', '-v7.3');
                    
                    phys_cell{ss, sca} = physio_reg;
                end % sca (scan)
            end % ss (sess)
        end % s (sn)
        varargout{1} = phys_cell;
    
    case 'GLM:mdtb:design_glm7' % GLM with each condition modelled as a regerssor. The instruction for each TASK is also modeled as a separate regressor
        %%% This case will calculate the design matrix with the instruction
        %%% period for each task separated and coming before the task.
        % Example: sc1_sc2_mdtb('GLM:mdtb:design_glm7', 'sn', [3]);
        sn            = returnSubjs; %% list of subjects
        experiment_num    = 1;           %% sc1 or sc2?
        ppmethod      = '';          %% 'stc' or ''? The default is set to ''
        deriv         = 1;           %% 'temp', 'temp_disp', or 'none'?
        glm           = 72;           %% the glm number      
        
        vararginoptions(varargin,{'sn', 'experiment_num', 'ppmethod', 'deriv', 'glm'});
        
        announceTime = 5;
                
        % load in task information
        C  = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc = getrow(C, C.StudyNum == experiment_num);
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% SPM and SPM_info files will be saved in glmDir
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                imDir  = 'imaging_data';
                prefix = 'r';
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                imDir  = 'imaging_data_stc';
                prefix = 'ra';
        end
        dircheck(glmDir)
        
        for s = sn
            fprintf('******************** creating SPM.mat file for %s ********************\n', subj_name{s});
            
            T = [];
            A = dload(fullfile(baseDir, experiment,'data', subj_name{s},sprintf('%s_%s.dat', experiment, subj_name{s})));
            A = getrow(A, A.runNum >= funcRunNum(1) & A.runNum <= funcRunNum(2));
            
            glmSubjDir = fullfile(glmDir, subj_name{s});
            dircheck(glmSubjDir);
            
            % starting to set up spm job parameters
            J.dir            = {glmSubjDir};
            J.timing.units   = 'secs';
            J.timing.RT      = 1.0;
            switch ppmethod
                case 'stc'
                    J.timing.fmri_t  = 48; %% there are 48 slices
                    J.timing.fmri_t0 = 24; %% set it to the middle slice
                case ''
                    J.timing.fmri_t  = 16; 
                    J.timing.fmri_t0 = 1; %% set it to the middle slice
            end
            
            % from Maedbh code: annoying but reorder behavioural runs slightly for 2
            % subjects...
            switch experiment
                case 'sc1'
                    if strcmp(subj_name{s},'s18')
                        %                 runTrue = [51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,60];
                        runTruestr = {'01','02','03','04','05','06','07','08','09','11','12','13','14','15','16','10'};
                    elseif strcmp(subj_name{s},'s21')
                        %                 runTrue = [51,52,53,54,55,56,57,58,59,60,61,63,64,65,66,62];
                        runTruestr = {'01','02','03','04','05','06','07','08','09', '10', '11','13','14','15','16','12'};
                    else
%                         runTrue = runB;
                        runTruestr = run{1};
                    end
                    runTrue = runB;
                case 'sc2'
                    runTrue    = runB;
                    runTruestr = run{2};
            end % for experiment 1 the order of the runs for two subjects is different
            
            for r = runLst
                P = getrow(A,A.runNum == runTrue(r)); %% contains the data (start time, end timem etc.) for behavioral runs without the data for instructions
                N = cell(1, numTRs - numDummys);   %% preallocating the variable that will contain the paths pointing to the image files
                for it = 1:(numTRs-numDummys)
                    N{it} = fullfile(baseDir, 'sc1', imDir, subj_name{s}, sprintf('%srun_%s.nii,%d', prefix, runTruestr{r}, it));
                end % scans (volume images)
                J.sess(r).scans = N; % path to scans
                
                %%% I intend to have separate regressors for each task
                %%% condition. Here, I first get the onset times for the
                %%% instruction that comes before a task condition.
                %%%%% for each task condition, first the instruction and
                %%%%% then the task itself
                ic0 = 2; %% task condition index (not including the first index which is the instruction when reading from Cc)
                ic  = 1; %% I call this the universal index :)
                for it = 1:length(P.taskNum) %it = 1:length(P.taskNum)-1 % "it" is task index (each task can have multiple conditions)
                    %%% P.taskNum is the number of tasks so it does not
                    %%% include the instruction. But it includes rest!
                    % Instructions first
                    ST  = find(strcmp(P.taskName,Cc.taskNames{ic0}));
                    instruct_onset = P.realStartTime(ST)-J.timing.RT*numDummys; %% get the instruction start time for the first task                   
%                     instruct_onset = P.realStartTime(it)-J.timing.RT*numDummys; %% get the instruction start time for the first task
                    J.sess(r).cond(ic).name     = Cc.condNames{1};
                    J.sess(r).cond(ic).onset    = instruct_onset(1); % correct start time for numDummys and announcetime included (not for instruct)
                    J.sess(r).cond(ic).duration = Cc.duration(1);  % duration of trials (+ fixation cross) we are modeling
                    J.sess(r).cond(ic).tmod     = 0;
                    J.sess(r).cond(ic).orth     = 0;
                    J.sess(r).cond(ic).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                    
                    S.SN    = s;
                    S.run   = r;
                    S.inst  = 1;
                    
                    S.instime         = instruct_onset;
                    S.taskName_after  = P.taskName(ST);
                    if ST>1
                        S.taskName_before = P.taskName(ST-1);
                    else
                        S.taskName_before = 'NaN';
                    end
                    S.instOrder = ST;
                    
                    S.task  = Cc.taskNum(1);
                    S.cond  = 0;
                    S.TN    = {Cc.condNames{1}};
                    S.sess  = sess(r);
                    S.deriv = 0;
                
                    T  = addstruct(T,S);
                    switch deriv
                        case 1
                            % adding the row representing the temporal derivative
                            S.SN    = s;
                            S.run   = r;
                            S.inst  = 1;
                            
                            S.instime         = instruct_onset;
                            S.taskName_after  = P.taskName(ST);
                            if ST>1
                                S.taskName_before = P.taskName(ST-1);
                            else
                                S.taskName_before = cellstr('NaN');
                            end
                            S.instOrder = ST;
                            
                            S.task  = Cc.taskNum(1);
                            S.cond  = 0;
                            S.TN    = {Cc.condNames{1}};
                            S.sess  = sess(r);
                            S.deriv = 1;

                            T  = addstruct(T,S);
                        case 2
                            % adding the rows for temporal and dispersion
                            % derivatives
                            S.SN    = s;
                            S.run   = r;
                            S.inst  = 1;
                            
                            S.instime         = instruct_onset;
                            S.taskName_after  = P.taskName(ST);
                            if ST>1
                                S.taskName_before = P.taskName(ST-1);
                            else
                                S.taskName_before = 'NaN';
                            end
                            S.instOrder = ST;
                            
                            S.task  = Cc.taskNum(1);
                            S.cond  = 0;
                            S.TN    = {Cc.condNames{1}};
                            S.sess  = sess(r);
                            S.deriv = 1;

                            T  = addstruct(T,S);
                            
                            S.SN    = s;
                            S.run   = r;
                            S.inst  = 1;
                            
                            S.instime         = instruct_onset;
                            S.taskName_after  = P.taskName(ST);
                            if ST>1
                                S.taskName_before = P.taskName(ST-1);
                            else
                                S.taskName_before = 'NaN';
                            end
                            S.instOrder = ST;
                            
                            S.task  = Cc.taskNum(1);
                            S.cond  = 0;
                            S.TN    = {Cc.condNames{1}};
                            S.sess  = sess(r);
                            S.deriv = 2;

                            T  = addstruct(T,S);
                    end

                    % Find Number of Conditions for this task (The conditions for this task have the same instruction)
                    numCond = length(find(Cc.taskNum == Cc.taskNum(ic0)));
                    
                    ic  = ic + 1; 
                    for cond=1:numCond 
                        
                        D  = dload(fullfile(baseDir, experiment,behavDir, subj_name{s},sprintf('%s_%s_%s.dat', experiment, subj_name{s}, Cc.taskNames{ic0})));
                        R  = getrow(D,D.runNum==runB(r)); % functional runs
                        ST = find(strcmp(P.taskName,Cc.taskNames{ic0}));

                        switch experiment % the onsets for sc1 and sc2 are determined differently!
                            case 'sc1' 
                                % determine trialType (ugly -- but no other way)
                                if isfield(R,'trialType')
                                    tt = (R.trialType==Cc.trialType(ic0));
                                else
                                    tt = Cc.trialType(ic0);
                                end
                                if strcmp(Cc.taskNames{ic0},'visualSearch')
                                    tt = (R.setSize==Cc.trialType(ic0));
                                elseif strcmp(Cc.taskNames{ic0},'nBack') || strcmp(Cc.taskNames{ic0},'nBackPic')
                                    tt = (R.respMade==Cc.trialType(ic0));
                                elseif strcmp(Cc.taskNames{ic0},'motorImagery') || strcmp(Cc.taskNames{ic0},'ToM'),
                                    tt = 1;
                                end
                                onset = [P.realStartTime(ST)+R.startTimeReal(tt)+announceTime-(J.timing.RT*numDummys)];
                            case 'sc2'
                                onset=[P.realStartTime(ST)+R.startTimeReal(R.condition==Cc.trialType(ic0))+announceTime-(J.timing.RT*numDummys)];
                        end % switch experiment

                        % loop through trial-types (ex. congruent or incongruent)
                        J.sess(r).cond(ic).name     = Cc.condNames{ic0};
                        J.sess(r).cond(ic).onset    = onset; % correct start time for numDummys and announcetime included (not for instruct)
                        J.sess(r).cond(ic).duration = Cc.duration(ic0);  % duration of trials (+ fixation cross) we are modeling
                        J.sess(r).cond(ic).tmod     = 0;
                        J.sess(r).cond(ic).orth     = 0;
                        J.sess(r).cond(ic).pmod     = struct('name', {}, 'param', {}, 'poly', {});

                        S.SN    = s;
                        S.run   = r;
                        S.inst  = 0; % instruction flag
                        
                        S.instime         = 0;
                        S.taskName_after  = 'none';
                        S.taskName_before = 'none';
                        
                        S.task  = Cc.taskNum(ic0);
                        S.cond  = Cc.condNum(ic0);
                        S.TN    = {Cc.condNames{ic0}};
                        S.sess  = sess(r);
                        S.deriv = 0;
                        
                        T  = addstruct(T,S); 
                        switch deriv
                            case 1 
                                % adding the row representing the temporal derivative
                                S.SN    = s;
                                S.run   = r;
                                S.inst  = 0;
                                
                                S.instime         = 0;
                                S.taskName_after  = 'none';
                                S.taskName_before = 'none';
                        
                                S.task  = Cc.taskNum(ic0);
                                S.cond  = Cc.condNum(ic0);
                                S.TN    = {Cc.condNames{ic0}};
                                S.sess  = sess(r);
                                S.deriv = 1;

                                T  = addstruct(T,S);
                            case 2 
                                % adding the rows for temporal and dispersion
                                % derivatives
                                S.SN    = s;
                                S.run   = r;
                                S.inst  = 0;
                                
                                S.instime         = 0;
                                S.taskName_after  = 'none';
                                S.taskName_before = 'none';
                                
                                S.task  = Cc.taskNum(ic0);
                                S.cond  = Cc.condNum(ic0);
                                S.TN    = {Cc.condNames{ic0}};
                                S.sess  = sess(r);
                                S.deriv = 1;

                                T  = addstruct(T,S);

                                S.SN    = s;
                                S.run   = r;
                                S.inst  = 1;
                                
                                S.instime         = 0;
                                S.taskName_after  = 'none';
                                S.taskName_before = 'none';
                                
                                S.task  = Cc.taskNum(ic0);
                                S.cond  = Cc.condNum(ic0);
                                S.TN    = {Cc.condNames{ic0}};
                                S.sess  = sess(r);
                                S.deriv = 2;

                                T  = addstruct(T,S);
                        end
                        
                        ic  = ic + 1;
                        ic0 = ic0 + 1;
                    end % cond
                end %taskNum
                J.sess(r).multi     = {''};
                J.sess(r).regress   = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf       = inf; 
            end % run
            J.fact = struct('name', {}, 'levels', {});
            switch deriv
                case 0
                    J.bases.hrf.derivs = [0 0];
                case 1
                    J.bases.hrf.derivs = [1 0];
                case 2
                    J.bases.hrf.derivs = [1 1];
            end                                    
            J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
            J.volt             = 1;
            J.global           = 'None';
            J.mask             = {fullfile(baseDir, 'sc1', imDir,subj_name{s},'rmask_noskull.nii,1')};
            J.mthresh          = 0.05;
            J.cvi_mask         = {fullfile(baseDir, 'sc1', imDir,subj_name{s},'rmask_gray.nii')};
            J.cvi              =  'fast';
            
            spm_rwls_run_fmri_spec(J);
            
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            fprintf('******************** glm_%d (SPM.mat) has been saved for %s ********************\n\n',glm, subj_name{s}); 
        end        
    case 'GLM:mdtb:design_glm8' % GLM with each task modeled as a 30 sec block regressor
        % Example: sc1_sc2_mdtb('GLM:mdtb:design_glm8', 'sn', [3]);
        sn            = returnSubjs; %% list of subjects
        experiment_num    = 1;           %% sc1 or sc2?
        ppmethod      = '';          %% 'stc' or ''? The default is set to ''
        deriv         = 1;           %% 0, 1, or 2 for no derivative, temporal, and temporal + dispersion?
        glm           = 8;          %% the glm number      
        
        vararginoptions(varargin,{'sn', 'experiment_num', 'ppmethod', 'deriv', 'glm'});
                
        % load in task information
        C     = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc    = getrow(C, C.StudyNum == experiment_num);
        Tasks = unique(Cc.taskNames,'rows','stable');                       % get the task names
        Tasks(strcmp(Tasks, 'Instruct') | strcmp(Tasks, 'Instruct2')) = []; % .dat file with all the info for the tasks does not have 'Instruct', so I'm eliminating it here!
        nTask      = unique(length(Tasks));                                 % how many tasks there are? for sc1: 18 (including rest) and sc2: 33 (including rest)

        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        announceTime = 5; % there is a 5 sec interval between instruction onset and task onset.
        
        %%% SPM and SPM_info files will be saved in glmDir
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                imDir  = 'imaging_data';
                prefix = 'r';
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                imDir  = 'imaging_data_stc';
                prefix = 'ra';
        end
        dircheck(glmDir)
        
        for s = sn 
            fprintf('******************** creating SPM.mat file for %s ********************\n', subj_name{s});
            
            S.SN = s;
            
            T = []; % T will be saved as SPM_info.mat
            A = dload(fullfile(baseDir, experiment, 'data', subj_name{s}, sprintf('%s_%s.dat', experiment, subj_name{s})));
            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
            
            glmSubjDir = fullfile(glmDir, subj_name{s});
            dircheck(glmSubjDir);
            
            % starting to set up spm job parameters
            J.dir            = {glmSubjDir};
            J.timing.units   = 'secs';
            J.timing.RT      = 1.0;
            switch ppmethod
                case 'stc'
                    J.timing.fmri_t  = 48; %% there are 48 slices
                    J.timing.fmri_t0 = 24; %% set it to the middle slice
                case ''
                    J.timing.fmri_t  = 16;
                    J.timing.fmri_t0 = 1; %% set it to the middle slice
            end
            
            % from Maedbh code: annoying but reorder behavioural runs slightly for 2
            % subjects...
            switch experiment
                case 'sc1'
                    if strcmp(subj_name{s},'s18')
                        %                 runTrue = [51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,60];
                        runTruestr = {'01','02','03','04','05','06','07','08','09','11','12','13','14','15','16','10'};
                    elseif strcmp(subj_name{s},'s21')
                        %                 runTrue = [51,52,53,54,55,56,57,58,59,60,61,63,64,65,66,62];
                        runTruestr = {'01','02','03','04','05','06','07','08','09', '10', '11','13','14','15','16','12'};
                    else
                        %                         runTrue = runB;
                        runTruestr = run{1};
                    end
                    runTrue = runB;
                case 'sc2'
                    runTrue    = runB;
                    runTruestr = run{2};
            end % for experiment 1 the order of the runs for two subjects is different
            
            % loop through runs
            for r = runLst
                S.run  = r;
                S.sess = sess(r);
                % get the task data for the run
                P = getrow(A,A.runNum == runTrue(r));
                
                % get the paths to the images for each run
                N = cell(1, numTRs - numDummys);   %% preallocating the variable that will contain the paths pointing to the image files
                for it = 1:(numTRs-numDummys)
                    N{it} = fullfile(baseDir, 'sc1', imDir, subj_name{s}, sprintf('%srun_%s.nii,%d', prefix, runTruestr{r}, it));
                end % scans (volume images)
                J.sess(r).scans = N; % path to scans
                
                % loop through tasks
                itt = 1; % this additional index will be used to fill in the fields for the J structure
                for it = 1:nTask 
                    % The order of tasks are different for each run, to
                    % have a common order for the tasks, I will be reading
                    % from the Cc file for all the runs and subjects
                    ST = find(strcmp(P.taskName,Tasks{it}));
                    for taskType = 1:2 % there are two taskTypes: instruction; not instructions
                        if taskType == 1 % instructions
                            % get the isntruction onset
                            instruct_onset = P.realStartTime(ST)- J.timing.RT*numDummys; %% get the instruction start time for the first task 
                            
                            % filling in the fields for SPM_info.mat
                            S.deriv     = 0;              % deriv is used to identify the derivative regressors
                            S.task      = 0;              % task Number
                            S.TN        = {'Instruct'};   % task name (TN)
                            S.inst      = 1;              % is it instruction (1) or not (0)?
                            S.instOrder = ST;             % instOrder is defined by the task that comes after the instruction
                            S.time      = instruct_onset; % instruction onset time
                            % Determine taskName_after and taskName_before
                            % this instruction
                            S.taskName_after  = P.taskName(ST);
                            if ST > 1
                                S.taskName_before = P.taskName(ST-1);
                            elseif ST == 1
                                S.taskName_before = {'NaN'};
                            end
                            T  = addstruct(T, S);
                            % the temporal or temporal + dispersion
                            % derivatives
                            switch deriv
                                case 1
                                    S.deriv = 1;
                                    T  = addstruct(T, S);
                                case 2
                                    S.deriv = 2;
                                    T  = addstruct(T, S);
                            end
                            
                            % filling in the fields for SPM.mat
                            J.sess(r).cond(itt).name     = 'Instruct';
                            J.sess(r).cond(itt).onset    = instruct_onset; % correct start time for numDummys and announcetime included (not for instruct)
                            J.sess(r).cond(itt).duration = 5;              % instructions last for 5 sec
                            J.sess(r).cond(itt).tmod     = 0;
                            J.sess(r).cond(itt).orth     = 0;
                            J.sess(r).cond(itt).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                            
                            itt = itt + 1;
                        elseif taskType == 2 % not instructions
                            % get the task onset (instruction onset + announceTime)
                            onset = P.realStartTime(ST) - J.timing.RT*numDummys +announceTime;
                            
                            % filling in the fields for SPM_info.mat
                            S.deriv     = 0;
                            S.task      = it;
                            S.TN        = {Tasks{it}};
                            S.inst      = 0;
                            S.instOrder = 0;
                            S.time      = onset;
                            S.taskName_after  = {'none'}; % taskName before and after are only defined for instructions
                            S.taskName_before = {'none'};
                            
                            T  = addstruct(T, S);
                            
                            switch deriv 
                                case 1
                                    S.deriv = 1;
                                    T  = addstruct(T, S);
                                case 0
                                    S.deriv = 2;
                                    T  = addstruct(T, S);
                            end
                            
                            % filling in the fields for SPM.mat
                            J.sess(r).cond(itt).name     = Tasks{it};
                            J.sess(r).cond(itt).onset    = onset;
                            J.sess(r).cond(itt).duration = 30;             % each task lasts for 30 sec
                            J.sess(r).cond(itt).tmod     = 0;
                            J.sess(r).cond(itt).orth     = 0;
                            J.sess(r).cond(itt).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                            
                            itt = itt + 1;
                        end % if it's instructions or not?
                    end % taskType (looping over instruction and non-instructions)
                end % it (tasks)
                J.sess(r).multi = {''};
                J.sess(r).regress = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf = inf;                                        % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
            end % r (runs)
            J.fact = struct('name', {}, 'levels', {});
            switch deriv
                case 0 % no derivatives
                    J.bases.hrf.derivs = [0 0];
                case 1 % temporal derivative
                    J.bases.hrf.derivs = [1 0];
                case 2 % temporal and dispersion derivative
                    J.bases.hrf.derivs = [1 1];
            end % including temporal and dispersion derivatives or not?                                    
            J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
            J.volt             = 1;
            J.global           = 'None';
            J.mask             = {fullfile(baseDir, 'sc1', imDir,subj_name{s},'rmask_noskull.nii,1')};
            J.mthresh          = 0.05;
            J.cvi_mask         = {fullfile(baseDir, 'sc1', imDir,subj_name{s},'rmask_gray.nii')};
            J.cvi              =  'fast';
            
            spm_rwls_run_fmri_spec(J);
            
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            fprintf('******************** glm_%d (SPM.mat) has been saved for %s ********************\n\n',glm, subj_name{s}); 
        end % s (sn)
    case 'GLM:mdtb:run_glm'
        %%% This is time consuming, but I run it to get the whitening
        %%% filter that will then be applied to the time series.
        % Example: sc1_sc2_mdtb('GLM:mdtb:run_glm', 'sn', [2])
        sn         = returnSubjs;   %% list of subjects
        glm        = 8;            %% The glm number :)
        experiment_num = 1;
        ppmethod   = '';            %% was the preprocessing done with stc included? Input 'stc' for pp with slice timing and 'no_stc' for pp without it
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'ppmethod'})
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'

        %%% setting the directory paths I need
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
        end
        
        for s = sn
            
            fprintf('******************** estimating glm%d parameters for %s ********************\n', glm, subj_name{s});
            
            glmSubjDir = fullfile(glmDir, subj_name{s});
            load(fullfile(glmSubjDir,'SPM.mat'));
            
            SPM.swd = glmSubjDir;
            spm_rwls_spm(SPM);
            
            fprintf('******************** glm%d parameters estimated for %s ********************\n\n', glm, subj_name{s});
        end %sn        
    case 'GLM:mdtb:contrast'
        %%% Calculating contrast images.
        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM).
        % This case is written so that it works with both GLM 7 and GLM 8.
        % Reminder: GLM 7 was written with each condition as a separate
        % regressor and a regressor for each of the instructions. GLM 8 was
        % written with each task modeled as a 30 sec block and instructions
        % modeled as a separate regressor.
        % Example1: sc1_sc2_mdtb('GLM:mdtb:contrast', 'sn', [2], 'glm', 9, 'which', 'task')
        % Example2: sc1_sc2_mdtb('GLM:mdtb:contrast', 'sn', [3], 'glm', 72, 'which', 'cond')
        
        sn         = returnSubjs;        %% list of subjects
        glm        = 7;           %% The glm number :)
        experiment_num = 1;
        ppmethod   = '';          %% was the preprocessing done with stc included? Input 'stc' for pp with slice timing and 'no_stc' for pp without it
        con_vs     = 'average_1'; %% set it to 'rest' or 'average' (depending on the contrast you want)
        which      = 'task';      %% it can be set to either cond or task. set it to 'task for GLM_8 and 'cond' for GLM_7
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'ppmethod', 'con_vs', 'which'})
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% setting directory paths I need
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
        end
        
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            T    = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            
            % t contrast for tasks
            ucondition = unique(T.(which));
            idx = 1;
            for tt = 1:length(ucondition) % 0 is "instruct" regressor
                switch con_vs
                    case 'rest_task' % contrast against rest
                        con                                  = zeros(1,size(SPM.xX.X,2));
                        con(:,logical((T.(which) == ucondition(tt)) .* (T.deriv == 0))) = 1;
                        tmp = zeros(1, size(SPM.xX.X, 2));
                        tmp(:, strcmp(T.TN, 'rest')) = 1;
                        sum_rest = sum(tmp);
                        tmp = tmp./sum_rest;
                        con                                  = con/abs(sum(con));
                        con(:, tmp~=0) = -1/sum_rest;
                    case 'average_1' % contrast against the average of the other tasks including the instructions
                        % New: Eva, Oct 2nd
                        con        = zeros(1,size(SPM.xX.X, 2));
                        con(1,logical((T.(which) == ucondition(tt)) .* (T.deriv == 0))) = 1./sum(logical((T.(which) == ucondition(tt)) .* (T.deriv == 0)));
                        con(1,logical((T.(which) ~= ucondition(tt)) .* (T.deriv == 0))) = -1./sum(logical((T.(which) ~= ucondition(tt)) .* (T.deriv == 0)));
                    case 'average_2' % contrast against the average of the other tasks not including the instructions
                        con        = zeros(1,size(SPM.xX.X, 2));
                        % TO TRY below - no instructions as contrasts
                        con(1,logical((T.(which) == ucondition(tt)) .* (T.deriv == 0))) = 1./sum(logical((T.(which) == ucondition(tt)) .* (T.deriv == 0)));
                        con(1,logical((T.(which) ~= ucondition(tt)) .* (T.deriv == 0) .* (T.inst == 0))) = -1./sum(logical((T.(which) ~= ucondition(tt)) .* (T.deriv == 0) .* (T.inst == 0)));
                    case 'rest'
                        con                                  = zeros(1,size(SPM.xX.X,2));
                        con(:,logical((T.(which) == ucondition(tt)) .* (T.deriv == 0))) = 1;
                        con                                  = con/abs(sum(con));
                end
                
                name = sprintf('%s-%s',char(unique(T.TN(T.(which) == ucondition(tt)))), con_vs);
                
                SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
                idx=idx+1;
            end % tt (conditions)
            SPM = spm_contrasts(SPM,1:length(SPM.xCon));
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glmDir, subj_name{s}, 'SPM_light.mat'), 'SPM');

            % rename contrast images and spmT images
             conName = {'con','spmT'};
            for i = 1:length(SPM.xCon)
                for n = 1:numel(conName)
                    oldName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%2.4d.nii',conName{n},i));
                    newName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                    movefile(oldName{i}, newName{i});
                end % conditions (n, conName: con and spmT)
            end % i (contrasts)
        end % sn 
    case 'GLM:mdtb:contrast_task'
        %%% Calculating contrast images for tasks.
        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM). The contrast for the task is created as
        % the average of the beta values for the conditions of a task
        % Example: sc1_sc2_mdtb('GLM:mdtb:contrast_task', 'sn', [3])
        
        sn         = returnSubjs;        %% list of subjects
        glm        = 72;           %% The glm number :)
        experiment_num = 1;
        ppmethod   = '';          %% was the preprocessing done with stc included? Input 'stc' for pp with slice timing and 'no_stc' for pp without it
        con_vs     = 'average_1'; %% set it to 'rest' or 'average_1' or 'average_2' (depending on the contrast you want)
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'ppmethod', 'con_vs'})
        
        % gt the task info
        C   = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc  = getrow(C,C.StudyNum == experiment_num);
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% setting directory paths I need
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
        end
        
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            T    = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            T_rd = getrow(T, T.deriv == 0);
            
            % t contrast for tasks
            utask = unique(T_rd.task);
            
            idx = 1;
            for tt = 1:length(utask) % 0 is "instruct" regressor
                switch con_vs
                    case 'rest' % contrast against rest
                        con                                = zeros(1,size(SPM.xX.X,2));
%                         con(:,T_rd.task == utask(tt))      = 1;
                        con(:,logical((T.task == utask(tt)) .* (T.deriv == 0))) = 1;
                        con                                = con/abs(sum(con));
                    case 'average' % contrast against the average of all the tasks or all the other tasks???
                        con                                = zeros(1,size(SPM.xX.X,2));
                        con(:,T_rd.task == task(tt))       = 1;
                        con                                = con/abs(sum(con));
                        con                                = bsxfun(@minus, con, 1/length(T.cond));
                    case 'average_1' % contrast vs the average of all the tasks
                        con        = zeros(1,size(SPM.xX.X, 2));
                        con(1,logical((T.task == utask(tt)) .* (T.deriv == 0))) = 1./sum(logical((T.task == utask(tt)) .* (T.deriv == 0)));
                        con(1,logical((T.task ~= utask(tt)) .* (T.deriv == 0))) = -1./sum(logical((T.task ~= utask(tt)) .* (T.deriv == 0)));
                    case 'average_2' % contrast vs average of all the tasks except for instructions
                        con        = zeros(1,size(SPM.xX.X, 2));
                        % TO TRY below - no instructions as contrasts
                        con(1,logical((T.task == utask(tt)) .* (T.deriv == 0))) = 1./sum(logical((T.task == utask(tt)) .* (T.deriv == 0)));
                        con(1,logical((T.task ~= utask(tt)) .* (T.deriv == 0) .* (T.inst == 0))) = -1./sum(logical((T.task ~= utask(tt)) .* (T.deriv == 0) .* (T.inst == 0)));
                end
                % fix the name!!!!!!!!!!!!!!!!
                name = sprintf('%s-%s_taskCon',char(unique(Cc.taskNames(Cc.taskNum == utask(tt)))), con_vs);
                
                SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
                idx=idx+1;
            end % tt (tasks)
            SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glmDir, subj_name{s}, 'SPM_light.mat'), 'SPM');

            % rename contrast images and spmT images
            conName = {'con','spmT'};
            for i = 1:length(SPM.xCon)
                for n = 1:numel(conName)
                    oldName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%2.4d.nii',conName{n},i));
                    newName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                    movefile(oldName{i}, newName{i});
                end % conditions (n, conName: con and spmT)
            end % i (contrasts)
        end % sn        
        
    case 'SURF:mdtb:map_con'
        % projects individual contrast map volume files for the conditions
        % to the workbench surface.
        % Example: sc1_sc2_mdtb('SURF:mdtb:map_con', 'sn', [3])
    
        sn         = returnSubjs;        %% list of subjects
        ppmethod   = '';          %% with or without stc
        atlas_res  = 32;          %% set it to 32 or 164
        experiment_num = 1;           %% enter 1 for sc1 and 2 for sc2
        glm        = 8;           %% glm number
        con_vs     = 'average_1'; %% set it to 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment_num', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                wbDir  = 'surfaceWB';
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                wbDir  = 'surfaceWB_stc';
        end
        glmSurfDir = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm));
        dircheck(glmSurfDir);
        
        for s = sn
            fprintf('******************** start mapping contrasts to surface for %s ********************\n', subj_name{s});
            subjSurfDir = fullfile(baseDir, 'sc1', wbDir, 'data', subj_name{s});
            dircheck(fullfile(glmSurfDir, subj_name{s})); %% directory to save the contrast maps
            
            T  = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            Td = getrow(T, T.deriv == 0); %% get the non-derivative regressor names
            conNames = unique(Td.TN);
            for h = 1:2 % two hemispheres
                white   = fullfile(subjSurfDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                pial    = fullfile(subjSurfDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                C1      = gifti(white);
                C2      = gifti(pial);
                for ic = 1:length(conNames)
                    
                    conMapName      = sprintf('con_%s-%s', conNames{ic}, con_vs);
                    images{1}       = fullfile(glmDir, subj_name{s}, sprintf('%s.nii', conMapName));
                    column_name{1}  = sprintf('con_%s-%s.nii', conNames{ic}, con_vs);
                    outfile         = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.con_%s-%s.func.gii', ...
                        subj_name{s}, hemI{h}, conNames{ic}, con_vs));
                    G               = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names', column_name, ...
                        'anatomicalStruct',hemName{h});
                    save(G, outfile);
                    fprintf('******************** mapping to surface for %s hemi %s contrast %s done********************\n', subj_name{s}, ...
                        hemI{h}, conMapName);
                end % ic (condition/contrast)   
            end % hemi
        end % sn    
    case 'SURF:mdtb:map_con_task'
        % projects individual contrast map volume files for tasks to
        % WorkBench surface
        % Example: sc1_sc2_mdtb('SURF:mdtb:map_con_task', 'sn', [3])
    
        sn         = returnSubjs;        %% list of subjects
        ppmethod   = '';          %% with or without stc
        atlas_res  = 32;          %% set it to 32 or 164
        experiment_num = 1;           %% enter 1 for sc1 and 2 for sc2
        glm        = 72;           %% glm number
        con_vs     = 'average_1'; %% set it to 'rest' or 'average_1' or 'average_2'
        
        % gt the task info
        C   = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc  = getrow(C,C.StudyNum == experiment_num);
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment_num', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                wbDir  = 'surfaceWB';
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                wbDir  = 'surfaceWB_stc';
        end
        glmSurfDir = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm));
        dircheck(glmSurfDir);
        
        for s = sn
            fprintf('******************** start mapping contrasts to surface for %s ********************\n', subj_name{s});
            subjSurfDir = fullfile(baseDir, 'sc1', wbDir, 'data', subj_name{s});
            dircheck(fullfile(glmSurfDir, subj_name{s})); %% directory to save the contrast maps
            
            % get the task names
            conNames = unique(Cc.taskNames);
            
            for h = 1:2 % two hemispheres
                white   = fullfile(subjSurfDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                pial    = fullfile(subjSurfDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                C1      = gifti(white);
                C2      = gifti(pial);
                for ic = 1:length(conNames)
%                     if ~ strcmp(conNames{ic}, 'rest')
                        conMapName      = sprintf('con_%s-%s_taskCon', conNames{ic}, con_vs);
                        images{1}       = fullfile(glmDir, subj_name{s}, sprintf('%s.nii', conMapName));
                        column_name{1}  = sprintf('con_%s-%s_taskCon.nii', conNames{ic}, con_vs);
                        
                        outfile    = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.con_%s-%s_taskCon.func.gii', subj_name{s}, hemI{h}, conNames{ic}, con_vs));
                        G          = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names', column_name, ...
                            'anatomicalStruct',hemName{h});
                        save(G, outfile);
                        fprintf('******************** mapping to surface for %s hemi %s contrast %s done********************\n', subj_name{s}, ...
                            hemI{h}, conMapName);
%                     end
                    
                end % ic (condition/contrast)   
            end % hemi
        end % sn
    case 'SURF:mdtb:groupmap_con'
        % creates group average contrast maps for task contrasts
        % Example: sc1_sc2_mdtb('SURF:mdtb:groupmap_con', 'sn', [3])
    
        sn         = returnSubjs;   %% list of subjects
        ppmethod   = '';     %% with or without stc
        atlas_res  = 32;     %% set it to 32 or 164
        experiment_num = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 8;      %% glm number
        replaceNaN = 1;      %% replacing NaNs
        con_vs     = 'average_1'; %% contrast was calculated against 'rest' or 'average'        
        smooth     = 1;      %% add smoothing
        kernel     = 1;      %% for smoothing
        which      = 'task'; %% 'task' for glm8 and 'cond' for glm7
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment_num', 'glm', 'replaceNaN', 'con_vs', 'smooth', 'kernel', 'which'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        switch which
            case 'task' % task for glm8
                conNames = unique(Cc.taskNames);
            case 'cond' % condition for glm7
                conNames = unique(Cc.condNames);
        end %% do you want the group maps for tasks or conditions
        
        experiment = sprintf('sc%d', experiment_num);
        
        switch ppmethod
            case ''
                wbDir  = 'surfaceWB';
            case 'stc'
                wbDir  = 'surfaceWB_stc';
        end
        
        % go to the directory where fs_LR atlas is.
        groupSurfDir     = fullfile(baseDir, experiment, wbDir, 'data', sprintf('group%dk', atlas_res));
        groupSurfDir_glm = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), sprintf('group%dk', atlas_res));
        dircheck(groupSurfDir_glm);
        cd(groupSurfDir);
        
        for h = 1:2 % two hemispheres
            for cc = 1:length(conNames)
                for s = 1:length(sn)
                    %%% make the group metric file for each contrasts
%                     infilenames{s}   = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                    infilenames{s}   = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_%s-%s.func.gii', subj_name{sn(s)}, hemI{h}, conNames{cc}, con_vs));
                    
                    if smooth
                        surfFile    = fullfile(groupSurfDir,sprintf('fs_LR.32k.%s.inflated.surf.gii',hemI{h}));
                        surf_smooth(infilenames{s},'surf',surfFile,'kernel',kernel); % smooth outfilenames - it will prefix an 's'
                    end
%                     s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                    s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s.func.gii', subj_name{sn(s)}, hemI{h}, conNames{cc}, con_vs));
                
                end % sn
%                 outfilenames    = fullfile(groupSurfDir_glm,sprintf('%s.con_%s-%s_taskCon.func.gii',hemI{h},taskNames{cc}, con_vs));
%                 summaryname     = fullfile(groupSurfDir_glm,sprintf('%s.group.con_%s-%s_taskCon.func.gii',hemI{h},taskNames{cc}, con_vs));
                
                outfilenames    = fullfile(groupSurfDir_glm,sprintf('%s.con_%s-%s.func.gii',hemI{h},conNames{cc}, con_vs));
                summaryname     = fullfile(groupSurfDir_glm,sprintf('%s.group.con_%s-%s.func.gii',hemI{h},conNames{cc}, con_vs));
                
                surf_groupGiftis(infilenames, 'outfilenames', {outfilenames}, 'groupsummary', summaryname, 'replaceNaNs', replaceNaN);
                if smooth % also save the smoothed versions
%                     s_outfilenames    = fullfile(groupSurfDir_glm,sprintf('s%s.con_%s-%s_taskCon.func.gii', hemI{h},taskNames{cc}, con_vs));
                    s_outfilenames    = fullfile(groupSurfDir_glm,sprintf('s%s.con_%s-%s.func.gii', hemI{h},conNames{cc}, con_vs));
%                     s_summaryname     = fullfile(groupSurfDir_glm,sprintf('s%s.group.con_%s-%s_taskCon.func.gii', hemI{h},taskNames{cc}, con_vs));
                    s_summaryname     = fullfile(groupSurfDir_glm,sprintf('s%s.group.con_%s-%s.func.gii', hemI{h},conNames{cc}, con_vs));
                    surf_groupGiftis(s_infilenames, 'outfilenames', {s_outfilenames}, 'groupsummary', s_summaryname, 'replaceNaNs', replaceNaN);
                end
                
                % delet the smoothed contrasts for each subject
                for s = 1:length(sn)
                    delete(fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, conNames{cc}, con_vs)));
                end
                
                fprintf('******************** group average contrast for %s vs %s is created! ********************\n\n', conNames{cc}, con_vs);
            end % contrasts(ic)
        end % hemi(h)
    case 'SURF:mdtb:groupmap_con_task'
        % creates group average contrast maps for task contrasts
        % Example: sc1_sc2_mdtb('SURF:mdtb:groupmap_con_task', 'sn', [3])
    
        sn         = returnSubjs;   %% list of subjects
        ppmethod   = '';     %% with or without stc
        atlas_res  = 32;     %% set it to 32 or 164
        experiment_num = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 8;      %% glm number
        replaceNaN = 1;      %% replacing NaNs
        con_vs     = 'rest'; %% contrast was calculated against 'rest' or 'average'        
        smooth     = 1;      %% add smoothing
        kernel     = 1;      %% for smoothing
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment_num', 'glm', 'replaceNaN', 'con_vs', 'smooth', 'kernel'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        taskNames = unique(Cc.taskNames);
        
        experiment = sprintf('sc%d', experiment_num);
        
        switch ppmethod
            case ''
                wbDir  = 'surfaceWB';
            case 'stc'
                wbDir  = 'surfaceWB_stc';
        end
        
        % go to the directory where fs_LR atlas is.
        groupSurfDir     = fullfile(baseDir, experiment, wbDir, 'data', sprintf('group%dk', atlas_res));
        groupSurfDir_glm = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), sprintf('group%dk', atlas_res));
        dircheck(groupSurfDir_glm);
        cd(groupSurfDir);
        
        for h = 1:2 % two hemispheres
            for cc = 1:length(taskNames)
                for s = 1:length(sn)
                    %%% make the group metric file for each contrasts
%                     infilenames{s}   = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                    infilenames{s}   = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_%s-%s.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                    
                    if smooth
                        surfFile    = fullfile(groupSurfDir,sprintf('fs_LR.32k.%s.inflated.surf.gii',hemI{h}));
                        surf_smooth(infilenames{s},'surf',surfFile,'kernel',kernel); % smooth outfilenames - it will prefix an 's'
                    end
%                     s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                    s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                
                end % sn
%                 outfilenames    = fullfile(groupSurfDir_glm,sprintf('%s.con_%s-%s_taskCon.func.gii',hemI{h},taskNames{cc}, con_vs));
%                 summaryname     = fullfile(groupSurfDir_glm,sprintf('%s.group.con_%s-%s_taskCon.func.gii',hemI{h},taskNames{cc}, con_vs));
                
                outfilenames    = fullfile(groupSurfDir_glm,sprintf('%s.con_%s-%s.func.gii',hemI{h},taskNames{cc}, con_vs));
                summaryname     = fullfile(groupSurfDir_glm,sprintf('%s.group.con_%s-%s.func.gii',hemI{h},taskNames{cc}, con_vs));
                
                surf_groupGiftis(infilenames, 'outfilenames', {outfilenames}, 'groupsummary', summaryname, 'replaceNaNs', replaceNaN);
                if smooth % also save the smoothed versions
%                     s_outfilenames    = fullfile(groupSurfDir_glm,sprintf('s%s.con_%s-%s_taskCon.func.gii', hemI{h},taskNames{cc}, con_vs));
                    s_outfilenames    = fullfile(groupSurfDir_glm,sprintf('s%s.con_%s-%s.func.gii', hemI{h},taskNames{cc}, con_vs));
%                     s_summaryname     = fullfile(groupSurfDir_glm,sprintf('s%s.group.con_%s-%s_taskCon.func.gii', hemI{h},taskNames{cc}, con_vs));
                    s_summaryname     = fullfile(groupSurfDir_glm,sprintf('s%s.group.con_%s-%s.func.gii', hemI{h},taskNames{cc}, con_vs));
                    surf_groupGiftis(s_infilenames, 'outfilenames', {s_outfilenames}, 'groupsummary', s_summaryname, 'replaceNaNs', replaceNaN);
                end
                
                % delet the smoothed contrasts for each subject
                for s = 1:length(sn)
                    delete(fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs)));
                end
                
                fprintf('******************** group average contrast for %s vs %s is created! ********************\n\n', taskNames{cc}, con_vs);
            end % contrasts(ic)
        end % hemi(h)
        
    case 'SUIT:mdtb:suit_parcel2native'
        % maps each atlas of the suit into individual space
        % Example: sc1_sc2_mdtb('SUIT:mdtb:suit_parcel2native', 'sn', [3])
        
        sn         = returnSubjs;
        parcelType = 'Buckner_7';                         %% set it to Buckner_7 or Buckner_17
        parcelDir  = fullfile(suitToolDir, 'atlasesSUIT'); %% directory where the nifti image is stored
        experiment_num = 1;
        
        vararginoptions(varargin, {'sn', 'parcelType', 'parcelDir', 'experiment_num', 'ppmethod'});
        
        experiment = sprintf('sc%d', experiment_num);
                
        spm('defaults','fmri');
        spm_jobman('initcfg');

        for s = sn
            J.Affine       = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'Affine_c_anatomical_seg1.mat'); %% the affine transformation already calculated using suit_normalize_dartel
            J.flowfield    = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'u_a_c_anatomical_seg1.nii');    %% the flow field already calculated
            J.resample{1}  = fullfile(parcelDir, sprintf('%sNetworks.nii', parcelType));                                          %% Buckner parcellation is suit space
            J.ref          = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'maskbrainSUITGrey.nii');        %% the reference image 
            
            suit_reslice_dartel_inv(J); %% creates sprintf('wi%sNetworks.nii', parcelType)
            
            fprintf('******************** %s transformation to native space for %s done! ********************\n\n', parcelType, subj_name{s}); 
        end % sn
    case 'SUIT:mdtb:reslice'
        % this case is used to reslice volumes into suit space
        % before you run all the cases for the group map, you have to run
        % this case to map all the contrast maps to suit space.
        % Example: sc1_sc2_mdtb('SUIT:mdtb:reslice', 'sn', [3])
        
        sn         = returnSubjs;                   %% list of subjects
        ppmethod   = '';                     %% with or without stc
        experiment_num = 1;                      %% enter 1 for sc1 and 2 for sc2
        type       = 'con';                  %% enter the image you want to reslice to suit space
        glm        = 7;                      %% glm number
        mask       = 'cereb_prob_corr_grey'; %% the cerebellar mask to be used:'cereb_prob_corr_grey' or 'cereb_prob_corr' or 'dentate_mask'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'experiment_num', 'glm', 'type', 'mask'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        switch ppmethod
            case ''
                glmDir        = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                imageDir      = fullfile(baseDir, experiment, 'imaging_data');
                glmSuitDir    = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
            case 'stc'
                glmDir        = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                imageDir      = fullfile(baseDir, experiment, 'imaging_data_stc');
                glmSuitDir    = fullfile(baseDir, experiment, suitDir, sprintf('glm%d_stc', glm));
        end
        dircheck(glmSuitDir);
        
        for s = sn
            dircheck(fullfile(glmSuitDir, subj_name{s}));
            glmSubjDir = fullfile(glmDir, subj_name{s});
            outDir     = fullfile(glmSuitDir, subj_name{s});
            
            job.subj.affineTr  = {fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'Affine_c_anatomical_seg1.mat')};
            job.subj.flowfield = {fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s},'u_a_c_anatomical_seg1.nii')};
            job.subj.mask      = {fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, sprintf('%s.nii', mask))};
            job.vox            = [2 2 2];
            
            fprintf('******************** using dartel to reslice %s images to suit space for %s ********************\n', type, subj_name{s});
            
            switch type
                case 'beta'
                    images = 'beta_';
                    source = dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                    
                    job.subj.resample = {source.name};
                    suit_reslice_dartel(job);
                case 'con'  
                    % reslice all the contrasts (vs rest and vs average)
                    images = 'con';
                    source = dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                    
                    job.subj.resample = {source.name};
                    suit_reslice_dartel(job);
                case 'spmT'
                    images = 'spmT_';
                    source = dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                    
                    job.subj.resample = {source.name};
                    suit_reslice_dartel(job);
                case 'ResMS'
                    images = 'ResMS';
                    source = dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                    
                    job.subj.resample = {source.name};
                    suit_reslice_dartel(job);
                case 'cerebellarGrey'
                    source = dir(fullfile(baseDir, experiment, suitDir,'anatomicals',subj_name{s},'c1anatomical.nii')); % image to be resliced
                    cd(fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}));
                    
                    job.subj.resample = {source.name};
                    suit_reslice_dartel(job);
                case 'time_series'
                    dircheck(fullfile(imageDir, subj_name{s}));
                    imageSubjDir = fullfile(imageDir, subj_name{s});
                    outDir       = fullfile(baseDir, experiment, suitDir, 'imaging_data_resliced', subj_name{s});
                    dircheck(outDir);
                    switch ppmethod
                        case ''
                            images    = 'rrun';
                        case 'stc'
                            images    = 'rarun';
                    end % switch ppmethod
                    
                    % there are 598 time points (images) that need to be
                    % resliced to the suit space.
                    for rr = 1:length(runLst)
                        for itp = 1:numTRs - numDummys
                            filenames{itp} = fullfile(imageSubjDir, sprintf('%s_%0.2d.nii,%d', images, rr, itp));
                        end % itp (time points)
                        job.subj.resample = filenames';
                        suit_reslice_dartel(job);
                    end % rr (runs)
            end % switch type
            
            if ~strcmp(type,'cerebellarGrey')
                if strcmp(type, 'time_series')
                    source = fullfile(imageSubjDir, '*wd*');
                    
                else
                    source = fullfile(glmSubjDir,'*wd*');
                end
                dircheck(fullfile(outDir));
                
                destination = fullfile(baseDir, experiment, suitDir, sprintf('glm%d',glm), subj_name{s});
                movefile(source, destination);
            end
            fprintf('******************** %s for glm %d have been resliced into suit space for %s ********************\n\n', type, glm, subj_name{s});
        end % sn
%         % plotting the flatmap
%         V = spm_vol(fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), sprintf('wd%s_Instruct-rest.nii', type)));
%         D = suit_map2surf(V,'stats','nanmean');
%         suit_plotflatmap(D, 'cmap', colormap(jet(256)), 'cscale', [-3, 3]);
%         caxis([-3, 3]);
%         colorbar;
    case 'SUIT:mdtb:groupmap_con_cond'
        % creates group average for the condition contrast maps.
        % you need to reslice all the images to suit space before running
        % this case
        % Example: sc1_sc2_mdtb('SUIT:mdtb:groupmap_con_cond', 'sn', [3]);
        
        sn         = returnSubjs;                   %% list of subjects
        ppmethod   = '';                     %% with or without stc
        experiment_num = 1;                      %% enter 1 for sc1 and 2 for sc2
        type       = 'con';                  %% enter the image you want to reslice to suit space
        glm        = 7;                      %% glm number
        con_vs     = 'rest';                 %% is the contrast calculated vs 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'experiment_num', 'glm', 'type', 'mask'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        conNames = unique(Cc.condNames);
        
        experiment = sprintf('sc%d', experiment_num);
        
        switch ppmethod
            case ''
                glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
                glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), 'group');
            case 'stc'
                glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d_stc', glm));
                glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d_stc', glm), 'group');
        end
        dircheck(glmSuitDir);
        dircheck(glmSuitGroupDir);
        
        % preallocating!
        maps = []; %% the structure that will have alll?! the info I will need (hopefully :))
        for cc = 1:length(conNames)
            for s = 1:length(sn)                
                infilename{s} = fullfile(glmSuitDir, subj_name{sn(s)}, sprintf('wd%s_%s-%s.nii', type, conNames{cc}, con_vs));
                
                % load in individual cerebellar maps
                V_tmp = spm_vol(infilename{s});
                X     = spm_read_vols(V_tmp);
                X_vec = X(:);
                
                maps_tmp.conNames           = cellstr(conNames{cc});
                maps_tmp.conBaseName        = cellstr(con_vs);
                maps_tmp.subjectMaps{1}     = X;
                maps_tmp.subjectMapsVec     = X_vec';
                maps_tmp.subjectMapNames{1} = infilename{s};
                maps_tmp.glm                = 7;
                maps_tmp.sn                 = sn(s);
                maps_tmp.vol                = V_tmp;
                
                maps  = addstruct(maps, maps_tmp);
            end % sn
            outfilename = fullfile(glmSuitGroupDir, sprintf('group_%s_%s-%s.nii', type, conNames{cc}, con_vs));
            opt.dmtx    = 1;
            % calculate the average across subjects
            spm_imcalc(infilename', outfilename, 'nanmean(X)', opt);
            
            % saving the gifti file for group map
            V = spm_vol(fullfile(glmSuitGroupDir, sprintf('group_%s_%s-%s.nii', type, conNames{cc}, con_vs)));
            D = suit_map2surf(V,'stats','nanmean');
            
            G = surf_makeFuncGifti(D, 'anatomicalStruct', 'Cerebellum', 'columnNames', {sprintf('group_%s_%s-%s', type, conNames{cc}, con_vs)});
            
            save(G, fullfile(glmSuitGroupDir, sprintf('Cereb.group.con_%s-%s.func.gii', conNames{cc}, con_vs)));
            save(fullfile(glmSuitGroupDir, sprintf('indMaps_%s_%s-vs-%s.mat', type, conNames{cc}, con_vs)), 'maps', '-v7.3');
            fprintf('******************** %s group average for %s vs %s is created! ********************\n\n', type, conNames{cc}, con_vs);
            
        end % contrasts (cc)
        save(fullfile(glmSuitGroupDir, sprintf('indMaps_%s-vs-%s.mat', type, con_vs)), 'maps', '-v7.3');
    case 'SUIT:mdtb:groupmap_con_task'
        % creates group map for task contrasts
        % Example: sc1_sc2_mdtb('SUIT:mdtb:groupmap_con_task', 'sn', [3]);
        
        sn         = returnSubjs;                   %% list of subjects
        ppmethod   = '';                     %% with or without stc
        experiment_num = 1;                      %% enter 1 for sc1 and 2 for sc2
        type       = 'con';                  %% enter the image you want to reslice to suit space
        glm        = 7;                      %% glm number
        con_vs     = 'rest';                 %% is the contrast calculated vs 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'experiment_num', 'glm', 'type', 'mask'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        taskNames = unique(Cc.taskNames);
        
        experiment = sprintf('sc%d', experiment_num);
        
        switch ppmethod
            case ''
                glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
                glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), 'group');
            case 'stc'
                glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d_stc', glm));
                glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d_stc', glm), 'group');
        end
        dircheck(glmSuitDir);
        dircheck(glmSuitGroupDir);
        
        % preallocating!
        maps = []; %% the structure that will have alll?! the info I will need (hopefully :))
        for cc = 1:length(taskNames)
            for s = 1:length(sn)                
                infilename{s} = fullfile(glmSuitDir, subj_name{sn(s)}, sprintf('wd%s_%s-%s_taskCon.nii', type, taskNames{cc}, con_vs));
                
                % load in individual cerebellar maps
                V_tmp = spm_vol(infilename{s});
                X     = spm_read_vols(V_tmp);
                X_vec = X(:);
                
                maps_tmp.conNames           = cellstr(taskNames{cc});
                maps_tmp.conBaseName        = cellstr(con_vs);
                maps_tmp.subjectMaps{1}     = X;
                maps_tmp.subjectMapsVec     = X_vec';
                maps_tmp.subjectMapNames{1} = infilename{s};
                maps_tmp.glm                = 7;
                maps_tmp.sn                 = sn(s);
                maps_tmp.vol                = V_tmp;
                
                maps  = addstruct(maps, maps_tmp);
            end % sn
            outfilename = fullfile(glmSuitGroupDir, sprintf('group_%s_%s-%s_taskCon.nii', type, taskNames{cc}, con_vs));
            opt.dmtx    = 1;
            % calculate the average across subjects
            spm_imcalc(infilename', outfilename, 'nanmean(X)', opt);
            
            % saving the gifti file for group map
            V = spm_vol(fullfile(glmSuitGroupDir, sprintf('group_%s_%s-%s_taskCon.nii', type, taskNames{cc}, con_vs)));
            D = suit_map2surf(V,'stats','nanmean');
            
            G = surf_makeFuncGifti(D, 'anatomicalStruct', 'Cerebellum', 'columnNames', {sprintf('group_%s_%s-%s', type, taskNames{cc}, con_vs)});
            
            save(G, fullfile(glmSuitGroupDir, sprintf('Cereb.group.con_%s-%s_taskCon.func.gii', taskNames{cc}, con_vs)));
            save(fullfile(glmSuitGroupDir, sprintf('indMaps_%s_%s-vs-%s.mat', type, taskNames{cc}, con_vs)), 'maps', '-v7.3');
            fprintf('******************** %s group average for %s vs %s is created! ********************\n\n', type, taskNames{cc}, con_vs);
        end % contrasts (cc)
        save(fullfile(glmSuitGroupDir, sprintf('indMaps_%s-vs-%s.mat', type, con_vs)), 'maps', '-v7.3');

    case 'Visualize:mdtb:suit'
        % takes in a vector (or map) for a subject and transfer it to suit space and plot it
        % on the flatmap
        % Example: sc1_sc2_mdtb('Visualize:mdtb:suit', corrmap);

        subj       = 2;           % the suit anatomical data for the subject is used in the mapping 
        experiment_num = 1;
        data       = varargin{1}; % the vector map you want to transfer to suit flatmap
        
        vararginoptions(varargin, {'sn', 'experiment_num'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % Determine the voxels we want to resample in SUIT space
        V = spm_vol(fullfile(baseDir, experiment,suitDir,'anatomicals','cerebellarGreySUIT.nii'));
        X = spm_read_vols(V);
        
        grey_threshold = 0.1; % gray matter threshold
        
        linIn1     = find(X > grey_threshold);
        [i1,j1,k1] = ind2sub(V.dim,linIn1');
        [x1,y1,z1] = spmj_affine_transform(i1, j1, k1, V.mat);
        
        %%% Determine voxel locations from the original ROI
        load(fullfile(baseDir, experiment, regDir, 'data', subj_name{subj}, 'regions_cerebellum_grey.mat')); % 'regions' are defined in 'ROI_define'
        
        Vmask      = spm_vol(fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{subj}, 'maskbrainSUITGrey.nii')); %% cerebellum grey matter mask
        [i3,j3,k3] = spmj_affine_transform(R{1}.data(:,1),R{1}.data(:,2),R{1}.data(:,3),inv(Vmask.mat));
        linIn3     = sub2ind(Vmask.dim, round(i3), round(j3), round(k3));
        
        %%% transform SUIT coords into anatomical space of the individual
        flowfield    = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{subj}, 'u_a_c_anatomical_seg1.nii');
        affine       = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{subj}, 'Affine_c_anatomical_seg1.mat');
        [Def, Aff]   = spmdefs_get_dartel(flowfield, affine);
        [x2, y2, z2] = spmdefs_transform(Def, Aff, x1, y1, z1);
        [i2, j2, k2] = spmj_affine_transform(x2, y2, z2, inv(Vmask.mat));
        
        %%% resample the weights into SUIT space for the lp vector
        
        Vout             = Vmask;
        Vout.dat         = zeros(Vout.dim);
        Vout.dat(linIn3) = data;
        Vout.dt          = [64 0];
        Vout.pinfo       = [1 0 0]';
        
        DataSUIT(ip,:) = spm_sample_vol(Vout,i2,j2,k2,1);
        
        V.dat   = zeros(V.dim);
        Vres(ip) = V;
        Vres(ip).dat(linIn1) = DataSUIT(ip,:);  % Offset by one to account for 1 being medial wall
        Vres(ip).fname       = sprintf('data_%2.2d.nii',ip);
        Vres(ip).pinfo       = [1 0 0]';
        
        myV = Vres(ip);
        
        %%% Now map the Lp vector to surface-based representation
        D = suit_map2surf(myV,'stats','nanmean');
        
        figure; suit_plotflatmap(D , 'cmap', colormap(jet(256)), 'cscale', [min(D(:)), max(D(:))]);
        caxis([min(D(:)), max(D(:))]);
        colorbar;
        title(sprintf('cerebellar regions with high corr with cortex parcel %0.2d for %s', ip, subj_name{subj}));
        
        keyboard;
end
end

% Local functions
function dircheck(dir)
if ~exist(dir,'dir');
    warning('%s doesn''t exist. Creating one now. You''re welcome! \n',dir);
    mkdir(dir);

end
end

function extractCMRRPhysio(varargin)
% -------------------------------------------------------------------------
% extractCMRRPhysio.m
% -------------------------------------------------------------------------
% Extract physiological log files from encoded "_PHYSIO" DICOM file
% generated by CMRR MB sequences (>=R015, >=VD13A)
%   E. Auerbach, CMRR, 2016
%
% Usage:
%    extractCMRRPhysio(DICOM_filename, [output_path]);
%
% This function expects to find a single encoded "_PHYSIO" DICOM file
% generated by the CMRR C2P sequences >=R015. It will extract and write
% individual log files (*_ECG.log, *_RESP.log, *_PULS.log, *_EXT.log,
% *_Info.log) compatible with the CMRR C2P sequences >=R013. Only log
% files with nonzero traces will be written.
%
% Inputs:
%    DICOM_filename = 'XXX.dcm'
%    output_path    = '/path/to/output/' (optional; if not specified, will
%                                         use path of input file)

% say hello
fprintf('\nextractCMRRPhysio: E. Auerbach, CMRR, 2016\n\n');

if (nargin < 1) || (nargin > 2)
    error('Invalid number of inputs.');
end
fn = varargin{1};
if (nargin == 2)
    outpath = varargin{2};
else
    [outpath, ~, ~] = fileparts(fn);
end

% first, verify our input is a DICOM file
if (2 ~= exist(fn,'file'))
    error('%s not found!', fn);
end
if (~isdicom(fn))
    error('%s not a DICOM file!', fn);
end

% read in the DICOM
fprintf('Attempting to read CMRR Physio DICOM format file...\n');
warning('off','images:dicominfo:attrWithSameName');
dcmInfo = dicominfo(fn);
if (~isempty(dcmInfo) && isfield(dcmInfo,'ImageType') && strcmp(dcmInfo.ImageType,'ORIGINAL\PRIMARY\RAWDATA\PHYSIO') ...
        && isfield(dcmInfo,'Private_7fe1_10xx_Creator') && strcmp(deblank(char(dcmInfo.Private_7fe1_10xx_Creator)),'SIEMENS CSA NON-IMAGE'))
    np = size(dcmInfo.Private_7fe1_1010,1);
    rows = dcmInfo.AcquisitionNumber;
    columns = np/rows;
    numFiles = columns/1024;
    if (rem(np,rows) || rem(columns,1024)), error('Invalid image size (%dx%d)!', columns, rows); end
    dcmData = reshape(dcmInfo.Private_7fe1_1010,[],numFiles)';
    % encoded DICOM format: columns = 1024*numFiles
    %                       first row: uint32 datalen, uint32 filenamelen, char[filenamelen] filename
    %                       remaining rows: char[datalen] data
    [~,~,endian] = computer;
    needswap = ~strcmp(endian,'L');
    for idx=1:numFiles
        datalen = typecast(dcmData(idx,1:4),'uint32');
        if needswap, datalen = swapbytes(datalen); end
        filenamelen = typecast(dcmData(idx,5:8),'uint32');
        if needswap, filenamelen = swapbytes(filenamelen); end
        filename = char(dcmData(idx,9:9+filenamelen-1));
        logData = dcmData(idx,1025:1025+datalen-1);
        outfn = fullfile(outpath, filename);
        fprintf('  Writing: %s\n', outfn);
        fp = fopen(outfn,'w');
        fwrite(fp, char(logData));
        fclose(fp);
    end
    fprintf('\nDone!\n');
else
    error('%s is not a valid encoded physio DICOM format file!', fn);
end
end