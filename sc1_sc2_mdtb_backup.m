function [ varargout ] = sc1_sc2_mdtb_backup( what, varargin )
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
lpDir        = 'Lp';                    %% Lps will be saved here
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

%%% ROIs
corticalParcels    = {'yeo_7WB', 'yeo_17WB', 'tesselsWB', 'desikan', 'cortex_grey', 'cortex'};
cerebellarParcels = {'Buckner_7', 'Buckner_17', 'cerebellum_grey', 'cerebellum_MDTB'};

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
structs = {'cortex', 'cerebellum'};
hemName = {'CortexLeft', 'CortexRight'};

% The total number of conditions in experiment 1 and 2, including separate
% instruction regerssor for each task
nCond_sc1 = 28 + 16;
nCond_sc2 = 31 + 16;
nTask     = 16; %% there are 16 tasks and hence 16 instructions

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
        
        sn   = [24, 25, 26, 27, 28];
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
                    
                    if ss == 2
                        sca_ind = sca + 8;
                    elseif ss == 1
                        sca_ind = sca;
                    end
                    save(fullfile(outDir, sprintf('%s_sess%d_scan%d_R.mat', subj_name{s}, ss, sca_ind)), 'R', '-v7.3')
                    close all
%                     keyboard;
                    fprintf('\n********** Physio regressors extracted for %s sess %d scan %d **********\n', subj_name{s}, ss, sca_ind);
                end % sca (scan)
            end % ss (sess)
        end % s (sn)
    case 'PHYS:mdtb:lin_reg' % regresses HRV and RVT on the regressor(s) for instructions
        %%% uses one subject
        % Example: sc1_sc2_mdtb('PHYS:mdtb:lin_reg')
        sn             = [24, 25, 26, 27, 28];
        experiment_num = 1; 
        glm            = 7;
        sess           = 1:2;
        scan           = 1:8;
        
        vararginoptions(varargin, {'sn', 'experiment_num', 'glm', 'sess', 'scan'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % GLM directory
        glmDir    = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        PhysioDir = fullfile(baseDir, 'Physio');
        
        for s = sn
            
            % load in the SPM file (This could take a while)
%             load(fullfile(glmDir, subj_name{s}, 'SPM.mat'));
            load(fullfile(PhysioDir, 'X.mat'));
            
            % load in SPM_info file
            T      = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            
            phys_cell = cell(length(sess), length(scan)); % preallocating the cell mat that will have all the linear reg results for a subject
            for ss = sess
                for sca = scan
                    
                    X_run1     = X(1:598, 1:end - 16);                                        % discarding the intercepts
                    Xuse       = X_run1;
                    
                    if ss == 2
                        sca_ind = sca + 8;
                    elseif ss == 1
                        sca_ind = sca;
                    end
                    
                    ind   = ((T.sess == ss) & (T.run == sca_ind) & (T.inst == 1) & (T.deriv == 0) ); % get the indices for run1, instructions, non-derivatives
                    X_ind = Xuse(:, ind);                                                          
                    X_reg = [zeros(3, 16); X_ind];                                              % My design matrix, adding zeros for dummies!
                    % load in HRV and RVT regressors
                    load(fullfile(PhysioDir, subj_name{s}, sprintf('%s_sess%d_scan%d_physio.mat', subj_name{s}, ss, sca_ind)));
                    
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
%         C  = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        C  = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM copy.txt'));
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
        % Example: sc1_sc2_mdtb('GLM:mdtb:design_glm8', 'sn', a);
        sn            = returnSubjs; %% list of subjects
        experiment_num    = 1;           %% sc1 or sc2?
        ppmethod      = '';          %% 'stc' or ''? The default is set to ''
        deriv         = 1;           %% 0, 1, or 2 for no derivative, temporal, and temporal + dispersion?
        glm           = 8;           %% the glm number      
        
        vararginoptions(varargin,{'sn', 'experiment_num', 'ppmethod', 'deriv', 'glm'});
                
        % load in task information
%         C  = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        C     = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM copy.txt'));
        Cc    = getrow(C, C.StudyNum == experiment_num);
        Tasks = unique(Cc.taskNames,'rows','stable'); % get the task names
        Tasks(strcmp(Tasks, 'Instruct')) = [];        % .dat file with all the info for the tasks does not have 'Instruct', so I'm eliminating it here!
        nTask = unique(length(Tasks));   % how many tasks there are? for sc1: 18 (including rest) and sc2: 33 (including rest)

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
                            if ST>1
                                S.taskName_before = P.taskName(ST-1);
                            else
                                S.taskName_before = 'NaN';
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
        glm        = 82;            %% The glm number :)
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
        % Example1: sc1_sc2_mdtb('GLM:mdtb:contrast', 'sn', [2, 4, 6, 8, 9, 10, 12, 14, 15], 'glm', 8, 'which', 'task')
        % Example2: sc1_sc2_mdtb('GLM:mdtb:contrast', 'sn', [3], 'glm', 72, 'which', 'cond')
        
        sn         = returnSubjs;        %% list of subjects
        glm        = 7;           %% The glm number :)
        experiment_num = 1;
        con_vs     = 'average_1'; %% set it to 'rest' or 'average' (depending on the contrast you want)
        which      = 'task';      %% it can be set to either cond or task. set it to 'task for GLM_8 and 'cond' for GLM_7
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'ppmethod', 'con_vs', 'which'})
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            T    = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            
            % t contrast for tasks
%             ucondition = unique(T_rd.(cond));
            ucondition = unique(T.(which));
            idx = 1;
            for tt = 1:length(ucondition) % 0 is "instruct" regressor
                switch con_vs
                    case 'rest_task' % contrast against rest
                        con                                  = zeros(1,size(SPM.xX.X,2));
%                         con(:,T_rd.(cond) == ucondition(tt)) = 1;
                        con(:,logical((T.(which) == ucondition(tt)) .* (T.deriv == 0))) = 1;
                        tmp = zeros(1, size(SPM.xX.X, 2));
                        tmp(:, strcmp(T.TN, 'rest')) = 1;
                        sum_rest = sum(tmp);
                        tmp = tmp./sum_rest;
                        con                                  = con/abs(sum(con));
                        con(:, tmp~=0) = -1/sum_rest;
                    case 'average_1' % contrast against the average of the other tasks including the instructions
%                         con                                  = zeros(1,size(SPM.xX.X,2));
%                         con(:,logical((T.(which) == ucondition(tt)) .* (T.deriv == 0))) = 1;
%                         con                                  = con/abs(sum(con));
                        %%% if you want to compare against the average of
                        %%% all the tasks including the one you are
                        %%% interested in:
%                         con                                  = bsxfun(@minus, con, 1/length(T.(which)));
                        %%% if you want to compare against the average of
                        %%% all the tasks other than the task you are
                        %%% interested in.
%                         con(con == 0)     = -1;
%                         con(T.deriv == 1) = 0;
%                         tmp            = con == -1;
%                         tmp_count      = sum(tmp); %% how many beta values other than the ones we are interested in
%                         con(con == -1) = con(con == -1)/tmp_count; 
                        % New: Eva, Oct 2nd
                        con        = zeros(1,size(SPM.xX.X, 2));
                        con(1,logical((T.(which) == ucondition(tt)) .* (T.deriv == 0))) = 1./sum(logical((T.(which) == ucondition(tt)) .* (T.deriv == 0)));
                        con(1,logical((T.(which) ~= ucondition(tt)) .* (T.deriv == 0))) = -1./sum(logical((T.(which) ~= ucondition(tt)) .* (T.deriv == 0)));
                        
%                         ncondition = ucondition(ucondition ~= ucondition(tt));
%                         con        = zeros(length(ncondition),
%                         size(SPM.xX.X, 2)); % Eva: shouldn't this be
%                         ncondition-1?
%                         for j = 1:length(ncondition)
%                             k = find(T.task == ucondition(tt) & T.deriv == 0);
%                             l = find(T.task == ncondition(j) & T.deriv == 0);
% 
%                             con(j, k) = 1;
%                             tmp_conin = con(j, k);
%                             tmp_con   = con(j, :);
%                             tmp_conin_sum = sum(tmp_conin);
%                             con(j, k) = con(j, k)./tmp_conin_sum;
% 
%                             tmp_con(l) = -1;
%                             tmp_con(tmp_con == 1) = 0;
%                             tmp_con_sum = -1*sum(tmp_con);
% 
%                             con(j, l) = -1;
%                             con(j, l) = con(j, l)./tmp_con_sum;
%                        end % j
                      %  if (any(abs(sum(con,2))>0.00001))
                      %      keyboard;
                      %  end;
                    case 'average_2' % contrast against the average of the other tasks not including the instructions
                        con        = zeros(1,size(SPM.xX.X, 2));
                        % TO TRY below - no instructions as contrasts
                        con(1,logical((T.(which) == ucondition(tt)) .* (T.deriv == 0))) = 1./sum(logical((T.(which) == ucondition(tt)) .* (T.deriv == 0)));
                        con(1,logical((T.(which) ~= ucondition(tt)) .* (T.deriv == 0) .* (T.inst == 0))) = -1./sum(logical((T.(which) ~= ucondition(tt)) .* (T.deriv == 0) .* (T.inst == 0)));
                    case 'rest'
                        con                                  = zeros(1,size(SPM.xX.X,2));
%                         con(:,T_rd.(cond) == ucondition(tt)) = 1;
                        con(:,logical((T.(which) == ucondition(tt)) .* (T.deriv == 0))) = 1;
                        con                                  = con/abs(sum(con));
                end
                switch which
                    case 'cond'
                        name = sprintf('%s-%s',char(unique(T.TN(T.(which) == ucondition(tt)))), con_vs);
                    case 'task'
%                         name = sprintf('%s-%s',char(unique(T.TN(T.it == ucondition(tt)))), con_vs);
                        name = sprintf('%s-%s',char(unique(T.TN(T.task == ucondition(tt)))), con_vs);
                end
                
                 SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
                %SPM.xCon(idx) = spm_FcUtil('Set',name, 'F', 'c',con',SPM.xX.xKXs);
                idx=idx+1;
            end % tt (conditions)
            SPM = spm_contrasts(SPM,1:length(SPM.xCon));
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glmDir, subj_name{s}, 'SPM_light.mat'), 'SPM');

            % rename contrast images and spmT images
             conName = {'con','spmT'};
%            conName = {'ess','spmF'};
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
        con_vs     = 'average_1'; %% set it to 'rest' or 'average_1' or 'average_2' (depending on the contrast you want)
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'con_vs'})
        
        % gt the task info
        C   = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM copy.txt'));
        Cc  = getrow(C,C.StudyNum == experiment_num);
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
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
    case 'GLM:mdtb:contrast_utransitions_names'
        %%% Calculating contrast images for unique transitions. The names
        %%% of the task before and task after are appended to the contrast
        %%% name.
        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM)
        % Example: sc1_sc2_mdtb('GLM:mdtb:contrast_utransitions_names', 'sn', [3])
        
        sn         = returnSubjs;   %% list of subjects
        glm        = 7;      %% The glm number :)
        experiment_num = 1;
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num'})
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        load(fullfile(glmDir, 'all_trans_task_names.mat'));
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            T    = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            
            idx = 1;
            for tt = 1:size(all_trans, 1)% 0 is "instruct" regressor
                con = zeros(1,size(SPM.xX.X,2));
                
                con(:,strcmp(T.taskName_before, all_trans(tt, 1)) + strcmp(T.taskName_after, all_trans(tt, 2)) == 2) = 1; % contrast against rest
                
                con  = con/abs(sum(con));
                name = sprintf('transition_%s_%s-rest', all_trans{tt, 1}, all_trans{tt, 2});
                
                SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
                idx=idx+1;
            end % tt (conditions)

            SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glmDir, subj_name{s}, 'SPM_light.mat'), 'SPM');

            % rename contrast images and spmT images
            conName = {'con','spmT'};
            for i = 1:size(all_trans, 1)
                for n = 1:numel(conName)
                    oldName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%2.4d.nii',conName{n},i));
                    newName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                    movefile(oldName{i}, newName{i});
                end % conditions (n, conName: con and spmT)
            end % i (contrasts)
        end % sn
    case 'GLM:mdtb:contrast_transitions_id'
        %%% Calculating contrast images for all transitions. There will be
        %%% 256 (16*16) transitions for sc1
        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM)
        % Example: sc1_sc2_mdtb('GLM:mdtb:contrast_transitions_id', 'sn', [3])
        
        sn             = returnSubjs;   %% list of subjects
        glm            = 7;      %% The glm number :)
        experiment_num = 1;
        con_vs         = 'rest'; %% contrasts were calculated vs 'rest' or 'average'
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num'})
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        load(fullfile(glmDir, 'all_trans_task_names.mat'));
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            T    = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            T_rd = getrow(T, T.deriv == 0);
            
            % load in trans_info
            load(fullfile(glmDir, 'trans_info.mat'));
            
            idx = 1;
            for tt = 1:max(t.run)
                for it = (unique(t.instOrder_run(t.run==tt)))'
                    
                    con = zeros(1,size(SPM.xX.X,2));
                    con(:,T_rd.run == tt & T_rd.instOrder == it & T_rd.inst == 1) = 1;
                     if con==0
                         keyboard;
                     end
                    con  = con/abs(sum(con));
                    name = sprintf('transition_%d-%s', idx, con_vs);
                    
                    SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
                    idx=idx+1;
                end
            end % tt run

            SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glmDir, subj_name{s}, 'SPM_light.mat'), 'SPM');

            % rename contrast images and spmT images
            conName = {'con','spmT'};
            for i = 1:length(t.instOrder_all)
                for n = 1:numel(conName)
                    oldName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%2.4d.nii',conName{n},i));
                    newName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                    movefile(oldName{i}, newName{i});
                end % conditions (n, conName: con and spmT)
            end % i (contrasts)
        end % sn
    case 'GLM:mdtb:transition_info'
        % gets the task transition info. Uses one subject to get the
        % transition info for all the subjects, assuming that task
        % transitions are the same across all the subjects.
        % Example: sc1_sc2_mdtb_backup('GLM:mdtb:transition_info')
        
        sn         = 3;  %% the subject you want to use to get the transition info
        experiment_num = 1;
        glm        = 8;
        
        vararginoptions(varargin, {'experiment_num', 'glm'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directory paths
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        T = load(fullfile(glmDir, subj_name{sn}, 'SPM_info'));
        
        if isfield(T, 'deriv') % derivatives were included in the glm
            t = getrow(T, T.deriv == 0 & T.inst == 1); % getting all the non-derivative regressors and regressors corresponding to instructions
            t = rmfield(t,{'SN','TN', 'deriv', 'task'});         % removing unnecessary fields
        elseif ~isfield(T, 'deriv') % derivatives were not included in the glm
            t = getrow(T, T.inst == 1);
            t = rmfield(t,{'SN','TN','task'});         % removing unnecessary fields
        end % checking if derivatives were included in the glm

        t.instOrder_all = (1:length(t.inst))';
        t.instOrder_run = t.instOrder;
        t               = rmfield(t, {'instOrder'});
        
        save(fullfile(glmDir, 'trans_info'), 't', '-v7.3');
        
    case 'ROI:mdtb:define'
        % defines the cortical ROIs using atlases' workbench gifti files for each subject
        % For the cerebellar parcels, you may need to run
        % 'SUIT:mdtb:suit_parcel2native' first!
        % Example: sc1_sc2_mdtb('ROI:mdtb:define', 'sn', [3])
        
        sn         = returnSubjs;
        atlas_res  = 32;        %% atlas resolution set to 32 or 164
        nNodes     = 162;       %% options: 162, 362, 642, 1002, 1442
        experiment_num = 1;
        glm        = 7;
        parcelType = 'yeo_7WB'; %% set it to 'tesselsWB', 'yeo_7WB', or 'yeo_17WB', 'Buckner_7', 'Buckner_17' (type of the parcels you want to use), 'cortex_cole'
        
        vararginoptions(varargin, {'sn', 'atlas', 'nNodes', 'experiment_num', 'glm', 'type'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directory paths
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                
        surfDir = fullfile(baseDir, experiment, sprintf('fs_LR_%d', atlas_res)); %% the directory where the atlas is located
        for s = sn
            fprintf('******************** Defining ROIs for %s ********************\n', subj_name{s})
            R = []; %% initializing the R variable which will contain the data for the subject's ROI
            
            mask = fullfile(glmDir, subj_name{s}, 'mask.nii');
            
            dircheck(fullfile(baseDir, experiment,regDir,'data',subj_name{s}));
            
            if ismember(parcelType, corticalParcels)
                idx       = 1;  %% parcel index
                subjWbDir = fullfile(baseDir, experiment, wbDir, 'data',subj_name{s});
                
                for h = 1:2 % two hemispheres
                    switch parcelType
                        case 'tesselsWB'
                            G          = gifti(fullfile(surfDir,sprintf('Icosahedron-%d.%dk.%s.label.gii', nNodes, atlas_res, hemI{h})));
                        case 'yeo_17WB'
                            G          = gifti(fullfile(surfDir, sprintf('Yeo_JNeurophysiol11_17Networks.%dk.%s.label.gii', atlas_res, hemI{h})));
                        case 'yeo_7WB'
                            G          = gifti(fullfile(surfDir, sprintf('Yeo_JNeurophysiol11_7Networks.%dk.%s.label.gii', atlas_res, hemI{h})));
                        case 'cortex_cole'
                            G          = gifti(fullfile(surfDir, sprintf('%s.%dk.%s.label.gii', parcelType, atlaas_res, hemI{h})));
                    end % switch type
                    
                    nReg = numel(unique(G.cdata))-1; % exclude 0 - medial wall
                    for r=1:nReg
                        R{idx}.type     = 'surf_nodes_wb';
                        R{idx}.location = find(G.cdata(:,1) == r);
                        R{idx}.white    = fullfile(subjWbDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                        R{idx}.pial     = fullfile(subjWbDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                        R{idx}.linedef  = [5, 0, 1];
                        R{idx}.image    = mask;
                        R{idx}.name     = sprintf('%s.parcel-%0.2d', hemI{h}, r);
                        
                        idx = idx + 1;
                    end
                    R = region_calcregions(R);
                end % h (hemi)
            elseif ismember(parcelType, cerebellarParcels)
                    parcelDir  = fullfile(suitToolDir, 'atlasesSUIT'); %% directory where the nifti image is stored
                    
                    switch parcelType
                        case 'Buckner_7'
                            CV          = spm_vol(fullfile(parcelDir, sprintf('%sNetworks.nii', parcelType)));
                        case 'Buckner_17'
                            CV          = spm_vol(fullfile(parcelDir, sprintf('%sNetworks.nii', parcelType)));
                        case 'Cerebellum_cole'
                            CV          = spm_vol(fullfile(parcelDir, sprintf('%s.nii', parcelType)));
                    end % switch parcelType for cerebellum
                    CX          = spm_read_vols(CV);
                    nparcel     = length(unique(CX)) - 1; % 0 will also be discarded
                    
                    for r = 1:nparcel
                        file = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, sprintf('iw_%sNetworks_u_a_c_anatomical_seg1.nii', parcelType));
                        R{r}.type  = 'roi_image';
                        R{r}.file  = file;
                        R{r}.name  = sprintf('cereb_parcel-%0.2d', r);
                        R{r}.value = r;
                        R{r}.image = mask;
                    end % i(regions/parcels)
                    R = region_calcregions(R);
            end % cortex or cerebellum
            fprintf('\n******************** %s Rois were defined for %s ********************\n\n', parcelType, subj_name{s});
            save(fullfile(baseDir, experiment, regDir, 'data', subj_name{s}, sprintf('regions_%s.mat',parcelType)),'R', '-v7.3');
        end % sn (subject)
    case 'ROI:mdtb:empty_parcel'
        % for each subject identifies the ROIs that are empty!
        % This case output will then be used to identify the ROIs that are
        % non-empty across all the subjects
        % Example: sc1_sc2_mdtb('ROI:mdtb:empty_parcel', 'sn', 3)
        
        sn           = returnSubjs;
        parcellation = 'Buckner_17'; %% the parcellation you want to inspect
        experiment_num   = 1;
        
        vararginoptions(varargin, {'sn', 'parcellation', 'experiment_num'})
        
        experiment = sprintf('sc%d', experiment_num);
        
        empty_roi = {};
        for s = sn
            % load in the region mat file
            load(fullfile(baseDir, experiment, regDir, 'data', subj_name{s}, sprintf('regions_%s.mat', parcellation)));
            
            % R is a structure variable and the field that needs to be
            % checked is the data field. So for each parcel (region) you
            % will have to check the data field and if that's empty, then
            % you can return the ROI number for which the data field is
            % empty alongside the id for the subject.
            for ip = 1:size(R, 2)
                data_tmp = R{ip}.data;
                if isempty(data_tmp)
                    empty_roi{s} = ip;
                end
            end % ip (parcel)
        end % sn
        varargout{1} = empty_roi;
    case 'ROI:mdtb:make_nifti'
        % makes nifti file for each ROI. Useful for checking.
        % Example: sc1_sc2_mdtb('ROI:mdtb:make_nifti', 'sn', [3]);
        
        sn         = returnSubjs;      %% subject for whom you want to create the ROI nifti files. (I always use s03 :))         
        parcelType = 'Buckner_7';
        experiment_num = 1;
        
        vararginoptions(varargin,{'sn', 'parcelType', 'experiment_num'});
        
        experiment = sprintf('sc%d', experiment_num);
           
        for s = sn
            % load ROI definition
            load(fullfile(baseDir, experiment, regDir, 'data', subj_name{s}, sprintf('regions_%s.mat', parcelType)));
            
            if ismember(parcelType, corticalParcels) % the parcelType entered is from cortex
                mask = fullfile(baseDir, experiment, regDir, 'data', subj_name{s}, 'cortical_mask_grey.nii');
            elseif ismember(parcelType, cerebellarParcels) % the parcelType entered is from cerebellum
                    mask = fullfile(baseDir, experiment, regDir, 'data', subj_name{s}, 'maskbrainSUITGrey.nii');
            end % a cortical ROI or a cerebellar ROI?
            
            % loop over rois
            for roi = 1:size(R,2)
                % Save region file as nifti
                cd(fullfile(baseDir, experiment, regDir, 'data', subj_name{s}));
                region_saveasimg(R{roi}, mask);      
            end % roi(region/parcel)
        end % sn (subject)
    case 'ROI:mdtb:beta_unn'
        % univariately normalizes beta values and saves both the normalized
        % and non-normalized beta values in a new structure. The case for
        % time series extraction only saves the non-normalized beta values.
        % for sc1, there are 88 beta values per subject. The normalization
        % method is 'runwise'.
        % Example: sc1_sc2_mdtb_backup('ROI:mdtb:beta_unn', 'sn', [2, 3, 4, 6, 8, 9, 10, 12, 14, 15, 17, 18, 19, 20]);
        
        sn         = returnSubjs;
        ppmethod   = '';
        experiment = 2;
        roi        = 'yeo_17WB'; %% other options are 'Buckner_17', 'yeo_7WB', and 'yeo_17WB'
        glm        = 8;
        oparcel    = 1;           %% the parcel that will be omited due to the fact that it was empty in some of the subjects
        %%% run 'ROI:mdtb:empty_parcel' to get oparcel. 
        discardp   = 1;           %% set this flag to 1 if you want to omit a parcel and 0 otherwise.
        
        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'glm', 'roi', 'oparcel', 'discardp'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm));
            case 'stc'
                glmDir = fullfile(baseDir, experiment, 'GLM_firstlevel_%d_stc', glm);
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d_stc', glm));
        end
        dircheck(betaDir); % beta values for each roi is saved here

        if ismember(roi, corticalParcels) % the roi is from cortex
                structure_cc = 'cortex';
        elseif ismember(roi, cerebellarParcels) % the roi is from cerebellum
                    structure_cc = 'cerebellum';
        end % a cortical roi or a cerebellar roi?
        
        for s = sn
            fprintf('******************** calculating univariate normalized betas for %s ********************\n', subj_name{s});
            B = [];
            
            % load in SPM_info to get the info for the tasks and
            % condditions
            T           = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            T_rd        = getrow(T, T.deriv == 0); %% removing the derivative regressors
            nregressors = length(T_rd.cond)/(length(runLst)); % number of regressors (deriv and non-deriv) in a run
            indd        = T.deriv == 0;
            
            % load in spm file
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'));
            
            % load in ROI file
            load(fullfile(baseDir, 'sc1', regDir, 'data', subj_name{s}, sprintf('regions_%s.mat', roi)));
            
            % Get the raw data files
            V=SPM.xY.VY;
            
            Y = region_getdata(V,R);  % Data is N x P
            
            for r=1:numel(R) % R is the output 'regions' structure from 'ROI_define'
                % Get betas (univariately prewhitened)
                fprintf('get the beta values for region %0.2d\n', r);
                
                switch discardp
                    case 0 % don't want to discard any parcel
                        [~,resMS,~,beta_hat,~,~] = rsa.spm.noiseNormalizeBeta(Y{r},SPM,'normmode','runwise');
                    case 1 % want to discard the parcel that is empty and is identified by oparcel
                        if r ~= oparcel
                            [~,resMS,~,beta_hat,~,~] = rsa.spm.noiseNormalizeBeta(Y{r},SPM,'normmode','runwise');
                        else
                            resMS    = [];
                            beta_hat = [];
                        end
                end % switch discardp
                
                % option to include empty betasUW matrix (no betas for some
                % ROIs for some subjects)
                if isempty(beta_hat)
%                     betasUW = zeros(size(beta_hat,1),1);
                    
                    B_tmp.betasNW{1, 1} = 0;
                    B_tmp.betasUW{1, 1} = 0;
                    B_tmp.mbetasNW      = zeros(1, size(SPM.xX.X, 2));
                    B_tmp.mbetasUW      = zeros(1, size(SPM.xX.X, 2));
%                     B_tmp.rmbetasUW     
%                     B_tmp.rmbetasNW
                else
                    betasUW = bsxfun(@rdivide,beta_hat,sqrt(resMS)); % univariate noise normalisation
                    B_tmp.betasNW{1, 1} = beta_hat;
                    B_tmp.betasUW{1, 1} = betasUW;
                    B_tmp.mbetasNW      = mean(beta_hat, 2)';
                    B_tmp.mbetasUW      = mean(betasUW, 2)';
                    
                    % calculate the mean beta across runs
                    %%% for all the voxels in the parcel
                    %%%% remove the intercepts
                    beta_hat_ri = beta_hat(1:end - length(runLst), :);
                    betasUW_ri  = betasUW(1:end - length(runLst), :);
                    
                    beta_hat_ri_rd = beta_hat_ri(indd, :)';
                    betasUW_ri_rd  = betasUW_ri(indd, :)';
                    
                    nvox = size(beta_hat, 2); % number of voxels in a parcel.
                    beta_hat_ri_reshaped = reshape(beta_hat_ri_rd, [nvox, nregressors, length(runLst)]);
                    betasUW_ri_reshaped  = reshape(betasUW_ri_rd, [nvox, nregressors, length(runLst)]);
                    
                    mr_beta_hat_ri_reshaped = mean(beta_hat_ri_reshaped, 3);
                    mr_betasUW_ri_reshaped  = mean(betasUW_ri_reshaped, 3);
                    
                    B_tmp.mrbetasNW{1, 1} = mr_beta_hat_ri_reshaped;
                    B_tmp.mrbetasUW{1, 1} = mr_betasUW_ri_reshaped;
                    
                    %%% for the mean beta in a parcel
                    B_tmp.mrmbetasNW = mean(mr_beta_hat_ri_reshaped, 1);
                    B_tmp.mrmbetasUW = mean(mr_betasUW_ri_reshaped, 1);
                end
                
                B_tmp.region_num    = r;
                B_tmp.region_name   = R{r}.name;
                B_tmp.SN            = s;
                B_tmp.structure_cc  = structure_cc; %% is the structure from cortex or cerebellum?
                
                B = addstruct(B, B_tmp);
            end % R (regions)
            
            % save the betas
            dircheck(fullfile(betaDir, subj_name{s}));
            save(fullfile(betaDir, subj_name{s}, sprintf('beta_regions_%s.mat', roi)), 'B', '-v7.3');
            fprintf('******************** betas calculated! ********************\n\n');
        end % sn 
    case 'ROI:mdtb:extract_ts'
        %%% This case will be used to get the raw, adjusted, predicted, and
        %%% residual time series as well as the beta values. Before running
        %%% this case, make sure that you have ran 'GLM:mdtb:study1_glm4'
        %%% to get the design matrix that you want.
        %%% You have to run it for each roi separately
        % Example: sc1_sc2_mdtb('ROI:mdtb:extract_ts', 'sn', [3])
        sn         = returnSubjs;       %% subject list
        ppmethod   = '';         %% Did you do the slice timing correction as a preprocessing step?
        experiment = 1;
        glm        = 7;          %% the GLM number
        stat       = 'mean'; %% this option is used to get the time series for the roi. It can be 'not_mean' or 'mean'
        condition  = [];         %% the regressors of interest used in getting the predicted time series
        
        roi = 'tesselsWB_162'; %% Enter the name of the ROI mat file for which you want to get the time series 
        % to the analysis for cortex and the cerebellum the ROI files can
        % be: 'regions_cerebellum_grey' and 'regions_162_tessellation_hem',
        % or 'regions_tesselsWB_162'

        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'glm', 'roi', 'stat'});
        
        % load in task information
        if isempty(condition)
            taskName = 'all';
        else
            C   = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
            Cc  = getrow(C,C.StudyNum == experiment);
            Ccc = getrow(Cc, Cc.condNum == condition(1));

            taskName = Ccc.taskNames;
        end
        
        experiment = sprintf('sc%d', experiment); %% experiment number is converted to 'sc1' or 'sc2'

        switch ppmethod
            case 'stc' 
                glmDir = sprintf('GLM_firstlevel_%d_stc', glm);
                tsDir  = fullfile(baseDir, experiment, sprintf('TS_GLM_%d_stc', glm));
            case ''
                glmDir = sprintf('GLM_firstlevel_%d', glm);
                tsDir = fullfile(baseDir, experiment, sprintf('TS_GLM_%d', glm));
        end
        dircheck(tsDir)
        
        %%% the number of conditions is different between studies
        switch experiment
            case 'sc1'
                % for sc1, there are 28 conditions (not including the instruction)
                % There are 16 unique tasks and each task has its own
                % instruction period, so there will be 16 unique
                % instruction regressors. The total number of
                % regressors are then calculated as 28 + 16. If we also
                % include the temporal derivative of each regressor, then
                % there will be (28+16)*2 regressor in total.
                nCond = nCond_sc1; 
            case 'sc2'
                % for sc2, there are 31 conditions (not including the instruction)
                % There are 16 unique tasks and each task has its own
                % instruction period, so there will be 16 unique
                % instruction regressors. The total number of
                % regressors are then calculated as 31 + 16. If we also
                % include the temporal derivative of each regressor, then
                % there will be (31+16)*2 regressor in total.
                nCond = nCond_sc2;
        end
        
        %%% creating the vector that will contain the regressor numbers for
        %%% the conditions of interest.
        %%%% if the condition vector is not empty then you have to create a
        %%%% vector that contains the numbers of conditions for all the
        %%%% runs.
        if ~isempty(condition)
            regressors = []; %% I'm going to set it to be empty for now
            for rn = 1:length(runLst)
                for nc = 1:length(condition)
                    regressors = [regressors, condition(nc) + (rn - 1)*nCond];
                end
            end
        end 
        
        ts_allsubs = cell(1, length(sn)); %% the cell variable containing the time series for each subject.
        for s = sn
            fprintf('******************** getting the time series for %s ********************\n', subj_name{s});
            %%%% load in the SPM.mat file ...
            load(fullfile(baseDir, experiment, glmDir, subj_name{s}, 'SPM.mat'));
            
            %%%% load in the region file for the cerebellum
            load(fullfile(baseDir, experiment, regDir, 'data', subj_name{s}, sprintf('regions_%s.mat', roi)));
            
            %%%% get the evoked response (not my code!)
            [y_raw, y_adj, y_hat, y_res, B, y_filt] = region_getts(SPM, R, 'stats', stat, 'reg_interest', condition);
            
            %%%% create a structure variable that will contain all the time
            %%%% series for the cerebellum
            ts.y_raw  = y_raw;
            ts.y_filt = y_filt;
            ts.y_adj  = y_adj;
            ts.y_hat  = y_hat;
            ts.y_res  = y_res;
            ts.beta   = B;
            
            ts.ppmethod   = ppmethod;
            ts.glm        = glm;
            ts.roi        = roi;
            ts.regressors = condition;
            
            ts_allsubs{s} = ts;

            %%% Save the time series for the cerebellum
            dircheck(fullfile(tsDir, subj_name{s}));
            save(fullfile(tsDir, subj_name{s}, sprintf('Ts_%s_regions_%s.mat', taskName, roi)), 'ts', '-v7.3');

            fprintf('******************** predicted, raw, and adjusted time series for cerebellum and cortex saved for %s ********************\n\n', subj_name{s})
        end % sn
        varargout{1} = ts_allsubs;    
    case 'ROI:mdtb:extract_ts_pw'
        % extract and prewhiten the time series.
        % it will be like the "filtering" thing I did in the first code for
        % time delay analysis.
        % I am using parts of the code from
        % sc1_sc2_connectivity('TS_get_ts'). It uses the imaging files and
        % ResMS file that has already been calculated and stored in glm
        % directory. Basically, it's:
        % V = SPM.xY.VY;
        % Y = region_getdata(V, R{r});
        % Vres = spm_vol(fullfile(glmDir, subj_name{s}, 'ResMS.nii'));
        % resMS = region_getdata(Vres, R{r});
        % Y = bsxfun(@rdivide, Y, sqrt(resMS));
    case 'ROI:mdtb:average_ts'
        %%% calculating the average time series across runs and saving it
        %%% in a structure variable that will be used in further analysis
        %%% steps.
        % Example: sc1_sc2_mdtb('ROI:mdtb:average_ts', 'sn', a, 'roi', 'regions_tesselsWB_162')
        sn         = returnSubjs;       %% subject list
        ppmethod   = '';         %% Did you do the slice timing correction as a preprocessing step?
        experiment = 1;
        glm        = 7;          %% the GLM number
        condition  = [];         %% the regressors of interest used in getting the predicted time series      
        
        roi = 'cerebellum_grey'; %% Enter the name of the ROI mat file for which you want to get the time series 
        % to the analysis for cortex and the cerebellum the ROI files can
        % be: 'regions_cerebellum_grey' and 'regions_162_tessellation_hem',
        % 'regions_tesselsWB_162'

        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'glm', 'roi'});
        
        if isempty(condition)
            taskName = 'all';
        else
            % load in task information
            C  = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
            Cc = getrow(C,C.StudyNum == experiment);
            Ccc = getrow(Cc, Cc.condNum == condition(1));

            taskName = Ccc.taskNames;
        end
        
        experiment = sprintf('sc%d', experiment); %% experiment number is converted to 'sc1' or 'sc2'
        
        switch ppmethod
            case 'stc' 
                tsDir  = fullfile(baseDir, experiment, sprintf('TS_GLM_%d_stc', glm));
            case ''
                tsDir = fullfile(baseDir, experiment, sprintf('TS_GLM_%d', glm));
        end

        avg_ts_allsubs = cell(1, length(sn)); %% the cell variable containing the time series for each subject.
        for s = sn
            fprintf('******************** calculating the average time series across runs for %s ********************\n', subj_name{s});
            % load in the time series structure variable
            load(fullfile(tsDir, subj_name{s}, sprintf('Ts_%s_regions_%s.mat', taskName, roi)));
            yfilt_tmp = (ts.y_filt)';
            yhat_tmp  = (ts.y_hat)';
            yraw_tmp  = (ts.y_raw)';
            
            % separating the time series for each run. runs will be the
            % third dimension in the matrix
            yfilt = reshape(yfilt_tmp, [size(yfilt_tmp, 1), ntp, length(runLst)]);
            yhat = reshape(yhat_tmp, [size(yhat_tmp, 1), ntp, length(runLst)]);
            yraw = reshape(yraw_tmp, [size(yraw_tmp, 1), ntp, length(runLst)]);
            
            % calculate the average time series across 16 runs
            avg_yfilt = mean(yfilt, 3);
            avg_yhat  = mean(yhat, 3);
            avg_yraw  = mean(yraw, 3);
            
            avgts.y_filt = avg_yfilt';
            avgts.y_hat  = avg_yhat';
            avgts.y_raw  = avg_yraw';
            
            avgts.ppmethod   = ppmethod;
            avgts.glm        = glm;
            avgts.roi        = roi;
            avgts.regressors = condition;
            
            avg_ts_allsubs{s} = avgts;
            
            % save the structure variable
            save(fullfile(tsDir, subj_name{s}, sprintf('avgTs_%s_%s.mat', taskName, roi)), 'avgts', '-v7.3');
            fprintf('******************** average time series across runs for %s calculated! ********************\n\n', subj_name{s});
        end % sn
        
        varargout{1} = avg_ts_allsubs;
        
    case 'SURF:mdtb:map_con'
        % projects individual contrast map volume files for the conditions
        % to the workbench surface.
        % Example: sc1_sc2_mdtb('SURF:mdtb:map_con', 'sn', [3])
    
        sn         = returnSubjs;        %% list of subjects
        ppmethod   = '';          %% with or without stc
        atlas_res  = 32;          %% set it to 32 or 164
        experiment = 1;           %% enter 1 for sc1 and 2 for sc2
        glm        = 8;           %% glm number
        con_vs     = 'average_1'; %% set it to 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment);
        
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
    case 'SURF:mdtb:map_con_new'
        % projects individual contrast map volume files for the conditions
        % to the workbench surface.
        % Example: sc1_sc2_mdtb('SURF:mdtb:map_con', 'sn', [3])
    
        sn         = returnSubjs;        %% list of subjects
        ppmethod   = '';          %% with or without stc
        atlas_res  = 32;          %% set it to 32 or 164
        experiment = 1;           %% enter 1 for sc1 and 2 for sc2
        glm        = 9;           %% glm number
%         con_vs     = 'average_1'; %% set it to 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment);
        
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
            subjSurfDir = fullfile(baseDir, experiment, wbDir, 'data', subj_name{s});
            dircheck(fullfile(glmSurfDir, subj_name{s})); %% directory to save the contrast maps
            
            T  = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            Td = getrow(T, T.deriv == 0); %% get the non-derivative regressor names
            conNames = unique(Td.TN);
            for h = 1:2 % two hemispheres
                white   = fullfile(subjSurfDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                pial    = fullfile(subjSurfDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                C1      = gifti(white);
                C2      = gifti(pial);
                for se = [1, 2]
                    % map ResMS
                    ResMSname = fullfile(glmDir, subj_name{s},'ResMS.nii');
                    ResMS     = surf_vol2surf(C1.vertices,C2.vertices,ResMSname);
                    Data      = zeros(size(ResMS.cdata,1),length(conNames));
                    for ic = 1:length(conNames)
                        betaNum = find(T.cond==ic & T.sess==se);
                        for j=1:length(betaNum)
                            filenames{j}=fullfile(glmDir,sprintf('beta_%04d.nii',betaNum(j)));
                        end
                        Beta       = surf_vol2surf(C1.vertices,C2.vertices,filenames);
                        Data(:,ic) = mean(Beta.cdata,2)./sqrt(ResMS.cdata);
                    end % ic (condition/contrast)
                    fprintf('\n');
                    Data(:, length(conNames) + 1) = 0;
                    
                    Data = bsxfun(@minus,Data,mean(Data,2));
                    
                    G       = surf_makeFuncGifti(Data,'anatomicalStruct',hemName{h},'columnNames',[T.condNames(1:end);'Rest']);
                    outfile = fullfile(subjSurfDir,sprintf('%s.%s.%s.wcon.sess%d.%s.func.gii',subj_name{s},hemName{h},experiment,se,atlas_res));
                    save(G,outfile);
                    
                end % se
                 
            end % hemi
        end % sn
    case 'SURF:mdtb:map_con_task'
        % projects individual contrast map volume files for tasks to
        % WorkBench surface
        % Example: sc1_sc2_mdtb('SURF:mdtb:map_con_task', 'sn', [3])
    
        sn         = returnSubjs;        %% list of subjects
        ppmethod   = '';          %% with or without stc
        atlas_res  = 32;          %% set it to 32 or 164
        experiment = 1;           %% enter 1 for sc1 and 2 for sc2
        glm        = 72;           %% glm number
        con_vs     = 'average_1'; %% set it to 'rest' or 'average_1' or 'average_2'
        
        % gt the task info
        C   = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM copy.txt'));
        Cc  = getrow(C,C.StudyNum == experiment);
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment);
        
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
    case 'SURF:mdtb:map_con_utransitions_names'
        % projects individual contrast map volume files for the unique
        % transitions with the task names to workbench surface
        % Example: sc1_sc2_mdtb('SURF:mdtb:map_con_utransitions_names', 'sn', [3])
    
        sn         = returnSubjs;   %% list of subjects
        ppmethod   = '';     %% with or without stc
        atlas_res  = 32;     %% set it to 32 or 164
        experiment = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 7;      %% glm number
        con_vs     = 'rest'; %% set it to 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                wbDir  = 'surfaceWB';
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                wbDir  = 'surfaceWB_stc';
        end
        load(fullfile(glmDir, 'all_trans_task_names.mat'));
        glmSurfDir = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm));
        dircheck(glmSurfDir);
        
        for s = sn
            fprintf('******************** start mapping contrasts to surface for %s ********************\n', subj_name{s});
            subjSurfDir = fullfile(baseDir, experiment, wbDir, 'data', subj_name{s});
            dircheck(fullfile(glmSurfDir, subj_name{s})); %% directory to save the contrast maps
            
            for tt = 1:size(all_trans, 1)
                for h = 1:2 % two hemispheres
                    white   = fullfile(subjSurfDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                    pial    = fullfile(subjSurfDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                    C1      = gifti(white);
                    C2      = gifti(pial);
                    conMapName      = sprintf('con_transition_%s_%s-%s', all_trans{tt, 1}, all_trans{tt, 2}, con_vs);
                    images{1}       = fullfile(glmDir, subj_name{s}, sprintf('%s.nii', conMapName));
                    column_name{1}  = sprintf('con_transition_%s_%s-%s.nii', all_trans{tt, 1}, all_trans{tt, 2}, con_vs);

                    outfile    = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.con_transition_%s_%s-%s.func.gii', subj_name{s}, hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                    G          = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names', column_name, ...
                                    'anatomicalStruct',hemName{h});
                    save(G, outfile);
                    fprintf('******************** mapping to surface for %s hemi %s contrast %s done********************\n', subj_name{s}, ...
                        hemI{h}, conMapName);
                end % hemi
            end % tt (transitions)
        end % sn    
    case 'SURF:mdtb:map_con_transitions_id'
        % projects individual contrast map volume files for all the
        % transitions (256 for sc1) to workbench surface
        % Example: sc1_sc2_mdtb('SURF:mdtb:map_con_transitions_id', 'sn', [3])
    
        sn         = returnSubjs;   %% list of subjects
        ppmethod   = '';     %% with or without stc
        atlas_res  = 32;     %% set it to 32 or 164
        experiment = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 7;      %% glm number
        con_vs     = 'rest'; %% set it to 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                wbDir  = 'surfaceWB';
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                wbDir  = 'surfaceWB_stc';
        end
        load(fullfile(glmDir, 'trans_info.mat'));
        glmSurfDir = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm));
        dircheck(glmSurfDir);
        
        for s = sn
            fprintf('******************** start mapping contrasts to surface for %s ********************\n', subj_name{s});
            subjSurfDir = fullfile(baseDir, experiment, wbDir, 'data', subj_name{s});
            dircheck(fullfile(glmSurfDir, subj_name{s})); %% directory to save the contrast maps
            
            for tt = 1:length(t.instOrder_all)
                for h = 1:2 % two hemispheres
                    white   = fullfile(subjSurfDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                    pial    = fullfile(subjSurfDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                    C1      = gifti(white);
                    C2      = gifti(pial);
                    conMapName      = sprintf('con_transition_%d-%s', t.instOrder_all(tt), con_vs);
                    images{1}       = fullfile(glmDir, subj_name{s}, sprintf('%s.nii', conMapName));
                    column_name{1}  = sprintf('con_transition_%d-%s.nii', t.instOrder_all(tt), con_vs);

                    outfile    = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.con_transition_%d-%s.func.gii', subj_name{s}, hemI{h}, t.instOrder_all(tt), con_vs));
                    G          = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names', column_name, ...
                                    'anatomicalStruct',hemName{h});
                    save(G, outfile);
                    fprintf('******************** mapping to surface for %s hemi %s contrast %s done********************\n', subj_name{s}, ...
                        hemI{h}, conMapName);
                end % hemi
            end % tt (transitions)
        end % sn
    case 'SURF:mdtb:groupmap_con'
        % creates group average contrast maps for task contrasts
        % Example: sc1_sc2_mdtb('SURF:mdtb:groupmap_con', 'sn', [3])
    
        sn         = returnSubjs;   %% list of subjects
        ppmethod   = '';     %% with or without stc
        atlas_res  = 32;     %% set it to 32 or 164
        experiment = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 8;      %% glm number
        replaceNaN = 1;      %% replacing NaNs
        con_vs     = 'average_1'; %% contrast was calculated against 'rest' or 'average'        
        smooth     = 1;      %% add smoothing
        kernel     = 1;      %% for smoothing
        which      = 'task'; %% 'task' for glm8 and 'cond' for glm7
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment', 'glm', 'replaceNaN', 'con_vs', 'smooth', 'kernel', 'which'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc       = getrow(C, C.StudyNum == experiment);
        switch which
            case 'task' % task for glm8
                conNames = unique(Cc.taskNames);
            case 'cond' % condition for glm7
                conNames = unique(Cc.condNames);
        end %% do you want the group maps for tasks or conditions
        
        experiment = sprintf('sc%d', experiment);
        
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
        experiment = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 8;      %% glm number
        replaceNaN = 1;      %% replacing NaNs
        con_vs     = 'rest'; %% contrast was calculated against 'rest' or 'average'        
        smooth     = 1;      %% add smoothing
        kernel     = 1;      %% for smoothing
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment', 'glm', 'replaceNaN', 'con_vs', 'smooth', 'kernel'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc       = getrow(C, C.StudyNum == experiment);
        taskNames = unique(Cc.taskNames);
        
        experiment = sprintf('sc%d', experiment);
        
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
    case 'SURF:mdtb:groupmap_con_utransitions_names'
        % creates group average contrast maps for unique transitions with
        % task names
        % Example: sc1_sc2_mdtb('SURF:mdtb:groupmap_con_utransitions_names', 'sn', [3])
    
        sn         = returnSubjs;   %% list of subjects
        ppmethod   = '';     %% with or without stc
        atlas_res  = 32;     %% set it to 32 or 164
        experiment = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 7;      %% glm number
        replaceNaN = 1;      %% replacing NaNs
        smooth     = 1;      %% add smoothing
        kernel     = 1;      %% for smoothing
        con_vs     = 'rest'; %% contrast was calculated against 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment', 'glm', 'replaceNaN', 'convs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                wbDir  = 'surfaceWB';
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                wbDir  = 'surfaceWB_stc';
        end
        load(fullfile(glmDir, 'all_trans_task_names.mat'));

        % go to the directory where the group data is 
        groupSurfDir = fullfile(baseDir, experiment, wbDir, 'data', sprintf('group%dk', atlas_res));
        cd(groupSurfDir);
        
        for tt = 1:size(all_trans, 1)
            for h = 1:2 % two hemispheres
                for s = 1:length(sn)
                    %%% make the group metric file for each contrasts
                    infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_transition_%s_%s-%s.func.gii', subj_name{s}, hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                    if smooth
                        surfFile    = fullfile(groupSurfDir,sprintf('fs_LR.32k.%s.inflated.surf.gii',hemI{h}));
                        surf_smooth(infilenames{s},'surf',surfFile,'kernel',kernel); % smooth outfilenames - it will prefix an 's'
                    end
                    s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_transition_%s_%s-%s.func.gii', subj_name{s}, hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                end % sn
                outfilenames    = fullfile(groupSurfDir,sprintf('%s.con_transition_%s_%s-%s.func.gii', hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                summaryname     = fullfile(groupSurfDir,sprintf('%s.group.con_transition_%s_%s-%s.func.gii', hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                
                surf_groupGiftis(infilenames, 'outfilenames', {outfilenames}, 'groupsummary', summaryname, 'replaceNaNs', replaceNaN);
                if smooth % also save the smoothed versions
                    s_outfilenames    = fullfile(groupSurfDir,sprintf('s%s.con_transition_%s_%s-%s.func.gii', hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                    s_summaryname     = fullfile(groupSurfDir,sprintf('s%s.group.con_transition_%s_%s-%s.func.gii', hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                    surf_groupGiftis(s_infilenames, 'outfilenames', {s_outfilenames}, 'groupsummary', s_summaryname, 'replaceNaNs', replaceNaN);
                end
                
                % delet the smoothed contrasts for each subject
                for s = 1:length(sn)
                    delete(fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_transition_%s_%s-%s.func.gii', subj_name{s}, hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs)));
                end
                
                fprintf('******************** group average contrast transition from %s to %s vs %s is created! ********************\n\n', all_trans{tt, 1}, all_trans{tt, 2}, con_vs);
            end % hemi(h)
        end % tt (transitions)
    case 'SURF:mdtb:groupmap_con_transitions_id'
        % creates group average contrast maps for all the transitions
        % Example: sc1_sc2_mdtb('SURF:mdtb:groupmap_con_transitions_id', 'sn', [3])
    
        sn         = returnSubjs;   %% list of subjects
        ppmethod   = '';     %% with or without stc
        atlas_res  = 32;     %% set it to 32 or 164
        experiment = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 7;      %% glm number
        replaceNaN = 1;      %% replacing NaNs
        smooth     = 1;      %% add smoothing
        kernel     = 1;      %% for smoothing
        con_vs     = 'rest'; %% contrast was calculated against 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment', 'glm', 'replaceNaN', 'convs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                wbDir  = 'surfaceWB';
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                wbDir  = 'surfaceWB_stc';
        end
        load(fullfile(glmDir, 'trans_info.mat'));

        % go to the directory where the group data is 
        groupSurfDir = fullfile(baseDir, experiment, wbDir, 'data', sprintf('group%dk', atlas_res));
        cd(groupSurfDir);
        
        for tt = 1:length(t.instOrder_all)
            for h = 1:2 % two hemispheres
                for s = 1:length(sn)
                    %%% make the group metric file for each contrasts
                    infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_transition_%d-%s.func.gii', subj_name{sn(s)}, hemI{h}, t.instOrder_all(tt), con_vs));
                    if smooth
                        surfFile    = fullfile(groupSurfDir,sprintf('fs_LR.32k.%s.inflated.surf.gii',hemI{h}));
                        surf_smooth(infilenames{s},'surf',surfFile,'kernel',kernel); % smooth outfilenames - it will prefix an 's'
                    end
                    s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_transition_%d-%s.func.gii', subj_name{sn(s)}, hemI{h}, t.instOrder_all(tt), con_vs));
                end % sn
                outfilenames    = fullfile(groupSurfDir,sprintf('%s.con_transition_%d-%s.func.gii', hemI{h}, t.instOrder_all(tt), con_vs));
                summaryname     = fullfile(groupSurfDir,sprintf('%s.group.con_transition_%d-%s.func.gii', hemI{h}, t.instOrder_all(tt), con_vs));
                
                surf_groupGiftis(infilenames, 'outfilenames', {outfilenames}, 'groupsummary', summaryname, 'replaceNaNs', replaceNaN);
                if smooth % also save the smoothed versions
                    s_outfilenames    = fullfile(groupSurfDir,sprintf('s%s.con_transition_%d-%s.func.gii', hemI{h}, t.instOrder_all(tt), con_vs));
                    s_summaryname     = fullfile(groupSurfDir,sprintf('s%s.group.con_transition_%d-%s.func.gii', hemI{h}, t.instOrder_all(tt), con_vs));
                    surf_groupGiftis(s_infilenames, 'outfilenames', {s_outfilenames}, 'groupsummary', s_summaryname, 'replaceNaNs', replaceNaN);
                end % if smooth
                
                % delet the smoothed contrasts for each subject
                for s = 1:length(sn)
                    delete(fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_transition_%d-%s.func.gii', subj_name{sn(s)}, hemI{h}, t.instOrder_all(tt), con_vs)));
                end
                
                fprintf('******************** group average contrast transition for %d vs %s is created! ********************\n\n', t.instOrder_all(tt), con_vs);
            end % hemi(h)
        end % tt (transitions)    
    case 'SURF:mdtb:noiseCeiling_get_utransitions_names_data'
        % Gets the data for calculating the noise ceiling
        % Uses the case for noise ceiling calculations
        % Creates surface gifti files for noise ceilings and saves them in
        % the workbench group directory
        % Example: sc1_sc2_mdtb('SURF:mdtb:noiseCeiling_get_utransitions_names_data')
        
        ppmethod   = '';     %% with or without stc
        experiment = 1;      %% sc1 or sc2;
        glm        = 7;
        atlas_res  = 32;     %% set it to 32 or 64
        smooth     = 1;      %% use the data with smoothing or without smoothing
        kernel     = 1;      %% the kernel used to smooth the data
        con_vs     = 'rest'; %% contrasts were calculated vs 'rest' or 'average'
        
        vararginoptions(varargin, {'ppmethod', 'experiment', 'glm', 'smooth', 'con_vs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                wbDir  = 'surfaceWB';
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                wbDir  = 'surfaceWB_stc';
        end
        
        % load in all_transMat (contains task_before and task_after names)
        load(fullfile(glmDir, 'all_trans_task_names.mat'));
        
        nTrans = size(all_trans, 1);
        
        % Load in data
        for h = 1:2 % two hemispheres
            % group gifti for each task transition
            for trans=1:nTrans
                switch smooth % use data with or without smoothing
                    case 1
                        T = gifti(fullfile(baseDir, experiment, wbDir, 'data',...
                            sprintf('group%dk', atlas_res), sprintf('s%s.con_transition_%s_%s-%s.func.gii', hemI{h}, all_trans{trans, 1}, all_trans{trans, 2}, con_vs)));
                    case 0
                        T = gifti(fullfile(baseDir, experiment, wbDir, 'data',...
                            sprintf('group%dk', atlas_res), sprintf('%s.con_transition_%s_%s-%s.func.gii', hemI{h}, all_trans{trans, 1}, all_trans{trans, 2}, con_vs)));

                end % switch smooth
                groupData(:,:,trans) = T.cdata;
            end % trans
            
            % calculate low and high noise ceilings for hemI{h}
            [noise_low, noise_high] = sc1_sc2_mdtb('SURF:mdtb:noiseCeiling_calculate', groupData);
            
            % creates surface gifti files for the noise ceilings
            G_high = surf_makeFuncGifti(noise_high, 'anatomicalStruct', hemName{h}, 'columnNames', {'noiseCeiling_high'});
            G_low  = surf_makeFuncGifti(noise_low, 'anatomicalStruct', hemName{h}, 'columnNames', {'noiseCeiling_low'});
            
            % save the gifti files
            switch smooth
                case 1 % with smoothing 
                    save(G_high, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('s%s.transition_names-%s_noiseCeiling_high.func.gii', hemI{h}, con_vs)));
                    save(G_low, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('s%s.transition_names-%s_noiseCeiling_low.func.gii', hemI{h}, con_vs)));
                case 0 % wihtout smoothing
                    save(G_high, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('%s.transition_names-%s_noiseCeiling_high.func.gii', hemI{h}, con_vs)));
                    save(G_low, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('%s.transition_names-%s_noiseCeiling_low.func.gii', hemI{h}, con_vs)));
            end % switch smooth
            
            fprintf('******************** low and high noise ceilings created for %s ********************', hemI{h});
        end % h
    case 'SURF:mdtb:noiseCeiling_get_transitions_id_data'
        % Gets the data for calculating the noise ceiling
        % Uses the case for noise ceiling calculations
        % Creates surface gifti files for noise ceilings and saves them in
        % the workbench group directory
        % Example: sc1_sc2_mdtb('SURF:mdtb:noiseCeiling_get_transitions_id_data')
        
        ppmethod   = '';     %% with or without stc
        experiment = 1;      %% sc1 or sc2;
        glm        = 7;
        atlas_res  = 32;     %% set it to 32 or 64
        con_vs     = 'rest'; %% contrast vs 'rest' or 'average'
        smooth     = 1;      %% use the smoothed data?
        kernel     = 1;      %% smoothing kernel used
        
        vararginoptions(varargin, {'ppmethod', 'experiment', 'glm', 'atlas_res', 'con_vs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                wbDir  = 'surfaceWB';
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                wbDir  = 'surfaceWB_stc';
        end
        
        % load in all_transMat (contains task_before and task_after names)
        load(fullfile(glmDir, 'trans_info.mat'));
        
        nTrans = length(t.instOrder_all);
        
        % Load in data
        for h = 1:2 % two hemispheres
            % group gifti for each task transition
            for trans=1:nTrans
                switch smooth % use the data with or without smoothing
                    case 1 % with smoothing
                        T = gifti(fullfile(baseDir, experiment, wbDir, 'data', sprintf('group%dk', atlas_res),...
                            sprintf('s%s.con_transition_%d-%s.func.gii', hemI{h}, t.instOrder_all(trans), con_vs)));
                    case 0 % without smoothing
                        T = gifti(fullfile(baseDir, experiment, wbDir, 'data', sprintf('group%dk', atlas_res),...
                            sprintf('%s.con_transition_%d-%s.func.gii', hemI{h}, t.instOrder_all(trans), con_vs)));
                end % switch smooth
                groupData(:,:,trans) = T.cdata;
            end % trans
            
            % calculate low and high noise ceilings for hemI{h}
            [noise_low, noise_high] = sc1_sc2_mdtb('SURF:mdtb:noiseCeiling_calculate', groupData);
            
            % creates surface gifti files for the noise ceilings
            G_high = surf_makeFuncGifti(noise_high, 'anatomicalStruct', hemName{h}, 'columnNames', {'noiseCeiling_high'});
            G_low  = surf_makeFuncGifti(noise_low, 'anatomicalStruct', hemName{h}, 'columnNames', {'noiseCeiling_low'});
            
            % save the gifti files
            switch smooth % save giftis for data with or without smoothing
                case 1 % with smoothing
                    save(G_high, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('s%s.transition_id-%s_noiseCeiling_high.func.gii', hemI{h}, con_vs)));
                    save(G_low, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('s%s.transition_id-%s_noiseCeiling_low.func.gii', hemI{h}, con_vs)));
                case 0 % without smoothing
                    save(G_high, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('%s.transition_id-%s_noiseCeiling_high.func.gii', hemI{h}, con_vs)));
                    save(G_low, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('%s.transition_id-%s_noiseCeiling_low.func.gii', hemI{h}, con_vs)));
            end % switch smooth
            fprintf('******************** low and high noise ceilings created for %s hemi ********************\n\n', hemI{h});
        end % h
    case 'SURF:mdtb:noiseCeiling_calculate'
        % Calculates and returns the noise ceilings
        % Example: sc1_sc2_mdtb('SURF:mdtb:noise_ceiling', groupData)
        
        Data = varargin{1}; %% group data for calculation of noise ceilings
        
        nNode  = size(Data,1);
        nSubj  = size(Data, 2);
        
        % preallocating the corr mats
        r_low  = zeros(nNode,nSubj);
        r_high = zeros(nNode,nSubj);

        % perform noise ceilings
        for s = 1:nSubj
            data1      = squeeze(Data(:, s, :));
            data2_low  = squeeze(nanmean(Data(:, ~ismember(1:nSubj, s), :), 2));
            data2_high = squeeze(nanmean(Data(:, :, :), 2));
%             % remove means
%             d1 = bsxfun(@minus,data1,mean(data1,2));
%             d2_l = bsxfun(@minus,data2_low,mean(data2_low,2));
%             d2_h = bsxfun(@minus,data2_high,mean(data2_high,2));
%             Pearson correlation formula
%             var1        = d1.^2;
%             var2_low    = d2_l.^2;
%             var2_high   = d2_h.^2;
%             cv_l        = var1.*var2_low;
%             cv_h        = var1.*var2_high;
%             r_low(:,i)  = sum(cv_l,2)./sqrt(sum(var1,2).*sum(var2_low,2));
%             r_high(:,i) = sum(cv_h,2)./sqrt(sum(var1,2).*sum(var2_high,2)); 
%             
            % This is not optimum!    
            for n = 1:nNode
                r_low(n, s)  = corr(data1(n, :)', data2_low(n, :)');
                r_high(n, s) = corr(data1(n, :)', data2_high(n, :)');
            end % n
        end % s (subjects)

        % overall noise ceiling across subjects
        noise_low   = nanmean(r_low,2);
        noise_high  = nanmean(r_high,2);
        
        % plot
        figure
        subplot(121);
        histogram(noise_low);
        subplot(122);
        histogram(noise_high);
        
        % output noise_low, noise_high
        varargout{1} = noise_low;
        varargout{2} = noise_high;
    case 'SURF:contrasts_single_gifti'
        % takes all the gifti files for the contrasts for each subject and
        % merge them all together in a single file. It can do it both on
        % individual level and group level.
        % Example: sc1_sc2_mdtb_backup('SURF:contrasts_single_gifti', 'sn', 3)
        
        sn          = returnSubjs;    %%
        glm         = 8;              %%
        experiment  = 1;              %%
        ppmethod    = '';             %%
        which       = 'task';         %%
        igroup      = 0;              %% do it for the group average files or for the individual subjects?
        con_vs      = 'average_2'; %% choose 'average_1', 'average_2', or 'rest', or 'rest_taskCon' (if glm7 and task)
%         atlas_res   = 32;           %% atlas resolution can either be 32 or 64
        replaceNaNs = 1;              %% set to 1 or 0 
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment', 'ppmethod', 'which', 'con_vs', 'atlas_res', 'replaceNaNs'})
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment);
        switch which
            case 'task' % task for glm8
                conNames = unique(Cc.taskNames);
            case 'cond' % condition for glm7
                conNames = unique(Cc.condNames);
        end %% do you want the group maps for tasks or conditions
        
        experiment = sprintf('sc%d', experiment);
        
        % in 'sc1_sc2_taskConds.txt' file, instruct is not coded as a
        % task/condition name. So I will have to add that to the list of
        % names
        conNames = ['Instruct'; conNames];
        
        switch ppmethod
            case ''
                wbDir = fullfile(baseDir, experiment, 'surfaceWB');
            case 'stc'
                wbDir = fullfile(baseDir, experiment, 'surfaceWB_stc');
        end
        
        dircheck(wbDir);

        switch igroup % do it for group data or individual data
            case 1 % do it for the group data which is already smoothed?
            case 0 % do it for each subject separately
                for s = sn
                    dircheck(fullfile(wbDir, sprintf('glm%d', glm), subj_name{s}))
%                     subWbDir = fullfile(wbDir, subj_name{s});
                    % group all the contrast surface maps into a single
                    % file
                    %%% creating a single file for each hemisphere
                    for h = 1:2
                        for cc = 1:length(conNames)
                            infilenames{cc}   = fullfile(wbDir, sprintf('glm%d', glm), subj_name{s},sprintf('%s.%s.con_%s-%s.func.gii',...
                                subj_name{s}, hemI{h}, conNames{cc}, con_vs));
                            columnName{cc} = sprintf('%s-%s', conNames{cc}, con_vs);
                        end % cc (condition)
                        cd(fullfile(wbDir, sprintf('glm%d', glm), subj_name{s}));
                        outfilename = sprintf('%s.%s.con_%s-%s.func.gii', subj_name{s}, hemI{h}, which, con_vs);
                        surf_groupGiftis(infilenames, 'outfilenames', {outfilename}, 'outcolnames', columnName, 'replaceNaNs', replaceNaNs);
                        fprintf('a single gifti file for contrasts for %s hemi successfully created for %s\n', hemI{h}, subj_name{s})
                    end % h (hemi)  
                end % s (sn)
        end % switch igroup
        
    case 'SUIT:mdtb:suit_parcel2native'
        % maps each atlas of the suit into individual space
        % Example: sc1_sc2_mdtb('SUIT:mdtb:suit_parcel2native', 'sn', [3])
        
        sn         = returnSubjs;
        parcelType = 'Buckner_7';                         %% set it to Buckner_7 or Buckner_17
        parcelDir  = fullfile(suitToolDir, 'atlasesSUIT'); %% directory where the nifti image is stored
        experiment = 1;
        
        vararginoptions(varargin, {'sn', 'parcelType', 'parcelDir', 'experiment', 'ppmethod'});
        
        experiment = sprintf('sc%d', experiment);
                
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
        experiment = 1;                      %% enter 1 for sc1 and 2 for sc2
        type       = 'con';                  %% enter the image you want to reslice to suit space
        glm        = 7;                      %% glm number
        mask       = 'cereb_prob_corr_grey'; %% the cerebellar mask to be used:'cereb_prob_corr_grey' or 'cereb_prob_corr' or 'dentate_mask'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'experiment', 'glm', 'type', 'mask'});
        
        experiment = sprintf('sc%d', experiment);
        
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
        experiment = 1;                      %% enter 1 for sc1 and 2 for sc2
        type       = 'con';                  %% enter the image you want to reslice to suit space
        glm        = 7;                      %% glm number
        con_vs     = 'rest';                 %% is the contrast calculated vs 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'experiment', 'glm', 'type', 'mask'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc       = getrow(C, C.StudyNum == experiment);
        conNames = unique(Cc.condNames);
        
        experiment = sprintf('sc%d', experiment);
        
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
        experiment = 1;                      %% enter 1 for sc1 and 2 for sc2
        type       = 'con';                  %% enter the image you want to reslice to suit space
        glm        = 7;                      %% glm number
        con_vs     = 'rest';                 %% is the contrast calculated vs 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'experiment', 'glm', 'type', 'mask'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc       = getrow(C, C.StudyNum == experiment);
        taskNames = unique(Cc.taskNames);
        
        experiment = sprintf('sc%d', experiment);
        
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
    case 'SUIT:mdtb:groupmap_con_utransitions_names'
        % creates group map for unique transitions with task names
        % Example: sc1_sc2_mdtb('SUIT:mdtb:groupmap_con_utransitions_names', 'sn', [3])
        
        sn         = returnSubjs;   %% list of subjects
        ppmethod   = '';     %% with or without stc
        experiment = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 7;      %% glm number
        con_vs     = 'rest'; %% contrast was calculated against 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment', 'glm', 'replaceNaN', 'convs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                glmDir          = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm)); % where the all_trans matrix is stored
                glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
                glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), 'group');
            case 'stc'
                glmDir          = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm)); % where the all_trans matrix is stored
                glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d_stc', glm));
                glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d_stc', glm), 'group');
        end
        dircheck(glmSuitDir);
        dircheck(glmSuitGroupDir);
        
        % load in all_trans.mat for the task names before and after
        load(fullfile(glmDir, 'all_trans_task_names.mat'));
        
        % preallocating!
        maps = []; %% the structure that will have alll?! the info I will need (hopefully :))
        for tt = 1:size(all_trans, 1)
            for s = 1:length(sn)
                infilename{s} = fullfile(glmSuitDir, subj_name{sn(s)}, sprintf('wdcon_transition_%s_%s-%s', all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                
                % load in individual cerebellar maps
                V_tmp = spm_vol(infilename{s});
                X     = spm_read_vols(V_tmp);
                X_vec = X(:);
                
                maps_tmp.conNames           = cellstr(sprintf('%s_%s', all_trans{tt, 1}, all_trans{tt, 2}));
                maps_tmp.conBaseName        = cellstr(con_vs);
                maps_tmp.subjectMaps{1}     = X;
                maps_tmp.subjectMapsVec     = X_vec';
                maps_tmp.subjectMapNames{1} = infilename{s};
                maps_tmp.glm                = 7;
                maps_tmp.sn                 = sn(s);
                maps_tmp.vol                = V_tmp;
                
                maps  = addstruct(maps, maps_tmp);
            end % sn
            outfilename = fullfile(glmSuitGroupDir, sprintf('group_con_transition_%s_%s-%s.nii', all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
            opt.dmtx    = 1;
            % calculate the average across subjects
            spm_imcalc(infilename', outfilename, 'nanmean(X)', opt);
            
            % saving the gifti file for group map
            V = spm_vol(fullfile(glmSuitGroupDir, sprintf('group_con_transition_%s_%s-%s.nii', all_trans{tt, 1}, all_trans{tt, 2}, con_vs)));
            D = suit_map2surf(V,'stats','nanmean');
            
            G = surf_makeFuncGifti(D, 'anatomicalStruct', 'Cerebellum', 'columnNames', {sprintf('group_con_transition_%s_%s-%s', all_trans{tt, 1}, all_trans{tt, 2}, con_vs)});
            
            save(G, fullfile(glmSuitGroupDir, sprintf('Cereb.group.con_transition_%s_%s-%s.func.gii', all_trans{tt, 1}, all_trans{tt, 2}, con_vs)));
            save(fullfile(glmSuitGroupDir, sprintf('indMaps_con_transitions_%s_%s-vs-%s.mat', all_trans{tt, 1}, all_trans{tt, 2}, con_vs)), 'maps', '-v7.3');
            fprintf('******************** con group average for transitions from %s to %s vs %s is created! ********************\n\n', all_trans{tt, 1}, all_trans{tt, 2}, con_vs);
        end % tt 
    case 'SUIT:mdtb:groupmap_con_transitions_id'
        % creates group map for all the transitions (256 for sc1)
        % creates group map for unique transitions with task names
        % Example: sc1_sc2_mdtb('SUIT:mdtb:groupmap_con_transitions_id', 'sn', [3])
        
        sn         = returnSubjs;   %% list of subjects
        ppmethod   = '';     %% with or without stc
        experiment = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 7;      %% glm number
        con_vs     = 'rest'; %% contrast was calculated against 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'ppmethod', 'atlas_res', 'experiment', 'glm', 'replaceNaN', 'convs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                glmDir          = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm)); % where the all_trans matrix is stored
                glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
                glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), 'group');
            case 'stc'
                glmDir          = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm)); % where the all_trans matrix is stored
                glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d_stc', glm));
                glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d_stc', glm), 'group');
        end
        dircheck(glmSuitDir);
        dircheck(glmSuitGroupDir);
        
        % load in all_trans.mat for the task names before and after
        load(fullfile(glmDir, 'trans_info.mat'));
        
        % preallocating!
        for tt = 1:length(t.instOrder_all),             
            maps = []; %% the structure that will have alll?! the info I will need (hopefully :))
            for s = 1:length(sn)
                infilename{s} = fullfile(glmSuitDir, subj_name{sn(s)}, sprintf('wdcon_transition_%d-%s.nii', t.instOrder_all(tt), con_vs));
                
                % load in individual cerebellar maps
                V_tmp = spm_vol(infilename{s});
                X     = spm_read_vols(V_tmp);
                
                X_vec = X(:);
                
                maps_tmp.conNames           = cellstr(sprintf('transition_%d', t.instOrder_all(tt)));
                maps_tmp.conBaseName        = cellstr(con_vs);
                maps_tmp.subjectMaps{1}     = X;
                maps_tmp.subjectMapsVec     = X_vec';
                maps_tmp.subjectMapNames{1} = infilename{s};
                maps_tmp.glm                = 7;
                maps_tmp.sn                 = sn(s);
                maps_tmp.vol                = V_tmp;
                
                maps  = addstruct(maps, maps_tmp);
            end % sn            
            outfilename = fullfile(glmSuitGroupDir, sprintf('group_con_transition_%d-%s.nii', t.instOrder_all(tt), con_vs));
            opt.dmtx    = 1;
            % calculate the average across subjects
            spm_imcalc(infilename', outfilename, 'nanmean(X)', opt);
            
            % saving the gifti file for group map
            % saving the gifti file for group map
            V = spm_vol(fullfile(glmSuitGroupDir, sprintf('group_con_transition_%d-%s.nii', t.instOrder_all(tt), con_vs)));
            D = suit_map2surf(V,'stats','nanmean');
            
            G = surf_makeFuncGifti(D, 'anatomicalStruct', 'Cerebellum', 'columnNames', {sprintf('group_con_transition_%d-%s', t.instOrder_all(tt), con_vs)});
            
            save(G, fullfile(glmSuitGroupDir, sprintf('Cereb.group.con_transition_%d-%s.func.gii', t.instOrder_all(tt), con_vs)));
            
            save(fullfile(glmSuitGroupDir, sprintf('indMaps_con_transitions_id_%d-vs-%s.mat', t.instOrder_all(tt), con_vs)), 'maps', '-v7.3');
            
            fprintf('******************** con group average for transition %d vs %s is created! ********************\n\n', t.instOrder_all(tt), con_vs);
        end % tt    
    case 'SUIT:mdtb:noiseCeiling_get_transitions_id_data'
        % calculating noise ceilings for the cerebellum
        % Example: sc1_sc2_mdtb('SUIT:mdtb:noiseCeiling_get_transitions_id_data')
        
        ppmethod   = '';     %% with or without stc
        experiment = 1;      %% sc1 or sc2;
        glm        = 7;
        con_vs     = 'rest'; %% contrast vs 'rest' or 'average'
%         smooth     = 1;      %% use the smoothed data?
%         kernel     = 1;      %% smoothing kernel used
        
        vararginoptions(varargin, {'ppmethod', 'experiment', 'glm', 'atlas_res', 'con_vs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                suitGroupDir  = fullfile(baseDir, experiment, 'suit', sprintf('glm%d', glm), 'group');
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                suitGroupDir  = fullfile(baseDir, experiment, 'suit_stc', sprintf('glm%d', glm), 'group');
        end
        
        % load in all_transMat (contains task_before and task_after names)
        load(fullfile(glmDir, 'trans_info.mat'));
        
        nTrans = length(t.instOrder_all);
        
        % Load in data
        % suitmaps for each task transition
        for trans=1:nTrans
            fprintf('load in group data for transition %0.3d\n', trans)
            load(fullfile(suitGroupDir, sprintf('indMaps_con_transitions_id_%d-vs-%s', trans, con_vs)));
            
%             tmp = maps.subjectMaps; % this is a temp variable that will have the data for all the subjects.
%             for s = 1:length(maps.sn)
%                 groupData(:, :, :, s, trans) = tmp{s};
%             end
            groupData2(:, :, trans) = maps.subjectMapsVec';
        end % trans
        
        % calculate low and high noise ceilings for hemI{h}
%         [noise_low, noise_high] = sc1_sc2_mdtb('SUIT:mdtb:noiseCeiling_calculate', groupData);
        [noise_low, noise_high] = sc1_sc2_mdtb('SURF:mdtb:noiseCeiling_calculate', groupData2);
        
        myTmp     = maps.vol;
        V_low     = myTmp(1);
        V_low     = rmfield(V_low, 'fname');
        V_low.dat = noise_low;
        
        D_low = suit_map2surf(V_low,'stats','nanmean');
        
        V_high     = myTmp(1);
        V_high     = rmfield(V_high, 'fname');
        V_high.dat = noise_high;
        
        D_high = suit_map2surf(V_high,'stats','nanmean');
        
        
        % creates surface gifti files for the noise ceilings
        G_high = surf_makeFuncGifti(D_high, 'anatomicalStruct', 'Cerebellum', 'columnNames', {'noiseCeiling_high'});
        G_low  = surf_makeFuncGifti(D_low, 'anatomicalStruct', 'Cerebellum', 'columnNames', {'noiseCeiling_low'});
        
        % save the gifti files
        switch smooth % save giftis for data with or without smoothing
            case 0 % without smoothing
                save(G_high, fullfile(suitGroupDir, sprintf('Cereb.transition_id-%s_noiseCeiling_high.func.gii', con_vs)));
                save(G_low, fullfile(suitGroupDir, sprintf('Cereb.transition_id-%s_noiseCeiling_low.func.gii', con_vs)));
        end % switch smooth
        fprintf('******************** low and high noise ceilings created for %s hemi ********************\n\n', hemI{h});
    case 'SUIT:mdtb:noiseCeiling_calculate'
        % Calculates and returns the noise ceilings
        % Example: sc1_sc2_mdtb('SURF:mdtb:noise_ceiling', groupData)
        
        Data = varargin{1}; %% group data for calculation of noise ceilings
        
%         [nx, ny, nz]  = size(Data{});
        nSubj  = size(Data, 4);
        nx = size(Data, 1);
        ny = size(Data, 2);
        nz = size(Data, 3);
        
        % preallocating the corr mats
        r_low  = zeros(nx, ny, nz, nSubj);
        r_high = zeros(nx, ny, nz, nSubj);

        % perform noise ceilings
        for s = 1:nSubj
            data1      = squeeze(Data(:, :, :, s, :));
            data2_low  = squeeze(nanmean(Data(:, :, :, ~ismember(1:nSubj, s), :), 4));
            data2_high = squeeze(nanmean(Data(:, :, :, :, :), 4));
%             % remove means
%             d1 = bsxfun(@minus,data1,mean(data1,2));
%             d2_l = bsxfun(@minus,data2_low,mean(data2_low,2));
%             d2_h = bsxfun(@minus,data2_high,mean(data2_high,2));
%             Pearson correlation formula
%             var1        = d1.^2;
%             var2_low    = d2_l.^2;
%             var2_high   = d2_h.^2;
%             cv_l        = var1.*var2_low;
%             cv_h        = var1.*var2_high;
%             r_low(:,i)  = sum(cv_l,2)./sqrt(sum(var1,2).*sum(var2_low,2));
%             r_high(:,i) = sum(cv_h,2)./sqrt(sum(var1,2).*sum(var2_high,2)); 
%             
            % This is not optimum!    
            for ix = 1:nx
                for iy = 1:ny
                    for iz = 1:nz
                        r_low(ix, iy, iz, s)  = corr(squeeze(data1(ix, iy, iz, :)), squeeze(data2_low(ix, iy, iz, :)));
                        r_high(ix, iy, iz, s) = corr(squeeze(data1(ix, iy, iz, :)), squeeze(data2_high(ix, iy, iz, :)));
                    end % iz
                end % iy
            end % ix
        end % s (subjects)

        % overall noise ceiling across subjects
        noise_low   = nanmean(r_low,4);
        noise_high  = nanmean(r_high,4);
        
        % plot
        figure
        subplot(121);
        histogram(noise_low(:));
        subplot(122);
        histogram(noise_high(:));
        
        % output noise_low, noise_high
        varargout{1} = noise_low;
        varargout{2} = noise_high;
        
    case 'SUIT:cole_parcellation'
        % this case is used to transfer the cerebellar parcellation in the
        % Cole study to SUIT space.
        % Example: sc1_sc2_mdtb('SUIT:cole_parcellation') 
        
    case 'LAG:mdtb:lccf'
        %%% calculating the lagged cross covariance function for the time
        %%% series averaged across runs.
        % Example: sc1_sc2_mdtb('LAG:mdtb:lccf', 'sn', [3], 'roi1', 'cerebellum_grey', 'roi2', 'cerebellum_grey', 'tsType', 'raw')
        sn         = returnSubjs;       %% subject list
        ppmethod   = '';         %% Did you do the slice timing correction as a preprocessing step?
        experiment = 1;
        glm        = 7;          %% the GLM number
        condition  = [];         %% the regressors of interest used in getting the predicted time series      
        tsType     = 'filtered';      %% set this option to raw or predicted
        maxShift   = 4;          %% maximum shift you want to use to calculate the lagged cross covariance funtion
        
        roi1 = '162_tessellation_hem';
        roi2 = 'cerebellum_grey';      %% roi2 can be the same as roi1 if you want to calculate time delays in one roi 
        
        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'glm', 'roi1', 'roi2', 'maxShift', 'tsType', 'condition'});
        
        if isempty(condition)
            taskName = 'all';
        else
            % load in task information
            C   = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
            Cc  = getrow(C,C.StudyNum == experiment);
            Ccc = getrow(Cc, Cc.condNum == condition(1));

            taskName = Ccc.taskNames;
        end
        
        experiment = sprintf('sc%d', experiment); %% experiment number is converted to 'sc1' or 'sc2'
        
        switch ppmethod
            case 'stc' 
                tsDir  = fullfile(baseDir, experiment, sprintf('TS_GLM_%d_stc', glm));
                lccDir = fullfile(baseDir, experiment, sprintf('LCC_GLM_%d_stc', glm));
            case ''
                tsDir  = fullfile(baseDir, experiment, sprintf('TS_GLM_%d', glm));
                lccDir = fullfile(baseDir, experiment, sprintf('LCC_GLM_%d', glm));
        end
        dircheck(lccDir)
        
        lcc_allsubs = cell(1, length(sn));
        for s = sn
            fprintf('******************** calculating lagged cross covariance for %s ********************\n', subj_name{s});
            % load in the average time series across runs
            roi1_ts = load(fullfile(tsDir, subj_name{s}, sprintf('avgTs_%s_regions_%s.mat', taskName, roi1)));
            roi1_ts = roi1_ts.avgts;
            roi2_ts = load(fullfile(tsDir, subj_name{s}, sprintf('avgTs_%s_regions_%s.mat', taskName, roi2)));
            roi2_ts = roi2_ts.avgts;
            
            switch tsType %% do the analysis with the raw time series or predicted time series?
                case 'raw'
                    roi1_ts = roi1_ts.y_raw;
                    roi2_ts = roi2_ts.y_raw;
                case 'predicted'
                    roi1_ts = roi1_ts.y_hat;
                    roi2_ts = roi2_ts.y_hat;
                case 'filtered'
                    roi1_ts = roi1_ts.y_filt;
                    roi2_ts = roi2_ts.y_filt;
            end
            
            [lcc, lags] = lcc_vectorized(roi1_ts, roi2_ts, maxShift);
            lcc_allsubs{s} = lcc;
            
            lcc_struct.lcc  = lcc;
            lcc_struct.lags = lags;
            
            lcc_struct.ppmethod   = ppmethod;
            lcc_struct.glm        = glm;
            lcc_struct.roi1       = roi1;
            lcc_struct.roi2       = roi2;
            lcc_struct.regressors = condition;
            lcc_struct.tsType     = tsType;
            
            dircheck(fullfile(lccDir, subj_name{s}));
            save(fullfile(lccDir, subj_name{s}, sprintf('lcc_%s_regions_%s_regions_%s_%s.mat', taskName, roi1, roi2, tsType)), 'lcc_struct', '-v7.3');
            fprintf('******************** lagged cross covariance calculated for %s ********************\n\n', subj_name{s});    
        end
        varargout{1} = lcc_allsubs;
    case 'LAG:mdtb:td'
        %%% estimate the time delays using parabolic interpolation
        % Example: sc1_sc2_mdtb('LAG:mdtb:td', 'sn', [3])
        sn         = returnSubjs;       %% subject list
        ppmethod   = '';         %% Did you do the slice timing correction as a preprocessing step?
        experiment = 1;
        glm        = 7;          %% the GLM number
        condition  = [];         %% the regressors of interest used in getting the predicted time series      
        tsType     = 'raw';      %% set this option to raw or predicted
        maxLag     = 5;          %% Lag values higher than this will be set to NaN (noise) 
        Tr         = 1;          %% Tr value is set to 1 sec as the default value
        
        roi1 = '162_tessellation_hem';
        roi2 = 'cerebellum_grey';      %% roi2 can be the same as roi1 if you want to calculate time delays in one roi 
        
        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'glm', 'roi1', 'roi2', 'maxLag', 'tsType', 'Tr', 'condition'});
        
        if isempty(condition)
            taskName = 'all';
        else
            % load in task information
            C  = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
            Cc = getrow(C,C.StudyNum == experiment);
            Ccc = getrow(Cc, Cc.condNum == condition(1));

            taskName = Ccc.taskNames;
        end
        
        experiment = sprintf('sc%d', experiment); %% experiment number is converted to 'sc1' or 'sc2'
        
        switch ppmethod
            case 'stc' 
                lccDir  = fullfile(baseDir, experiment, sprintf('LCC_GLM_%d_stc', glm));
                tdDir = fullfile(baseDir, experiment, sprintf('TD_GLM_%d_stc', glm));
            case ''
                lccDir  = fullfile(baseDir, experiment, sprintf('LCC_GLM_%d', glm));
                tdDir = fullfile(baseDir, experiment, sprintf('TD_GLM_%d', glm));
        end
        dircheck(tdDir)
        
        td_allsubs = cell(1, length(sn));
        for s = sn
            fprintf('******************** Estimate TD matrix for %s ********************\n', subj_name{s});
            % load in the lcc structure variable
            lcc_tmp = load(fullfile(lccDir, subj_name{s}, sprintf('lcc_%s_regions_%s_regions_%s_%s.mat', taskName, roi1, roi2, tsType)));
            lcc_tmp = lcc_tmp.lcc_struct;
            
            % estimate the time delay using parabolic interpolation
            [Td, Pc, P, I] = td_estimate(lcc_tmp.lcc, lcc_tmp.lags, maxLag, Tr);
            
            Pcz = fisher_z(Pc);
            Pz  = fisher_z(P);
            
            td_struct.td  = Td;
            td_struct.pc  = Pcz;
            td_struct.p   = Pz;
            td_struct.ind = I;
            
            td_struct.ppmethod   = ppmethod;
            td_struct.glm        = glm;
            td_struct.roi1       = roi1;
            td_struct.roi2       = roi2;
            td_struct.regressors = condition;
            td_struct.tsType     = tsType;
            
            td_allsubs{s} = td_struct;
            
            dircheck(fullfile(tdDir, subj_name{s}));
            save(fullfile(tdDir, subj_name{s}, sprintf('td_%s_regions_%s_regions_%s_%s.mat', taskName, roi1, roi2, tsType)), 'td_struct', '-v7.3');
            fprintf('******************** TD estimated for %s ********************\n\n', subj_name{s});
        end
        varargout{1} = td_allsubs;
    case 'LAG:mdtb:lp'
        % calculating the lag projection map 
        % Example: sc1_sc2_mdtb('LAG:mdtb:lp', 'sn', [3])
        sn         = returnSubjs;       %% subject list
        ppmethod   = '';         %% Did you do the slice timing correction as a preprocessing step?
        experiment = 1;
        glm        = 7;          %% the GLM number
        condition  = [];         %% the regressors of interest used in getting the predicted time series      
        tsType     = 'raw';      %% set this option to raw or predicted
        lpmethod   = 'max';      %% the method used to project the TD matrix into 1D. Options are: 'max', 'marek', 'percentile'
        dim        = 1;          %% the dimension along which the TD will be projected
        
        roi1 = '162_tessellation_hem';
        roi2 = 'cerebellum_grey';      %% roi2 can be the same as roi1 if you want to calculate time delays in one roi 
        
        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'glm', 'roi1', 'roi2', 'maxLag', 'tsType', 'Tr', 'condition', 'lpmethod'});
        
        if isempty(condition)
            taskName = 'all';
        else
            % load in task information
            C  = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
            Cc = getrow(C,C.StudyNum == experiment);
            Ccc = getrow(Cc, Cc.condNum == condition(1));

            taskName = Ccc.taskNames;
        end
        
        experiment = sprintf('sc%d', experiment); %% experiment number is converted to 'sc1' or 'sc2'
        
        switch ppmethod
            case 'stc' 
                lpDir  = fullfile(baseDir, experiment, sprintf('LP_GLM_%d_stc', glm));
                tdDir = fullfile(baseDir, experiment, sprintf('TD_GLM_%d_stc', glm));
            case ''
                lpDir  = fullfile(baseDir, experiment, sprintf('LP_GLM_%d', glm));
                tdDir = fullfile(baseDir, experiment, sprintf('TD_GLM_%d', glm));
        end
        dircheck(lpDir)
        
        lp_allsubs = cell(1, length(sn));
        for s = sn
            fprintf('******************** Calculating Lp(%s) for %s ********************\n', lpmethod, subj_name{s});
            % load in td_struct
            td_struct_tmp = load(fullfile(tdDir, subj_name{s}, sprintf('td_%s_regions_%s_regions_%s_%s.mat', taskName, roi1, roi2, tsType)));
            td_struct_tmp = td_struct_tmp.td_struct;
            
            % calculate the projection
            [ lp, pp ] = lag_projection(td_struct_tmp, dim, lpmethod);
            
            lp_struct.lp  = lp;
            lp_struct.pp  = pp;
            
            lp_struct.lpmethod   = lpmethod;
            lp_struct.ppmethod   = ppmethod;
            lp_struct.glm        = glm;
            lp_struct.roi1       = roi1;
            lp_struct.roi2       = roi2;
            lp_struct.regressors = condition;
            lp_struct.tsType     = tsType;
            
            lp_allsubs{s} = lp_struct;
            
            dircheck(fullfile(lpDir, subj_name{s}));
            save(fullfile(lpDir, subj_name{s}, sprintf('lp_%s_regions_%s_regions_%s_%s_%s.mat', taskName, roi1, roi2, tsType, lpmethod)), 'lp_struct', '-v7.3');
            
            fprintf('******************** Lp(%s) calculated for %s ********************\n\n', lpmethod, subj_name{s});
        end
        varargout{1} = lp_allsubs;
    case 'LAG:mdtb:suit:map_cerebellum'
        % creates the flatmap of lag projection for the cerebellum 
        % Example:sc1_sc2_mdtb('LAG:mdtb:suit:map_cerebellum', 'sn', [3])
        
        sn         = returnSubjs;       %% subject list
        ppmethod   = '';         %% Did you do the slice timing correction as a preprocessing step?
        experiment = 1;
        glm        = 7;          %% the GLM number
        condition  = [];         %% the regressors of interest used in getting the predicted time series      
        tsType     = 'raw';      %% set this option to raw or predicted
        lpmethod   = 'max';      %% the method used to project the TD matrix into 1D. Options are: 'max', 'marek', 'percentile'
        
        roi1 = '162_tessellation_hem';
        roi2 = 'cerebellum_grey';      %% roi2 can be the same as roi1 if you want to calculate time delays in one roi 
        
        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'glm', 'roi1', 'roi2', 'maxLag', 'tsType', 'Tr', 'condition', 'lpmethod'});
        
        if isempty(condition)
            taskName = 'all';
        else
            % load in task information
            C  = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
            Cc = getrow(C,C.StudyNum == experiment);
            Ccc = getrow(Cc, Cc.condNum == condition(1));

            taskName = Ccc.taskNames;
        end
        
        experiment = sprintf('sc%d', experiment); %% experiment number is converted to 'sc1' or 'sc2'
        
        switch ppmethod
            case 'stc' 
                lpDir  = fullfile(baseDir, experiment, sprintf('LP_GLM_%d_stc', glm));
                fmDir = fullfile(baseDir, experiment, sprintf('FM_GLM_%d_stc', glm));
            case ''
                lpDir  = fullfile(baseDir, experiment, sprintf('LP_GLM_%d', glm));
                fmDir = fullfile(baseDir, experiment, sprintf('FM_GLM_%d', glm));
        end
        dircheck(fmDir)
        
        fm_allsubs = cell(1, length(sn));
        for s = sn
            % load in the lp structure variable
            fprintf('******************** Creating the input to suit_plotflatmap for %s ********************\n\n', subj_name{s});
            
            % Determine the voxels we want to resample in SUIT space
            V = spm_vol(fullfile(baseDir, experiment,suitDir,'anatomicals','cerebellarGreySUIT.nii'));
            X = spm_read_vols(V);

            grey_threshold = 0.1; % gray matter threshold

            linIn1     = find(X > grey_threshold);
            [i1,j1,k1] = ind2sub(V.dim,linIn1');
            [x1,y1,z1] = spmj_affine_transform(i1, j1, k1, V.mat);
            
            
            %%% load in the lp and pp
            load(fullfile(lpDir, subj_name{s}, sprintf('lp_%s_regions_%s_regions_%s_%s_%s.mat', taskName, roi1, roi2, tsType, lpmethod)));
            data_map_lp = lp_struct.lp; %% I want to creat a flatmap for the Lp I got using the average of lccf across runs.
            data_map_pp = lp_struct.pp; %% I want to creat a flatmap for the Lp I got using the average of lccf across runs. 
            
            %%% Determine voxel locations from the original ROI
            load(fullfile(baseDir, experiment, regDir, 'data', subj_name{s}, sprintf('regions_%s.mat', roi2))); % 'regions' are defined in 'ROI_define'

            Vmask      = spm_vol(fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'maskbrainSUITGrey.nii')); %% cerebellum grey matter mask
            [i3,j3,k3] = spmj_affine_transform(R{1}.data(:,1),R{1}.data(:,2),R{1}.data(:,3),inv(Vmask.mat));
            linIn3     = sub2ind(Vmask.dim, round(i3), round(j3), round(k3));

            %%% transform SUIT coords into anatomical space of the individual
            flowfield    = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'u_a_c_anatomical_seg1.nii');
            affine       = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'Affine_c_anatomical_seg1.mat');
            [Def, Aff]   = spmdefs_get_dartel(flowfield, affine);
            [x2, y2, z2] = spmdefs_transform(Def, Aff, x1, y1, z1);
            [i2, j2, k2] = spmj_affine_transform(x2, y2, z2, inv(Vmask.mat));

            %%% resample the weights into SUIT space for the lp vector
            for r = 1 : size(data_map_lp,1)
                Vout_lp             = Vmask;
                Vout_lp.dat         = zeros(Vout_lp.dim);
                Vout_lp.dat(linIn3) = data_map_lp(r,:);
                Vout_lp.dt          = [64 0];
                Vout_lp.pinfo       = [1 0 0]';

                DataSUIT(r,:) = spm_sample_vol(Vout_lp,i2,j2,k2,1);

                V.dat   = zeros(V.dim);
                Vres_lp(r) = V;
                Vres_lp(r).dat(linIn1) = DataSUIT(r,:);  % Offset by one to account for 1 being medial wall
                Vres_lp(r).fname       = sprintf('data_%2.2d.nii',r);
                Vres_lp(r).pinfo       = [1 0 0]';
            end;

            %%% Now map the Lp vector to surface-based representation
            D_lp = suit_map2surf(Vres_lp,'stats','nanmean');
            
            %%% resample the weights into SUIT space for the pp vector
            for r = 1 : size(data_map_pp,1)
                Vout_pp             = Vmask;
                Vout_pp.dat         = zeros(Vout_pp.dim);
                Vout_pp.dat(linIn3) = data_map_pp(r,:);
                Vout_pp.dt          = [64 0];
                Vout_pp.pinfo       = [1 0 0]';

                DataSUIT(r,:) = spm_sample_vol(Vout_pp,i2,j2,k2,1);

                V.dat   = zeros(V.dim);
                Vres_pp(r) = V;
                Vres_pp(r).dat(linIn1) = DataSUIT(r,:);  % Offset by one to account for 1 being medial wall
                Vres_pp(r).fname       = sprintf('data_%2.2d.nii',r);
                Vres_pp(r).pinfo       = [1 0 0]';
            end;

            %%% Now map the Lp vector to surface-based representation
            D_pp = suit_map2surf(Vres_pp,'stats','nanmean');
            
            fm_struct.lp      = D_lp;
            fm_struct.Vres_lp = Vres_lp;
            fm_struct.pp      = D_pp;
            fm_struct.Vres_pp = Vres_pp;
            
            fm_struct.lpmethod   = lpmethod;
            fm_struct.ppmethod   = ppmethod;
            fm_struct.glm        = glm;
            fm_struct.roi1       = roi1;
            fm_struct.roi2       = roi2;
            fm_struct.regressors = condition;
            fm_struct.tsType     = tsType;
            
            fm_allsubs{s} = fm_struct;
            
            dircheck(fullfile(fmDir, subj_name{s}));
            save(fullfile(fmDir, subj_name{s}, sprintf('fm_%s_regions_%s_regions_%s_%s_%s.mat', taskName, roi1, roi2, tsType, lpmethod)), 'fm_struct', '-v7.3');
        end
        
        varargout{1} = fm_allsubs;
    case 'LAG:mdtb:suit:plot_flatmap'
        % plots the flatmap of the cerebellar measure
        % Example: sc1_sc2_mdtb('LAG:mdtb:suit:plot_flatmap', 'sn', [3]);
        
        sn         = returnSubjs;       %% subject list
        ppmethod   = '';         %% Did you do the slice timing correction as a preprocessing step?
        experiment = 1;
        glm        = 7;          %% the GLM number
        condition  = [];         %% the regressors of interest used in getting the predicted time series      
        tsType     = 'raw';      %% set this option to raw or predicted
        lpmethod   = 'max';      %% the method used to project the TD matrix into 1D. Options are: 'max', 'marek', 'percentile'
        
        colorscale_lp = [-0.2, 0.2];      %% colorscale for the lp flatmap. This will be used in suit_plotflatmap and caxis (to change the colorbar)
        colorscale_pp = [-0.6, 0.6];      %% colorscale for the peak projection flatmap. This will be used in suit_plotflatmap and caxis (to change the colorbar)
        
        
        roi1 = '162_tessellation_hem';
        roi2 = 'cerebellum_grey';      %% roi2 can be the same as roi1 if you want to calculate time delays in one roi 
        
        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'glm', 'roi1', 'roi2', 'maxLag', 'tsType', 'Tr', 'condition', 'lpmethod', 'colorscale_lp', 'colorscale_pp'});
        
        if isempty(condition)
            taskName = 'all';
        else
            % load in task information
            C  = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
            Cc = getrow(C,C.StudyNum == experiment);
            Ccc = getrow(Cc, Cc.condNum == condition(1));

            taskName = Ccc.taskNames;
        end
        
        experiment = sprintf('sc%d', experiment); %% experiment number is converted to 'sc1' or 'sc2'
        
        switch ppmethod
            case 'stc' 
                fmFigDir  = fullfile(baseDir, experiment, sprintf('FM_FIG_GLM_%d_stc', glm));
                fmDir = fullfile(baseDir, experiment, sprintf('FM_GLM_%d_stc', glm));
            case ''
                fmFigDir  = fullfile(baseDir, experiment, sprintf('FM_FIG_GLM_%d', glm));
                fmDir = fullfile(baseDir, experiment, sprintf('FM_GLM_%d', glm));
        end
        dircheck(fmFigDir)
        %%% plot the lp map
        for s = sn
            dircheck(fullfile(fmFigDir, subj_name{s}));
            % load in the struture variable for the flatmap
            fprintf('******************** Plot the lp flatmap for Subject %s ********************\n\n', subj_name{s});
            load(fullfile(fmDir, subj_name{s}, sprintf('fm_%s_regions_%s_regions_%s_%s_%s.mat', taskName, roi1, roi2, tsType, lpmethod)));
            suit_plotflatmap(fm_struct.lp , 'cmap', colormap(jet(256)), 'cscale', colorscale_lp);
            caxis(colorscale_lp);
            colorbar;
            title(sprintf('fm_lp_%s_regions_%s_regions_%s_%s_%s.mat', taskName, roi1, roi2, tsType, lpmethod), 'Interpreter', 'none');
            h = gcf;
            savefig(h, fullfile(fmFigDir, subj_name{s}, sprintf('fm_lp_%s_regions_%s_regions_%s_%s_%s.fig', taskName, roi1, roi2, tsType, lpmethod)));
            close(h);
        end
        
        %%% plot the pp map
        for s = sn
            % load in the struture variable for the flatmap
            fprintf('******************** Plot the pp flatmap for Subject %s ********************\n\n', subj_name{s});
            load(fullfile(fmDir, subj_name{s}, sprintf('fm_%s_regions_%s_regions_%s_%s_%s.mat', taskName, roi1, roi2, tsType, lpmethod)));
            suit_plotflatmap(fm_struct.pp , 'cmap', colormap(jet(256)), 'cscale', colorscale_pp);
            caxis(colorscale_lp);
            colorbar;
            title(sprintf('fm_pp_%s_regions_%s_regions_%s_%s_%s.mat', taskName, roi1, roi2, tsType, lpmethod), 'Interpreter', 'none');
            h = gcf;
            savefig(h, fullfile(fmFigDir, subj_name{s}, sprintf('fm_pp_%s_regions_%s_regions_%s_%s_%s.fig', taskName, roi1, roi2, tsType, lpmethod)));
            close(h);
        end
    case 'LAG:mdtb:group_lp'
    case 'LAG:mdtb:tmap_lp'

    case 'EA:mdtb:beta_dataframe'
        % creates an "Evoked Activation" dataframe
        % Example: sc1_sc2_mdtb('EA:mdtb:beta_dataframe', 'sn', [2, 3, 4, 7, 8])
        
        sn         = returnSubjs;       %% subject list
        ppmethod   = '';         %% Did you do the slice timing correction as a preprocessing step?
        experiment = 2;
        glm        = 7;          %% the GLM number
%         condition  = [];         %% the regressors of interest used in getting the predicted time series
        normmode   = 'UW';       %% set this to 'UW' for univariately normalized betas or 'NW' for non-normalized betas
        oparcel    = 1;          %% this is the empty parcel identified in 'ROI:mdtb:empty_parcel'
        
        %%% The rois for which you want to scatterplot the beta values
        roi1 = 'yeo_17WB';   % the cortical ROI
        roi2 = 'Buckner_17'; % the cerebellar ROI
        %%%% => nstructure_cc = 2
        nstructure_cc = 2;

        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'glm', 'roi1', 'roi2', 'regSelect'});
        
        experiment = sprintf('sc%d', experiment); %% experiment number is converted to 'sc1' or 'sc2'
        
        switch ppmethod
            case ''
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm));
                glmDir  = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm)); % where SPM_info.mat is
            case 'stc'
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm));
                glmDir  = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm)); % where SPM_info.mat is
        end
        
        df = [];
        fprintf('******************** creating the dataframe ********************\n');
        for s = sn
            fprintf('getting the data for %s\n', subj_name{s});
            % load in task information
            T = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            
            df_tmp.sn = s;
            for scc = 1:nstructure_cc
                eval(['roi = ', sprintf('roi%d;', scc)]);
                
                % load in the beta file for the roi
                Bfile = load(fullfile(betaDir, subj_name{s}, sprintf('beta_regions_%s.mat', roi)));
                Bfile = Bfile.B;
                
                nparcel = length(unique(Bfile.region_num)); % number of parcels
                
                % which structure?
                if ismember(roi, corticalParcels)
                    df_tmp.structure = cellstr('cortex');
                elseif ismember(roi, cerebellarParcels)
                    df_tmp.structure = cellstr('cerebellum');
                end % 'cortex' or 'cerebellum'
                
                % load in the region beta file
                mbp     = Bfile.(sprintf('mbetas%s', normmode)); % mean across voxels of a parcel
                bp_all  = Bfile.(sprintf('betas%s', normmode));  % all the beta values in a parcel
                
                % omit the parcel that you don't want the data for :)
                %%% it is the empty parcel that is identified in 'ROI:mdtb:empty_parcel'
                ips = 1:nparcel;  
                if ismember(roi, corticalParcels)
                    nHemiParcel = nparcel/2;
                    ips([oparcel, oparcel + nHemiParcel]) = [];
                elseif ismember(roi, cerebellarParcels)
                    ips(oparcel) = [];
                end
                
                for ip = ips
                    % get the bp for ipth parcel
                    bp = bp_all{ip}';
                    
                    % if the structure is the cortex then there will be two
                    % hemispheres! otherwise, there is no hemi :))
                    if ismember(roi, corticalParcels)
                        nHemiParcel = nparcel/2;
                        if ip <= nHemiParcel % in the left hemi
                            df_tmp.region_num = ip;
                            df_tmp.hemi       = cellstr('L');
                        elseif ip > nHemiParcel % in the right hemi
                            df_tmp.region_num = ip - nHemiParcel;
                            df_tmp.hemi       = cellstr('R');
                        end % fill in the fields for hemi and region_num 
                    elseif ismember(roi, cerebellarParcels)
                        df_tmp.region_num = ip;
                        df_tmp.hemi       = cellstr('none');
                    end % if it's cortex there are two hemi
                    % separate left and right hemi???
                    for i = 0:length(unique(T.cond))- 1 % including the instructions
                        
                        % remove the intercepts. One intercept for each run!
                        mbp_ri = mbp(:, 1:end - length(runLst));
                        bp_ri  = bp(:, 1:end - length(runLst));
                        
                        % remove the derivatives!
                        T_rd      = getrow(T, T.deriv == 0); % remove the derivative regressors
                        indd      = T.deriv == 0;            % get the indices for the non-derivative regressors
                        mbp_ri_rd = mbp_ri(:, indd);         % remove the derivative regressors
                        bp_ri_rd  = bp_ri(:, indd);          
                        
                        % get the condition
                        T_rd_c = getrow(T_rd, T_rd.cond == i); % get the dataframe for condition i
                        indc   = T_rd.cond == i;               % get the indices for the condition
                        
                        mbp_ri_rd_c = mbp_ri_rd(ip, indc);
                        bp_ri_rd_c  = bp_ri_rd(:, indc);
                        
                        nvox        = size(bp_ri_rd_c, 1); % number of voxels in a parcel 
                        
                        % reshape so that the runs are on the 3rd dimension
                        nregressors       = length(T_rd_c.cond)/length(runLst); % number of regressors (deriv and non-deriv) in a run
                        mbp_ri_rd_c_rshpd = reshape(mbp_ri_rd_c, [1, nregressors, length(runLst)]);                        
                        bp_ri_rd_c_rshpd  = reshape(bp_ri_rd_c, [nvox, nregressors, length(runLst)]);
                        
                        % create separate condition names for each instruction.
                        % Instructions will be like instruct_01, instruct_02, ...
                        % for the "cond" field all of the instructions will be
                        % labeled as 0, but they will be in different rows (the first 16 rows in all of the fields)
                        %%% the first nTask rows are for the instructions.
                        if i == 0 % the indicator for the instructions is 0. Separating 16 instructions!
                            for k = 1:nTask % there are 16 instructions
                                
                                %%% calculate mean and std across runs
                                mr_mbp_ri_rd_c_rshpd = mean(mbp_ri_rd_c_rshpd(:, k, :), 3);
                                sr_mbp_ri_rd_c_rshpd = std(mbp_ri_rd_c_rshpd(:, k, :), 0, 3);
                                
                                mr_bp_ri_rd_c_rshpd = mean(bp_ri_rd_c_rshpd(:, k, :), 3);
                                sr_bp_ri_rd_c_rshpd = std(bp_ri_rd_c_rshpd(:, k, :), 0, 3);
                                                                
                                df_tmp.mrmparcel = mr_mbp_ri_rd_c_rshpd;
                                df_tmp.srmparcel = sr_mbp_ri_rd_c_rshpd;
                                
                                df_tmp.mrparcel{1, 1} = mr_bp_ri_rd_c_rshpd;
                                df_tmp.srparcel{1, 1} = sr_bp_ri_rd_c_rshpd;
                              
                                df_tmp.cName = cellstr(sprintf('%s_%2.2d', char(unique(T_rd_c.TN)), k));
                                df_tmp.cond  = i;
                                df_tmp.task  = 0; %% instructions are task 0
                                df = addstruct(df, df_tmp);
                            end
                        else % do the whole thing for all the other task conditions
                            
                            %%% calculate mean and std across runs
                            mr_mbp_ri_rd_c_rshpd = mean(mbp_ri_rd_c_rshpd, 3);
                            sr_mbp_ri_rd_c_rshpd = std(mbp_ri_rd_c_rshpd, 0, 3);
                            
                            mr_bp_ri_rd_c_rshpd = mean(bp_ri_rd_c_rshpd, 3);
                            sr_bp_ri_rd_c_rshpd = std(bp_ri_rd_c_rshpd, 0, 3);
                                                        
                            df_tmp.mrmparcel = mr_mbp_ri_rd_c_rshpd;
                            df_tmp.srmparcel = sr_mbp_ri_rd_c_rshpd;
                            
                            df_tmp.mrparcel{1, 1} = mr_bp_ri_rd_c_rshpd;
                            df_tmp.srparcel{1, 1} = sr_bp_ri_rd_c_rshpd;
                            
                            df_tmp.cName = unique(T_rd_c.TN);
                            df_tmp.cond  = i;
                            df_tmp.task  = unique(T_rd_c.task);
                            df = addstruct(df, df_tmp);
                        end
                    end % i (condition)
                    
                end % ip
            end % scc
        end % sn
        fprintf('******************** Beta dataframe created ********************\n\n')
        save(fullfile(betaDir, sprintf('Beta%sDf_GLM_%d_regions_%s_regions_%s.mat', normmode, glm, roi1, roi2)), 'df', '-v7.3');
        varargout{1} = df;
    case 'EA:mdtb:beta_plot'
        % uses the dataframe structure created in the previous case to plot
        % the cerebellar activation vs cortical activation.
        % Example:sc1_sc2_mdtb('EA:mdtb:beta_plot', 'sn', [2, 3, 4, 6, 7, 8, 9, 10, 12, 14, 15, 17, 18, 19, 20, 21, 22, 24, 25, 27, 28, 29, 30, 31])
        
        %%% use xyplot, scatterplot and tapply to reproduce cool figures
        
        sn         = returnSubjs;        %% list of the subjects
        ppmethod   = '';          %% use the data with or without slice timing correction
        roi1       = 'yeo_17WB';
        roi2       = 'Buckner_17';
        glm        = 7;
        experiment = 2;           %% sc1 or sc2?
        ifig       = 0;           %% set this index to 0 if you don't want the plot for individual subjects
        ihemi      = 0;           %% set this flag to 1 if you want the scatterplots for each hemisphere separately
        category   = 'cond';      %% set it to 'cond' or 'task'
        normmode   = 'UW';        %% with or without univariate noise normalization
        parcelNum  = 2;           %% parcelNum for which you want to do the plotting
        oparcel    = 1;
        xvar       = 'cortex';
        yvar       = 'cerebellum';
        
        plotvar = 'mrmparcel'; %% the variable I want to use for plotting
        
        vararginoptions(varargin, {'sn', 'ppmethod', 'roi1', 'roi2', 'glm', 'experiment', 'category', 'normmode', 'parcelNum', 'ifig', 'oparcel'});
        
        % load in the task information for the subjects
        D = dload(fullfile(baseDir, 'sc1_sc2_taskConds_GLM.txt'));
        D = getrow(D, D.StudyNum == experiment);
        
        experiment = sprintf('sc%d', experiment); %% experiment number is converted to 'sc1' or 'sc2'
        
        if parcelNum == oparcel
            error('!!!!!!!!!!!!!!!!!!!!requesting the scatter plot for an empty parcel!!!!!!!!!!!!!!!!!!!!');
        end
        
        switch ppmethod
            case ''
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm)); 
            case 'stc'
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d_stc', glm));
        end
        
        % load in the beta dataframe
        load(fullfile(betaDir, sprintf('Beta%sDf_GLM_%d_regions_%s_regions_%s.mat', normmode, glm, roi1, roi2)))
        
        switch ifig
            case 1 % the individual figure flag is set to 1
                switch ihemi
                    case 0 % you don't want the scatter plots for each hemi separately
                        for s = sn'
                            % create the beta vector for each ROI
                            d_sub = getrow(df, df.sn == s);
                            T_sub = tapply(d_sub, {'structure', 'region_num', 'cond'}, {'mrmparcel'});
                            x_axe = getrow(T_sub, strcmp(T_sub.structure, xvar));
                            x_axe2 = getrow(x_axe, x_axe.region_num == parcelNum);
                            
                            y_axe = getrow(T_sub, strcmp(T_sub.structure, yvar));
                            y_axe2 = getrow(y_axe, y_axe.region_num == parcelNum);
                            
                            figure;
                            scatterplot(x_axe2.(plotvar), y_axe2.(plotvar), 'markertype', 'v',...
                                'markercolor', [0.6, 0.2 0], 'markerfill', [0.6, 0.2, 0],...
                                'markersize', 6, 'label', D.taskNames);
                            
                            h           = lsline;
                            h.Color     = 'b';
                            h.LineWidth = 3;
                            
                            title(sprintf('cerebellum vs cortex for parcel %0.2d', parcelNum));
                            xlabel(xvar);
                            ylabel(yvar);
                        end % sn
                    case 1 % you want the scatter plots for each hemi separately
                        for s = sn
                            % create the beta vector for each ROI
                            d_sub = getrow(df, df.sn == s);
                            
                            % get the dataframe for each structure
                            xdf_sub = getrow(d_sub, strcmp(d_sub.structure, xvar));
                            ydf_sub = getrow(d_sub, strcmp(d_sub.structure, yvar));
                            
                            % get the dataframe for the parcels you want to do the plotting
                            if strcmp(xvar, 'cortex')
                                %%% xdf_L and xdf_R are the data structures for left and right
                                %%% hemispheres for that particular parcel. So the only
                                %%% variables remaining will be subjects and conditions.
                                %%% all you need to do is to calculate the average between left
                                %%% and right for each subject and each condition
                                xdf_sub_L = getrow(xdf_sub, strcmp(xdf_sub.hemi, 'L'));
                                xdf_sub_R = getrow(xdf_sub, strcmp(xdf_sub.hemi, 'R'));
                                
                                xdf_sub_L_p = getrow(xdf_sub_L, xdf_sub_L.region_num == parcelNum);
                                xdf_sub_R_p = getrow(xdf_sub_R, xdf_sub_R.region_num == parcelNum);
                                
                                %             xdf_L_s = getrow(xdf_L, strcmp(xdf_L.hemi, 'L'));
                                %             xdf_R_s = getrow(xdf_R, strcmp(xdf_R.hemi, 'R'));
                                
                                T_sub_L = tapply(xdf_sub_L_p, {category}, {plotvar});
                                T_sub_R = tapply(xdf_sub_R_p, {category}, {plotvar});
                                
                                x_L = T_sub_L.(plotvar);
                                x_R = T_sub_R.(plotvar);
                            end
                            ydf_sub_p = getrow(ydf_sub, ydf_sub.region_num == parcelNum);
                            
                            T_sub_y = tapply(ydf_sub_p, {category}, {plotvar});
                            
                            y = T_sub_y.(plotvar);
                            
                            figure;
                            subplot(121)
                            scatterplot(x_L, y, 'markertype', 'v', 'markercolor', [0.6, 0.2 0], 'markerfill', [0.6, 0.2, 0], 'markersize', 6, 'label', df.cName);
                            h           = lsline;
                            h.Color     = 'b';
                            h.LineWidth = 3;
                            
                            title(sprintf('vs left hemi for parcel %d in %s', parcelNum, subj_name{s}))
                            xlabel(xvar)
                            ylabel(yvar)
                            
                            subplot(122)
                            scatterplot(x_R, y, 'markertype', 'v', 'markercolor', [0.6, 0.2 0], 'markerfill', [0.6, 0.2, 0], 'markersize', 6, 'label', df.cName);
                            h           = lsline;
                            h.Color     = 'b';
                            h.LineWidth = 3;
                            
                            title(sprintf('vs right hemi for parcel %d in %s', parcelNum, subj_name{s}))
                            xlabel(xvar)
                            ylabel(yvar)
                        end % sn
                end %  switch ihemi
            case 0 % you don't want the scatter plots for each subject
                switch ihemi
                    case 0 % you don't want separate scatterplots for each cortical hemi
                        T_all = tapply(df, {'structure', 'region_num', 'cond'}, {'mrmparcel'});
                        x_axe = getrow(T_all, strcmp(T_all.structure, xvar));
                        x_axe2 = getrow(x_axe, x_axe.region_num == parcelNum);
                        
                        y_axe = getrow(T_all, strcmp(T_all.structure, yvar));
                        y_axe2 = getrow(y_axe, y_axe.region_num == parcelNum);
                        
                        figure;
                        scatterplot(x_axe2.(plotvar), y_axe2.(plotvar), 'markertype', 'v',...
                            'markercolor', [0.6, 0.2 0], 'markerfill', [0.6, 0.2, 0],...
                            'markersize', 6, 'label', D.taskNames);
                        
                        h           = lsline;
                        h.Color     = 'b';
                        h.LineWidth = 3;
                        
                        title(sprintf('cerebellum vs cortex for parcel %0.2d', parcelNum));
                        xlabel(xvar);
                        ylabel(yvar);
                    case 1 % you want separate scatterplots for each hemi
                        figure;
                        subplot(121)
                        %         xyplot(xdf_L_p.(plotvar), ydf_p.(plotvar), [], 'split', df.(category), 'errorbars', 'ellipse');
                        scatterplot(T_L.(plotvar), T_Y.(plotvar), 'markertype', 'v',...
                            'markercolor', [0.6, 0.2 0], 'markerfill', [0.6, 0.2, 0],...
                            'markersize', 6, 'label', D.taskNames);
                        h           = lsline;
                        h.Color     = 'b';
                        h.LineWidth = 3;
                        title(sprintf('vs left hemi for parcel %d', parcelNum))
                        xlabel(xvar)
                        ylabel(yvar)
                        
                        subplot(122)
                        scatterplot(T_R.(plotvar), T_Y.(plotvar), 'markertype', 'v',...
                            'markercolor', [0.6, 0.2 0], 'markerfill', [0.6, 0.2, 0],...
                            'markersize', 6, 'label', D.taskNames);
                        h           = lsline;
                        h.Color     = 'b';
                        h.LineWidth = 3;
                        title(sprintf('vs right hemi for parcel %d', parcelNum))
                        xlabel(xvar)
                        ylabel(yvar)
                end
                
        end % switch ifig
        
        
        % scatter plot for the whole structures
        figure;
        subplot(121)
        %         xyplot(xdf_L_p.(plotvar), ydf_p.(plotvar), [], 'split', df.(category), 'errorbars', 'ellipse');
        scatterplot(T_L_s.(plotvar), T_Y_s.(plotvar), 'markertype', 'v',...
            'markercolor', [0.6, 0.2 0], 'markerfill', [0.6, 0.2, 0],...
            'markersize', 6, 'label', D.taskNames);
        h           = lsline;
        h.Color     = 'b';
        h.LineWidth = 3;
        title(sprintf('vs left hemi'))
        xlabel(xvar)
        ylabel(yvar)
        
        subplot(122)
        scatterplot(T_R_s.(plotvar), T_Y_s.(plotvar), 'markertype', 'v',...
            'markercolor', [0.6, 0.2 0], 'markerfill', [0.6, 0.2, 0],...
            'markersize', 6, 'label', D.taskNames);
        h           = lsline;
        h.Color     = 'b';
        h.LineWidth = 3;
        title(sprintf('vs right hemi'))
        xlabel(xvar)
        ylabel(yvar)
    case 'EA:mdtb:beta_corr_net'
        % This case calculates the correlations between activations in
        % networks across structures (cerebellum and cortex)
        % example: sc1_sc2_mdtb('EA:mdtb:beta_corr_net', 'sn', [2, 3, 4, 6, 8, 9, 10, 12, 14, 15, 17, 18, 19, 20, 21, 22, 24, 25])
        
        sn         = returnSubjs;        %% subject list
        ppmethod   = '';          %% Did you do the slice timing correction as a preprocessing step?
        experiment = 1;
        glm        = 7;           %% the GLM number
        normmode   = 'UW';        %% set this to 'UW' for univariately normalized betas or 'NW' for non-normalized betas
        corrvar    = 'mrmparcel'; %% corr between mrmparcel will be calculated.
        ihemi      = 0;           %% set this flag to 1 if you want to calculate the corr for each hemisphere separately
%         oparcel    = 1;           %% oparcel is identified after running the case 'ROI:mdtb:empty_parcel'
        
        %%% The rois for which you want to scatterplot the beta values
        roi1 = 'yeo_17WB';   % the cortical ROI
        roi2 = 'Buckner_17'; % the cerebellar ROI

        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'glm', 'roi1', 'roi2', 'ihemi'});
        
        experiment = sprintf('sc%d', experiment); %% experiment number is converted to 'sc1' or 'sc2'
        
        switch ppmethod
            case ''
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm));
            case 'stc'
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm));
        end
        
        % load in the beta dataframe
        load(fullfile(betaDir, sprintf('Beta%sDf_GLM_%d_regions_%s_regions_%s.mat', normmode, glm, roi1, roi2)))
        
        % I'm going to use mrmparcel field of the dataframe. The other
        % option I first thought about was mrparcel. But calculating the
        % correlation using mrparcel doesn't make sense as there are
        % different numbers of voxels within matching and non-matching
        % networks of the cortex and the cerebellum.
        % The idea is to loop over the parcels and calculate the
        % correlations between parcels in the cerebellum with the parcels
        % in the cortex (which can be done hemisphere-wise).
        % First I calculate the average across subjects and then correlate
        % them. But I can also calculate the correlation matrix for
        % individual subjects. 
        % I expect to see high correlations between matching networks in
        % the cortex and the cerebellum. But I won't be surprised if there
        % are high correlations between non-matching networks as well.
        
%         ss = 1; % subject index
        for s = 1:length(sn)
            df_sub     = getrow(df, df.sn == sn(s));
            structures = unique(df_sub.structure);
            switch ihemi
                case 0 % you don't want to separate two hemispheres
                    T_sub = tapply(df_sub, {'structure', 'region_num', 'cond'}, {corrvar});
%                     scs      = unique(T_sub.structure);
                    x        = 0;
                    for scc = 1:length(structures)
                        T_sub_sc = getrow(T_sub, strcmp(T_sub.structure, structures{scc}));
                        parcels  = unique(T_sub_sc.region_num);
                        
                        for ip = 1:length(parcels)
                            T_sub_sc_p = getrow(T_sub_sc, T_sub_sc.region_num == parcels(ip));
                            mrmparcel_subs(:, ip + x*length(parcels), s) = bsxfun(@minus, T_sub_sc_p.(corrvar), mean(T_sub_sc_p.(corrvar)));
                        end % ip
                        x = x+1;
                    end % scc
                case 1 % you want to separate the two hemispheres
                    for scc = 1:length(structures)
                        df_sub_sc = getrow(df_sub, strcmp(df_sub.structure, structures{scc}));
                        hemis     = unique(df_sub_sc.hemi);
                        x         = length(hemis) - 1;
                        for ih = 1:length(hemis)
                            df_sub_sc_h = getrow(df_sub_sc, strcmp(df_sub_sc.hemi, hemis{ih}));
                            parcels     = unique(df_sub_sc_h.region_num);
                            
                            for ip = 1:length(parcels)
                                df_sub_sc_h_p = getrow(df_sub_sc_h, df_sub_sc_h.region_num == parcels(ip));
                                mrmparcel_subs(:, ip + x * length(parcels), s) = bsxfun(@minus, df_sub_sc_h_p.(corrvar), mean(df_sub_sc_h_p.(corrvar)));
                            end % ip
                            x = x+1;
                        end % ih
                    end % scc
            end % switch ihemi
            
%             ss = ss + 1;
        end % sn
        msubs = mean(mrmparcel_subs, 3);
        R     = corr(msubs);
        
        figure; imagesc(R); axis square; colorbar;
        xlabel('network');
        ylabel('network');
        title('inter-network correlations of beta values')
        
        varargout{1} = R;
        varargout{2} = mrmparcel_subs;        
    case 'EA:mdtb:beta_modelfit'
        % uses the dataframe created in beta_dataframe case to fit a line
        % to the relationship between the activations in the cerebellum and
        % cortex.
        % this can be done at the group level, or at the individual level.
        % at the individual level, it fits the line for each subject
        % separately. At the group level, it fits the line to the
        % activations averaged across subjects. 
        % It also has the option of doing the fitting for the whole
        % structures (cerebellum vs cortex) and also within matching
        % networks (buckner vs yeo).
        % after the fitting, the residuals will be calculated and saved!
        % Example: sc1_sc2_mdtb('EA:mdtb:beta_modelfit', 'sn', [2])
                
        sn         = returnSubjs;         %% list of the subjects
        ppmethod   = '';           %% use the data with or without slice timing correction
        roi1       = 'yeo_17WB';   %%
        roi2       = 'Buckner_17'; %%
        glm        = 7;            %% glm number
        experiment = 1;            %% sc1 or sc2?
        category   = 'cond';       %% set this option to 'task' or 'cond' 
        ihemi      = 1;            %% set this flag to 1 if you want to do the fitting using the activations averaged across the whole structures and 0 if you want to do it network- (parcel-) wise
        igroup     = 1;            %% set this to 1 if you want to fit the line to the group level data and to 0 if you want to fit the line to individual subjects
        normmode   = 'UW';
        pf         = 5;            %% the degree of the polynomial that you want to fit to the data
        visualize  = 0;            %% set this to 1 if you want to plot and 0 if you don't want the plot
        
        fitvar     = 'mrmparcel'; %% the variable I want to fit the line to
        
        vararginoptions(varargin, {'sn', 'ppmethod', 'roi1', 'roi2', 'glm', 'experiment', 'category', 'ihemi', 'igroup', 'fitvar', 'normmode', 'pf', 'visualize'});
        
        % load in the task information for the subjects
        D = dload(fullfile(baseDir, 'sc1_sc2_taskConds_GLM.txt'));
        D = getrow(D, D.StudyNum == experiment);
        
        experiment = sprintf('sc%d', experiment); %% experiment number is converted to 'sc1' or 'sc2'
        
        switch ppmethod
            case ''
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm)); 
            case 'stc'
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d_stc', glm));
        end
        
        % parcellations were done for which structures? usually, it's
        % cerebellum and cortex
        if ismember(roi1, corticalParcels)
            x = 'cortex';
%         elseif
        end % if it's cortex 
        
        if ismember(roi2, cerebellarParcels)
            y = 'cerebellum';
%         elseif
        end % if it's cerebellum
        
        % load in the beta dataframe
        load(fullfile(betaDir, sprintf('Beta%sDf_GLM_%d_regions_%s_regions_%s.mat', normmode, glm, roi1, roi2)))
        
        switch ihemi
            case 1 % you want to do the fitting network- and hemisphere- wise
                switch igroup % do the fitting at the individual or group level?
                    case 0 % individual level
                        T_all = tapply(df, {'sn', 'structure', 'hemi', 'region_num', category}, {fitvar});
                        dfFit_all = [];
                        dfFit_sub = [];
                        for s = sn
                            dfFit_tmp = [];
                            dfFit_tmp.sn = s;
                            
                            T_all_sub = getrow(T_all, T_all.sn == s);
                            yax_sub = getrow(T_all_sub, strcmp(T_all_sub.structure, y)); % cerebellum on y
                            xax_sub = getrow(T_all_sub, strcmp(T_all_sub.structure, x)); % cortex on x
                            
                            parcels = unique(T_all_sub.region_num);
                            hemis   = unique(xax_sub.hemi);
                            
                            for ip = parcels'
                                fprintf('fitting %d polynomial to data for %s for parcel %0.2d\n', pf, subj_name{s}, ip);
                                dfFit_tmp.region_num = ip;
                                
                                yax_sub_p = getrow(yax_sub, yax_sub.region_num == ip);
                                xax_sub_p = getrow(xax_sub, xax_sub.region_num == ip);
                                
                                switch visualize
                                    case 1
                                        ih = 1;
                                        figure;
                                end % if visualize

                                for h = hemis'
                                    dfFit_tmp.hemi = h;
                                    xax_sub_p_h  = getrow(xax_sub_p, strcmp(xax_sub_p.hemi, h));
                                    xvar_sub_p_h = xax_sub_p_h.(fitvar);
                                    yvar_sub_p   = yax_sub_p.(fitvar);
                                    
                                    % the fitting
                                    [p_sub_p_h, S_sub_p_h]        = polyfit(xvar_sub_p_h, yvar_sub_p, pf);
                                    [yfit_sub_p_h, delta_sub_p_h] = polyval(p_sub_p_h, xvar_sub_p_h, S_sub_p_h);
                                    R2_sub_p_h                    = (S_sub_p_h.normr/norm(yvar_sub_p - mean(yvar_sub_p)))^2;
                                    
                                    switch visualize
                                        case 1 % if you want to visualize
                                            subplot(sprintf('12%d', ih))
                                            scatterplot(xvar_sub_p_h, yvar_sub_p, 'markertype', 'v',...
                                                'markercolor', [0.6, 0.2 0], 'markerfill', [0.6, 0.2, 0],...
                                                'markersize', 6, 'label', D.taskNames);
                                            hf           = lsline;
                                            hf.Color     = 'b';
                                            hf.LineWidth = 3;
                                            hold on
                                            myxp = min(xvar_sub_p_h):0.01:max(xvar_sub_p_h);
                                            [myyp, ~ ] = polyval(p_sub_p_h, myxp, S_sub_p_h);
                                            plot(myxp, myyp, 'r--', 'LineWidth', 2);
                                            ih = ih + 1;
                                    end % visualize or not
                                    
                                    % create the dataframe
                                    dfFit_tmp.pcoef = p_sub_p_h;
                                    dfFit_tmp.S     = S_sub_p_h;
                                    dfFit_tmp.yfit  = yfit_sub_p_h';
                                    dfFit_tmp.yres  = (yvar_sub_p - yfit_sub_p_h)';
                                    dfFit_tmp.delta = delta_sub_p_h';
                                    dfFit_tmp.R2    = R2_sub_p_h;
                                    
                                    dfFit_sub = addstruct(dfFit_sub, dfFit_tmp);
                                end % h
                            end % ip 
                            dfFit_all = addstruct(dfFit_all, dfFit_sub);
                            % save the dataframe for each subject
                            dircheck(fullfile(betaDir, subj_name{s}));
                            save(fullfile(betaDir, subj_name{s}, sprintf('%s_Beta%spolyFit%d_region_%s_hemi_region_%s.mat', subj_name{s}, normmode, pf, roi1, roi2)), 'dfFit_sub', '-v7.3');
                        end % sn
                        save(fullfile(betaDir, sprintf('Beta%spolyFit%d_allsub_region_%s_hemi_region_%s.mat', normmode, pf, roi1, roi2)), 'dfFit_all', '-v7.3');
                    case 1 % group level
                        T_all   = tapply(df, {'structure', 'hemi', 'region_num', category}, {fitvar});

                        yax = getrow(T_all, strcmp(T_all.structure, y)); % cerebellum on y
                        xax = getrow(T_all, strcmp(T_all.structure, x)); % cortex on x
                        
                        parcels = unique(T_all.region_num); % parcels
                        hemis   = unique(xax.hemi);         % hemispheres
                        
                        dfFitAvg = []; % the dataframe containing the fit variables
                        for ip = parcels'
                            fprintf('fitting %d polynomial to data for parcel %0.2d\n', pf, ip);
                            dfFit_tmp.region_num = ip;
                            
                            yax_p = getrow(yax, yax.region_num == ip);
                            xax_p = getrow(xax, xax.region_num == ip);
                            
                            switch visualize
                                case 1
                                    ih = 1;
                                    figure;
                            end % if visualize
                            
                            for h = hemis'
                                dfFit_tmp.hemi = h;
                                
                                xax_p_h = getrow(xax_p, strcmp(xax_p.hemi, h));
                                yvar_p  = yax_p.(fitvar);
                                xvar_p_h  = xax_p_h.(fitvar);
                                
                                % the fitting!
                                [p_p_h, S_p_h]        = polyfit(xvar_p_h, yvar_p, pf);
                                [yfit_p_h, delta_p_h] = polyval(p_p_h, xvar_p_h, S_p_h);
                                R2_p_h                = (S_p_h.normr/norm(yvar_p - mean(yvar_p)))^2;
                                
                                switch visualize
                                    case 1 % if you want to visualize
                                        subplot(sprintf('12%d', ih))
                                        scatterplot(xvar_p_h, yvar_p, 'markertype', 'v',...
                                            'markercolor', [0.6, 0.2 0], 'markerfill', [0.6, 0.2, 0],...
                                            'markersize', 6, 'label', D.taskNames);
                                        hf           = lsline;
                                        hf.Color     = 'b';
                                        hf.LineWidth = 3;
                                        hold on
                                        myxp = min(xvar_p_h):0.01:max(xvar_p_h);
                                        [myyp, ~ ] = polyval(p_p_h, myxp, S_p_h);
                                        plot(myxp, myyp, 'r--', 'LineWidth', 2);
                                        ih = ih + 1;
                                end % visualize or not
                                
                                % create a dataframe
                                dfFit_tmp.pcoef = p_p_h;
                                dfFit_tmp.S     = S_p_h;     % the same structure variable returned by polyfit and used in polyval
                                dfFit_tmp.yfit  = yfit_p_h';
                                dfFit_tmp.yres  = (yvar_p - yfit_p_h)';
                                dfFit_tmp.delta = delta_p_h'; % delta is an estimate of the standard deviation of the error in predicting a future observation at x by p(x)
                                dfFit_tmp.R2    = R2_p_h;
                                
                                dfFitAvg = addstruct(dfFitAvg, dfFit_tmp);
                            end % h
                        end % ip
                        % save the dataframe
                        dircheck(betaDir);
                        save(fullfile(betaDir, sprintf('Beta%spolyFit%d_group_region_%s_hemi_region_%s.mat', normmode, pf, roi1, roi2)), 'dfFitAvg', '-v7.3');
                end
            case 0 % you want to do the fitting using the data averaged across the whole structures
                switch igroup % do the fitting at the individual or group level?
                    case 0 % individual level
                        T_all = tapply(df, {'sn', 'structure', 'region_num', category}, {fitvar});
                        
                        dfFit_all = [];
                        dfFit_sub = [];
                        for s = sn
                            dfFit_tmp = [];
                            dfFit_tmp.sn = s;
                            
                            T_all_sub = getrow(T_all, T_all.sn == s);
                            yax_sub = getrow(T_all_sub, strcmp(T_all_sub.structure, y)); % cerebellum on y
                            xax_sub = getrow(T_all_sub, strcmp(T_all_sub.structure, x)); % cortex on x
                            
                            parcels = unique(T_all_sub.region_num);
                            
                            for ip = parcels'
                                fprintf('fitting %d polynomial to data for %s for parcel %0.2d\n', pf, subj_name{s}, ip);
                                dfFit_tmp.region_num = ip;
                                
                                yax_sub_p = getrow(yax_sub, yax_sub.region_num == ip);
                                xax_sub_p = getrow(xax_sub, xax_sub.region_num == ip);
                                
                                yvar_sub_p = yax_sub_p.(fitvar);
                                xvar_sub_p = xax_sub_p.(fitvar);
                                
                                % the fitting
                                [p_sub_p, S_sub_p]        = polyfit(xvar_sub_p, yvar_sub_p, pf);
                                [yfit_sub_p, delta_sub_p] = polyval(p_sub_p, xvar_sub_p, S_sub_p);
                                R2_sub_p                  = (S_sub_p.normr/norm(yvar_sub_p - mean(yvar_sub_p)))^2;
                                
                                switch visualize
                                    case 1 % if you want to visualize
                                        figure;
                                        scatterplot(xvar_sub_p, yvar_sub_p, 'markertype', 'v',...
                                            'markercolor', [0.6, 0.2 0], 'markerfill', [0.6, 0.2, 0],...
                                            'markersize', 6, 'label', D.taskNames);
                                        hf           = lsline;
                                        hf.Color     = 'b';
                                        hf.LineWidth = 3;
                                        hold on
                                        myxp = min(xvar_sub_p):0.01:max(xvar_sub_p);
                                        [myyp, ~ ] = polyval(p_sub_p, myxp, S_sub_p);
                                        plot(myxp, myyp, 'r--', 'LineWidth', 2);
                                end % visualize or not
                                
                                % create the dataframe
                                dfFit_tmp.pcoef = p_sub_p;
                                dfFit_tmp.S     = S_sub_p;
                                dfFit_tmp.yfit  = yfit_sub_p';
                                dfFit_tmp.yres  = (yvar_sub_p - yfit_sub_p)';
                                dfFit_tmp.delta = delta_sub_p';
                                dfFit_tmp.R2    = R2_sub_p;
                                
                                dfFit_sub = addstruct(dfFit_sub, dfFit_tmp);
                            end % ip 
                            dfFit_all = addstruct(dfFit_all, dfFit_sub);
                            % save the dataframe for each subject
                            dircheck(fullfile(betaDir, subj_name{s}));
                            save(fullfile(betaDir, subj_name{s}, sprintf('%s_Beta%spolyFit%d_region_%s_region_%s.mat', subj_name{s}, normmode, pf, roi1, roi2)), 'dfFit_sub', '-v7.3');
                        end % sn
                        save(fullfile(betaDir, sprintf('Beta%spolyFit%d_allsub_region_%s_region_%s.mat', normmode, pf, roi1, roi2)), 'dfFit_all', '-v7.3');
                    case 1 % group level
                        T_all   = tapply(df, {'structure', 'region_num', category}, {fitvar});
                        parcels = unique(T_all.region_num);
                        
                        yax = getrow(T_all, strcmp(T_all.structure, y)); % cerebellum on y
                        xax = getrow(T_all, strcmp(T_all.structure, x)); % cortex on x
                        
                        dfFitAvg = []; % the dataframe containing the fit variables
                        for ip = parcels'
                            fprintf('fitting %d polynomial to data for parcel %0.2d\n', pf, ip);
                            dfFit_tmp.region_num = ip;
                            
                            yax_p = getrow(yax, yax.region_num == ip);
                            xax_p = getrow(xax, xax.region_num == ip);
                            
                            yvar_p = yax_p.(fitvar);
                            xvar_p = xax_p.(fitvar);
                            
                            % the fitting!
                            [p_p, S_p]        = polyfit(xvar_p, yvar_p, pf);
                            [yfit_p, delta_p] = polyval(p_p, xvar_p, S_p);
                            R2_p              = (S_p.normr/norm(yvar_p - mean(yvar_p)))^2;
                            
                            switch visualize
                                case 1 % if you want to visualize
                                    figure;
                                    scatterplot(xvar_p, yvar_p, 'markertype', 'v',...
                                        'markercolor', [0.6, 0.2 0], 'markerfill', [0.6, 0.2, 0],...
                                        'markersize', 6, 'label', D.taskNames);
                                    h           = lsline;
                                    h.Color     = 'b';
                                    h.LineWidth = 3;
                                    hold on
                                    myxp = min(xvar_p):0.01:max(xvar_p);
                                    [myyp, ~ ] = polyval(p_p, myxp, S_p);
                                    plot(myxp, myyp, 'r--', 'LineWidth', 2);
                            end % visualize or not
                            
                            % create a dataframe
                            dfFit_tmp.pcoef = p_p;
                            dfFit_tmp.S     = S_p;     % the same structure variable returned by polyfit and used in polyval
                            dfFit_tmp.yfit  = yfit_p';
                            dfFit_tmp.yres  = (yvar_p - yfit_p)';
                            dfFit_tmp.delta = delta_p'; % delta is an estimate of the standard deviation of the error in predicting a future observation at x by p(x)
                            dfFit_tmp.R2    = R2_p;
                            
                            dfFitAvg = addstruct(dfFitAvg, dfFit_tmp);
                        end % ip
                        % save the dataframe
                        dircheck(betaDir);
                        save(fullfile(betaDir, sprintf('Beta%spolyFit%d_group_region_%s_region_%s.mat', normmode, pf, roi1, roi2)), 'dfFitAvg', '-v7.3');
                end % switch igroup
        end % switch ihemi
        
        %%% code scratch for the ttest
%         myDf = getrow(dfFit_all, dfFit_all.region_num == 5)
%         [h, p, ci, stats] = ttest(myDf.yres);
    case 'EA:mdtb:beta_modelfit_condCV'
        % This case do the same linear fitting as before, but each time
        % does the fitting with one condition left out. Then it "tries" to
        % predict the linear regression estimate for the left out
        % condition.
        % the dataframe is the same as the previous one (polyfit) with "CV"
        % added to the name to show that the task cv was done.
        % Example: sc1_sc2_mdtb('EA:mdtb:beta_modelfit_condCV', 'sn', [2, 3, 4, 6, 7, 8, 9, 10])
        
        sn         = returnSubjs;
        experiment = 1;
        ppmethod   = '';
        glm        = 7;
        roi1       = 'yeo_17WB';
        roi2       = 'Buckner_17'; %
        normmode   = 'UW';
        igroup     = 0;            % do the whole thing on group level or individdual level data? default: individual level
        ihemi      = 0;            % do it for each hemi separately or for the whole cortex. default: whole structure
        fitvar     = 'mrmparcel';  % the variable for which you want to do the fitting
        pf         = 1;            % the polynomial fit degree
        category   = 'cond';
        
        vararginoptions(varargin, {'sn', 'experiment', 'ppmethod', 'glm', 'roi1', 'roi2', 'normmode', 'igroup', 'fitvar', 'pf', 'category'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm));
            case 'stc'
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d_stc', glm));
        end
        
        % parcellations were done for which structures? usually, it's
        % cerebellum and cortex
        if ismember(roi1, corticalParcels)
            x = 'cortex';
%         elseif
        end % if it's cortex 
        
        if ismember(roi2, cerebellarParcels)
            y = 'cerebellum';
%         elseif
        end % if it's cerebellum
        
        % load in the beta dataframe (loads 'df' into workspace).
        load(fullfile(betaDir, sprintf('Beta%sDf_GLM_%d_regions_%s_regions_%s.mat', normmode, glm, roi1, roi2)));
        
        switch ihemi
            case 1 % for each hemi separately
                switch igroup
                    case 0 % on individual level
                    case 1 % on group level
                        T_all = tapply(df, {'structure', 'hemi', 'region_num', category}, {fitvar});
                end % switch igroup
            case 0 % for the whole cortex
                switch igroup
                    case 0 % on individual level
                        T_all = tapply(df, {'sn', 'structure', 'region_num', category}, {fitvar});
                        df_all = [];
                        for s = sn
                            dfFit_sub = [];
                            dfFit_tmp = [];
                            dfFit_tmp.sn = s;
                            
                            T_all_sub = getrow(T_all, T_all.sn == s);
                            yax_sub = getrow(T_all_sub, strcmp(T_all_sub.structure, y)); % cerebellum on y
                            xax_sub = getrow(T_all_sub, strcmp(T_all_sub.structure, x)); % cortex on x
                            
                            parcels = unique(T_all_sub.region_num);
                            conds   = unique(T_all_sub.cond);
                            
                            % overal fit has already been calculated using
                            % the previous case. For plotting R2 and R2cv
                            % vs. model order, you may use the data from
                            % the previous case which has already been
                            % saved.
         
                            for ip = parcels' 
                                dfFit_tmp.region_num = ip;
                                
                                % get the data for parcel ip
                                yax_sub_p = getrow(yax_sub, yax_sub.region_num == ip);
                                xax_sub_p = getrow(xax_sub, xax_sub.region_num == ip);
                                
                                
                                icc = 1;
                                for ic = conds'
                                    dfFit_tmp.leave_c = ic; % the left out condition
                                    
                                    fprintf('fitting %d polynomial to data for %s for parcel %d, leaving out condition %d\n', pf, subj_name{s}, ip, ic);
                                    
                                    % get the data points for the left out
                                    % condition in parcel ip
                                    yax_sub_p_lc = getrow(yax_sub_p, yax_sub_p.cond == ic); % cerebellum on y axis
                                    xax_sub_p_lc = getrow(xax_sub_p, xax_sub_p.cond == ic); % cortex on x a xis

                                    % get the data for all the conditions
                                    % except for the left out condition
                                    yax_sub_p_nc = getrow(yax_sub_p, yax_sub_p.cond ~= ic);
                                    xax_sub_p_nc = getrow(xax_sub_p, xax_sub_p.cond ~= ic);
                                    
                                    % get the y and x variables for the
                                    % left out condition
                                    yvar_sub_p_lc = yax_sub_p_lc.(fitvar);
                                    xvar_sub_p_lc = xax_sub_p_lc.(fitvar);
                                    
                                    % get the y and x variables for all the
                                    % conditions except for the left out
                                    % one.
                                    yvar_sub_p_nc = yax_sub_p_nc.(fitvar);
                                    xvar_sub_p_nc = xax_sub_p_nc.(fitvar);
                                    
                                    % the fitting (using the datapoints for all the condition except for the left out condition)
                                    [p_sub_p_nc, S_sub_p_nc]        = polyfit(xvar_sub_p_nc, yvar_sub_p_nc, pf);
%                                     [yfit_sub_p_nc, delta_sub_p_nc] = polyval(p_sub_p_nc, xvar_sub_p_nc, S_sub_p_nc);
%                                     R2_sub_nc_p                     = (S_sub_p_nc.normr/norm(yvar_sub_p_nc - mean(yvar_sub_p_nc)))^2;
                                    
                                    % calculate the residuals for the left
                                    % out condition
                                    [yfit_sub_p_lc(icc), ~] = polyval(p_sub_p_nc, xvar_sub_p_lc, S_sub_p_nc);
                                    
%                                     % create the dataframe
%                                     dfFit_tmp.pcoef    = p_sub_p_nc;
%                                     dfFit_tmp.S        = S_sub_p_nc;
%                                     dfFit_tmp.yfit_lc  = yfit_sub_p_lc(icc);
%                                     dfFit_tmp.yres     = (yvar_sub_p_nc - yfit_sub_p_nc)';
%                                     dfFit_tmp.delta    = delta_sub_p_nc';
%                                     dfFit_tmp.R2_nc    = R2_sub_nc_p;
%                                     
%                                     dfFit_tmp.yres_lc = (yvar_sub_p_lc - yfit_sub_p_lc(icc))';
%                                     
%                                     dfFit_sub = addstruct(dfFit_sub, dfFit_tmp);
                                    icc = icc+1;
                                end % ic
                                
                                % RSScv
                                RSScv = sum((yax_sub_p.mrmparcel' - yfit_sub_p_lc).^2)
                                % TSS
                                TSS   = sum((yax_sub_p.mrmparcel' - mean(yax_sub_p.mrmparcel')).^2)
                                % R2cv
                                R2cv = 1 - (RSScv/TSS)
                                
                                dfFit_tmp.R2cv = R2cv;
                                
                                % calculate the overal fit and R2
                                yvar_sub_p = yax_sub_p.(fitvar);
                                xvar_sub_p = xax_sub_p.(fitvar);
                                
                                [p_sub_p, S_sub_p] = polyfit(xvar_sub_p, yvar_sub_p, pf);
                                [yfit_sub_p, ~]    = polyval(p_sub_p, xvar_sub_p, S_sub_p);
                                
                                RSS = sum((yax_sub_p.mrmparcel - yfit_sub_p).^2)
                                
                                R2           = 1 - (RSS/TSS)
                                dfFit_tmp.R2 = R2;
                                
                                dfFit_sub = addstruct(dfFit_sub, dfFit_tmp);
                            end % ip
                            % save the dataframe for each subject
                            dircheck(fullfile(betaDir, subj_name{s}));
                            save(fullfile(betaDir, subj_name{s}, sprintf('%s_Beta%spolyFit%d_condCV_region_%s_region_%s.mat', subj_name{s}, normmode, pf, roi1, roi2)), 'dfFit_sub', '-v7.3');
                            df_all = addstruct(df_all, dfFit_sub);
                        end % sn
                        dircheck(fullfile(betaDir));
                        save(fullfile(betaDir, sprintf('Beta%spolyfit%d_all_condCV_region_%s_region_%s.mat', normmode, pf, roi1, roi2)), 'df_all', '-v7.3');
                        
                        varargout{1} = df_all;
                    case 1 % on group level
                end % switch igroup
        end % switch ihemi
    case 'EA:mdtb:modelfit_res_ttest'
        % for each region, does the ttest to see whether the beta for each
        % condition is significantly different from zero. H0:the linear fit
        % residuals is equal to zero for that condition.
        % Example: sc1_sc2_mdtb('EA:mdtb:modelfit_res_ttest');
        
        ppmethod   = '';
        experiment = 1;
        roi1       = 'yeo_17WB';
        roi2       = 'Buckner_17';
        ihemi      = 0;    % doi for each hemi separately
        igroup     = 0;    % use the residuals from the linear fits done at the subject level (igroup = 0)
        normmode   = 'UW';
        pf         = 1;    % the degree of the fitted polynomial
        glm        = 7;    % glm number
        alpha      = 0.05; % pvalue significance level
        
        vararginoptions(varargin, {'ppmethod', 'experiment', 'roi1', 'roi2', 'ihemi', 'igroup', 'normmode', 'pf', 'glm'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm));
            case 'stc'
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d_stc', glm));
        end
        
        % load in the dataframe that have all the info for all the subjects
        load(fullfile(betaDir, sprintf('Beta%spolyfit%d_allsub_region_%s_region_%s.mat', normmode, pf, roi1, roi2)));
        
        tDf = [];
        parcels = unique(dfFit_all.region_num);
        for ip = parcels'
            tmpDf = getrow(dfFit_all, dfFit_all.region_num == ip);
            
            x                 = tmpDf.yres;
            [h, p, ci, stats] = ttest(x, 0, 'Alpha', alpha);
            
            tDf_tmp.region_num = ip;
            tDf_tmp.h        = h;
            tDf_tmp.p        = p;
            tDf_tmp.ci{1}    = ci;
            tDf_tmp.stats{1} = stats;
            tDf.alpha        = alpha;
            
            tDf = addstruct(tDf, tDf_tmp);
            
            % ttestDirect(dfFit_all.yres,[dfFit_all.sn],2, 'onesample')
            
        end % ip
        
        save(fullfile(betaDir, sprintf('ttest_%0.2d_beta_polyFit%d_region_%s_region_%s.mat', alpha, pf, roi1, roi2)), 'tDf', '-v7.3');
        
        % use ttestDirect from dataframe
        % you can also use ttest from dataframe or matlab's ttest    
    case 'EA:mdtb:model_selection'
        % fits different degrees of models to the data for each parcel and
        % then plots R2 and R2cv to see which model is the best for each
        % individual subject.
        % Example: sc1_sc2_mdtb('EA:mdtb:model_selection', 'sn', a)
        
        sn         = returnSubjs;
        glm        = 7;
        experiment = 1;
        ppmethod   = '';
        roi1       = 'yeo_17WB';
        roi2       = 'Buckner_17'; %
        ihemi      = 0;
        igroup     = 0;
        morder     = 9;            % the maximum model order you want to consider
        normmode   = 'UW';
        visualize  = 1;
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment', 'ppmethod', 'roi1', 'roi2', 'ihemi', 'igroup', 'morder', 'normmode', 'visualize'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_glm_%d', glm));
            case 'stc'
                baseDir = fullfile(baseDir, experiment, sprintf('Beta_glm_%d_stc', glm));
        end
            
        switch ihemi
            case 0
                switch igroup
                    case 0
                        modelFit_df = [];
                        modelFit_df_tmp = [];
                        
                        % loop over subjects
                        for s = sn
                            modelFit_df_tmp.sn = s;

                            for im = 1:morder
                                
                                modelFit_df_tmp.order = im;
                                % load in the dataframe created using cross
                                % validation
                                beta_modelCV = load(fullfile(betaDir, sprintf('Beta%spolyfit%d_all_condCV_region_%s_region_%s.mat', normmode, im, roi1, roi2)));
                                beta_modelCV = beta_modelCV.df_all;
                                
                                parcels = unique(beta_modelCV.region_num);
                           
                                beta_modelCV_sub = getrow(beta_modelCV, beta_modelCV.sn == s);

                                for ip = parcels'
                                    modelFit_df_tmp.region_num = ip;
                                    beta_modelCV_sub_p = getrow(beta_modelCV_sub, beta_modelCV_sub.region_num == ip);
                                    R2_sub_pm         = beta_modelCV_sub_p.R2;
                                    R2cv_sub_pm       = beta_modelCV_sub_p.R2cv;
                                    
                                    
                                    modelFit_df_tmp.R2     = R2_sub_pm;
                                    modelFit_df_tmp.R2cv   = R2cv_sub_pm;
                                    modelFit_df = addstruct(modelFit_df, modelFit_df_tmp);
                                end % ip
                            end % im (model order)
                        end % s (sn)
                        
                        switch visualize
                            case 1 % you want to do the visualization for the subject
                                for s = sn
                                    fprintf('Do the plotting for %s\n', subj_name{s});
                                    tmp_df = getrow(modelFit_df, modelFit_df.sn == s);
                                    
                                    for ip = parcels'
                                        tmp_df_p = getrow(tmp_df, tmp_df.region_num == ip);
                                        R2_tmp_sub_p = tmp_df_p.R2;
                                        R2cv_tmp_sub_p = tmp_df_p.R2cv;
                                        order = tmp_df_p.order;
                                        
                                        % the plotting
                                        figure; plot(order, R2_tmp_sub_p, 'b*')
                                        hold on
                                        plot(order, R2cv_tmp_sub_p, 'ro')
                                        legend('R2', 'R2cv');
                                        title(sprintf('R2 and R2cv for %s parcel %0.1d', subj_name{s}, ip))
                                        keyboard;
                                    end % ip
                                    keyboard; 
                                end % s(sn)
                            case 0 % don't want to do the visulization
                        end % switch visualize
                    case 1
                end % igroup
            case 1
                switch igroup
                    case 0
                    case 1
                end % igroup
        end % ihemi
        
    case 'CONN:mdtb:wta_cerebellum'
        % correlate each cerebellar voxel time series with cortical
        % parcels/vertices time series and "register" the highest
        % correlation. Then a map is created for the cerebellum with every
        % cerebellar voxel having the highest correlation value.
        % This is just to see how correlated cerebellar voxels are to the
        % cortex and which regions have the highest correlations.
        % it's a sorta winner-take-all :)
        % Example: sc1_sc2_mdtb('CONN:mdtb:wta_cerebellum', 'sn', [3]);
        
        sn         = returnSubjs;              %
        experiment = 1;                 %
        ppmethod   = '';                %
        glm        = 7;                 %
        roi1       = 'tesselsWB_162';   %
        roi2       = 'cerebellum_grey'; %
        prun       = 0;                 % set this flag to 1 if you want to do the thing per run :)
        visualize  = 1;                 % set this fag to 1 if you want to visualize the maps and 0 if you don't
        corrvar    = 'beta';              % set this to 'beta' if you want to use beta values and 'ts' if you want to use raw time series.
        normmode   = 'UW';
        
        vararginoptions(varargin, {'sn', 'experiment', 'ppmethod', 'glm', 'roi1', 'roi2', 'prun', 'visualize', 'corrvar', 'normmode'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                tsDir   = fullfile(baseDir, experiment, sprintf('TS_GLM_%d', glm));
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm));
            case 'stc'
                tsDir   = fullfile(baseDir, experiment, sprintf('TS_GLM_%d_stc', glm));
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d_stc', glm));
        end

        for s = sn
            % load in the data
            switch corrvar
                case 'ts' % use time series to calculate correlations
                    switch prun % do it on the time series averaged across runs or for each run separately?
                        case 0 % load in the time series averaged across runs
                            roi1_ts = load(fullfile(tsDir, subj_name{s}, sprintf('avgTS_all_regions_%s.mat', roi1)));
                            roi2_ts = load(fullfile(tsDir, subj_name{s}, sprintf('avgTS_all_regions_%s.mat', roi2)));
                            var2 = roi2_ts.avgts.y_raw;
                            var1 = roi1_ts.avgts.y_raw;
                        case 1 % load in the time series concatenated across runs
                            roi1_ts = load(fullfile(tsDir, subj_name{s}, sprintf('TS_all_regions_%s.mat', roi1)));
                            roi2_ts = load(fullfile(tsDir, subj_name{s}, sprintf('TS_all_regions_%s.mat', roi2)));
                            var2 = roi2_ts.avgts.y_raw;
                            var1 = roi1_ts.avgts.y_raw;
                    end % switch prun
                case 'beta' % use beta values to calculate correlations
                    roi1_b = load(fullfile(betaDir, subj_name{s}, sprintf('beta_regions_%s.mat', roi1)));
                    roi2_b = load(fullfile(betaDir, subj_name{s}, sprintf('beta_regions_%s.mat', roi2)));
                    
                    roi1_b = roi1_b.B;
                    roi2_b = roi2_b.B;
                    switch prun
                        case 0
                            var1 = roi1_b.(sprintf('mrbetas%s', normmode));
                            var2 = roi2_b.(sprintf('mrbetas%s', normmode));
                        case 1
                            var1 = roi1_b.(sprintf('betas%s', normmode));
                            var2 = roi2_b.(sprintf('betas%s', normmode));
                    end % switch prun
                    
            end % switch corrvar
            
            
            
%             % pick a random voxel from the cerebellum (Just for testing)
%             svox    = randi(size(var2, 2)); %% selects a random integer from 1 to number of voxels/parcels
%             seed    = var2zm(:, svox);
%             % subtract the mean before calculating correlation
%             var1zm = bsxfun(@minus, var1, mean(var1));
%             seedzm = bsxfun(@minus, seed, mean(seed));
%             
%             R = corr(seed, var1);
%             
%             max(R)
%             
%             keyboard;
            
            % create a map of the correlations with the cortex
            % a winner-take-all algorithm
            % loop over all the cerebellar voxels, calculate the
            % correlations with the cortical parcels and assign the winner
            % correlation to that voxel and simply create a map.
            % this map will be telling us which voxels are highly
            % correlated with the cortex
            Rho = zeros(size(var2, 2), 1); % preallocating the correlation matrix
            for ivox = 1:size(var2, 1)
                seed          = var2(ivox, :);
                    
                % subtract the mean before calculating correlation
                var1zm = bsxfun(@minus, var1, mean(var1));
                seedzm = bsxfun(@minus, seed, mean(seed));
                
                R_tmp         = corr(seedzm', var1zm');
                R_tmp_abs     = abs(R_tmp);
                [~, mxind]    = max(R_tmp_abs);
                Rho(ivox, 1)  = R_tmp(mxind);
            end % ivox
            
            keyboard;
            
            % creating flatmaps to visualize using suit
            switch visualize
                case 1
                     % Determine the voxels we want to resample in SUIT space
                     V = spm_vol(fullfile(baseDir, experiment,suitDir,'anatomicals','cerebellarGreySUIT.nii'));
                     X = spm_read_vols(V);
                     
                     grey_threshold = 0.1; % gray matter threshold
                     
                     linIn1     = find(X > grey_threshold);
                     [i1,j1,k1] = ind2sub(V.dim,linIn1');
                     [x1,y1,z1] = spmj_affine_transform(i1, j1, k1, V.mat);
                     
                     %%% Determine voxel locations from the original ROI
                     load(fullfile(baseDir, experiment, regDir, 'data', subj_name{s}, sprintf('regions_%s.mat', roi2))); % 'regions' are defined in 'ROI_define'
                     
                     Vmask      = spm_vol(fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'maskbrainSUITGrey.nii')); %% cerebellum grey matter mask
                     [i3,j3,k3] = spmj_affine_transform(R{1}.data(:,1),R{1}.data(:,2),R{1}.data(:,3),inv(Vmask.mat));
                     linIn3     = sub2ind(Vmask.dim, round(i3), round(j3), round(k3));
                     
                     %%% transform SUIT coords into anatomical space of the individual
                     flowfield    = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'u_a_c_anatomical_seg1.nii');
                     affine       = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'Affine_c_anatomical_seg1.mat');
                     [Def, Aff]   = spmdefs_get_dartel(flowfield, affine);
                     [x2, y2, z2] = spmdefs_transform(Def, Aff, x1, y1, z1);
                     [i2, j2, k2] = spmj_affine_transform(x2, y2, z2, inv(Vmask.mat));
                     
                     %%% resample the weights into SUIT space for the lp vector
                     for r = 1 : size(Rho,2)
                         Vout             = Vmask;
                         Vout.dat         = zeros(Vout.dim);
                         Vout.dat(linIn3) = Rho(:,r);
                         Vout.dt          = [64 0];
                         Vout.pinfo       = [1 0 0]';
                         
                         DataSUIT(r,:) = spm_sample_vol(Vout,i2,j2,k2,1);
                         
                         V.dat   = zeros(V.dim);
                         Vres(r) = V;
                         Vres(r).dat(linIn1) = DataSUIT(r,:);  % Offset by one to account for 1 being medial wall
                         Vres(r).fname       = sprintf('data_%2.2d.nii',r);
                         Vres(r).pinfo       = [1 0 0]';
                     end % r (number of regions)
                     
                     %%% Now map the Lp vector to surface-based representation
                     D = suit_map2surf(Vres,'stats','nanmean');
                     
                     figure; suit_plotflatmap(D , 'cmap', colormap(jet(256)), 'cscale', [min(D(:)), max(D(:))]);
                     caxis([min(D(:)), max(D(:))]);
                     colorbar;
                     title(sprintf('max corr with cortex parcels map for %s', subj_name{s}));
            end % switch visualize
            keyboard
        end % sn
    case 'CONN:mdtb:wta_cortex'
        % same as 'CONN:mdtb:wta_cerebellum' but for the cortex. So I'm
        % taking each cortical voxel/vertex, correlating it with all the
        % cerebellar voxel and then just "register" the highest correlation
        % value, assign it to that voxel/vertex and create a map of the
        % highest correlation values.
        % It will give a sense of the range of correlation values and which
        % cortical regions have the highest correlations with the
        % cerebellum.
        % Example: sc1_sc2_mdtb('CONN:mdtb:wta_cortex', 'sn', [3]);
    case 'CONN:mdtb:toy_vtoparcel'
        % This is a toy case just for testing and possibly visualization.
        % it takes one cerebellar voxel, correlate it with all the
        % voxels/vertices in a cortical parcel and create a cortical map 
        % Just to have a sense of correlation values 
        % Example: sc1_sc2_mdtb('CONN:mdtb:toy_vtoparcel', 'sn', [2])
        
        sn         = returnSubjs;
        ppmethod   = '';
        experiment = 1;
        glm        = 7;
        corrvar    = 'beta';
        prun       = 0;
        pparcel    = 0;                 % do it using the mean "var" for the parcel? set this to 1 if you want the mean var 
        roi1       = 'yeo_17WB';
        roi2       = 'cerebellum_grey';
        normmode   = 'UW';
        oparcel    = 1;                 % omit parcel number oparcel
        
        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'glm', 'corrvar', 'prun', 'roi1', 'roi2', 'normmode'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                tsDir   = fullfile(baseDir, experiment, sprintf('TS_GLM_%d', glm));
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm));
            case 'stc'
                tsDir   = fullfile(baseDir, experiment, sprintf('TS_GLM_%d_stc', glm));
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d_stc', glm));
        end
        
        for s = sn
            % load in the data
            switch corrvar
                case 'ts' % use the time series to calculate the correlation
                case 'beta' % use beta values to calculate the correlation
                    switch prun
                        case 1 % do it for each run separately
                        case 0 % do it with the data averaged across runs
                            B1 = load(fullfile(betaDir, subj_name{s}, sprintf('beta_regions_%s.mat', roi1)));
                            B2 = load(fullfile(betaDir, subj_name{s}, sprintf('beta_regions_%s.mat', roi2)));
                            B1 = B1.B;
                            B2 = B2.B;
                            switch pparcel
                                case 0 % do it for all the voxels in a parcel
                                    var1 = B1.(sprintf('mrbetas%s', normmode)); 
                                    var2 = B2.(sprintf('mrbetas%s', normmode));
                                case 1 % do it per parcel (using the mean)
                                    var1 = B1.(sprintf('mrmbetas%s', normmode)); 
                                    var2 = B2.(sprintf('mrmbetas%s', normmode));
                            end % switch pparcel
                    end % switch prun
            end % switch corrvar
            
            % right now the code is written so that it works with
            % cerebellum_grey and some cortical parcel. So the following
            % lines might change in future
            var2mat      = cell2mat(var2);
            var1cell     = var1;
            nCortParcels = length(var1cell);
            nCerebVoxels = size(var2mat, 1);
            RHO          = cell(nCortParcels, nCerebVoxels);
            for ivox = 1:nCerebVoxels
                for ip = 1:nCortParcels % looping over parcels
                    if ip ~= oparcel && ip ~= ((nCortParcels)/2)+ oparcel
                        var1_tmp = var1cell{ip};
                        varseed  = var2mat(ivox, :);
                        
                        % zero mean before correlation
                        varseed_zm = bsxfun(@minus, varseed, mean(varseed));
                        var1zm     = bsxfun(@minus, var1_tmp, mean(var1_tmp));
                        
                        % calculate the correlations
                        RHO{ip, ivox} = corr(varseed_zm', var1zm');
                    else
                        RHO{ip, ivox}  = NaN;
                    end % if ip is not oparcel
                end % ip (parcel)
            end % ivox
            
            % plotting the distribution of correlation coefficients with
            % the cerebellum for each cortical parcel
            for ip = 1: nCortParcels
                if ip ~= oparcel && ip ~= ((nCortParcels)/2)+ oparcel
                    tmp = cat(2, RHO{ip, :});
                    figure; histogram(tmp);
                    keyboard;
                end % if ip is not oparcel
            end % ip (parcel) AGAIN!
        end % sn
            
        keyboard
    case 'CONN:mdtb:toy_corr_overlap'
        % this case determines the overlap between the cerebellar maps of
        % 95th percentile correlation coef with cortical parcels:
        % loads in the cortical parcels data
        % calculate the correlation of all cerebellar voxels with the
        % voxels within each cortical parcel
        % a map of 95th percentile for each parcel
        % the overlap?
        % Example: sc1_sc2_mdtb('CONN:mdtb:toy_corr_overlap', 'sn', [2])
        
        sn         = returnSubjs;
        ppmethod   = '';
        experiment = 1;
        roi1       = 'yeo_17WB';
        roi2       = 'cerebellum_grey'; %
        corrvar    = 'beta';
        normmode   = 'UW';
        glm        = 7;
        coparcel   = 1;
        prun       = 0;
        fprc       = 0;                 % set this flag to 1 if you want to use the percentiles to create maps and 0 if you want to use "absolute" threshold  
        pcntile    = 95;                % the percentile
        thresh     = 0.7;               % the "absolute" threshold value
        
        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'roi1', 'roi2', 'corrvar', 'normmode', 'glm', 'coparcel', 'prun', 'fprc', 'th'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm));
            case 'stc'
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d_stc', glm));
        end
        
        for s = sn
            switch corrvar % use the time series or the beta values
                case 'ts' % use the time series
                case 'beta' % use the beta values
                    % load in the beta values for cortical ROI
                    B1 = load(fullfile(betaDir, subj_name{s}, sprintf('beta_regions_%s.mat', roi1)));
                    % load in the beta values for cerebellum grey
                    B2 = load(fullfile(betaDir, subj_name{s}, sprintf('beta_regions_%s.mat', roi2)));
                    
                    B1 = B1.B;
                    B2 = B2.B;
                    
                    nparcel = length(B1.region_num);
                    
                    switch prun 
                        case 0 % do it with the data averaged across runs
                            var1 = B1.(sprintf('mrbetas%s', normmode));
                            var2 = B2.(sprintf('mrbetas%s', normmode));
                            
                            % var2 is the file for cerebellum_grey but
                            % still, it's a cell array and I have to
                            % convert it to a matrix
                            var2mat = cell2mat(var2);
                        case 1 % do the whole thing for each run separately
                    end % prun
            end % switch corrvar
            
            % calculate the correlation of each voxel in each
            % cortical parcel with the cerebellum
            corr_map       = cell(1, nparcel); % this cell array will have the correlation maps for each parcel
            corr_map_label = cell(1, nparcel); % this cell array will have ones (or ip) for voxels where the correlation is in 95th percentile
            for ip = 1:nparcel
                if ip ~= coparcel && ip ~= coparcel + (nparcel/2)
                    var1_tmp     = var1{ip};
                    corr_map{ip} = corr(var1_tmp', var2mat');
                    switch fprc
                        case 1 % use percentile
                            threshold = prctile(corr_map{ip}(:), pcntile);
                        case 0 % use "absolute" threshold
                            threshold = thresh;
                    end % switch fprc (use the percentile or "absolute" threshold)
                    labeled      = bsxfun(@ge, corr_map{ip}, threshold);
                    labeled2     = sum(labeled);
                    y            = labeled2 ~=0;
                    
                    corr_map_label{ip} = y;
                else
                    corr_map{ip}       = NaN;
                    corr_map_label{ip} = NaN;
                end % if not coparcel (the parcel that I'm going to omit)
                
            end % ip (parcel)
            
            % loop over the parcels and project the map onto the
            % suit flatmap
            for ip = 1:nparcel
                if ip ~= coparcel && ip ~= (nparcel/2) + coparcel
                    % Determine the voxels we want to resample in SUIT space
                    V = spm_vol(fullfile(baseDir, experiment,suitDir,'anatomicals','cerebellarGreySUIT.nii'));
                    X = spm_read_vols(V);
                    
                    grey_threshold = 0.1; % gray matter threshold
                    
                    linIn1     = find(X > grey_threshold);
                    [i1,j1,k1] = ind2sub(V.dim,linIn1');
                    [x1,y1,z1] = spmj_affine_transform(i1, j1, k1, V.mat);
                    
                    %%% Determine voxel locations from the original ROI
                    load(fullfile(baseDir, experiment, regDir, 'data', subj_name{s}, sprintf('regions_%s.mat', roi2))); % 'regions' are defined in 'ROI_define'
                    
                    Vmask      = spm_vol(fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'maskbrainSUITGrey.nii')); %% cerebellum grey matter mask
                    [i3,j3,k3] = spmj_affine_transform(R{1}.data(:,1),R{1}.data(:,2),R{1}.data(:,3),inv(Vmask.mat));
                    linIn3     = sub2ind(Vmask.dim, round(i3), round(j3), round(k3));
                    
                    %%% transform SUIT coords into anatomical space of the individual
                    flowfield    = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'u_a_c_anatomical_seg1.nii');
                    affine       = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'Affine_c_anatomical_seg1.mat');
                    [Def, Aff]   = spmdefs_get_dartel(flowfield, affine);
                    [x2, y2, z2] = spmdefs_transform(Def, Aff, x1, y1, z1);
                    [i2, j2, k2] = spmj_affine_transform(x2, y2, z2, inv(Vmask.mat));
                    
                    %%% resample the weights into SUIT space for the lp vector
                    
                    Vout             = Vmask;
                    Vout.dat         = zeros(Vout.dim);
                    Vout.dat(linIn3) = corr_map_label{ip};
                    Vout.dt          = [64 0];
                    Vout.pinfo       = [1 0 0]';
                    
                    DataSUIT = spm_sample_vol(Vout,i2,j2,k2,1);
                    
                    V.dat                = zeros(V.dim);
                    Vres             = V;
                    Vres.dat(linIn1) = DataSUIT;  % Offset by one to account for 1 being medial wall
                    Vres.fname       = sprintf('data_%2.2d.nii',ip);
                    Vres.pinfo       = [1 0 0]';
                    
%                     myV = Vres(ip);
                    
                    %%% Now map the Lp vector to surface-based representation
                    D = suit_map2surf(Vres,'stats','nanmean');
                    
                    figure; suit_plotflatmap(D , 'cmap', colormap(jet(256)), 'cscale', [min(D(:)), max(D(:))]);
                    caxis([min(D(:)), max(D(:))]);
                    colorbar;
                    title(sprintf('cerebellar regions with high corr with cortex parcel %0.2d for %s', ip, subj_name{s}));
                    
                    keyboard;
                end % if not the parcel you want to omit (parcel 1 is usually omited)
            end % ip (parcels again)
        end % sn
        
    case 'SHINE:mdtb:ts_dataframe'
        % creates the dataframe which will then be used in applying shines
        % method. You have to run the reslicing case before you run this
        % case so that all subjects raw time series are in the suit space.
        % maybe you can use something like the case for mapping 1D maps to suit space! 
        % Example: sc1_sc2_mdtb('SHINE:mdtb:ts_dataframe', 'sn', [2, 3, 4, 6, 8, 9, 10])
        
        sn         = returnSubjs;
        ppmethod   = '';
        glm        = 7;
        experiment = 1;
        roi        = 'Buckner_17'; % the roi for which you want to do the pca thing
        nrun       = 16;                % there are 16 runs 
        nsess      = 2;                 % there are two sessions. * runs per sessions
        
        vararginoptions(varargin, {'sn', 'ppmethod', 'glm', 'experiment', 'roi'});
        
        experiment = sprintf('sc%d', experiment);
        
        ntp = 598; % there are 598 time points.
        
        switch ppmethod
            case ''
                tsDir = fullfile(baseDir, experiment, sprintf('TS_GLM_%d', glm));
            case 'stc'
                tsDir = fullfile(baseDir, experiment, sprintf('TS_GLM_%d_stc', glm));
        end % switch ppmethod
        
        ts_df = [];
        for s = sn
            fprintf('******************** get the time series data for %s ********************\n', subj_name{s});
            ts_df_tmp.sn = s;
            % load in the time series data for the roi
            load(fullfile(tsDir, subj_name{s}, sprintf('Ts_all_regions_%s', roi)));
            %%% separate the sessions and runs
            for ss = 1:nsess
                ts_df_tmp.sess = ss;
                rpsess = ((ss - 1)* (nrun/2)) + 1:(ss * nrun/2); % list of runs per session
                for rr = rpsess
                    % I will include two types of indices for runs:
                    % urun_num is the unified run number that runs from 1
                    % to 16. run_num is the run numbers per session that
                    % goes from 1 to 8.
                    ts_df_tmp.urun_num = rr;
                    if ss == 1
                        ts_df_tmp.run_num  = rr;
                    else
                        ts_df_tmp.run_num = rr - (nrun/2);
                    end % ss == 1 or ss == 2?
                    
                    % get the range of numbers for time points in the run
                    % rr
                    ir = (rr - 1)*ntp + 1;
                    fr = rr * ntp;
                    % get the time series (y_raw) and (y_hat) for each run
                    ry_res = ts.y_res(ir:fr, :);
                    ry_hat = ts.y_hat(ir:fr, :);
                    
                    ts_df_tmp.y_res = ry_res;
                    ts_df_tmp.y_hat = ry_hat;
                end % rr (runs)
            end % ss (sessions)
            ts_df = addstruct(ts_df, ts_df_tmp);
        end % sn
        
        keyboard;
      
    case 'Visualize:mdtb:suit'
        % takes in a vector (or map) for a subject and transfer it to suit space and plot it
        % on the flatmap
        % Example: sc1_sc2_mdtb('Visualize:mdtb:suit', corrmap);

        subj       = 2;           % the suit anatomical data for the subject is used in the mapping 
        experiment = 1;
        data       = varargin{1}; % the vector map you want to transfer to suit flatmap
        
        vararginoptions(varargin, {'sn', 'experiment'});
        
        experiment = sprintf('sc%d', experiment);
        
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

    case 'CHECK:mdtb:glm_design_check' % This case will be used to check regressors 
        % loads in the design matrix for the GLM and SPM_info and plots
        % each regressor.
        % Example: sc1_sc2_mdtb('CHECK:mdtb:glm_design_check')
        
        sn         = 3;  %% subject used to check the design matrix
        glm        = 8;  %% the glm I want to check
        ppmethod   = '';
        experiment = 1;
        
        vararginoptions(varargin, {'sn', 'glm', 'ppmethod', 'experiment'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
        end
        
        % load in the SPM file
        load(fullfile(glmDir, subj_name{sn}, 'SPM.mat'));
        X = SPM.xX.X;
        X = X(:, 1:end - 16); % removing the intercepts
        
        nreg = size(X, 2)/16; % number of regressors
        
        % load in SPM_info file
        T         = load(fullfile(glmDir, subj_name{sn}, 'SPM_info.mat'));
        taskNames = T.TN;
        
        for r = runLst
            fprintf('Checking regressors for glm%d in %s and run %0.2d\n', glm, subj_name{sn}, r);
            lb_tp = (r - 1)* ntp + 1;
            ub_tp = lb_tp + (ntp - 1);
            
            lb_b = (r - 1)*nreg + 1;
            ub_b = lb_b + (nreg - 1);
            
            X_run = X(lb_tp:ub_tp, lb_b:ub_b);
            
            for it = 1:2:length(unique(taskNames))*2
                figure;
                plot(X_run(:, it), 'r');
                title(sprintf('the regressor for %s in run %0.2d', taskNames{it}, r));
                keyboard;
            end % it (tasks)
        end % r (run)
        
        %%% Checking the convolution of boxcars with lengths of 5 and 30
        %%% with canonical hrf
        oneCycle_task = [zeros(1, 30), ones(1, 30), zeros(1, 30)];
        oneCycle_inst = [zeros(1, 5), ones(1, 5), zeros(1, 5)];
        
        TR  = 1;
        hrf = spm_hrf(TR);
        
        reg_task = conv(hrf, oneCycle_task);
        reg_inst = conv(hrf, oneCycle_inst);
        
        figure;
        subplot(121);
        plot(reg_task);
        subplot(122);
        plot(reg_inst);
        
        keyboard;
    case 'CHECK:mdtb:glm_design_s026'
        % checking s026
        % Example: sc1_sc2_mdtb('CHECK:mdtb:glm_design_s026');
        
        for r = 1:length(runLst)
            r
            D_26  = dload(fullfile(baseDir, 'sc1',behavDir, 's26', 'sc1_s26_GoNoGo.dat'));
            R_26  = getrow(D_26,D_26.runNum==runB(r))
            D_02  = dload(fullfile(baseDir, 'sc1',behavDir, 's02', 'sc1_s02_GoNoGo.dat'));
            R_02  = getrow(D_02,D_02.runNum==runB(r))
        end
    case 'CHECK:mdtb:modify_betaUnn'
        % I ran the beta unn case for s02 to s10 and then added a field to
        % the beta data structure. So for these subjects that field didn't
        % exist. Using this case, I'm adding the field to the beta data
        % structure. Use this case to add fields to the beta data structure
        % after modifying the case :)
        % Example: sc1_sc2_mdtb('CHECK:mdtb:modify_betaUnn', 'sn', [3])
        sn         = returnSubjs;
        ppmethod   = '';
        experiment = 1;
        roi        = 'yeo_17WB';
        glm        = 7;
        oparcel    = 1;
        
        vararginoptions(varargin, {'sn', 'ppmethod', 'experiment', 'glm', 'roi', 'oparcel'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm));
                glmDir  = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
            case 'stc'
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d_stc', glm));
                glmDir  = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
        end
        dircheck(betaDir); % beta values for each roi is saved here
        
%         if ismember(roi, corticalParcels) % the roi is from cortex
%             structure_cc = 'cortex';
%         elseif ismember(roi, cerebellarParcels) % the roi is from cerebellum
%             structure_cc = 'cerebellum';
%         end % a cortical roi or a cerebellar roi?

        for s = sn
            fprintf('******************** calculating univariate normalized betas for %s ********************\n', subj_name{s});
            
            % load in ROI file
            load(fullfile(baseDir, experiment, regDir, 'data', subj_name{s}, sprintf('regions_%s.mat', roi)));
            
            % load in SPM_info to get the info for the tasks and
            % condditions
            T           = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            T_rd        = getrow(T, T.deriv == 0); %% removing the derivative regressors
            nregressors = length(T_rd.cond)/length(runLst); % number of regressors (deriv and non-deriv) in a run
            indd        = T.deriv == 0;
            
            % load in the beta file
            load(fullfile(betaDir, subj_name{s}, sprintf('beta_regions_%s.mat', roi)));
            
            % removing the old fields :(
            % I ran the case for some of the ROIs and created the new fields
            % then I corrected the whole thing and ran the case again. Here,
            % I'm just removing the fields that were created before.
            if isfield(B, 'mrbetasNW') 
                B = rmfield(B, 'mrbetasNW');
            end
            if isfield(B, 'mrbetasUW') 
                B = rmfield(B, 'mrbetasUW');
            end
            if isfield(B, 'mrmbetasNW') 
                B = rmfield(B, 'mrmbetasNW');
            end
            if isfield(B, 'mrmbetasUW') 
                B = rmfield(B, 'mrmbetasUW');
            end
            
            B_tmp = [];
            for r=1:numel(R) % R is the output 'regions' structure from 'ROI_define'
% %                 B_tmp.region_name = R{r}.name;
%                 B_tmp.structure_cc = structure_cc;
%                 B                 = addstruct(B, B_tmp);
                
                
                % load in SPM_info to get the info for the tasks and
                % condditions
                T           = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
                T_rd        = getrow(T, T.deriv == 0); %% removing the derivative regressors
                nregressors = length(T_rd.cond)/length(runLst); % number of regressors (deriv and non-deriv) in a run
                indd        = T.deriv == 0;
                
                % calculate the mean beta across runs
                %%% for all the voxels in the parcel
                %%%% remove the intercepts
                
                if r ~= oparcel && r ~= oparcel + 17
                    beta_hat = B.betasNW{r}';
                    betasUW  = B.betasUW{r}';
                    beta_hat_ri = beta_hat(:, 1:end - length(runLst));
                    betasUW_ri  = betasUW(:, 1:end - length(runLst));
                    
                    beta_hat_ri_rd = beta_hat_ri(:, indd);
                    betasUW_ri_rd  = betasUW_ri(:, indd);
                    
                    nvox = size(B.betasNW{r}, 2); % number of voxels in a parcel.
                    beta_hat_ri_reshaped = reshape(beta_hat_ri_rd, [nvox, nregressors, length(runLst)]);
                    betasUW_ri_reshaped  = reshape(betasUW_ri_rd, [nvox, nregressors, length(runLst)]);
                    
                    mr_beta_hat_ri_reshaped = mean(beta_hat_ri_reshaped, 3);
                    mr_betasUW_ri_reshaped  = mean(betasUW_ri_reshaped, 3);
                    
                    mrm_beta_hat_ri_reshaped = mean(mr_beta_hat_ri_reshaped, 1);
                    mrm_betasUW_ri_reshaped  = mean(mr_betasUW_ri_reshaped, 1);
                else
                    mr_beta_hat_ri_reshaped = zeros(1, nregressors);
                    mr_betasUW_ri_reshaped  = zeros(1, nregressors);
                    
                    mrm_beta_hat_ri_reshaped = zeros(1, nregressors);
                    mrm_betasUW_ri_reshaped  = zeros(1, nregressors);
                end
                
                B_tmp.mrbetasNW{1, 1} = mr_beta_hat_ri_reshaped;
                B_tmp.mrbetasUW{1, 1} = mr_betasUW_ri_reshaped;
                
                %%% for the mean beta in a parcel
                B_tmp.mrmbetasNW = mrm_beta_hat_ri_reshaped;
                B_tmp.mrmbetasUW = mrm_betasUW_ri_reshaped;
                
                B = addstruct(B, B_tmp);
            end % R (regions)

            % save the betas
            dircheck(fullfile(betaDir, subj_name{s}));
            save(fullfile(betaDir, subj_name{s}, sprintf('beta_regions_%s.mat', roi)), 'B', '-v7.3');
            fprintf('******************** betas calculated! ********************\n\n');
        end % sn 
    case 'CHECK:mdtb:task_switch'
        % "transition matrix" for all the subjects is the same (I checked).
        % This case also saves a matrix which contains the task_before and
        % task_after names. This matrix will then be used to calculate the
        % contrasts for the task transitions.
        % Example: sc1_sc2_mdtb('CHECK:mdtb:task_switch', 'sn', [3]); 
        
        sn            = returnSubjs; %% list of subjects
        experiment    = 1;    %% sc1 or sc2?
        ppmethod      = '';   %% 'stc' or ''? The default is set to ''
        glm           = 7;    %% the glm number      
        
        vararginoptions(varargin,{'sn', 'experiment', 'ppmethod', 'deriv', 'glm'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod 
            case ''
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
            case 'stc'
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
        end
        
        % check whether transition matrix for all the subjects is the same
        
        
        for s = sn
            % load in SPM_info
            T = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            [i{s}, j{s}, k{s}] = pivottable(T.taskName_before,T.taskName_after,T.inst,'sum');
        end
        [ind_before, ind_after] = find(~isnan(i{3}));
        all_trans = [j{3}(ind_before), k{3}(ind_after)];
        R = i;
        % check whether transition matrix for all the subjects is the same
        %%% delete empty cells (empty cells are for subjects whose data is not present)
        R = R(~cellfun('isempty',R));
        j = j(~cellfun('isempty',j));
        k = k(~cellfun('isempty',k));
        for l = 1:length(R)
            for m = 1:length(R)
                tmp1 = R{l};
                tmp2 = R{m};
                % replace NaNs with zero
                tmp1(isnan(tmp1)) = 0;
                tmp2(isnan(tmp2)) = 0;
                if isequal(tmp1, tmp2)
                    shared_transition(l, m) = 1;
                else
                    shared_transition(l, m) = 0;
                end
            end % m
        end % l
        % then check the matrix "shared_transition" and look for zeros.
        % there are no zeros!
        % as they are all the same, then you can do the following steps on
        % one of them.
        

        save(fullfile(glmDir, 'all_trans_task_names.mat'), 'all_trans', '-v7.3');
        varargout{1} = i;
        varargout{2} = j;
        varargout{3} = k;
    case 'CHECK:mdtb:cp_cerebellar_masks'
        % just copies the maskbrainSUITGrey.nii (the subject cerebellar
        % mask) to the RegionOfInterest/data/subj_name{s}.
        % Example: sc1_sc2_mdtb('CHECK:mdtb:cp_cerebellar_masks', 'sn', 3);
        
        sn         = returnSubjs; 
        experiment = 1;
        
        vararginoptions(varargin, {'sn', 'experiment'});
        
        experiment = sprintf('sc%d', experiment);
        
        for s = sn
            source = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'maskbrainSUITGrey.nii');
            destDir = fullfile(baseDir, experiment, regDir, 'data', subj_name{s});
            
            copyfile(source, destDir);
            
            fprintf('******************** cerebellar grey mask for %s is copied ********************\n\n', subj_name{s});
        end % sn (subjects)
    case 'CHECK:mdtb:regions_beta'
        % checks whether region betas with different regions lead into the
        % same mean across the whole structure
        % Example: sc1_sc2_mdtb('CHECK:mdtb:regions_beta', 'sn', [2, 3, 4, 6, 8, 9, 10])
        
        sn         = returnSubjs;
        experiment = 1;
        ppmethod   = '';
        roi_old    = 'cerebellum_grey';
        roi_new    = 'Buckner_7';
        glm        = 7;
        
        vararginoptions(varargin, {'sn', 'experiment', 'ppmethod', 'roi_old', 'roi_new', 'glm'});
        
        experiment = sprintf('sc%d', experiment);
        
        switch ppmethod
            case ''
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d', glm)); 
            case 'stc'
                betaDir = fullfile(baseDir, experiment, sprintf('Beta_GLM_%d_stc', glm));
        end
        
        for s = sn
            % load in old
            B_old = load(fullfile(betaDir, subj_name{s}, sprintf('beta_regions_%s.mat', roi_old)));
            B_old = B_old.B;
            % load in new
            B_new = load(fullfile(betaDir, subj_name{s}, sprintf('beta_regions_%s.mat', roi_new)));
            B_new = B_new.B;
            
            % mbetasUW is the field you will be using
            mbetasUW_old = B_old.mbetasUW;
            mbetasUW_new = B_new.mbetasUW;
            
            % calculate the average across the whole structure
            if size(mbetasUW_old, 1) == 1
                m_whole_old = mbetasUW_old;
            else
                m_whole_old = mean(mbetasUW_old);
            end % nparcels = 1
            if size(mbetasUW_new, 1) == 1
                m_whole_new = mbetasUW_new;
            else
                m_whole_new = mean(mbetasUW_new);
            end % npparcels = 1
            
            % calculate the whole mean
            ms_old = mean(m_whole_old); 
            ms_new = mean(m_whole_new);
            
            % plot
            figure; scatterplot(m_whole_old', m_whole_new');
            x = 1:length(m_whole_old);
            y = x;
            hold on;
            plot(x, y);
            title(sprintf('%s, old mean = %d, new mean = %d', subj_name{s}, ms_old, ms_new));
        end % sn
    case 'CHECK:mdtb_rename_con_transitions'
        % contrasts were created but there was a mistake in the code and so
        % not all of them were renamed to con_transitions_. This case is
        % fixing that!
        % Example: sc1_sc2_mdtb('CHECK:mdtb_rename_con_transitions', 'sn', 3);
        
        sn = returnSubjs;
        con_vs = 'rest';
        
        vararginoptions(varargin, {'sn', 'con_vs'});
        
        experiment = 'sc1';
        glmDir     = fullfile(baseDir, experiment, 'GLM_firstlevel_7');
        
        for s = sn
            cd(fullfile(glmDir, subj_name{s}));
            
            conName = {'con','spmT'};
            for i = 1:176
                for n = 1:numel(conName)
                    oldName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_transition_%d.nii',conName{n},i));
                    newName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_transition_%d-%s.nii',conName{n}, i, con_vs));
                    movefile(oldName{i}, newName{i});
                end % conditions (n, conName: con and spmT)
            end % i (contrasts)
            
        end
    case 'CHECK:rename_physio'
        % Example: sc1_sc2_mdtb_backup('CHECK:rename_physio')
        
        sn = [24, 25, 26, 27, 28];
        PhysioDir = fullfile(baseDir, 'Physio');
        
        cd(PhysioDir)
        
        for s = sn
            for sess = 2
                for sca = 1:8
                    cd(fullfile(PhysioDir, subj_name{s}));
                    oldName = sprintf('%s_sess%d_scan%d_lin_reg.mat', subj_name{s}, sess, sca);
                    newName = sprintf('%s_sess%d_scan%d_lin_reg.mat', subj_name{s}, sess, sca+8);
                    movefile(oldName, newName);
                    
%                     oldName = sprintf('%s_sess%d_scan%d_R.mat', subj_name{s}, sess, sca);
%                     newName = sprintf('%s_sess%d_scan%d_R.mat', subj_name{s}, sess, sca+8);
%                     movefile(oldName, newName);
                    
                end
            end
        end
    case 'MDTB:move_files'
        % moving files to the server
        % Example: sc1_sc2_mdtb_backup('MDTB:move_files')
        
        sn         = returnSubjs; 
        ppmethod   = '';
        experiment = 1;
        glm        = 8;
        con_vs     = 'average_2';
        nTrans     = 272;
        copywhich  = 'transitions';
        serverDir  = '/Volumes/MotorControl/data/super_cerebellum_new';

        vararginoptions(varargin, {'sn', 'ppmethod', 'glm', 'experiment', 'con_vs', 'nTrans', 'copywhich'});
          
        switch copywhich
            case 'physio'        % copying physio files for the new glm
                
                sn             = [24, 25];
                experiment     = 'sc1';
                dicomDirServer = fullfile(serverDir, experiment, 'imaging_data_dicom');
                
                for s = sn
                    for sess = 1:2
                        source      = fullfile(dicomDirServer, sprintf('%s_%d', subj_name{s}, sess));
                        destination = fullfile(baseDir, experiment, 'imaging_data_dicom', sprintf('%s_%d', subj_name{s}, sess));
                        dircheck(fullfile(destination, 'physio'));
                        
                        sourceFile = fullfile(source, 'physio');
                        
                        [success(s, sess), Message{s, sess}, ~] = copyfile(sourceFile, fullfile(destination, 'physio'));
                        
                        if success(s, sess) == 1
                            fprintf('physio data for %s sess %d was copied successfully\n', subj_name{s}, sess);
                        else 
                            fprintf('physio data copying failed for %s sess %d\n', subj_name{s}, sess);
                        end
                    end % sess (session)
                end % s (sn)             
            case 'transitions'   % copying files for transitions
                
                experiment = sprintf('sc%d', experiment);
                % copy files from group32k to server
                for tt = 1:nTrans
                    for h = 1:2
                        source = fullfile(baseDir, experiment, 'surfaceWB', 'data', 'group32k');
                        destination = fullfile(serverDir, experiment, 'surfaceWB', sprintf('glm%d', glm), 'group32k');
                        dircheck(destination);
                        
%                         sourceFile = fullfile(source, sprintf('s%s.group.con_transition_%d-%s.func.gii', hemI{h}, tt, con_vs));
                        sourceFile = fullfile(source, sprintf('s%s.con_transition_%d-%s.func.gii', hemI{h}, tt, con_vs));
                        
                        [success(tt, h), Message{tt, h}, ~] = copyfile(sourceFile, destination);
                        if success(tt, h) == 1
                            fprintf('%s coppied to the server\n', sourceFile);
                        else
                            fprintf('copying %s to the server failed\n', sourceFile)
                        end
                    end % h
                end % tt
            case 'taskCon'       % copying files for task contrasts
                
                % load in task information
                C  = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
                Cc = getrow(C, C.StudyNum == experiment);
                
                experiment = sprintf('sc%d', experiment);
                                
                for cc = 1:length(Cc.taskNames)
                    source = fullfile(baseDir, experiment, 'surfaceWB', 'data', 'group32k');
                    destination = fullfile(serverDir, experiment, 'surfaceWB', sprintf('glm%d', glm), 'group32k');
                    dircheck(destination);
                    for h = 1:2
                        sourceFile = fullfile(source, sprintf('s%s.group.con_%s-%s_taskCon.func.gii', hemI{h}, Cc.taskNames{cc}, con_vs));
                        [success(cc, h), Message{cc, h}, ~] = copyfile(sourceFile, destination);
                        
                        if success(cc, h) == 1
                            fprintf('%s coppied to the server\n', sourceFile);
                        else
                            fprintf('copying %s to the server failed\n', sourceFile)
                        end
                        
                        sourceFile = fullfile(source, sprintf('s%s.con_%s-%s_taskCon.func.gii', hemI{h}, Cc.taskNames{cc}, con_vs));
                        [success(cc, h), Message{cc, h}, ~] = copyfile(sourceFile, destination);
                        
                        if success(cc, h) == 1
                            fprintf('%s coppied to the server\n', sourceFile);
                        else
                            fprintf('copying %s to the server failed\n', sourceFile)
                        end
                    end % h 
                end % cc
                
                keyboard
            case 'noise_ceiling' % copying files for noise ceilings
                
                experiment = sprintf('sc%d', experiment);
                % copy files from group32k to server
                for h = 1:2
                    source = fullfile(baseDir, experiment, 'surfaceWB', 'data', 'group32k');
                    destination = fullfile(serverDir, experiment, 'surfaceWB', sprintf('glm%d', glm), 'group32k');
                    dircheck(destination);
                    
                    %                         sourceFile = fullfile(source, sprintf('s%s.group.con_transition_%d-%s.func.gii', hemI{h}, tt, con_vs));
                    sourceFile = fullfile(source, sprintf('s%s.transition_id-%s_noiseCeiling_high.func.gii', hemI{h}, con_vs));
                    
                    [success(tt, h), Message{tt, h}, ~] = copyfile(sourceFile, destination);
                    if success(tt, h) == 1
                        fprintf('%s coppied to the server\n', sourceFile);
                    else
                        fprintf('copying %s to the server failed\n', sourceFile)
                    end
                    
                    sourceFile = fullfile(source, sprintf('s%s.transition_id-%s_noiseCeiling_low.func.gii', hemI{h}, con_vs));
                    
                    [success(tt, h), Message{tt, h}, ~] = copyfile(sourceFile, destination);
                    if success(tt, h) == 1
                        fprintf('%s coppied to the server\n', sourceFile);
                    else
                        fprintf('copying %s to the server failed\n', sourceFile)
                    end
                end % h
            case 'GLM'
                experiment  = sprintf('sc%d', experiment);
                destination = fullfile(serverDir, experiment, sprintf('GLM_firstlevel_%d', glm));

                switch ppmethod
                    case ''
                        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                    case 'stc'
                        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d_stc', glm));
                end
                dircheck(glmDir);
                
                for s = sn
                    sourceFolder = fullfile(glmDir, subj_name{s});
                    destFolder   = fullfile(destination, subj_name{s});
                    dircheck(destFolder);
                    cd(sourceFolder);
                    system(sprintf('find -name -exec cp  -R %s %s;', sourceFolder, destFolder));
%                     [success(s), message{s}, ~] = copyfile(sourceFolder, destFolder, 'f');
%                     system('find -name "*.m" -exec cp {} target \;')
                    
                    if success(s) == 1
                        fprintf('%s coppied to the server\n', sourceFolder);
                    else
                        fprintf('copying %s to the server failed\n', sourceFolder);
                        fprintf('%s\n', message{s});
                    end
                end % sn
            case 'region_cortex'
                experiment = sprintf('sc%d', experiment);
                sourceDir  = fullfile(serverDir, experiment, 'RegionOfInterest', 'data');
                destDir    = fullfile(baseDir, experiment, 'RegionOfInterest', 'data');
                
                for s = sn
                    sourceFile  = fullfile(sourceDir, subj_name{s}, 'regions_cortex.mat');
                    destDirSubj = fullfile(destDir, subj_name{s});
                    [success(s), Message{s}, ~] = copyfile(sourceFile, destDirSubj);
                end % sn      
        end
end
end

% Local functions
function dircheck(dir)
if ~exist(dir,'dir')
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

function [ lcc, lags ] = lcc_vectorized(roi1, roi2, shift)
    %%%% here's a vectorized version for calculation of the lagged
    %%%% cross covariance function between the time series from two
    %%%% ROIs. The code is written so that the ROIs do not have the
    %%%% same number of voxels.
    roi1 = zscore(roi1);
    roi2 = zscore(roi2);

    [~, nv1] = size(roi1); %% lr1 is the length of the time series and na is the number of voxels
    [~, nv2] = size(roi2); 

    lccf_calc = single(zeros(nv1, nv2, (2*shift)+1)); %% initializing the matrix that will contain the lagged correlation values

    lags = -shift:1:shift; %% vector containing the lag values

    for l = 1:length(lags)

        %%% introduce shifts
        %%%% one time series is kept fixed, its length will be shorter
        %%%% though cause it has to have the same length as the shifted
        %%%% one.
        if lags(l) >= 0
            roi1Lagged = roi1(1:end-abs(lags(l)), :);
            roi2Lagged = roi2(abs(lags(l))+1:end, :);

        else 
            roi1Lagged = roi1(abs(lags(l))+1:end,:);
            roi2Lagged = roi2(1:end-abs(lags(l)),:);

        end

        %%% zero-maen the time series before calculating cross-covariance    
        roi1zm = bsxfun(@minus, roi1Lagged, mean(roi1Lagged));
        roi2zm = bsxfun(@minus, roi2Lagged, mean(roi2Lagged));

        %%% calculate cross covariance.
        C = (roi1zm' * roi2zm)./(size(roi1zm, 1)-1);
        cc = C;


        %%% calculate cross correlation    
%             cc = bsxfun(@rdivide, C, std(roi1).*std(roi2));            

        lccf_calc(:, :, l) = cc;
    end
    lcc = lccf_calc;
end

function [Td, Pc, P, I] = td_estimate(lcc_td, lags, noiseLag, Tr)
    %%% Estimating the time delay using lagged cross covariance
    nV1_td_calc = size(lcc_td, 1);
    nV2_td_calc = size(lcc_td, 2);
    %%% using Anish Mitra's method of estimating the delay
    %%% value: parabolic interpolation
    %%%% initializing time delay matrix (TD) and peak correlation matrix (PC)
    tmpTd = zeros(nV1_td_calc, nV2_td_calc); %% will contain the estimated lag
    tmpPc = zeros(nV1_td_calc, nV2_td_calc); %% will contain the estimated peak lag
    tmpP  = zeros(nV1_td_calc, nV2_td_calc); %% will contain the peak of the lccf
    tmpI  = zeros(nV1_td_calc, nV2_td_calc); %% will contain the index of the peak of lccf

    %%%% lag estimation using parabolic interpolation
    for in = 1:nV1_td_calc
        for jn = 1:nV2_td_calc
            tmplcc = double(squeeze(lcc_td(in, jn, :)));

            peak_cov = [];
            peak_lag = [];
            index = [];
            % MAX = 5; %Maximal lag to be considered
            zero = find(lags==0); %Index for zero lag

            %Local Maximum
            if tmplcc(zero) > 0
                [tmpP(in, jn), tmpI(in, jn)] = max(tmplcc);
                index = [tmpI(in, jn)-1 tmpI(in, jn) tmpI(in, jn)+1]; %These are the three x values to be used for parabolic interpolation
            end

            %Local Minimum
            if tmplcc(zero) < 0
                [tmpP(in, jn), tmpI(in, jn)] = min(tmplcc);
                index = [tmpI(in, jn)-1 tmpI(in, jn) tmpI(in, jn)+1]; %These are the three x values to be used for parabolic interpolation
            end

            if isempty(index) == 1
                peak_lag = NaN;
                peak_cov = NaN;
            end

            if isempty(peak_lag)
                if sum(index > 9) ~= 0 || sum(index < 1) ~= 0
                    peak_lag = Inf;
                    peak_cov = Inf;

                else
                    ypoints = tmplcc(index);
                    xpoints = lags(index)';
                    X = [(xpoints.^2) xpoints ones(3,1)]; %Matrix of x-terms
                    constants = X\(ypoints);
                %     constants = mldivide(X, ypoints')
                    peak_lag = -.5*(constants(2)/constants(1)); %From the first derivative
                %     peak_lag = peak_lag*(1/k);
                    peak_cov = (constants(1)*peak_lag^2) + (constants(2)*peak_lag) + constants(3); %peak_cov = abs(peak_cov);
                    peak_lag = Tr * peak_lag; %Get the lag in seconds

                end
            end
            tmpTd(in, jn) = peak_lag;
            tmpPc(in, jn) = peak_cov;
        end
    end

    if abs(peak_lag) > noiseLag
        peak_lag = NaN;
        peak_cov = NaN;
        tmpTd(in, jn) = peak_lag;
        tmpPc(in, jn) = peak_cov;
    end
    
    Td = tmpTd;
    Pc = tmpPc;
    P  = tmpP; %% will contain the peak of the cross-covariance function that will be ued in evaluation.
    I  = tmpI; %% will contain the index corresponding to the peak of the lccf.
end

function [ lpi, ppi ] = lag_projection(td_lp, dim, lpmethod_lp)
    % calculate the lag projection given the structure variable containing
    % the time delay matrices.
    td_tmp = td_lp.td;
    p_tmp  = td_lp.p;
    
    % convert the inf values to NaN
    td_tmp(isinf(td_tmp)) = NaN;
    p_tmp(isinf(p_tmp))   = NaN;
    
    switch lpmethod_lp
        %%% for each case, I will create a matrix that will be
        %%% multiplied by the Td matrix (element-wise) to single out
        %%% the elements in the Td that I want to average across for
        %%% projection.
        case 'marek' %% This case will implement the method used in Marek's 2018 cerebellum paper
            %%%% marek calculated the column-wise average of the Td matrix, so I will multiply the Td matrix by a matrix of all-ones. 
            w_mat    = ones(size(td_tmp))/size(td_tmp, dim); 
        case 'max' %% get the time delay corresponding to maximum correlation coef (before intepolation)
            r = size(p_tmp, 1);
            ic = size(p_tmp, 2);

            w_mat      = zeros(r, ic);
            
            %%% if p_tmp is a square matrix, first you need to zero or NaN
            %%% all the diagonal elements
            if size(p_tmp, 1) == size(p_tmp, 2)
                p_tmp = p_tmp - eye(size(p_tmp));
            end
            
            [~, mxind] = max(abs(p_tmp), [], 'omitnan'); %% the max will be calculated ignoring the nan values.

            w_mat(sub2ind([r, ic], mxind, 1:size(w_mat, 2))) = 1; %% This matrix will cotain 1s in the place where P matrix was max.   
        case 'percentile' %% get the time delay values corresponding to correlation coefs in the percentile defined and average them.
            prc        = 95;
            
            %%% if p_tmp is a square matrix, first you need to zero or NaN
            %%% all the diagonal elements
            if size(p_tmp, 1) == size(p_tmp, 2)
                p_tmp = p_tmp - eye(size(p_tmp));
            end

            Iprc  = prctile(p_tmp, prc, 1);   %% get the 95th percentile of the data in each column of the P matrix
            w_mat = bsxfun(@ge, p_tmp, Iprc); %% will be 1 for every element that is higher than the percentile and zero otherwise.
            w_mat = bsxfun(@rdivide, w_mat, sum(w_mat));
        case 'wighted_average' %% weight the time delay values by the correlation coefs and then average them.
            w_mat      = p_tmp;
    end
    td_tmp     = td_tmp .* w_mat;        
    lpi = nansum(td_tmp, dim);

    switch lpmethod_lp
        case 'wighted_average'
            lpi = bsxfun(@rdivide, lpi, sum(lpi));
    end

    %%% projecting the peak correlation coefficient.
    %%% It will be a measure of how on average an ROI is correlated to
    %%% the rest of ROIs.
    ppi = nanmean(p_tmp);
end

function [Z] = fisher_z(R)
    %%% this case will calculate the fisher z transformation for the
    %%% correlation coefficients.
    %%% z? = .5[ln(1+r) ? ln(1-r)]
    Z = 0.5.*(log(R + ones(size(R))) - log(ones(size(R)) - R));
end