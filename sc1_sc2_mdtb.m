function [ varargout ] = sc1_sc2_mdtb( what, varargin )
%[VARARGOUT] = SC1_SC2(WHAT, VARARGIN) This functions does every analysis
%that can be done on the mdtb dataset.
%   case 'PHYS:mdtb:...' cases used to get the physiological signals
%   case 'GLM:mdtb:desing_glm7' design matrix with conditions modeled as
%   separate regressors, with instructions modeled as separate regressors
%   case 'GLM:mdtb:design_glm8' design matrix with tasks modeled as 30-sec
%   blocks with instructions modeled as separate regressors
%   case 'GLM:mdtb:contrast' case used to get the contrasts
%   case 'GLM:mdtb:contrast_tasks' same as the previous case but it will be
%   used with glm7 (the glm with conditions modeled as separate regressors)
%   case 'SURF:mdtb:map_con' case used to map contrasts to the surface
%   case 'SURF:mdtb:groupmap_con' case used to get the group map for the
%   contrasts.

numDummys = 3;   % number of dummy scans per run
numTRs    = 601; % number of scans per run

%%% setting path for the working directories
%  baseDir = '/Volumes/MotorControl/data/super_cerebellum_new/';
% baseDir = '/Users/ladan/Documents/Project-Cerebellum/Cerebellum_Data';
% baseDir = '/home/ladan/Documents/Data/Cerebellum-MDTB';
% baseDir = '/Users/jdiedrichsen/Data/super_cerebellum_new';
baseDir = '/Volumes/diedrichsen_data$/data/super_cerebellum';
baseDir = '/srv/diedrichsen/data/super_cerebellum/';

%%% setting directory names
behavDir     ='/data';                  %% behavioral data directory.
suitDir      = 'suit';                  %% directory where the anatomicals used in creating flatmaps are stored.
wbDir        = 'surfaceWB';
regDir       = 'RegionOfInterest';      %% The ROI directory 
encodeDir    = 'encoding';

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
    case 'PHYS:mdtb:get_log'    % creates log files for RESP and PULS
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
    case 'PHYS:mdtb:get_reg'    % creates a text file for PULS and RESP regressors
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
    case 'PHYS:mdtb:lin_reg'    % regresses HRV and RVT on the regressor(s) for instructions
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
            
            % construct the design matrix with each instructions modeled as
            % a separate regressor
            % all the design matrix are the same so crating it using one
            % session and run would be enough
            % load in the SPM file (This could take a while)
%             load(fullfile(glmDir, subj_name{s}, 'SPM.mat'));
            load(fullfile(PhysioDir, 'X.mat'));
            
            X_run1     = X(1:598, 1:end - 16);                                        % discarding the intercepts
            Xuse       = X_run1;
            
            % load in SPM_info file
            T      = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            ind    = ((T.sess == 1) & (T.run == 1) & (T.inst == 1) & (T.deriv == 0) ); % get the indices for run1, instructions, non-derivatives
            X_ind  = Xuse(:, ind);
            X_reg = [zeros(3, 16); X_ind];

            phys_cell = cell(length(sess), length(scan)); % preallocating the cell mat that will have all the linear reg results for a subject
            for ss = sess
                for sca = scan
                    if ss == 2
                        sca_ind = sca + 8;
                    elseif ss == 1
                        sca_ind = sca;
                    end
                    
                    % load in HRV and RVT regressors
                    load(fullfile(PhysioDir, subj_name{s}, sprintf('%s_sess%d_scan%d_physio.mat', subj_name{s}, ss, sca_ind)));
                    
                    HRV = physio.model.R(:, 1);
                    RVT = physio.model.R(:, 2);
                    
                    % Do OLS regression for HRV
                    physio_lm.HRV = fitlm(X_reg, HRV);
                    % Do OLS regression for RVT
                    physio_lm.RVT = fitlm(X_reg, RVT);
                    
                    % represent instructions with a single regressor:
                    Xreg_uni = sum(X_reg, 2);
                    physio_lm.HRV_uni = fitlm(Xreg_uni, HRV);
                    physio_lm.RVT_uni = fitlm(Xreg_uni, RVT);
                    
                    save(fullfile(PhysioDir, subj_name{s}, sprintf('%s_sess%d_scan%d_lin_reg.mat', subj_name{s}, ss, sca_ind)), 'physio_lm', '-v7.3');
                    
                    phys_cell{sca_ind} = physio_lm;
                end % sca (scan)
            end % ss (sess)
        end % s (sn)
        varargout{1} = phys_cell;
    case 'PHYS:mdtb:linReg_df'  % creates a dataframe which can be used for visualizations later
        % Example: sc1_sc2_mdtb('PHYS:mdtb:linReg_df')
        
        sn   = [24, 25, 26, 27, 28];
        scan = 1:16;
        
        vararginoptions(varargin, {'sn', 'scan'});
        
        PhysDir = fullfile(baseDir, 'Physio');
        b_reg = [];
        S = [];
        for s = sn
            S.sn = s;
            for sca = scan
                S.run = sca;
                % load in the .mat file for the subject
                a = dir(fullfile(PhysDir, subj_name{s}, sprintf('*scan%d_lin_reg.mat', sca)));
                load(fullfile(PhysDir, subj_name{s}, a.name));
                % for HRV
                hrv_model    = physio_lm.HRV;
                hrvCoefTable = hrv_model.Coefficients;
                S.hrv_b      = hrvCoefTable.Estimate';
                S.hrv_p      = hrvCoefTable.pValue';
                
                % for RVT
                rvt_model    = physio_lm.RVT;
                rvtCoefTable = rvt_model.Coefficients;
                S.rvt_b      = rvtCoefTable.Estimate';
                S.rvt_p      = rvtCoefTable.pValue';
                
                % for HRV using a single regressor for all the instructions
                hrv_uni_model   = physio_lm.HRV_uni;
                hrvUniCoefTable = hrv_uni_model.Coefficients;
                S.hrvUni_b      = hrvUniCoefTable.Estimate';
                S.hrvUni_p      = hrvUniCoefTable.pValue';
                
                % for RVT using a single regressor for all the instructions
                rvt_uni_model   = physio_lm.RVT_uni;
                rvtUniCoefTable = rvt_uni_model.Coefficients;
                S.rvtUni_b      = rvtUniCoefTable.Estimate';
                S.rvtUni_p      = rvtUniCoefTable.pValue';
                
                b_reg = addstruct(b_reg, S);
            end % sca (scan)
        end % s (sn)
        save(fullfile(PhysDir, 'linReg_model.mat'), 'b_reg', '-v7.3')
        varargout{1} = b_reg;
    case 'PHYS:mdtb:linReg_viz' % visualizing the results of linear regression
        % two different types of plots are created.
        % Example: sc1_sc2_mdtb('PHYS:mdtb:linReg_viz')
        
        sn   = [24, 25, 26, 27, 28];
        
        vararginoptions(varargin, {'sn'})
        
        PhysDir    = fullfile(baseDir, 'Physio');
        PhysFigDir = fullfile(PhysDir, 'figure');
        dircheck(PhysFigDir);
        % load in the dataframe
        load(fullfile(PhysDir, 'linReg_model.mat'));
        
        params = {'hrv', 'rvt'};
        
        nSub_phys = 5; %% I only have data for 5 subjects

        % plot type 1: beta vs run for each subject
        for s = sn 
            sub_df = getrow(b_reg, b_reg.sn == s);
            dircheck(fullfile(PhysFigDir, subj_name{s}));
            
            for ireg = 2:17 % the first one is the intercept                
                for physParam = 1:2 % there are two physio param for which I'm doing the linear reg model
                    var = sub_df.(sprintf('%s_b', params{physParam}))(:, ireg);
                    figure; plot(var, '-o', 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', 'b',...
                        'MarkerEdgeColor', 'b');
                    xlabel('run')
                    ylabel('Beta coefficient');
                    hold on
                    % plot a horizental line showing the y = 0
                    y = zeros(1, 16); %% there are 16 runs
                    plot(y, 'r', 'LineWidth', 1.5);
                    
                    %%% find significant betas
                    a = sub_df.(sprintf('%s_p', params{physParam}))(:, ireg) < 0.05;
                    x = sub_df.(sprintf('%s_b', params{physParam}))(:, ireg);
                    x(~a) = NaN;
                    
                    hold on
                    plot(x, 'ro', 'MarkerSize', 10);
                    
                    hold on
                    y = y + max(var) + 0.05;
                    y(~a) = NaN;
                    
                    hold on
                    plot(y, 'k*', 'MarkerSize', 10);
                    
                    title(sprintf('Beta %d vs run for %s', ireg, subj_name{s}));
                    
                    h = gcf;
                    
                    savefig(h, fullfile(PhysFigDir, subj_name{s}, sprintf('%s_%s_Beta%d_vs_run.fig', subj_name{s}, params{physParam}, ireg)));
                    saveas(h, fullfile(PhysFigDir, subj_name{s}, sprintf('%s_%s_Beta%d_vs_run.png', subj_name{s}, params{physParam}, ireg)));
                    close all;
                end % physParam
            end % ireg (regressors)
        end % s (sn) 
        
        % plot type 2: each beta averaged across runs vs subjects
        for ireg = 2:17
            for physParam = 1:2
                var_name = sprintf('%s_b', params{physParam});
                b_table  = tapply(b_reg, {'sn'}, {var_name});
                
                var = b_table.(var_name)(:, ireg);
                figure; plot(var, '-o', 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', 'b', ...
                    'MarkerEdgeColor', 'b');
                xlabel('subject');
                ylabel('Beta coefficient');
                
                hold on;
                y = zeros(1, nSub_phys); %% I only have data for 5 subjects
                plot(y, 'r', 'LineWidth', 1.5);
                
                title(sprintf('averaged Beta %d vs subject', ireg));
                
                h = gcf;
                    
                savefig(h, fullfile(PhysFigDir, sprintf('%s_averaged_Beta%d_vs_run.fig', params{physParam}, ireg)));
                saveas(h, fullfile(PhysFigDir, sprintf('%s_averaged_Beta%d_vs_run.png', params{physParam}, ireg)));
            end % physParams
        end % ireg (regressors)

    case 'GLM:mdtb:design_glm7' % GLM with each condition modelled as a regressor. The instruction for each TASK is also modeled as a separate regressor
        %%% This case will calculate the design matrix with the instruction
        %%% period for each task separated and coming before the task.
        % Example: sc1_sc2_mdtb('GLM:mdtb:design_glm7', 'sn', [20]);
        
        sn             = returnSubjs; %% list of subjects
        experiment_num = 1;           %% sc1 or sc2?
        glm            = 7;           %% the glm number      
        
        vararginoptions(varargin,{'sn', 'experiment_num', 'glm'});
        
        announceTime = 5;
                
        % load in task information
        C  = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM7.txt'));
        Cc = getrow(C, C.StudyNum == experiment_num);
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% SPM and SPM_info files will be saved in glmDir
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        imDir  = 'imaging_data';
        prefix = 'r';
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
            J.timing.fmri_t  = 16; 
            J.timing.fmri_t0 = 1; %% set it to the middle slice
            
            % from Maedbh code: annoying but reorder behavioural runs slightly for 2
            % subjects...
            switch experiment
                case 'sc1'
                    if strcmp(subj_name{s},'s18')
                        %                 runTrue = [51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,60];
                        runTruestr = {'01','02','03','04','05','06','07','08','09','16','10','11','12','13','14','15'};
                    elseif strcmp(subj_name{s},'s21')
                        %                 runTrue = [51,52,53,54,55,56,57,58,59,60,61,63,64,65,66,62];
                        runTruestr = {'01','02','03','04','05','06','07','08','09','10','16','11','12','13','14','15'};
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
                    J.sess(r).cond(ic).name     = 'Instruct';
                    J.sess(r).cond(ic).onset    = instruct_onset(1); % correct start time for numDummys and announcetime included (not for instruct)
                    J.sess(r).cond(ic).duration = 5;  % duration of trials (+ fixation cross) we are modeling
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
                    S.CN    = {Cc.condNames{1}};
                    S.sess  = sess(r);
                
                    T  = addstruct(T,S);

                    % Find Number of Conditions for this task (The conditions for this task have the same instruction)
                    numCond = length(find(Cc.taskNum == Cc.taskNum(ic0)));
                    
                    ic  = ic + 1; 
                    for cond=1:numCond                         
                        D  = dload(fullfile(baseDir, experiment,behavDir, subj_name{s},sprintf('%s_%s_%s.dat', experiment, subj_name{s}, Cc.taskNames{ic0})));
                        R  = getrow(D,D.runNum==runB(r)); % functional runs
                        ST = find(strcmp(P.taskName,Cc.taskNames{ic0}));
                        

                        switch(experiment) % the onsets for sc1 and sc2 are determined differently!
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
                                elseif strcmp(Cc.taskNames{ic0},'motorImagery') || strcmp(Cc.taskNames{ic0},'ToM')
                                    tt = 1;
                                end
                            case 'sc2'
                                  tt = (R.condition==Cc.trialType(ic0));
                                  if strcmp(Cc.taskNames{ic0},'CPRO')
                                    tt = 1;
                                  elseif strcmp(Cc.taskNames{ic0},'ToM2')
                                    tt = 1;
                                  end
                        end % switch experiment
                        
                        switch experiment
                            case 'sc1'
                                if strcmp(Cc.taskNames{ic0},'intervalTiming')
%                                 if strcmp(Cc.taskNames{ic0},'arithmetic') || strcmp(Cc.taskNames{ic0},'intervalTiming') || strcmp(Cc.taskNames{ic0},'checkerBoard')
                                    onset = [P.realStartTime(ST)+R.startTimeReal+announceTime-(J.timing.RT*numDummys)];
                                else
                                    onset = [P.realStartTime(ST)+R.startTimeReal(tt)+announceTime-(J.timing.RT*numDummys)];
                                end
                            case 'sc2'
                                onset = [P.realStartTime(ST)+R.startTimeReal(tt)+announceTime-(J.timing.RT*numDummys)];
                        end
                        

                        % loop through trial-types (ex. congruent or incongruent)
                        J.sess(r).cond(ic).name     = Cc.condNames{ic0};
                        J.sess(r).cond(ic).onset    = onset; % correct start time for numDummys and announcetime included (not for instruct)
                        J.sess(r).cond(ic).duration = ones(1,length(onset))*Cc.duration(ic0);  % duration of trials (+ fixation cross) we are modeling
                        J.sess(r).cond(ic).tmod     = 0;
                        J.sess(r).cond(ic).orth     = 0;
                        J.sess(r).cond(ic).pmod     = struct('name', {}, 'param', {}, 'poly', {});

                        S.SN    = s;
                        S.run   = r;
                        S.inst  = 0; % instruction flag
                        
                        S.instime         = 0;
                        S.taskName_after  = {'none'};
                        S.taskName_before = {'none'};
                        
                        S.task  = Cc.taskNum(ic0);
                        S.cond  = Cc.condNum(ic0);
                        S.CN    = {Cc.condNames{ic0}};
                        S.TN    = {Cc.taskNames{ic0}};
                        S.sess  = sess(r);
                        
                        T  = addstruct(T,S); 
                        
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
            J.bases.hrf.derivs = [0 0];                                  
            J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
            J.volt             = 1;
            J.global           = 'None';
            J.mask             = {fullfile(baseDir, 'sc1', imDir,subj_name{s},'rmask_noskull.nii,1')};
            J.mthresh          = 0.05;
            J.cvi_mask         = {fullfile(baseDir, 'sc1', imDir,subj_name{s},'rmask_gray.nii')};
            J.cvi              =  'fast';
            
%             spm_rwls_run_fmri_spec(J);
            
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            fprintf('******************** glm_%d (SPM.mat) has been saved for %s ********************\n\n',glm, subj_name{s}); 
        end        
    case 'GLM:mdtb:design_glm8' % GLM with each task modeled as a 30 sec block regressor
        % Example: sc1_sc2_mdtb('GLM:mdtb:design_glm8', 'sn', [3]);
        sn             = returnSubjs; %% list of subjects
        experiment_num = 1;           %% sc1 or sc2?
        glm            = 8;           %% the glm number      
        
        vararginoptions(varargin,{'sn', 'experiment_num', 'glm'});
                
        % load in task information
        C     = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));

        Cc    = getrow(C, C.StudyNum == experiment_num);
        Tasks = unique(Cc.taskNames,'rows','stable');                       % get the task names
        nTask      = unique(length(Tasks));                                 % how many tasks there are? for sc1: 18 (including rest) and sc2: 33 (including rest)

        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        announceTime = 5; % there is a 5 sec interval between instruction onset and task onset.
        
        %%% SPM and SPM_info files will be saved in glmDir
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        imDir  = 'imaging_data';
        prefix = 'r';
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
            J.timing.fmri_t  = 16;
            J.timing.fmri_t0 = 1; %% set it to the middle slice
            
            % from Maedbh code: annoying but reorder behavioural runs slightly for 2
            % subjects...
            switch experiment
                case 'sc1'
                    if strcmp(subj_name{s},'s18')
                        %                 runTrue = [51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,60];
                        runTruestr = {'01','02','03','04','05','06','07','08','09','16','10','11','12','13','14','15'};
                    elseif strcmp(subj_name{s},'s21')
                        %                 runTrue = [51,52,53,54,55,56,57,58,59,60,61,63,64,65,66,62];
                        runTruestr = {'01','02','03','04','05','06','07','08','09','10','16','11','12','13','14','15'};
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
                            S.task      = it;
                            S.TN        = {Tasks{it}};
                            S.inst      = 0;
                            S.instOrder = 0;
                            S.time      = onset;
                            S.taskName_after  = {'none'}; % taskName before and after are only defined for instructions
                            S.taskName_before = {'none'};
                            
                            T  = addstruct(T, S);
                            
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
            J.bases.hrf.derivs = [0 0];                                  
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
        sn             = returnSubjs;   %% list of subjects
        glm            = 8;         %% The glm number :)
        experiment_num = 1;
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num'})
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'

        %%% setting the directory paths I need
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        for s = sn
            
            fprintf('******************** estimating glm%d parameters for %s ********************\n', glm, subj_name{s});
            
            glmSubjDir = fullfile(glmDir, subj_name{s});
            load(fullfile(glmSubjDir,'SPM.mat'));
            
            SPM.swd = glmSubjDir;
            spm_rwls_spm(SPM);
            
            fprintf('******************** glm%d parameters estimated for %s ********************\n\n', glm, subj_name{s});
        end %sn        
    
    case 'mdtb:uw'
        % univariately prewhitens the betas
        % Example: sc1_sc2_mdtb('mdtb:uw', 'sn', 2, 'experiment_num', 1, 'glm', 8)
        
        sn             = returnSubjs;
        experiment_num = 1;
        glm            = 8;
        type           = 'beta';        %% it can be set to 'beta' or 'con';
        
        vararginoptions(varargin, {'sn', 'experiment_num', 'glm'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        for s = sn
            
            % load in SPM_info file for task information
            T = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            nBeta = length(T.TN) + 16; % there are 16 intercepts for each run
            % load in the ResMS image
            VresMS = spm_vol(fullfile(glmDir,subj_name{s},'ResMS.nii'));
            ResMS  = spm_read_vols(VresMS);
            ResMS(ResMS==0)=NaN;
            % prewhiten each beta value
            %%% load in each beta image
            for ib = 1: nBeta 
                Vb = spm_vol(fullfile(glmDir, subj_name{s}, sprintf('beta_%04.f.nii', ib)));
                B  = spm_read_vols(Vb);
                
                % divide it by the ResMs
                UWb = bsxfun(@rdivide, B, sqrt(ResMS)); 
                
                % save the new image file
                spm_write_vol(Vb,UWb)
            end % ib (number of betas)
            
        end % s (sn)
    case 'GLM:mdtb:contrast'
        %%% Calculating contrast images.
        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM).
        % This case is written so that it works with both GLM 7 and GLM 8.
        % Reminder: GLM 7 was written with each condition as a separate
        % regressor and a regressor for each of the instructions. GLM 8 was
        % written with each task modeled as a 30 sec block and instructions
        % modeled as a separate regressor.
        % Example1: sc1_sc2_mdtb('GLM:mdtb:contrast', 'sn', [17, 18], 'glm', 8, 'which', 'task')
        % Example2: sc1_sc2_mdtb('GLM:mdtb:contrast', 'sn', [3], 'glm', 72, 'which', 'cond')
        
        sn             = returnSubjs;    %% list of subjects
        glm            = 8;              %% The glm number :)
        experiment_num = 1;
        con_vs         = 'average_4'; %% set it to 'rest' or 'average' (depending on the contrast you want)
        which          = 'cond';      %% it can be set to either cond or task. set it to 'task for GLM_8 and 'cond' for GLM_7
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'con_vs', 'which'})
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
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
                        con(:,logical(T.(which) == ucondition(tt))) = 1;
                        tmp = zeros(1, size(SPM.xX.X, 2));
                        tmp(:, strcmp(T.TN, 'rest')) = 1;
                        sum_rest = sum(tmp);
                        tmp      = tmp./sum_rest;
                        con      = con/abs(sum(con));
                        con      = con - tmp;
                    case 'average_1' % contrast against the average of the other tasks including the instructions
                        % New: Eva, Oct 2nd
                        con                     = zeros(1,size(SPM.xX.X, 2));
                        con(1,logical(T.(which) == ucondition(tt))) = 1./sum(logical(T.(which) == ucondition(tt)));
                        con(1,logical(T.(which) ~= ucondition(tt))) = -1./sum(logical(T.(which) ~= ucondition(tt)));
                    case 'average_2' % contrast against the average of the other tasks not including the instructions
                        con        = zeros(1,size(SPM.xX.X, 2));
                        % TO TRY below - no instructions as contrasts
                        con(1,logical(T.(which) == ucondition(tt))) = 1./sum(logical(T.(which) == ucondition(tt)));
                        con(1,logical((T.(which) ~= ucondition(tt)) .* (T.inst == 0))) = -1./sum(logical((T.(which) ~= ucondition(tt)) .* (T.inst == 0)));
                    case 'average_4' % contrast against the average of all the tasks (not including the instructions)
                        con        = zeros(1,size(SPM.xX.X, 2));
                        % TO TRY below - no instructions as contrasts
                        con(1,logical(T.(which)      == ucondition(tt))) = 1./sum(logical(T.(which) == ucondition(tt)));                        
                        con(1, logical(T.inst == 0)) = con(1, logical(T.inst == 0)) - 1./sum(T.inst == 0);                        
                    case 'rest'
                        con                     = zeros(1,size(SPM.xX.X,2));
                        con(:,logical(T.(which) == ucondition(tt))) = 1;
                        con                     = con/abs(sum(con));
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
        
        sn             = returnSubjs;  %% list of subjects
        glm            = 7;            %% The glm number :)
        experiment_num = 1;
        con_vs         = 'rest';  %% set it to 'rest' or 'average_1' or 'average_2' (depending on the contrast you want)
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'con_vs'})
        
        % gt the task info
        C   = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc  = getrow(C,C.StudyNum == experiment_num);
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
%             SPM  = rmfield(SPM,'xCon');
            T    = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            
            % t contrast for tasks
            utask = unique(T.task);
            
            idx = 1;
            for tt = 1:length(utask) % 0 is "instruct" regressor
                switch con_vs
                    case 'rest' % contrast against rest
                        con                  = zeros(1,size(SPM.xX.X,2));
                        con(:,logical(T.task == utask(tt))) = 1;
                        con                  = con/abs(sum(con));
                    case 'average_4' % contrast against the average of all the tasks or all the other tasks???
                        con        = zeros(1,size(SPM.xX.X, 2));
                        % TO TRY below - no instructions as contrasts
                        con(1,logical(T.task      == utask(tt))) = 1./sum(logical(T.task == utask(tt)));                        
                        con(1, logical(T.inst == 0)) = con(1, logical(T.inst == 0)) - 1./sum(T.inst == 0);   
                    case 'average_1' % contrast vs the average of all the tasks
                        con        = zeros(1,size(SPM.xX.X, 2));
                        con(1,logical(T.task == utask(tt))) = 1./sum(logical(T.task == utask(tt)));
                        con(1,logical(T.task ~= utask(tt))) = -1./sum(logical(T.task ~= utask(tt)));
                    case 'average_2' % contrast vs average of all the tasks except for instructions
                        con        = zeros(1,size(SPM.xX.X, 2));
                        % TO TRY below - no instructions as contrasts
                        con(1,logical(T.task == utask(tt))) = 1./sum(logical(T.task == utask(tt)));
                        con(1,logical(T.task ~= utask(tt) .* (T.inst == 0))) = -1./sum(logical((T.task ~= utask(tt)) .* (T.inst == 0)));
                    case 'average_6' % contrast vs average of all tasks except for instruction for glm4 
                        % SPM_info.mat file for glm4 doesn't have the inst
                        % field.
                        con = zeros(1, size(SPM.xX.X, 2));
                        con(1,logical(T.task == utask(tt))) = 1./sum(logical(T.task == utask(tt)));
                        con(1,logical(T.task ~= utask(tt) .* (T.task ~= 0))) = -1./sum(logical((T.task ~= utask(tt)) .* (T.task ~= 0)));
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
            for i = 1: 17%length(SPM.xCon)
                for n = 1:numel(conName)
                    oldName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%2.4d.nii',conName{n},i));
                    newName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                    movefile(oldName{i}, newName{i});
                end % conditions (n, conName: con and spmT)
            end % i (contrasts)
        end % sn        
    case 'GLM:mdtb:contrast_F'
        % Calculating contrast images for overall F-contrast between
        % tasks / conditions         
        sn             = returnSubjs;    %% list of subjects
        glm            = 8;              %% The glm number :)
        experiment_num = 1;
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num'})
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        for s = sn
            glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm),subj_name{s});
 
            fprintf('******************** calculating F contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, 'SPM.mat'))
            T    = load(fullfile(glmDir, 'SPM_info.mat'));
            
            % F contrast across all tasks
            numTasks = max(T.task); 
            con = zeros(numTasks,size(SPM.xX.X,2));
            for i=1:numTasks
                   con(i,T.task==i & T.inst==0)=1-1/numTasks;
                   con(i,T.task~=i & T.inst==0)=-1/numTasks;
            end
            
            SPM.xCon(1) = spm_FcUtil('Set','AllTask', 'F', 'c',con',SPM.xX.xKXs);
            SPM.swd = glmDir;
            spm_contrasts(SPM,1:length(SPM.xCon));
        end % sn 
    case 'GLM:mdtb:contrast_F_summary'
            % Calculating contrast images for overall F-contrast between
        % tasks / conditions         
        sn             = returnSubjs;    %% list of subjects
        glm            = 8;              %% The glm number :)
        experiment_num = 1;
        D=[]; 
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num'})
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        for s = sn
            glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm),subj_name{s});
            regDir = fullfile(baseDir,'sc1','RegionOfInterest','data',subj_name{s});
            
            fprintf('******************** doing calc %s ********************\n', subj_name{s});
            R{1}=region('image',fullfile(regDir,'cortical_mask_grey_corr.nii'),0.5);
            R{2}=region('image',fullfile(regDir,'regions_cerebellum_suit.nii'),0.5);
            R = region_calcregions(R); 
            
            V= spm_vol(fullfile(glmDir,'spmF_0001.nii')); 
            data = region_getdata(V,R); 
            T.sn =s;
            T.numzero = [sum(data{1}==0), sum(data{2}==0)];
            data{1}=data{1}(data{1}>0);
            data{2}=data{2}(data{2}>0);
            T.avrgF = [mean(data{1}),mean(data{2})]; 
            T.prcF = [prctile(data{1},95),prctile(data{2},95)];
            D=addstruct(D,T); 
        end % sn 
        subplot(1,2,1);
        myboxplot([],D.avrgF); 
        title('Average F-value');
        set(gca,'Ylim',[0 max(D.avrgF(:))+0.2]);
        drawline(1,'dir','horz')
        
        subplot(1,2,2);
        myboxplot([],D.percF)
        set(gca,'Ylim',[0 max(D.percF(:)+1])
        title('95% F-value');
        drawline(finv(0.95,16,1000),'dir','horz')
        varargout={D};  
    case 'SURF:mdtb:map_beta'
        % Maps the betas and the ResMS for a specific GLM to the cortex
        % saves the beta values as gifti files - in the 
        % sc1/surfaceWB/glmX/sxx
        sn             = returnSubjs;
        experiment_num = 1;
        glm            = 7;
        atlas_res      = 32;
        
        vararginoptions(varargin, {'sn', 'experiment_num', 'glm', 'atlas_res'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting glm and surfaceWB directory
        glmDir      = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        numRuns = length(runLst);
        
        for s = sn
            Y = [];
            fprintf('******************** start mapping betas to surface for %s ********************\n', subj_name{s});
            subjSurfDir = fullfile(baseDir, 'sc1', wbDir, subj_name{s});  % Surface subject dir for anatomical information 
            subjSurfGLMDir = fullfile(baseDir, 'sc1', wbDir, sprintf('glm%d',glm),subj_name{s}); % Surface subject dir for functional data for a specific GLM 
            dircheck(subjSurfGLMDir); 
            T = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            
            for h = 1:2 % two hemispheres
                white   = fullfile(subjSurfDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                pial    = fullfile(subjSurfDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                C1      = gifti(white);
                C2      = gifti(pial);
                
                % map ResMS to the surface
                fprintf('******************** mapping ResMS to surface for %s hemi %s ********************\n', hemI{h}, subj_name{s});
                ResMsImage{1}    = fullfile(glmDir, subj_name{s}, 'ResMS.nii');
                ResMs_colName{1} = 'ResMS.nii';
                ResMs_outfile    = fullfile(subjSurfGLMDir, sprintf('%s.%s.%s.ResMS.func.gii', subj_name{s}, hemI{h}, experiment));
                
                G_ResMs = surf_vol2surf(C1.vertices,C2.vertices,ResMsImage,'column_names', ResMs_colName, ...
                        'anatomicalStruct',hemName{h});
                save(G_ResMs, ResMs_outfile);
                
                % Start mapping betas to surface
                filenames = dir(fullfile(glmDir,subj_name{s},'beta*'));
                outfile   = fullfile(subjSurfGLMDir, sprintf('%s.%s.%s.beta.func.gii',subj_name{s}, hemI{h}, experiment));
                
                % Get all the beta files 
                for t = 1:length(filenames)
                    fileList{t}    = {filenames(t).name};
                    column_name{t} = {filenames(t).name};
                end 
                for f=1:length(fileList)
                    images(f) = spm_vol(fullfile(glmDir,subj_name{s},fileList{f}));
                end 
                
                % Map to the surface and save 
                G = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names', column_name, ...
                        'anatomicalStruct',hemName{h});
                save(G, outfile);

                % write out new structure ('Y_info')
                % Y.data = bsxfun(@rdivide, G.cdata', sqrt(G_ResMs.cdata')); % UW betas
                % Y.data(end-numRuns+1:end, :)= []; % deleting the intercepts;
                % Y = addstruct(Y, T);
                % Y.nonZeroInd=B.index';
                
                % outName = fullfile(encodeSubjDir,sprintf('Y_info_glm%d_cortex_%s.mat',glm,hemI{h}));
                % save(outName,'Y','-v7.3');
                % fprintf('cortical vertices: (Y data) computed for %s \n',subj_name{s});
                % clear B R Y
            end % hemi
        end % sn
    case 'SURF:mdtb:map_con'
        % projects individual contrast map volume files for the conditions
        % to the workbench surface.
        % Example: sc1_sc2_mdtb('SURF:mdtb:map_con', 'sn', 2, 'glm', 8, 'experiment_num', 1)
    
        sn             = returnSubjs; %% list of subjects
        atlas_res      = 32;          %% set it to 32 or 164
        experiment_num = 1;           %% enter 1 for sc1 and 2 for sc2
        glm            = 7;           %% glm number
        con_vs         = 'rest'; %% set it to 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting glm and surfaceWB directory
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        glmSurfDir = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm));
        dircheck(glmSurfDir);
        
        for s = sn
            fprintf('******************** start mapping contrasts to surface for %s ********************\n', subj_name{s});
            subjSurfDir = fullfile(baseDir, 'sc1', wbDir, 'data', subj_name{s});
            dircheck(fullfile(glmSurfDir, subj_name{s})); %% directory to save the contrast maps
            
            T  = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            conNames = unique(T.TN, 'stable');
            for h = 1:2 % two hemispheres
                white   = fullfile(subjSurfDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                pial    = fullfile(subjSurfDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                C1      = gifti(white);
                C2      = gifti(pial);
                for ic = 1:length(conNames)
                    
                    conMapName      = sprintf('con_%s-%s', conNames{ic}, con_vs);
                    images{1}       = fullfile(glmDir, subj_name{s}, sprintf('%s.nii', conMapName));
                    column_name{1}  = sprintf('con_%s-%s.nii', conNames{ic}, con_vs);
                    outfile         = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.con_%s-%s.%dk.func.gii', ...
                        subj_name{s}, hemI{h}, conNames{ic}, con_vs, atlas_res));
                    G               = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names', column_name, ...
                        'anatomicalStruct',hemName{h});
                    save(G, outfile);
                    fprintf('******************** mapping to surface for %s hemi %s contrast %s done********************\n', subj_name{s}, ...
                        hemI{h}, conMapName);
                end % ic (condition/contrast)   
            end % hemi
        end % sn    
    case 'SURF:mdtb:map_con_UW'
        % maps the contrasts and ResMS to the surface and also univariately
        % noise normalizes the contrasts using ResMS;
        % Example: sc1_sc2_mdtb('SURF:mdtb:map_con_UW', 'experiment_num', 2, 'glm', 7)
        
        sn             = returnSubjs; %% list of subjects
        atlas_res      = 32;          %% set it to 32 or 164
        experiment_num = 1;           %% enter 1 for sc1 and 2 for sc2
        glm            = 7;           %% glm number
        con_vs         = 'average_4'; %% set it to 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting glm and surfaceWB directory
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        glmSurfDir = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm));
        dircheck(glmSurfDir);
        
        for s = sn
            fprintf('******************** start mapping contrasts to surface for %s ********************\n', subj_name{s});
            subjSurfDir = fullfile(baseDir, 'sc1', wbDir, 'data', subj_name{s});
            dircheck(fullfile(glmSurfDir, subj_name{s})); %% directory to save the contrast maps
            
            T        = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            conNames = unique(T.TN, 'stable');
            for h = 1:2 % two hemispheres
                white   = fullfile(subjSurfDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                pial    = fullfile(subjSurfDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                C1      = gifti(white);
                C2      = gifti(pial);
                
                % map the ResMS to surface
                fprintf('******************** mapping ResMS to surface for %s hemi %s ********************\n', hemI{h}, subj_name{s});
                ResMsImage{1}    = fullfile(glmDir, subj_name{s}, 'ResMS.nii');
                ResMs_colName{1} = 'ResMS.nii';
                ResMs_outfile    = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.ResMS.func.gii', subj_name{s}, hemI{h}));
                
                G_ResMs = surf_vol2surf(C1.vertices,C2.vertices,ResMsImage,'column_names', ResMs_colName, ...
                        'anatomicalStruct',hemName{h});
                save(G_ResMs, ResMs_outfile);
                for ic = 1:length(conNames)
                    conMapName      = sprintf('con_%s-%s', conNames{ic}, con_vs);
                    images{1}       = fullfile(glmDir, subj_name{s}, sprintf('%s.nii', conMapName));
                    column_name{1}  = sprintf('con_%s-%s.nii', conNames{ic}, con_vs);
                    outfile         = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.con_%s-%s.%dk.func.gii', ...
                        subj_name{s}, hemI{h}, conNames{ic}, con_vs, atlas_res));
                    G               = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names', column_name, ...
                        'anatomicalStruct',hemName{h});
                    save(G, outfile);
                    % univariate prewhitening
                    Data(:,ic) = G.cdata ./ sqrt(G_ResMs.cdata);
                    fprintf('******************** mapping to surface for %s hemi %s contrast %s done********************\n', subj_name{s}, ...
                        hemI{h}, conMapName);
                    
                    G_ind       = surf_makeFuncGifti(Data(:,ic),'anatomicalStruct',hemName{h},'columnNames', conNames);
                    outfile_ind = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.wcon_%s-%s.%dk.func.gii', ...
                                        subj_name{s}, hemI{h}, conNames{ic}, con_vs, atlas_res));
                    save(G_ind, outfile_ind);
                end % ic (condition/contrast) 
%                 Data    = bsxfun(@minus,Data,mean(Data,2));
                G_single       = surf_makeFuncGifti(Data,'anatomicalStruct',hemName{h},'columnNames', conNames);
                outfile_single = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.wcon-%s.%dk.func.gii', ...
                    subj_name{s}, hemI{h}, con_vs, atlas_res));
                save(G_single, outfile_single);
            end % hemi
        end % sn 
    case 'SURF:mdtb:map_con_task'
        % projects individual contrast map volume files for tasks to
        % WorkBench surface
        % Example: sc1_sc2_mdtb('SURF:mdtb:map_con_task', 'sn', [3])
    
        sn             = returnSubjs; %% list of subjects
        atlas_res      = 32;          %% set it to 32 or 164
        experiment_num = 1;           %% enter 1 for sc1 and 2 for sc2
        glm            = 4;           %% glm number
        con_vs         = 'average_6'; %% set it to 'rest' or 'average_1' or 'average_2'
        
        % gt the task info
        C   = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc  = getrow(C,C.StudyNum == experiment_num);
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting glmDir and surfaceWB directory
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        glmSurfDir = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm));
        dircheck(glmSurfDir);
        
        for s = sn
            fprintf('******************** start mapping contrasts to surface for %s ********************\n', subj_name{s});
            subjSurfDir = fullfile(baseDir, 'sc1', wbDir, 'data', subj_name{s});
            dircheck(fullfile(glmSurfDir, subj_name{s})); %% directory to save the contrast maps
            
            % get the task names
            conNames = unique(Cc.taskNames, 'stable');
            
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
        % creates group average contrast maps for task contrasts (from Eva)
        % Example: sc1_sc2_mdtb('SURF:mdtb:groupmap_con', 'experiment_num', 2, 'glm', 7, 'which', 'cond')
    
        sn             = returnSubjs; %% list of subjects
        atlas_res      = 32;          %% set it to 32 or 164
        experiment_num = 1;           %% enter 1 for sc1 and 2 for sc2
        glm            = 4;           %% glm number
        replaceNaN     = 1;           %% replacing NaNs
        con_vs         = 'average_4'; %% contrast was calculated against 'rest' or 'average'        
        smooth         = 1;           %% add smoothing
        kernel         = 1;           %% for smoothing
        which          = 'task';      %% 'task' for glm8 and 'cond' for glm7
        normmode       = 'UW';        %% can be set to either 'UW' or 'NW';
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', ...
            'replaceNaN', 'con_vs', 'smooth', 'kernel', 'which', 'normmode'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        switch which
            case 'task' % task for glm8
                conNames = unique(Cc.taskNames, 'stable');
            case 'cond' % condition for glm7
                conNames = unique(Cc.condNames, 'stable');
        end %% do you want the group maps for tasks or conditions
        
        % in 'sc1_sc2_taskConds.txt' file, instruct is not coded as a
        % task/condition name. So I will have to add that to the list of
        % names
        conNames = ['Instruct'; conNames];
        
        experiment = sprintf('sc%d', experiment_num);
        
        % go to the directory where fs_LR atlas is.
        atlasDir         = fullfile(baseDir, sprintf('fs_LR_%d', atlas_res));
        groupSurfDir_glm = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), sprintf('group%dk', atlas_res));
        dircheck(groupSurfDir_glm);
        cd(atlasDir);
        
        for h = 1:2 % two hemispheres
            for cc = 1:length(conNames)
                for s = 1:length(sn)
                    %%% make the group metric file for each contrasts
                    switch normmode
                        case 'UW' % with noise normalization
                            infilenames{s}   = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), ...
                                subj_name{sn(s)},sprintf('%s.%s.wcon_%s-%s.%dk.func.gii', subj_name{sn(s)}, hemI{h}, conNames{cc}, con_vs, atlas_res));
                        case 'NW' % no noise normalization
                            infilenames{s}   = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), ...
                                subj_name{sn(s)},sprintf('%s.%s.con_%s-%s.%dk.func.gii', subj_name{sn(s)}, hemI{h}, conNames{cc}, con_vs, atlas_res));
                    end

                    if smooth
                        surfFile    = fullfile(atlasDir,sprintf('fs_LR.32k.%s.inflated.surf.gii',hemI{h}));
                        surf_smooth(infilenames{s},'surf',surfFile,'kernel',kernel); % smooth outfilenames - it will prefix an 's'
                    end
                    switch normmode
                        case 'UW' % with univariate noise normalization
                            s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.wcon_%s-%s.%dk.func.gii', ...
                                subj_name{sn(s)}, hemI{h}, conNames{cc}, con_vs, atlas_res));
                        case 'NW' % no univariate noise normalization
                            s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s.%dk.func.gii', ...
                                subj_name{sn(s)}, hemI{h}, conNames{cc}, con_vs, atlas_res));
                    end
                end % sn
                
                switch normmode
                    case 'UW' % with noise normalization
                        outfilenames    = fullfile(groupSurfDir_glm,sprintf('%s.wcon_%s-%s.%dk.func.gii',hemI{h},conNames{cc}, con_vs, atlas_res));
                        summaryname     = fullfile(groupSurfDir_glm,sprintf('%s.group.wcon_%s-%s.%dk.func.gii',hemI{h},conNames{cc}, con_vs, atlas_res));
                    case 'NW' % no noise normalization
                        outfilenames    = fullfile(groupSurfDir_glm,sprintf('%s.con_%s-%s.%dk.func.gii',hemI{h},conNames{cc}, con_vs, atlas_res));
                        summaryname     = fullfile(groupSurfDir_glm,sprintf('%s.group.con_%s-%s.%dk.func.gii',hemI{h},conNames{cc}, con_vs, atlas_res));
                end

                surf_groupGiftis(infilenames, 'outfilenames', {outfilenames}, 'groupsummary', summaryname, 'replaceNaNs', replaceNaN);
                if smooth % also save the smoothed versions
                    switch normmode
                        case 'UW' % univariate noise normalization
                            s_outfilenames    = fullfile(groupSurfDir_glm,sprintf('s%s.wcon_%s-%s.%dk.func.gii', hemI{h},conNames{cc}, con_vs, atlas_res));
                            s_summaryname     = fullfile(groupSurfDir_glm,sprintf('s%s.group.wcon_%s-%s.%dk.func.gii', hemI{h},conNames{cc}, con_vs, atlas_res));
                        case 'NW' % no noise normalization
                            s_outfilenames    = fullfile(groupSurfDir_glm,sprintf('s%s.con_%s-%s.func.gii', hemI{h},conNames{cc}, con_vs));
                            s_summaryname     = fullfile(groupSurfDir_glm,sprintf('s%s.group.con_%s-%s.func.gii', hemI{h},conNames{cc}, con_vs));
                    end
                    surf_groupGiftis(s_infilenames, 'outfilenames', {s_outfilenames}, 'groupsummary', s_summaryname, 'replaceNaNs', replaceNaN);
                end
                
%                 % delet the smoothed contrasts for each subject
%                 for s = 1:length(sn)
%                     delete(fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, conNames{cc}, con_vs)));
%                 end
                
                fprintf('******************** group average contrast for %s vs %s is created! ********************\n\n', conNames{cc}, con_vs);
            end % contrasts(ic)
        end % hemi(h)
    case 'SURF:mdtb:groupmap_con_task'
        % creates group average contrast maps for task contrasts
        % Example: sc1_sc2_mdtb('SURF:mdtb:groupmap_con_task', 'sn', [3])
    
        sn             = returnSubjs;   %% list of subjects
        atlas_res      = 32;     %% set it to 32 or 164
        experiment_num = 1;      %% enter 1 for sc1 and 2 for sc2
        glm            = 4;      %% glm number
        replaceNaN     = 1;      %% replacing NaNs
        con_vs         = 'rest'; %% contrast was calculated against 'rest' or 'average'        
        smooth         = 1;      %% add smoothing
        kernel         = 1;      %% for smoothing
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'replaceNaN', 'con_vs', 'smooth', 'kernel'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        taskNames = unique(Cc.taskNames, 'stable');
        
        experiment = sprintf('sc%d', experiment_num);
        
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
                    infilenames{s}   = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                    
                    if smooth
                        surfFile    = fullfile(groupSurfDir,sprintf('fs_LR.32k.%s.inflated.surf.gii',hemI{h}));
                        surf_smooth(infilenames{s},'surf',surfFile,'kernel',kernel); % smooth outfilenames - it will prefix an 's'
                    end
%                     s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                    s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                
                end % sn
%                 outfilenames    = fullfile(groupSurfDir_glm,sprintf('%s.con_%s-%s_taskCon.func.gii',hemI{h},taskNames{cc}, con_vs));
%                 summaryname     = fullfile(groupSurfDir_glm,sprintf('%s.group.con_%s-%s_taskCon.func.gii',hemI{h},taskNames{cc}, con_vs));
                
                outfilenames    = fullfile(groupSurfDir_glm,sprintf('%s.con_%s-%s_taskCon.func.gii',hemI{h},taskNames{cc}, con_vs));
                summaryname     = fullfile(groupSurfDir_glm,sprintf('%s.group.con_%s-%s_taskCon.func.gii',hemI{h},taskNames{cc}, con_vs));
                
                surf_groupGiftis(infilenames, 'outfilenames', {outfilenames}, 'groupsummary', summaryname, 'replaceNaNs', replaceNaN);
                if smooth % also save the smoothed versions
%                     s_outfilenames    = fullfile(groupSurfDir_glm,sprintf('s%s.con_%s-%s_taskCon.func.gii', hemI{h},taskNames{cc}, con_vs));
                    s_outfilenames    = fullfile(groupSurfDir_glm,sprintf('s%s.con_%s-%s_taskCon.func.gii', hemI{h},taskNames{cc}, con_vs));
%                     s_summaryname     = fullfile(groupSurfDir_glm,sprintf('s%s.group.con_%s-%s_taskCon.func.gii', hemI{h},taskNames{cc}, con_vs));
                    s_summaryname     = fullfile(groupSurfDir_glm,sprintf('s%s.group.con_%s-%s_taskCon.func.gii', hemI{h},taskNames{cc}, con_vs));
                    surf_groupGiftis(s_infilenames, 'outfilenames', {s_outfilenames}, 'groupsummary', s_summaryname, 'replaceNaNs', replaceNaN);
                end
                
                % delet the smoothed contrasts for each subject
                for s = 1:length(sn)
                    delete(fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs)));
                end
                
                fprintf('******************** group average contrast for %s vs %s is created! ********************\n\n', taskNames{cc}, con_vs);
            end % contrasts(ic)
        end % hemi(h)
    case 'SURF:mdtb:groupmap_con_groupGiftis'
        % takes in all the giftis for contrasts and make a single gifti 
        % Example: sc1_sc2_mdtb('SURF:mdtb:groupmap_con_groupGiftis', 'experiment_num', 2, 'glm', 7, 'which', 'cond', 'con_vs', 'average_4')
        
        experiment_num = 1;
        glm            = 8;
        con_vs         = 'average_4';
        atlas_res      = 32;
        which          = 'task';
        replaceNaN     = 1;           %% replacing NaNs

        vararginoptions(varargin, {'experiment_num', 'glm', 'con_vs', 'atlas_res', 'which', 'replaceNaN'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        switch which
            case 'task' % task for glm8
                conNames = unique(Cc.taskNames);
            case 'cond' % condition for glm7
                conNames = unique(Cc.condNames);
        end %% do you want the group maps for tasks or conditions
                
        % in 'sc1_sc2_taskConds.txt' file, instruct is not coded as a
        % task/condition name. So I will have to add that to the list of
        % names
        conNames = ['Instruct'; conNames];
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmSurfGroupDir = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), sprintf('group%dk', atlas_res));
        
        %%% creating a single file for each hemisphere
        for h = 1:2
            for cc = 1:length(conNames)
                infilenames{cc}   = fullfile(glmSurfGroupDir,sprintf('s%s.group.wcon_%s-%s.%dk.func.gii', hemI{h}, conNames{cc}, con_vs, atlas_res));
                columnName{cc} = sprintf('%s-%s', conNames{cc}, con_vs);
            end % cc (condition)
            cd(fullfile(glmSurfGroupDir));
            outfilename = sprintf('s%s.group.wcon_%s-%s.%dk.func.gii', hemI{h}, which, con_vs, atlas_res);
            surf_groupGiftis(infilenames, 'outfilenames', {outfilename}, 'outcolnames', columnName, 'replaceNaNs', replaceNaN);
            fprintf('a single gifti file for contrasts for %s hemi successfully created for\n', hemI{h})
        end % h (hemi)      
    case 'SURF:groupmap_meancond'
        % computes the mean across all conditions (relative to rest)
        glm = 4;
        experiment_num = 1;      %% enter 1 for sc1 and 2 for sc2
        atlas_res = 32;     %% set it to 32 or 164
        sn=returnSubjs;
        vararginoptions(varargin,{'glm'});
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        condNames = unique(Cc.condNames, 'stable');       
        experiment = sprintf('sc%d',experiment_num);
        
       groupSurfDir = fullfile(baseDir,experiment,'surfaceWB',sprintf('glm%d',glm),sprintf('group%dk',atlas_res));
       infilenames = cell(1,length(sn));
       for h=1:2 % hemisphere
           fprintf('Working on %s subjects:\n',hemName{h});
           for s=1:length(sn)
               fprintf('%d.',s);
               subjDir = fullfile(baseDir,experiment,'surfaceWB',sprintf('glm%d',glm),subj_name{sn(s)});
               for i=1:size(condNames,1)
                   G = gifti(fullfile(subjDir,sprintf('%s.%s.con_%s-rest.32k.func.gii',subj_name{sn(s)},hemI{h},condNames{i})));
                   allTasks(:,i) = G.cdata(:,1);
               end
               colName = {sprintf('conditionMean_%s',subj_name{sn(s)})};
               G_mean = surf_makeFuncGifti(mean(allTasks,2),'anatomicalStruct',hemName{h},'columnNames',colName);
               infilenames{s} = fullfile(subjDir,sprintf('%s.%s.con_meanConditions-rest.func.gii',subj_name{sn(s)},hemI{h}));
               save(G_mean,infilenames{s});
           end
           outfilenames    = fullfile(groupSurfDir,sprintf('%s.con_meanConditions-rest.func.gii',hemI{h}));
           summaryname     = fullfile(groupSurfDir,sprintf('%s.group.con_meanConditions-rest.func.gii',hemI{h}));    
           surf_groupGiftis(infilenames,'outfilenames',{outfilenames},'groupsummary',summaryname,'replaceNaNs',1);  
           fprintf('...done all.\n');
       end
    
    case 'SUIT:mdtb:suit_parcel2native'
        % maps each atlas of the suit into individual space
        % Example: sc1_sc2_mdtb('SUIT:mdtb:suit_parcel2native', 'sn', [3])
        
        sn             = returnSubjs;
        parcelType     = 'Buckner_7';                          %% set it to Buckner_7 or Buckner_17
        parcelDir      = fullfile(suitToolDir, 'atlasesSUIT'); %% directory where the nifti image is stored
        experiment_num = 1;
        
        vararginoptions(varargin, {'sn', 'parcelType', 'parcelDir', 'experiment_num'});
        
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
        % Example: sc1_sc2_mdtb('SUIT:mdtb:reslice', 'experiment_num', 2, 'glm', 7, 'type', 'beta')
        
        sn             = returnSubjs;            %% list of subjects
        experiment_num = 2;                      %% enter 1 for sc1 and 2 for sc2
        type           = 'beta';                  %% enter the image you want to reslice to suit space
        glm            = 7;                      %% glm number
        mask           = 'cereb_prob_corr_grey'; %% the cerebellar mask to be used:'cereb_prob_corr_grey' or 'cereb_prob_corr' or 'dentate_mask'
        con_vs         = 'average_4';            %% option for the contrasts and spmT
        
        vararginoptions(varargin,{'sn', 'experiment_num', 'glm', 'type', 'mask', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmDir        = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        imageDir      = fullfile(baseDir, experiment, 'imaging_data');
        glmSuitDir    = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
        dircheck(glmSuitDir);
        
        for s = sn
            dircheck(fullfile(glmSuitDir, subj_name{s}));
            glmSubjDir = fullfile(glmDir, subj_name{s});
            outDir     = fullfile(glmSuitDir, subj_name{s});
            
            job.subj.affineTr  = {fullfile(baseDir, 'sc1', suitDir, 'anatomicals', subj_name{s}, 'Affine_c_anatomical_seg1.mat')};
            job.subj.flowfield = {fullfile(baseDir, 'sc1', suitDir, 'anatomicals', subj_name{s},'u_a_c_anatomical_seg1.nii')};
            job.subj.mask      = {fullfile(baseDir, 'sc1', suitDir, 'anatomicals', subj_name{s}, sprintf('%s.nii', mask))};
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
                    source = dir(fullfile(glmSubjDir,sprintf('*%s*%s.nii',images, con_vs))); % images to be resliced
                    cd(glmSubjDir);
                    
                    job.subj.resample = {source.name};
                    suit_reslice_dartel(job);
                case 'spmT'
                    images = 'spmT_';
                    source = dir(fullfile(glmSubjDir,sprintf('*%s*%s.nii',images, con_vs))); % images to be resliced
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
                    
                    images    = 'rrun';
                    
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
    case 'SUIT:mdtb:UW_beta'
        % univariate prewhitening of beta values in the cerebellum
        % Example: sc1_sc2_mdtb('SUIT:mdtb:UW_beta', 'sn', 2, 'experiment_num', 1, 'glm', 8, 'data', 'grey')
        sn             = returnSubjs;
        experiment_num = 1;
        glm            = 8;
        data           = 'grey'; % 'grey_white' or 'grey' or 'grey_nan'
        
        vararginoptions(varargin, {'sn', 'experiment_num', 'glm', 'data'});
        
        experiment = sprintf('sc%d', experiment_num);

        % setting directories
        glmDirSuit = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
        glmDir     = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        numRuns = length(runLst);
        
        for s = sn
            Y = []; % new structure with prewhitened betas
            
            VresMS = spm_vol(fullfile(glmDirSuit, subj_name{s},'wdResMS.nii'));
            ResMS  = spm_read_vols(VresMS);
            ResMS(ResMS==0) = NaN;
            
            % Load over all grey matter mask
            if strcmp(data,'grey')
                V = spm_vol(fullfile(baseDir, experiment, 'suit','anatomicals',subj_name{s},'wdc1anatomical.nii')); % call from sc1
            else
                V = spm_vol(fullfile(baseDir, experiment, 'suit','anatomicals','cerebellarGreySUIT.nii')); % call from sc1
            end
            
            X = spm_read_vols(V);
            % Check if V.mat is the the same as wdResMS!!!
            grey_threshold = 0.1; % grey matter threshold
            indx           = find(X > grey_threshold);
            [i,j,k]        = ind2sub(size(X),indx');
            
            encodeSubjDir = fullfile(baseDir, experiment,encodeDir,sprintf('glm%d',glm),subj_name{s}); dircheck(encodeSubjDir);
            glmSubjDir    = fullfile(glmDir, subj_name{s});
            T             = load(fullfile(glmSubjDir,'SPM_info.mat'));
            
            switch data
                case 'grey'
                    % univariately pre-whiten cerebellar voxels
                    nam={};
                    for b = 1:length(T.SN)+16 % also prewhitening the intercepts
                        nam{1}  = fullfile(glmDirSuit,subj_name{s},sprintf('wdbeta_%2.4d.nii',b));
                        V       = spm_vol(nam{1});
                        B1(b,:) = spm_sample_vol(V,i,j,k,0);
                        B1(b,:) = bsxfun(@rdivide,B1(b,:),sqrt(ResMS(indx)')); % univariate noise normalisation
                    end % b
                    for b = 1:length(T.SN)+16
                        Yy         = zeros(1,V.dim(1)*V.dim(2)*V.dim(3));
                        Yy(1,indx) = B1(b,:);
                        Yy         = reshape(Yy,[V.dim(1),V.dim(2),V.dim(3)]);
                        Yy(Yy==0)  = NaN;
                        idx        = find(~isnan(Yy));
                        Yy         = Yy(:);
                        Bb(b,:)    = Yy(idx,:);
                    end % b
                    clear B1 indx
                    B1   = Bb;
                    indx = idx;
                case 'grey_nan'
                    load(fullfile(glmDirSuit, subj_name{s},'wdBetas_UW.mat'));
                    for b = 1:size(D,1)
                        dat     = D(b,:,:,:);
                        Vi.dat  = reshape(dat,[V.dim(1),V.dim(2),V.dim(3)]);
                        B1(b,:) = spm_sample_vol(Vi,i,j,k,0);
                    end % b
            end % switch data
            
            % write out new structure ('Y_info')
            Y.data = B1;
            Y.data(end-numRuns+1:end,:) = []; % deleting the intercepts
            Y.identity   = indicatorMatrix('identity_p',T.task);
            Y.nonZeroInd = repmat(indx',size(B1,1),1);
            Y.nonZeroInd(end-numRuns+1:end,:) = [];
            
            Y = addstruct(Y, T);
            
            outName = fullfile(encodeSubjDir,sprintf('Y_info_glm%d_%s.mat', glm, data));
            save(outName,'Y','-v7.3');
            fprintf('cerebellar voxels (%s) computed for %s \n', data, subj_name{s});
            clear B1 idx Bb indx
        end % s (sn)
    case 'SUIT:mdtb:map_flatmap'
        % maps betas, contrasts, spmT, and ResMS to flat map.
        % It also maps the univariaately prewhitened betas to the flatmap
        % first run with 'type', 'ResMS' to get the flatmap for the ResMS
        % and then us 'normmode', 'UW' to get the univariately prewhitened
        % betas and contrasts in suit space.
        % Example:sc1_sc2_mdtb('SUIT:mdtb:map_flatmap', 'type', 'beta')
        
        sn             = returnSubjs;
        glm            = 7;
        experiment_num = 2;
        type           = 'beta';      % can be set to 'beta', 'con', and 'ResMS'
        which          = 'cond';      % can be set to 'task' or 'cond'
        con_vs         = 'average_4'; % baseline for the contrast
        normmode       = 'UW';        % can be set to 'UW' or 'NW'
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'type', 'which', 'con_vs', 'normmode'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmSuitDir = fullfile(baseDir, experiment, 'suit', sprintf('glm%d', glm));
        glmDir     = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        for s = sn
            glmSuitSubjDir = fullfile(glmSuitDir, subj_name{s});
            glmSubjDir     = fullfile(glmDir, subj_name{s});

            T = load(fullfile(glmSubjDir,'SPM_info.mat'));
            switch type
                case 'beta'
                    filenames={};
                    for i=1:length(T.run) + 16 % also doing the intercepts
                        filenames{i} = fullfile(glmSuitSubjDir,sprintf('wdbeta_%4.4d.nii',i));
                        colNames{i}  = sprintf('wdbeta_%4.4d',i);
                    end % i
%                     filenames{i+1} = fullfile(glmDir,'ResMS.nii');
                    D = suit_map2surf(filenames,'stats','nanmean');
                    
                    switch normmode
                        case 'UW' % Prewhitening betas
                            %%% load in the ResMS gifti file
                            ResMS_gifti = gifti(fullfile(glmSuitSubjDir,sprintf('%s.Cereb.%s.ResMS.func.gii', subj_name{s}, experiment)));
                            ResMS_data = ResMS_gifti.cdata;
                            
                            wD = bsxfun(@rdivide, D, sqrt(ResMS_data));
                            woutfile = fullfile(glmSuitSubjDir, sprintf('%s.Cereb.%s.wbeta.func.gii',subj_name{s},experiment));
                            
                            wG = surf_makeFuncGifti(wD, 'anatomicalStruct', 'Cerebellum', 'columnNames', colNames); 
                            save(wG,woutfile);
                        case 'NW'
                            outfile = fullfile(glmSuitSubjDir, sprintf('%s.Cereb.%s.beta.func.gii',subj_name{s},experiment)); 
                            G = surf_makeFuncGifti(D, 'anatomicalStruct', 'Cerebellum', 'columnNames', colNames);
                            save(G,outfile);
                    end
                    fprintf('mapped %s %s \n',subj_name{s},experiment);
                case 'ResMS'
                    filenames{1} = fullfile(glmSuitSubjDir,'wdResMS.nii');
                    colNames{1}  = 'wdResMS';
%                     filenames{i+1} = fullfile(glmDir,'ResMS.nii');
                    outfile = fullfile(glmSuitSubjDir, sprintf('%s.Cereb.%s.ResMS.func.gii',subj_name{s},experiment));                    
                    D = suit_map2surf(filenames,'stats','nanmean');
                    
                    G = surf_makeFuncGifti(D, 'anatomicalStruct', 'Cerebellum', 'columnNames', colNames);                    
                    save(G,outfile);
                    
                    fprintf('mapped %s %s \n',subj_name{s},experiment);
                case 'con'
                    switch which
                        case 'task'
                            CN = 'TN';
                        case 'cond'
                            CN = 'CN';
                    end
                    condNames = unique(T.(CN), 'stable');
                    
                    filenames={};
                    for i=1:length(condNames)
                        filenames{i} = fullfile(glmSuitSubjDir,sprintf('wdcon_%s-%s.nii', condNames{i}, con_vs));
                        colNames{i}  = sprintf('%s-%s', condNames{i}, con_vs);
                    end % i
%                     filenames{i+1} = fullfile(glmDir,'ResMS.nii');
                    D = suit_map2surf(filenames,'stats','nanmean');
                    switch normmode
                        case 'UW' % Prewhitening betas
                            %%% load in the ResMS gifti file
                            ResMS_gifti = gifti(fullfile(glmSuitSubjDir,sprintf('%s.Cereb.%s.ResMS.func.gii', subj_name{s}, experiment)));
                            ResMS_data = ResMS_gifti.cdata;
                            
                            wD = bsxfun(@rdivide, D, sqrt(ResMS_data));
                            woutfile = fullfile(glmSuitSubjDir, sprintf('%s.Cereb.%s.wcon-%s.func.gii',subj_name{s},experiment, con_vs));
                            
                            wG = surf_makeFuncGifti(wD, 'anatomicalStruct', 'Cerebellum', 'columnNames', colNames);
                            save(wG,woutfile);
                        case 'NW'
                            outfile = fullfile(glmSuitSubjDir, sprintf('%s.Cereb.%s.con-%s.func.gii',subj_name{s},experiment, con_vs));
                            G = surf_makeFuncGifti(D, 'anatomicalStruct', 'Cerebellum', 'columnNames', colNames);
                            save(G,outfile);
                    end
                    
                    fprintf('mapped %s %s \n',subj_name{s},experiment);                   
            end
        end % s (sn) 
    case 'SUIT:mdtb:groupmap_con'
        % creates group average for the condition contrast maps.
        % you need to reslice all the images to suit space before running
        % this case
        % Example: sc1_sc2_mdtb('SUIT:mdtb:groupmap_con', 'sn', [2, 3, 4, 6, 8, 9, 10, 12, 14, 15]);
        
        sn             = returnSubjs;        %% list of subjects
        experiment_num = 2;                  %% enter 1 for sc1 and 2 for sc2
        type           = 'con';              %% enter the image you want to reslice to suit space
        glm            = 7;                  %% glm number
        con_vs         = 'average_4';        %% is the contrast calculated vs 'rest' or 'average'
        which          = 'cond';             %% you may choose 'cond' or 'task'
        
        vararginoptions(varargin,{'sn', 'experiment_num', 'glm', 'type', 'which', 'con_vs'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        switch which
            case 'task' % task for glm8
                conNames = unique(Cc.taskNames, 'stable');
            case 'cond' % condition for glm7
                conNames = unique(Cc.condNames, 'stable');
        end %% do you want the group maps for tasks or conditions
        
        % in 'sc1_sc2_taskConds.txt' file, instruct is not coded as a
        % task/condition name. So I will have to add that to the list of
        % names
        conNames = ['Instruct'; conNames];
        
        experiment = sprintf('sc%d', experiment_num);
        
        % Setting directories
        glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
        glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), 'group');
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
            
            save(G, fullfile(glmSuitGroupDir, sprintf('Cereb.con.con_%s-%s.func.gii', conNames{cc}, con_vs)));
            
            save(G, fullfile(glmSuitGroupDir, sprintf('Cereb.group.con_%s-%s.func.gii', conNames{cc}, con_vs)));
            save(fullfile(glmSuitGroupDir, sprintf('indMaps_%s_%s-vs-%s.mat', type, conNames{cc}, con_vs)), 'maps', '-v7.3');
            fprintf('******************** %s group average for %s vs %s is created! ********************\n\n', type, conNames{cc}, con_vs);
            
        end % contrasts (cc)
        save(fullfile(glmSuitGroupDir, sprintf('indMaps_%s-vs-%s.mat', type, con_vs)), 'maps', '-v7.3');
    case 'SUIT:mdtb:groupmap_con2'
        % creates group average for the condition contrast maps.
        % you need to reslice all the images to suit space before running
        % this case
        % Example: sc1_sc2_mdtb('SUIT:mdtb:groupmap_con2', 'sn', [2, 3, 4, 6, 8, 9, 10, 12, 14, 15]);
        
        sn             = returnSubjs;        %% list of subjects
        experiment_num = 2;                  %% enter 1 for sc1 and 2 for sc2
        type           = 'con';              %% enter the image you want to reslice to suit space
        glm            = 7;                  %% glm number
        con_vs         = 'average_4';        %% is the contrast calculated vs 'rest' or 'average'
        which          = 'cond';             %% you may choose 'cond' or 'task'
        
        vararginoptions(varargin,{'sn', 'experiment_num', 'glm', 'type', 'which', 'con_vs'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        switch which
            case 'task' % task for glm8
                conNames = unique(Cc.taskNames, 'stable');
            case 'cond' % condition for glm7
                conNames = unique(Cc.condNames, 'stable');
        end %% do you want the group maps for tasks or conditions
        
        % in 'sc1_sc2_taskConds.txt' file, instruct is not coded as a
        % task/condition name. So I will have to add that to the list of
        % names
        conNames = ['Instruct'; conNames];
        
        experiment = sprintf('sc%d', experiment_num);
        
        % Setting directories
        glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
        glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), 'group');
        dircheck(glmSuitDir);
        dircheck(glmSuitGroupDir);
        
        for s = 1:length(sn)
            tmp = gifti(fullfile(glmSuitDir, subj_name{sn(s)},...
                sprintf('%s.Cereb.%s.wcon-%s.func.gii', subj_name{sn(s)}, experiment, con_vs)));
            G(:, :, s) = tmp.cdata;
        end % sn
        groupG = nanmean(G, 3); % calculating group average for each condition/task
        
        G = surf_makeFuncGifti(groupG, 'anatomicalStruct', 'Cerebellum', 'columnNames', conNames);
        
        save(G, fullfile(glmSuitGroupDir, sprintf('Cereb.group.%s.wcon-%s.func.gii', experiment, con_vs)));
        fprintf('******************** %s group average for contrastsvs %s is created! ********************\n\n', type, con_vs);
    case 'SUIT:mdtb:groupmap_con_task'
        % creates group map for task contrasts
        % Example: sc1_sc2_mdtb('SUIT:mdtb:groupmap_con_task', 'sn', [3]);
        
        sn             = returnSubjs;                   %% list of subjects
        experiment_num = 1;                      %% enter 1 for sc1 and 2 for sc2
        type           = 'con';                  %% enter the image you want to reslice to suit space
        glm            = 7;                      %% glm number
        con_vs         = 'rest';                 %% is the contrast calculated vs 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'experiment_num', 'glm', 'type', 'mask'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        taskNames = unique(Cc.taskNames, 'stable');
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
        glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), 'group');
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
    case 'SUIT:mdtb:groupmap_con_groupGiftis'
        % Example:sc1_sc2_mdtb('SUIT:mdtb:groupmap_con_groupGiftis')
        
        experiment_num = 1;
        glm            = 7;
        con_vs         = 'average_4';
        which          = 'cond';
        replaceNaN     = 1;           %% replacing NaNs

        vararginoptions(varargin, {'experiment_num', 'glm', 'con_vs', 'which', 'replaceNaN'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        switch which
            case 'task' % task for glm8
                conNames = unique(Cc.taskNames);
            case 'cond' % condition for glm7
                conNames = unique(Cc.condNames);
        end %% do you want the group maps for tasks or conditions
                
        % in 'sc1_sc2_taskConds.txt' file, instruct is not coded as a
        % task/condition name. So I will have to add that to the list of
        % names
        conNames = ['Instruct'; conNames];
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), 'group');
        
        %%% creating a single file for each hemisphere
        for cc = 1:length(conNames)
            infilenames{cc}   = fullfile(glmSuitGroupDir,sprintf('Cereb.group.con_%s-%s.func.gii', conNames{cc}, con_vs));
            columnName{cc} = sprintf('%s-%s', conNames{cc}, con_vs);
        end % cc (condition)
        dircheck(glmSuitGroupDir)
        cd(fullfile(glmSuitGroupDir));
        outfilename = sprintf('Cereb.group.con_%s-%s.func.gii', which, con_vs);
        surf_groupGiftis(infilenames, 'outfilenames', {outfilename}, 'outcolnames', columnName, 'replaceNaNs', replaceNaN);
        fprintf('a single gifti file for contrasts for cerebellum successfully created\n')
        
    case 'Summ:mdtb:beta_dataframe'
        % creates a dataframe that has all the task information and the
        % univariately prewhitened beta values in cortex and the
        % cerebellum.
        % Example: sc1_sc2_mdtb('Summ:mdtb:beta_dataframe', 'sn', 2)
        
        sn             = returnSubjs;
        experiment_num = 1;
        glm            = 8;
        data_cereb     = 'grey';
        centre         = 'average'; % how do you want to center the betas ('average', 'rest')
        
        vararginoptions(varargin, {'sn', 'experiment_num', 'glm', 'centre'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        encodingDir = fullfile(baseDir, experiment, encodeDir, sprintf('glm%d', glm));

        hemis_num  = [1, 2];     % looping through two hemispheres (L = 1, R = 2, also cereb is set to 0)
        hemis      = {'L', 'R'};
        structures_num = [1, 2]; % looping through two structures (cortex = 1, cereb = 2)
        
        df = []; % dataframe that has all the variables for doing the scatter plots
        for s = sn 
            fprintf('getting data for %s\n', subj_name{s})
            temp.sn = s;
            for sc = 1:length(structures_num)
                temp.structure = structures_num(sc);
                if structures_num(sc) == 1 % cortex
                    for ih = 1:length(hemis)
                        temp.hemi = hemis_num(ih);
                        load(fullfile(encodingDir, subj_name{s}, sprintf('Y_info_glm%d_cortex_%s.mat', glm, hemis{ih})));
                        % remove unnecessary fields
                        y_rf = rmfield(Y, {'time', 'taskName_before', 'taskName_after', 'sess', 'TN', 'SN'});
                        
                        % get the baseline for betas
                        switch centre % centre against what?
                            case 'average' % centre against the average of all the tasks except for instructions
                                % remove the instructions
                                y_rf_ri = getrow(y_rf, y_rf.task ~= 0);
                                
                                % group the data with tapply
                                t_run = tapply(y_rf_ri, {'task', 'inst', 'instOrder'}, {'data'});
                                baseline = mean(t_run.data);
                            case 'rest'    % centre against the rest (for glm8 it is explicitly moedled)
                        end % switch centre
                        t_runa  = tapply(y_rf, {'task', 'inst', 'instOrder'}, {'data'});
                        uwBeta  = t_runa.data;
                        centred = bsxfun(@minus, uwBeta, baseline);
                        
                        % get the average beta value for structures{sc}
                        avgBeta = mean(centred, 2);
                        temp.mBeta = avgBeta';
%                         clear Y y_rf y_rf_ri t_run t_runa uwBeta centred avgBeta
                        df = addstruct(df, temp);
                    end % ih (hemispheres)
                else % the structure is the cerebellum and you don't need to loop through hemispheres
                    temp.hemi = 0;
                    load(fullfile(encodingDir, subj_name{s}, sprintf('Y_info_glm%d_%s.mat', glm, data_cereb)));
                    % remove unnecessary fields
                    y_rf = rmfield(Y, {'time', 'taskName_before', 'taskName_after', 'sess', 'TN', 'SN', ...
                        'identity', 'nonZeroInd'});
                    
                    % get the baseline for betas
                    switch centre % centre against what?
                        case 'average' % centre against the average of all the tasks except for instructions
                            % remove the instructions
                            y_rf_ri = getrow(y_rf, y_rf.task ~= 0);
                            
                            % group the data with tapply
                            t_run = tapply(y_rf_ri, {'task', 'inst', 'instOrder'}, {'data'});
                            baseline = nanmean(t_run.data);
                        case 'rest'    % centre against the rest (for glm8 it is explicitly moedled)
                    end % switch centre
                    t_runa  = tapply(y_rf, {'task', 'inst', 'instOrder'}, {'data'});
                    uwBeta  = t_runa.data;
                    centred = bsxfun(@minus, uwBeta, baseline);
                    
                    % get the average beta value for structures{sc}
                    avgBeta = mean(centred, 2);
                    temp.mBeta = avgBeta';
                    %                         clear Y y_rf y_rf_ri t_run t_runa uwBeta centred avgBeta
                    df = addstruct(df, temp);
                end % if the structure is cortex then you need to loop through hemispheres
            end % sc (structures)
        end % s (sn)
        keyboard

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
        
    case 'Houskeeping:renameSPM'
        % rename SPM directories
        % Example: sc1_sc2_mdtb('Houskeeping:renameSPM', 'experiment_num', 2, 'glm', 8)
        
        sn             = returnSubjs;
        experiment_num = 2;
        glm            = 7;
        
        vararginoptions(varargin, {'sn', 'experiment_num', 'glm'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        glmDir     = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        imagingDir = fullfile(baseDir, 'sc1/imaging_data'); %% all the imaging files are in sc1
        
        for s = sn
            fprintf('******************** changing directories for %s ********************\n', subj_name{s});
            newGLMDir   = fullfile(glmDir,subj_name{s});
            newRawDir   = fullfile(imagingDir,subj_name{s});
            
            % load SPM file
            load(fullfile(newGLMDir, 'SPM.mat'));
            
            SPM         = spmj_move_rawdata(SPM,newRawDir);
            SPM.swd     = fullfile(newGLMDir);
            save(fullfile(newGLMDir,'SPM.mat'),'SPM','-v7.3');
            varargout{1} = SPM;
        end % s (sn)
    case 'Houskeeping:move_files'
        % moving files to the server
        % Example: sc1_sc2_mdtb('Houskeeping:move_files', 'copywhich', 'SURF_files')
        
        sn             = returnSubjs;
        experiment_num = 1;
        glm            = 8;
        con_vs         = 'rest_task';
        nTrans         = 272;
        copywhich      = 'SURF_files';
        serverDir      = '/Volumes/MotorControl/data/super_cerebellum_new';
        atlas_res      = 32;
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'con_vs', 'nTrans', 'copywhich',...
            'atlas_res', 'serverDir'});
        
        experiment = sprintf('sc%d', experiment_num);
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
                
                % setting directories
                glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
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
                        % fprintf('%s\n', message{s});
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
            case 'contrast_spmT' % copying contrast and spmT files
                glmDir_local  = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                glmDir_server = fullfile(serverDir, experiment, sprintf('GLM_firstlevel_%d', glm));
                for s = sn
                    subj_local  = fullfile(glmDir_local, subj_name{s});
                    subj_server = fullfile(glmDir_server, subj_name{s});
                    
                    cop_files = dir(fullfile(subj_local, sprintf('*%s.nii', con_vs)));
                    
                    for icf = 1:length(cop_files)
                        
                        sourceFile = fullfile(subj_local, cop_files(icf).name);
                        destination = subj_server;
                        
                        [success(s, icf), Message{s, icf}, ~] = copyfile(sourceFile, destination);
                        if success(s, icf) == 1
                            fprintf('%s coppied to the server\n', cop_files(icf).name);
                        else
                            fprintf('copying %s to the server failed\n', cop_files(icf).name)
                        end
                    end % icf (files to be copied)
                    
                end % s (sn)
            case 'SURF_files' % copies contrast map and betas to the surface
                
                sourceDir = fullfile(baseDir, experiment, 'surfaceWB', sprintf('glm%d', glm));
                
                for s = sn
                    for h = 1:2
                        source = fullfile(sourceDir, subj_name{s});
                        destination = fullfile(serverDir, experiment, 'surfaceWB', sprintf('glm%d', glm), subj_name{s});
                        dircheck(destination);
                        
                        %                         sourceFile = fullfile(source, sprintf('s%s.group.con_transition_%d-%s.func.gii', hemI{h}, tt, con_vs));
                        sourceFile = fullfile(source, sprintf('%s.%s.wcon-%s.func.gii', subj_name{s}, hemI{h}, con_vs));
                        
                        [success(s, h), Message{s, h}, ~] = copyfile(sourceFile, destination);
                        if success(s, h) == 1
                            fprintf('%s coppied to the server\n', sourceFile);
                        else
                            fprintf('copying %s to the server failed\n', sourceFile)
                        end
                        
                        clear sourceFile success Message
                        
                        sourceFile = fullfile(source, sprintf('%s.%s.%s.beta.%dk.func.gii', subj_name{s}, hemI{h}, experiment, atlas_res));
                        
                        [success(s, h), Message{s, h}, ~] = copyfile(sourceFile, destination);
                        if success(s, h) == 1
                            fprintf('%s coppied to the server\n', sourceFile);
                        else
                            fprintf('copying %s to the server failed\n', sourceFile)
                        end
                        
                    end % h
                end % s (sn)
                
                clear source destination sourceFile
                % copying group files
                source = fullfile(sourceDir, sprintf('group%dk', atlas_res));
                destination = fullfile(serverDir, experiment, 'surfaceWB', sprintf('glm%d', glm), sprintf('group%dk', atlas_res));
                
                for h = 1:2
                    sourceFile = fullfile(source, sprintf('s%s.group.wcon_task-%s.%dk.func.gii', hemI{h}, con_vs, atlas_res));
                    
                    [successg(h), Messageg{h}, ~] = copyfile(sourceFile, destination);
                    if successg(h) == 1
                        fprintf('%s coppied to the server\n', sourceFile);
                    else
                        fprintf('copying %s to the server failed\n', sourceFile)
                    end
                end % h(hemi)    
        end
        
    case 'CHECK:mdtb:glms'
        % checking the regressors for all of the tasks and conditions in
        % all runs.
        % Example: sc1_sc2_mdtb('CHECK:mdtb:glms')
        sn             = 30;
        glm            = 7;
        run_num        = 5;
        experiment_num = 1;
        which          = 'CN';
        what_plot      = 'irun'; % can be set to 'irun'
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'run_num', 'what_plot'})
        
        experiment = sprintf('sc%d', experiment_num);
        
        glmSubjDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm), subj_name{sn});
        
        if glm == 4
            which = 'TN';
        elseif glm == 8
            which = 'TN';
        elseif glm == 7
            which = 'CN';
        end % glm4 info file is created differently
        % load in SPM
        load(fullfile(glmSubjDir, 'SPM.mat'));
        X = SPM.xX.X;
        % discarding the intercepts
        X = X(:, 1:end - 16);
        
        % load in SPM_info.mat
        T = load(fullfile(glmSubjDir, 'SPM_info.mat'));
        
        condNames = unique(T.(which), 'stable');
        
        switch what_plot
            case 'ireg' % plots a regressor across all the runs
                for i = 1:length(condNames)
                    tindx = strcmp(T.(which), condNames{i});
                    reg = X(:, tindx);
                    
                    figure;
                    plot(reg);
                    title(sprintf('regressor for glm %d %s %s %s', glm, which, condNames{i}, subj_name{sn}));
                end % i (conditions)
            case 'irun' % for each run, plots all the regressors
                for ir = runLst
                    a = (ir - 1) * 598 + 1;
                    b = ir * 598;
                    rindx = T.run == ir;
                    regs  = X(a:b, rindx);
                    
                    figure; 
                    plot(regs);
                    title(sprintf('all the regressors for run %d glm %d %s %s', ir, glm, which, subj_name{sn}));
                end % ir (runs)
        end
        keyboard; 
    case 'CHECK:mdtb:boxcars'
        % Checking the onsets and durations of the tasks and conditions
        % identified using the previous case:
        % stroop cong, stroop incong, Obj2back, Obj0back, Verbal2Back, Verbal0Back
        % IntervalTiming, Happy faces, Sad faces, Pleasant scenes,
        % Unpleasant scenes, go, NoGo.
        % Example: sc1_sc2_mdtb('CHECK:mdtb:data')
        
        sn             = 20; % subject that will be used for checking
        experiment_num = 1;
%         glm            = 7;
        all            = 0; % set this to 0 if you don't want to plot all the conditions
        
        vararginoptions(varargin, {'sn', 'experiment_num', 'glm', 'all'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        announceTime = 5;
        J.timing.RT  = 1.0;
        % load in task information
        C  = dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM7.txt'));
        Cc = getrow(C, C.StudyNum == experiment_num);
        
        switch all
            case 1
                condList = Cc.condNames;
                condList(1) = [];
            case 0
%                 condList = {'NoGo', 'Go', 'UnpleasantScenes', 'PleasantScenes', 'SadFaces', 'HappyFaces', ...
%                     'IntervalTiming', 'Object0Back', 'Object2Back', 'StroopIncon', 'StroopCon', ...
%                     'Verbal0Back', 'Verbal2Back', 'CheckerBoard'};
                condList = {'Math', 'DigitJudgement', 'IntervalTiming', 'CheckerBoard'};
        end
        
        % setting the directories
        bDir     = fullfile(baseDir, experiment, 'data');
        bSubjDir = fullfile(bDir, subj_name{sn});
        
        % get the onsets and durations
        A = dload(fullfile(bSubjDir,sprintf('%s_%s.dat', experiment, subj_name{sn})));
        A = getrow(A, A.runNum >= funcRunNum(1) & A.runNum <= funcRunNum(2));
        
        for r = 1
            for ic = 1:length(condList)
                figure;
                fprintf('condition is %s\n', condList{ic})
                P  = getrow(A,A.runNum == runB(r));
                TN = unique(Cc.taskNames(strcmp(Cc.condNames,condList{ic})));
                ST = find(strcmp(P.taskName,Cc.taskNames(strcmp(Cc.condNames,condList{ic}))));
                D  = dload(fullfile(bSubjDir,sprintf('%s_%s_%s.dat', experiment, subj_name{sn}, TN{1})));
                R  = getrow(D,D.runNum==runB(r)); % functional runs
                
                if isfield(R,'trialType')
                    tt = (R.trialType==Cc.trialType(strcmp(Cc.condNames,condList{ic})));
                else
                    tt = Cc.trialType(strcmp(Cc.condNames,condList{ic}));
                end
                if strcmp(Cc.taskNames(strcmp(Cc.condNames,condList{ic})),'visualSearch')
                    tt = (R.setSize==Cc.trialType(strcmp(Cc.condNames,condList{ic})));
                elseif strcmp(Cc.taskNames(strcmp(Cc.condNames,condList{ic})),'nBack') || strcmp(Cc.taskNames(strcmp(Cc.condNames,condList{ic})),'nBackPic')
                    tt = (R.respMade==Cc.trialType(strcmp(Cc.condNames,condList{ic})));
%                 elseif strcmp(Cc.taskNames(strcmp(Cc.condNames,condList{ic})),'motorImagery') || strcmp(Cc.taskNames(strcmp(Cc.condNames,condList{ic})),'ToM') || strcmp(Cc.taskNames(strcmp(Cc.condNames,condList{ic})),'checkerBoard')
                elseif strcmp(Cc.taskNames(strcmp(Cc.condNames,condList{ic})),'motorImagery') || strcmp(Cc.taskNames(strcmp(Cc.condNames,condList{ic})),'ToM')
                    tt = 1;
                end
                
                % Calculate the onset 
%                 onset = ([P.realStartTime(ST)+R.startTimeReal(tt)+announceTime-(J.timing.RT*numDummys)])';
                onset = ([P.realStartTime(ST)+R.startTimeReal+announceTime-(J.timing.RT*numDummys)])';
                % get the duration
                duration = Cc.duration(strcmp(Cc.condNames,condList{ic}));
%                 fprintf('duration for %s is %d\n', condList{ic}, duration);
                for io = 1:length(onset)
                    t = onset(io):(onset(io) + duration);
                    y = heaviside(t);
                    
                    subplot(211)
                    plot(t, y, 'LineWidth', 1.5); 
                    ylim([0, 1.2])
                    hold on
                    title(sprintf('boxcar for %s, duration %d', condList{ic}, duration))
                    
%                     subplot(312)
%                     stem(onset', zeros(1, length(onset)));
%                     title(sprintf('the onsets for %s', condList{ic}));

                end % io
                if isfield(R, 'trialType')
                        subplot(212)
                        stem(R.trialType)
                        title(sprintf('trial types for %s', condList{ic}));
                end % trying to figure out if there's sth wrong with the trial types
                
                
%                 subplot(211)
%                 plot(boxcar(min(trial_times)-5:end));
                
                
%                 canonical_hrf  = spm_hrf(1); % where 2 represents the TR (in seconds)
%                 stick_function = zeros(1,100);
%                 trial_times=round(onset)+ duration;
%                 stick_function(trial_times) = 1;
%                 subplot(211)
%                 plot(stick_function);
%                 xlim([min(onset)- 50, max(onset)+50]);
%                 title(sprintf('boxcars for %s', condList{ic}));
%                 
%                 regs = conv(stick_function, canonical_hrf);
%                 subplot(212)
%                 plot(regs);
%                 xlim([min(onset)- 50, max(onset)+50]);
%                 title(sprintf('regressor for %s', condList{ic}));
%                 keyboard;
                
            end % ic (conditions)
        end % r (runs)  
    case 'CHECK:mdtb'
        % figuring out which regressor is not correct. 
        % Example: sc1_sc2_mdtb('CHECK:mdtb', 'sn', 2);
        
        sn = 3; 
        experiment_num = 1;
        glm = 7;
        r   = 3;
        vararginoptions(varargin, {'sn', 'experiment_num', 'glm', 'r'});
        
        experiment = sprintf('sc%d', experiment_num);
        % load in taskinfo
        C  = dload(fullfile(baseDir, 'sc1_sc2_taskConds_GLM7.txt'));
        Cc = getrow(C, C.StudyNum == experiment_num);
        
        glmSubjDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm), subj_name{sn});
                
        load(fullfile(glmSubjDir, 'SPM.mat'));
        T = load(fullfile(glmSubjDir, 'SPM_info.mat'));
        
        X = SPM.xX.X(:, 1:end-16); % discarding the intercepts
        
        T_run = getrow(T, T.run == 2);
        indr  = T.run == r;
        
        a = (r - 1)* 598 + 1;
        b = r*598;
        X_run = X(a:b, indr); % get the regressors for run 1
        
        figure; % plotting all the regressors in one plot
        % look at the printed statements on the matlab's command window
        for ic = 1:29
            fprintf('plotting the reg for %s\n', T_run.TN{ic})
            plot(X_run(:, ic)); hold on            
        end
        figure; % adding the regressors for each task and plotting them
        % look at the printing statements on the command window
        tasks = unique(T_run.task, 'stable');
        for it = 1:length(tasks)
            name = Cc.taskNames(T_run.task == tasks(it));
            fprintf('task is %s\n', name{1})
            indt = T_run.task == tasks(it);
            A = X_run(:, indt);
            a = sum(A, 2);
            plot(a); hold on;
            keyboard;
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