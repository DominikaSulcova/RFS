%% RFS vs. laser - master script
% author:   Dominika Sulcova
%           MSH - Medical School Hamburg
% created:  2024   
% 
% This is a master script for a study comparing psychophisics and
% neurphysiological response to painful stimulation using a gold standard
% YAP-laser stimulation and a novel radiofrequency stimulator (RFS)
% develpped by the engineering faculty of the University of Oro Verde (AR).
% 
% Collected data: 
%   1) evoked potentials 
%           --> 2 stimulators: laser (LEP), RFS (RFEP)
%           --> 2 intensities: threshold (low), 1.5 * threshold (high)
%           --> 2 areas: left hand (left), right hand (right)
%   2) psychophysical evaluation of painful stimulation
%           --> threshold values (determined by pain-related descriptors)
%           --> intensity ratings acquired at single trial level 
%   3) resting state EEG (RS) 
%           --> 1.5 mins eyes open + 1.5 mins eyes closed
%           --> 2 repetitions: before / after stimulation
% 
% The first part of the script runs during the experimental session and
% allows to:    1) encode subject data
%               2) perform threshold tracking 
%               3) calculate the pain threshold based on psychophysical
%               descriptors --> imported from PsychoPy as .csv
% ==> output:   structure RFSxLASER_info
%               within the output file RFSxLASER_output.mat
% 
% The second part of the script performs EEG pre-processing of single
% subject data with follwing steps:
%   1)
% ==> output:   info, parameters -->    RFSxLASER_info 
%               measured variables -->  RFSxLASER_measures
%               processed data -->  	RFSxLASER_data
%               within the output file RFSxLASER_output.mat
% 
% The third part of the script provides group average visualization and
% exports the data in appropriate formates for the stat analysis in R.

%% ===================== PART 1: experimental session =====================
% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');    % MATLAB toolboxes
folder.raw = uigetdir(pwd, 'Coose the input folder');           % raw data --> at MSH, this should be the study folder at the V drive
folder.output = uigetdir(pwd, 'Choose the OneDrive folder');    % output folder --> One Drive: figures, loutput file, exports 
cd(folder.output)

% output
study = 'RFSxLASER';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);
figure_counter = 1;

% sound
load handel.mat
haleluja = y; clear y Fs

%% 1) fill in subject & session info 
% ----- section input -----
params.laser.intensity.low = 1.25;
params.laser.intensity.high = 1.75;
params.RFS.intensity.low = 15;
params.RFS.intensity.high = 35;
params.laser.pulse = 3;
params.laser.diameter = 5;
params.RFS.pulse = 10;
params.RFS.probe = 5;
% -------------------------

% load the info structure
if exist(output_file) == 2
    output_vars = who('-file', sprintf('%s', output_file));
    if ismember('RFSxLASER_info', output_vars)
        load(output_file, 'RFSxLASER_info')
    else
        RFSxLASER_info = struct;
        save(output_file, 'RFSxLASER_info','-append')
    end
else
    RFSxLASER_info = struct;
    save(output_file, 'RFSxLASER_info')
end

% get subject & session info
prompt = {'date:', 'subject:', 'age:', 'male:', 'handedness:', 'height (cm):', 'weight (kg):', 'right arm (cm):', 'left arm (cm):'};
dlgtitle = 'Subject session information';
dims = [1 35];
definput = {date, 'S000', '18', '1', '1', '175', '70', '60', '60'};
session_info = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput

% identify subject's index
subject_idx = str2num(session_info{2}(end-1:end));

% fill in the metadata structure
RFSxLASER_info(subject_idx).date = session_info{1};
RFSxLASER_info(subject_idx).ID = session_info{2};
RFSxLASER_info(subject_idx).age = str2num(session_info{3});
RFSxLASER_info(subject_idx).male = str2num(session_info{4});
RFSxLASER_info(subject_idx).handedness = str2num(session_info{5});
RFSxLASER_info(subject_idx).body.height = str2num(session_info{6});
RFSxLASER_info(subject_idx).body.weight = str2num(session_info{7});
RFSxLASER_info(subject_idx).body.arm_right = str2num(session_info{8});
RFSxLASER_info(subject_idx).body.arm_left = str2num(session_info{9});
RFSxLASER_info(subject_idx).stimulation.laser.intensity.low = params.laser.intensity.low;
RFSxLASER_info(subject_idx).stimulation.laser.intensity.high = params.laser.intensity.high;
RFSxLASER_info(subject_idx).stimulation.laser.pulse = params.laser.pulse;
RFSxLASER_info(subject_idx).stimulation.laser.diameter = params.laser.diameter;
RFSxLASER_info(subject_idx).stimulation.RFS.intensity.low = params.RFS.intensity.low;
RFSxLASER_info(subject_idx).stimulation.RFS.intensity.high = params.RFS.intensity.high;
RFSxLASER_info(subject_idx).stimulation.RFS.pulse = params.RFS.pulse;
RFSxLASER_info(subject_idx).stimulation.RFS.probe = params.RFS.probe;
clear session_info

% save to the output file
save(output_file, 'RFSxLASER_info', '-append');
clear params output_vars

%% 2) run threshold tracking 
% ----- section input -----
params.stimulus = {'laser', 'RFS'};
params.side = {'right', 'left'};
% ------------------------- 

% ask for condition order
prompt = {'starting stimulus:', 'starting side:'};
dlgtitle = sprintf('thresholding');
dims = [1 35];
definput = {'laser-RFS', 'right-left'};
thresholding_order = (inputdlg(prompt,dlgtitle,dims,definput));
clear prompt dlgtitle dims definput
if ~strcmp(params.stimulus{1}, thresholding_order{1})
    params.stimulus = params.stimulus([2, 1]);
end
if ~strcmp(params.side{1}, thresholding_order{2})
    params.side = params.side([2, 1]);
end
clear thresholding_order

% run for both stimuli
for s = 1:length(params.stimulus)
    for a = 1:length(params.side)
        % ask for numbers to be plotted big
        plot_fig = true;

        % press ENTER to start timing the session 
        fprintf('Threshold tracking: %s - %s hand\n', params.stimulus{s}, params.side{a});
        ok = 0; 
        while ok == 0    
            key = input('Press ENTER to start', 's');
            if isempty(key)  
                ok = 1;  
            end
        end    
        fprintf('Let''s go!\n');
        clear ok key  
    
        % initiate variables
        intensity = struct;
        trial_counter = 1;
        change_tracker = [];
        flip_counter = 0;
        continue_trials = true;
        
        % threshold tracking - loop through trials
        while continue_trials
            % encode current intensity
            prompt = {'intensity:'};
            dlgtitle = sprintf('trial n. %d', trial_counter);
            dims = [1 35];
            if strcmp(params.stimulus{s}, 'laser') && trial_counter == 1
                definput = {'0.75'};
            elseif strcmp(params.stimulus{s}, 'RFS') && trial_counter == 1
                definput = {'1'};
            else
                definput = {num2str(intensity(trial_counter - 1).intensity)};
            end
            trial.intensity = (inputdlg(prompt,dlgtitle,dims,definput));
            clear prompt dlgtitle dims definput
            
            % determine time
            trial.time = sprintf('%s', datetime);
            
            % append to the output vector
            intensity(trial_counter).trial = trial_counter;
            intensity(trial_counter).time = trial.time(end-7:end);
            intensity_input = true;
            while intensity_input
                if length(trial.intensity) > 0 && ~isempty(trial.intensity{:})
                    intensity(trial_counter).intensity = str2num(trial.intensity{:});
                    intensity_input = false;
                else
                    % ask to continue
                    answer = questdlg('No intensity was provided. Continue thresholding?', 'next step',...
                    'Continue', 'Stop', 'Continue');
                    if strcmp(answer, 'Continue')
                        prompt = {'intensity:'};
                        dlgtitle = sprintf('trial n. %d', trial_counter);
                        dims = [1 35];
                        if strcmp(params.stimulus{s}, 'laser') && trial_counter == 1
                            definput = {'0.75'};
                        elseif strcmp(params.stimulus{s}, 'RFS') && trial_counter == 1
                            definput = {'1'};
                        else
                            definput = {num2str(intensity(trial_counter - 1).intensity)};
                        end
                        trial.intensity = (inputdlg(prompt,dlgtitle,dims,definput));
                        clear prompt dlgtitle dims definput
                    else
                        % remove the last row from structure 'intensity'
                        intensity(end) = [];
                        trial_counter = trial_counter - 1;

                        % leave the loop
                        plot_fig = false;
                        intensity_input = false;                        
                    end
                end
            end
            clear answer
               
            % open the figure
            if plot_fig
                screen_size = get(0, 'ScreenSize');
                fig_intensity = figure('Position', [1, screen_size(4)/2, screen_size(3) / 2, screen_size(4) / 2], ...
                              'MenuBar', 'none', ...
                              'NumberTitle', 'off', ...
                              'Color', [1 1 1]); 
                text(0.5, 0.5, num2str(intensity(trial_counter).intensity), ...
                    'Color', 'r', ...  
                    'FontSize', 150, ...  
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle');                
                axis off;  
                
                % check if the figure should close
                key = input('Press ENTER to continue', 's');
                if isempty(key)
                    close(fig_intensity);  
                end
                clear key
    
                % check for flips
                if trial_counter > 1
                    change_tracker(end + 1) = intensity(trial_counter).intensity - intensity(trial_counter-1).intensity;
                end
                if trial_counter > 2
                    if change_tracker(end - 1) > 0 && change_tracker(end) < 0
                        flip_counter = flip_counter + 1;
                    elseif change_tracker(end - 1) < 0 && change_tracker(end) > 0
                        flip_counter = flip_counter + 1;
                    end
                end
                if flip_counter == 5
                    fprintf('Painful sensation was reached 3 times: as soon as the subject stops rating pain, thresholding is finished.\n')
                end
            end
          
            % ask to continue
            answer = questdlg('Continue with this thresholding block?', 'next step',...
                'Continue', 'Finish', 'Continue');
            if strcmp(answer, 'Continue')
                % update trial counter
                trial_counter = trial_counter + 1;
            else
                % leave the loop
                continue_trials = false;
            end
        end
        clear answer
            
        % append to info structure
        statement = sprintf('RFSxLASER_info(subject_idx).stimulation.%s.threshold.%s.measures = intensity;', params.stimulus{s}, params.side{a});
        eval(statement)
        if strcmp(params.stimulus{s}, 'laser')
            RFSxLASER_info(subject_idx).stimulation.laser.threshold.unit = 'J';
        elseif strcmp(params.stimulus{s}, 'RFS')
            RFSxLASER_info(subject_idx).stimulation.RFS.threshold.unit = '%MSO';
        end
    end
end

% save to the output file
save(output_file, 'RFSxLASER_info', '-append');
clear params s a intensity trial_counter continue_trials trial statement open screen_size ...
    intensity_input fig_intensity plot_fig change_tracker flip_counter

%% 3) determine thresholds 
% ----- section input -----
params.stimulus = {'laser', 'RFS'};
params.side = {'right', 'left'};
params.labels = {'nothing' 'warm' 'touch' 'pricking' 'burning' 'electric'};
params.values = [0 1 1 2 2 3];
params.folder = 'PsychoPy_threshold';
% -------------------------

% ask for subject number, if not defined
if ~exist('subject_idx')
    prompt = {'subject number:'};
    dlgtitle = 'subject';
    dims = [1 40];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    subject_idx = str2num(input{1,1});
end
clear prompt dlgtitle dims definput input

% load the info structure, if not loaded
if ~exist('RFSxLASER_info')
    load(output_file, 'RFSxLASER_info')
end 

% append the descriptors
RFSxLASER_info(subject_idx).stimulation.descriptors.labels = params.labels;
RFSxLASER_info(subject_idx).stimulation.descriptors.values = params.values;

% run for both stimuli and both hands
for s = 1:length(params.stimulus)
    for a = 1:length(params.side)
        % provide update
        fprintf('%s pain threshold - %s hand:\n', params.stimulus{s}, params.side{a})

        % determine abbreviation
        if strcmp(params.side{a}, 'right')
            params.area = 'RH';
        elseif strcmp(params.side{a}, 'left')
            params.area = 'LH';
        end
    
        % import PsychoPy .csv file
        fprintf('importing PsychoPy data...\n')
        file2import = dir(sprintf('%s\\%s\\%s\\*_%s*%s*.csv', folder.raw, RFSxLASER_info(subject_idx).ID, params.folder, params.stimulus{s}, params.area));
        % file2import = dir(sprintf('%s\\*%s*_%s*%s*.csv', folder.raw, RFSxLASER_info(subject_idx).ID, params.stimulus{s}, params.area));
        if size(file2import, 1) ~= 1
            error('ERROR: %d files attributed %s thresholding on the %s hand were found!\n', size(file2import, 1), params.stimulus{s}, params.side{a})
        end
        option = detectImportOptions(sprintf('%s\\%s', file2import.folder, file2import.name));
        option.SelectedVariableNames = {'trial', 'fixation_stopped', 'descriptor', 'pain'};
        ratings_table = readtable(sprintf('%s\\%s', file2import.folder, file2import.name), option);
        clear file2import option

        % get rid of NaNs
        table_idx = logical([]);
        for b = 1:height(ratings_table)
            if isnan(ratings_table.fixation_stopped(b))
                table_idx(b) = true;
            else
                table_idx(b) = false;
            end
        end
        ratings_table(table_idx, :) = [];
    
        % adjust ratings format
        fprintf('extracting trials...\n')
        ratings = struct([]);
        for b = 1:height(ratings_table)
            % if it is a valid trial
            if ~isnan(ratings_table{b, 1})
                % encode trial
                ratings(b).trial = ratings_table{b, 1};
                
                % encode timestamp
                ratings(b).time_psychopy = ratings_table{b, 2};
                
                % encode descriptors
                desc = ratings_table{b, 3};
                desc = regexprep(desc{:}, '[[]'',]' , '');
                ratings(b).descriptors = split(desc, ' ');
                
                % encode pain intensity
                int = ratings_table{b, 4};
                int = regexprep(int{:}, '[[]'',]' , '');
                ratings(b).pain = split(int, ' ');
            end
        end
    
        % check if the trial number matches
        statement = sprintf('n_trials = length(RFSxLASER_info(subject_idx).stimulation.%s.threshold.%s.measures);', params.stimulus{s}, params.side{a});
        eval(statement)
        if n_trials ~= length(ratings)
            error('ERROR: PsychoPy trials do not match MATLAB trials!\n')
        end
    
        % clasify trials and append to output structure
        fprintf('%d valid trials found\nclassifying ...\n', length(ratings)) 
        visual = struct;
        for c = 1:length(ratings)            
            % determine level of evoked sensation based on descriptors
            for d = 1:length(ratings(c).descriptors)
                desc_int = [];
                index = find(strcmp(ratings(c).descriptors{d}, params.labels));
                if ~isempty(index)
                    desc_int(d) = params.values(index);
                else
                    error('ERROR: Trial %d contains an unknown descriptor.\n', d);
                end                
            end
            ratings(c).sensation = max(desc_int);
            
            % prepare visualization input  
            statement = sprintf('visual.unit = RFSxLASER_info(subject_idx).stimulation.%s.threshold.unit;', params.stimulus{s});
            eval(statement)
            statement = sprintf('visual.intensity(c) = RFSxLASER_info(subject_idx).stimulation.%s.threshold.%s.measures(c).intensity;', ...
                params.stimulus{s}, params.side{a});
            eval(statement)
            if ratings(c).sensation == 0
                visual.color(c, :) = [0.6510    0.6510    0.6510];
            elseif ratings(c).sensation == 1
                visual.color(c, :) = [0    0.4471    0.7412];
            elseif ratings(c).sensation == 2
                visual.color(c, :) = [0.8000    0.0157    0.0157];
            elseif ratings(c).sensation == 3
                visual.color(c, :) = [0 0 0];
            end  
        end
    
        % identify transitions in sensation = thresholds
        fprintf('identifying thresholds...\n') 
        threshold.perception = [];
        threshold.pain = [];
        for e = 2:length(ratings)
            if ratings(e).sensation == 1 & ratings(e-1).sensation == 0
                threshold.perception(end+1) = mean(visual.intensity([e-1, e]));
            elseif ratings(e).sensation == 0 & ratings(e-1).sensation == 1
                threshold.perception(end+1) = mean(visual.intensity([e-1, e]));
            elseif ratings(e).sensation == 2 & ratings(e-1).sensation == 1
                threshold.pain(end+1) = mean(visual.intensity([e-1, e]));
            elseif ratings(e).sensation == 1 & ratings(e-1).sensation == 2
                threshold.pain(end+1) = mean(visual.intensity([e-1, e]));
            end
        end
    
        % identify transitions in pain rating = subjective pain threshold
        threshold.subjective = [];
        for f = 2:length(ratings)
            if length(ratings(f).pain) == 1
                if strcmp(ratings(f).pain{1}, 'yes') & strcmp(ratings(f-1).pain{1}, 'no')
                    threshold.subjective(end+1) = mean(visual.intensity([f-1, f]));
                elseif strcmp(ratings(f).pain{1}, 'no') & strcmp(ratings(f-1).pain{1}, 'yes')
                    threshold.subjective(end+1) = mean(visual.intensity([f-1, f]));
                end
            else
                error('ERROR: Trial %d contains invalid pain rating.\n', f);
            end
        end

        % encode thresholds to the output structure 
        statement = sprintf('RFSxLASER_info(subject_idx).stimulation.%s.threshold.%s.perception = round(mean(threshold.perception, 2));', ...
            params.stimulus{s}, params.side{a});
        eval(statement)
        statement = sprintf('RFSxLASER_info(subject_idx).stimulation.%s.threshold.%s.pain = round(mean(threshold.pain), 2);', ...
            params.stimulus{s}, params.side{a});
        eval(statement)
        statement = sprintf('RFSxLASER_info(subject_idx).stimulation.%s.threshold.%s.subjective = round(mean(threshold.subjective), 2);', ...
            params.stimulus{s}, params.side{a});
        eval(statement)
    
        % launch the figure
        visual.x = 1:length(ratings);
        fig = figure(figure_counter);
        set(fig, 'units','normalized','outerposition',[0 0 1 1])
        plot_thresholds(visual, threshold)
   
        % adjust labels
        if strcmp(params.stimulus{a}, 'laser')
            ylabel('laser energy (J)')
        elseif strcmp(params.stimulus{a}, 'RFS')
            ylabel('RFS intensity (%MSO)')
        end
        title(sprintf('%s: %s threshold tracking - %s hand', RFSxLASER_info(subject_idx).ID, params.stimulus{s}, params.side{a}))
           
        % save figure and update counter
        saveas(fig, sprintf('%s\\figures\\%s_threshold_%s_%s.png', folder.output, RFSxLASER_info(subject_idx).ID, params.stimulus{s}, params.side{a}))
        figure_counter = figure_counter + 1;
    
        % append ratings to the info structure
        statement = sprintf('measures_all = RFSxLASER_info(subject_idx).stimulation.%s.threshold.%s.measures;', params.stimulus{s}, params.side{a});
        eval(statement)
        fields = fieldnames(ratings);  
        for g = 1:length(ratings)            
            for h = 2:length(fields)
                measures_all(g).(fields{h}) = ratings(g).(fields{h});  
            end
        end
        statement = sprintf('RFSxLASER_info(subject_idx).stimulation.%s.threshold.%s.measures = measures_all;', params.stimulus{s}, params.side{a});
        eval(statement)

        % save to the output file
        save(output_file, 'RFSxLASER_info', '-append');
    end
end

% ask if the subject is done
answer = questdlg(sprintf('Do you want to continue to subject %d?', subject_idx + 1), 'Continue to the next subject?', 'YES', 'NO', 'NO'); 
switch answer
    case 'NO'
    case 'YES'
    	subject_idx = subject_idx + 1;
end 
clear params a b c d e f g h s file2import ratings_table desc int index desc_int threshold measures_all fields visual fig n_trials statement ratings 

%% ===================== PART 2: single subject data processing ================
% directories
if ~exist('folder') 
    folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');    % MATLAB toolboxes
    folder.raw = uigetdir(pwd, 'Coose the input folder');           % raw data --> at MSH, this should be the study folder at the V drive
    folder.output = uigetdir(pwd, 'Choose the OneDrive folder');    % output folder --> figures, loutput file, exports 
end
folder.processed = uigetdir(pwd, 'Choose the data folder');         % processed data --> wherever you want to store the voluminous EEG data
cd(folder.output)

% output
study = 'RFSxLASER';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);
if ~exist('figure_counter')  
    figure_counter = 1;
end

% current subject 
prompt = {'subject number:'};
dlgtitle = 'subject';
dims = [1 40];
definput = {''};
input = inputdlg(prompt,dlgtitle,dims,definput);
subject_idx = str2num(input{1,1});
clear prompt dlgtitle dims definput input

%% 1) import data for letswave, preview
% ----- section input -----
params.data = {'RS', 'LEP', 'RFS'};
params.preview = true;
params.folder = 'EEG';
params.data_n = 4;
params.event_n = 30;
params.downsample = 5;
params.eventcode = {'L  1', 'R  1'};
params.epoch = [-0.3 1];
params.suffix = 'preview';
params.eoi = 'Cz';
params.ylim = [-15 15];
% -------------------------
fprintf('section 1: data import\n')

% update info structure
load(output_file, 'RFSxLASER_info');

% add letswave 6 to the top of search path
addpath(genpath([folder.toolbox '\letswave 6']));
   
% cycle though datasets and import in letswave format
fprintf('loading:\n')
for a = 1:length(params.data)
    % identify the appropriate files
    file2import = dir(sprintf('%s\\%s\\%s\\*%s*%s*.vhdr', folder.raw, RFSxLASER_info(subject_idx).ID, params.folder, study, params.data{a}));

    % remove average files if necessary
    file2rmv = logical([]);
    for b = 1:length(file2import)
        if contains(file2import(b).name, 'avg') 
            file2rmv(b) = true;
        else
            file2rmv(b) = false;
        end
    end
    file2import(file2rmv) = [];
    
    % check that the number of files matches
    if size(file2import, 1) ~= params.data_n
        fprintf('WARNING: incorrect number (%d) of %s recordings was found!\n', size(file2import, 1), params.data{a})
    end

    % cycle through datasets
    for c = 1:length(file2import)
        % identify the filename
        filename = sprintf('%s\\%s', file2import(c).folder, file2import(c).name);
        
        % prepare the dataset name
        dataname = replace(file2import(c).name, ' ', '');
        underscores = strfind(dataname, '_');
        dataname = dataname(underscores(1) + 1 : underscores(end) - 1);
        dataname = replace(dataname, '_', ' ');
        if contains(dataname, 'LEP')
            dataname = replace(dataname, 'LEP', 'laser');
        end

        % provide update
        fprintf('%s ...\n', dataname)
        
        % encode the filename to metadata
        block = regexp(file2import(c).name, 'b(\d+)', 'tokens');
        block = str2num(block{1}{1});
        RFSxLASER_info(subject_idx).dataset(block - 1).block = block;
        RFSxLASER_info(subject_idx).dataset(block - 1).name = file2import(c).name;
    
        % import the dataset
        [dataset(subject_idx).raw(block - 1).header, dataset(subject_idx).raw(block - 1).data, ~] = RLW_import_VHDR(filename);
    
        % rename in the header
        dataset(subject_idx).raw(block - 1).header.name = dataname;
    end
end
fprintf('done.\n')

% save info to the output file
save(output_file, 'RFSxLASER_info', '-append');

% ask if the raw data should be saved in letswave format
answer = questdlg(sprintf('Should I save the raw data of subject %d in the letswave format?', subject_idx), 'Save the data?', 'YES', 'NO', 'NO'); 
switch answer
    case 'NO'
    case 'YES'
        % save raw data to the output folder
	    for d = 1:length(dataset(subject_idx).raw)
            data = dataset(subject_idx).raw(d).data;
            header = dataset(subject_idx).raw(d).header;
            save(sprintf('%s.mat', header.name), 'data');
            save(sprintf('%s.lw6', header.name), 'header');
        end
end 

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));
   
if params.preview
    % ask if the data should be roughly pre-processed and saved for preview
    answer = questdlg(sprintf('Do you want to save pre-processed data for a preview in letswave?'), 'Save preview data?', 'YES', 'NO', 'YES'); 
    switch answer
        case 'NO'
            preview_save = 0;
        case 'YES'
    	    preview_save = 1;
    end 
    
    % pre-process ERPs for preview and save for letswave if required
    fprintf('Pre-processing ERPs for preview:\n')
    for d = 1:length(dataset(subject_idx).raw)
        if ~isempty(strfind(dataset(subject_idx).raw(d).header.name, 'RS'))
        else 
            % provide update
            fprintf('%s ...\n', dataset(subject_idx).raw(d).header.name)
    
            % verify the number and type of triggers
            events = {dataset(subject_idx).raw(d).header.events.code};
            event_n = 0;
            for e = 1:length(events)
                if strcmp(events{e}, params.eventcode{1}) || strcmp(events{e}, params.eventcode{2})
                    event_n = event_n + 1;
                end
            end
            if event_n ~= params.event_n
                fprintf('ATTENTION: incorrect number of events (%d) were found!\n', event_n)
            end
    
            % select the data
            lwdata.header = dataset(subject_idx).raw(d).header;
            lwdata.data = dataset(subject_idx).raw(d).data;
    
            % downsample
            option = struct('x_dsratio', params.downsample, 'suffix', '', 'is_save',0);
            lwdata = FLW_downsample.get_lwdata(lwdata, option);
    
            % segment
            if contains(dataset(subject_idx).raw(d).header.name, 'laser')
                event_code = params.eventcode{1};
            elseif contains(dataset(subject_idx).raw(d).header.name, 'RFS')
                event_code = params.eventcode{2};
            end
            option = struct('event_labels', {event_code}, 'x_start', params.epoch(1), 'x_end', params.epoch(2), ...
                'x_duration', params.epoch(2)-params.epoch(1), 'suffix', '', 'is_save', 0);
            lwdata = FLW_segmentation.get_lwdata(lwdata, option);
        
            % remove DC + linear detrend
            option = struct('linear_detrend', 1, 'suffix', params.suffix, 'is_save', preview_save);
            lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
    
            % append to the dataset
            dataset(subject_idx).preview(d).header = lwdata.header;
            dataset(subject_idx).preview(d).data = lwdata.data;
        end    
    end
    dataset(subject_idx).preview([1, 2]) = [];       % remove empty rows
    
    % add fieldtrip to the top of search path
    addpath(genpath([folder.toolbox '\fieldtrip']));
    
    % determine electrode to plot
    eoi = find(strcmp({dataset(subject_idx).preview(1).header.chanlocs.labels}, params.eoi));
    
    % plot preview
    fig = figure(figure_counter);
    screen_size = get(0, 'ScreenSize');
    set(fig, 'Position', screen_size);
    figure_counter = figure_counter + 1;
    for f = 1:length(dataset(subject_idx).preview)
        % select data
        data_visual = squeeze(dataset(subject_idx).preview(f).data);
    
        % prepare fieldtrip input 
        data = [];
        data.label = {dataset(subject_idx).preview(f).header.chanlocs.labels};
        data.fsample = 1/dataset(subject_idx).preview(f).header.xstep;
        data.time = {params.epoch(1) : dataset(subject_idx).preview(f).header.xstep : params.epoch(2) - dataset(subject_idx).preview(f).header.xstep};  
        data.trial = {squeeze(mean(data_visual, 1))}; 
    
        % define the layout 
        cfg = [];
        cfg.layout = 'acticap-64ch-standard2';
        cfg.viewmode = 'vertical'; 
        cfg.showlabels = 'yes';  
        cfg.channel = 'all';  
        cfg.ylim = params.ylim;
    
        % plot all channels
        figure(figure_counter)
        ft_multiplotER(cfg, data);
        set(gcf, 'Position', screen_size);
        set(gcf, 'Name', dataset(subject_idx).preview(f).header.name);
    
        % determine subplot
        if contains(dataset(subject_idx).preview(f).header.name, 'laser')
            a = 0;
            fig_name = 'laser';
        else
            a = 4;
            fig_name = 'RFS';
        end
        if contains(dataset(subject_idx).preview(f).header.name, 'left')
            b = 0;
            fig_name = [fig_name ' left'];
        else
            b = 2;
            fig_name = [fig_name ' right'];
        end
        if contains(dataset(subject_idx).preview(f).header.name, 'high')
            c = 1;
            fig_name = [fig_name ' high'];
        else
            c = 2;
            fig_name = [fig_name ' low'];
        end
    
        % plot 
        figure(fig)
        subplot(2, 4, a + b + c)
        for d = 1:size(data_visual, 1)
            plot(data.time{1}, squeeze(data_visual(d, eoi, :)), "Color", [0.6471    0.8784    0.9804], "LineWidth", 1.1)
            hold on
        end
        plot(data.time{1}, squeeze(mean(data_visual(:, eoi, :), 1)), "Color", [0.0941    0.5059    0.7804], "LineWidth", 3)
        line([0, 0], params.ylim, 'Color', 'black', 'LineWidth', 2, 'LineStyle', '--')   
        title(fig_name)
        ylim(params.ylim)
        xlim([data.time{1}(1) data.time{1}(end)])
        set(gca, 'fontsize', 14)
        ylabel('amplitude (\muV)')
        xlabel('time (s)')
        set(gca, 'Layer', 'Top')
        set(gca, 'YDir', 'reverse');
    
        % update figure counter 
        figure_counter = figure_counter + 1;
    end
    
    % adjust and save output figure
    figure(fig)
    sgtitle(sprintf('%s: ERP preview - %s electrode', RFSxLASER_info(subject_idx).ID, params.eoi))
    saveas(fig, sprintf('%s\\figures\\%s_preview.png', folder.output, RFSxLASER_info(subject_idx).ID))
end

% clear and continue
fprintf('\nsection 1 finished.\n\n')
clear params a b c d e f block file2import file2rmv filename dataname underscores events event_n event_code lwdata ...
    data data_visual cfg eoi fig fig_name option screen_size answer preview_save

%% 2) pre-process all data
% ----- section input -----
params.suffix = {'dc' 'bandpass' 'notch' 'ds' 'reref' 'ep' 'dc'};
params.eventcode = {'L  1', 'R  1'};
params.eventcode_new = {'laser', 'RFS'};
params.interpolate = [-0.002 0.02];
params.shift = 0.002;
params.bandpass = [0.1 80];
params.downsample = 2;
params.interp_chans = 6;
params.event_n = 30;
params.epoch = [-0.3 1];
% -------------------------
fprintf('section 2: data pre-processing\n')

% update info structure
load(output_file, 'RFSxLASER_info');
cd(folder.processed)

% add letswave 6 to the top of search path
addpath(genpath([folder.toolbox '\letswave 6']));

% launch check figure
fig = figure(figure_counter);
set(fig, 'Position', [50, 30, 750, 750])
visual = struct;
visual.x = -0.01 : dataset(subject_idx).raw(1).header.xstep : 0.03;

% interpolate RFS artifact and shift
fprintf('*********** subject %d: interpolating RF artifact ***********\n', subject_idx)
d_rfs = 1;      % dataset counter
for d = 1:length(dataset(subject_idx).raw)
    if contains(dataset(subject_idx).raw(d).header.name, 'RFS')  
        % provide update
        fprintf('%s ...', dataset(subject_idx).raw(d).header.name(10:end))

        % interpolate signal around the RF artifact --> cubic interpolation
        [dataset(subject_idx).interpolated(d_rfs).header, dataset(subject_idx).interpolated(d_rfs).data, ~] = RLW_suppress_artifact_event(dataset(subject_idx).raw(d).header, dataset(subject_idx).raw(d).data,... 
            'xstart', params.interpolate(1), 'xend', params.interpolate(2), 'event_code', params.eventcode{2}, 'interp_method', 'pchip');

        % shift data by trigger duration
        shift_samples = round(params.shift / dataset(subject_idx).interpolated(d_rfs).header.xstep); 
        for c = 1:size(dataset(subject_idx).interpolated(d_rfs).data, 2)
            dataset(subject_idx).interpolated(d_rfs).data(1, c, 1, 1, 1, 1:end-shift_samples) = dataset(subject_idx).interpolated(d_rfs).data(1, c, 1, 1, 1, 1 + shift_samples:end);
        end

        % extract data for visualization --> Cz electrode, first trigger
        trigger = round(dataset(subject_idx).interpolated(d_rfs).header.events(2).latency / dataset(subject_idx).interpolated(d_rfs).header.xstep);
        visual.y(1, :) = dataset(subject_idx).raw(d).data(1, 23, 1, 1, 1, trigger - 25 : trigger + 75);                   % raw data
        visual.y(2, :) = dataset(subject_idx).interpolated(d_rfs).data(1, 23, 1, 1, 1, trigger - 25 : trigger + 75);      % iterpolated data

        % plot into the check figure
        figure(fig)
        subplot(2, 2, d_rfs)
        for a = 1:size(visual.y, 1)
            plot(visual.x, visual.y(a, :),  "LineWidth", 1.2)
            hold on
        end
        visual.ylim = get(gca, 'YLim');
        line([0, 0], visual.ylim, 'Color', 'black', 'LineWidth', 2, 'LineStyle', '--')   
        title(dataset(subject_idx).interpolated(d_rfs).header.name(6:end))        
        set(gca, 'fontsize', 10)
        ylabel('amplitude (\muV)')
        xlabel('time (s)')
        set(gca, 'Layer', 'Top')
        set(gca, 'YDir', 'reverse');
        if d_rfs == 1
            sgtitle(sprintf('subject %d RF artifact interpolation', subject_idx))
        end
    
        % update info structure
        if d_rfs == 1
            RFSxLASER_info(subject_idx).preprocessing(1).process = sprintf('RF artifact interpolated'); 
            RFSxLASER_info(subject_idx).preprocessing(1).params.method = 'pchip';
            RFSxLASER_info(subject_idx).preprocessing(1).params.limits = params.interpolate;
            RFSxLASER_info(subject_idx).preprocessing(1).suffix = [];
            RFSxLASER_info(subject_idx).preprocessing(1).date = sprintf('%s', date);
            RFSxLASER_info(subject_idx).preprocessing(2).process = sprintf('RF data shifted to correct for trigger delay');
            RFSxLASER_info(subject_idx).preprocessing(2).params.shift = params.shift;
            RFSxLASER_info(subject_idx).preprocessing(2).suffix = [];
            RFSxLASER_info(subject_idx).preprocessing(2).date = sprintf('%s', date);
        end
        
        % update counter
        d_rfs = d_rfs + 1;
    end
end
fprintf('done.\n')
fprintf('\n')

% save and update figure counter 
figure(fig)
saveas(fig, sprintf('%s\\figures\\%s_RFartifact.png', folder.output, RFSxLASER_info(subject_idx).ID))
figure_counter = figure_counter + 1;

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% first step of pre-processing
fprintf('*********** subject %d: first pre-processing ***********\n', subject_idx)
d_rfs = 1;      % RFS counter
for d = 1:length(dataset(subject_idx).raw)
    % provide update
    fprintf('%s:\n', dataset(subject_idx).raw(d).header.name(6:end))

    % select the data for pre-processing
    if contains(dataset(subject_idx).raw(d).header.name, 'RFS')
        % load the interpolated data
        lwdata.header = dataset(subject_idx).interpolated(d_rfs).header;
        lwdata.data = dataset(subject_idx).interpolated(d_rfs).data;

        % update the counter
        d_rfs = d_rfs + 1;
    else
        lwdata.header = dataset(subject_idx).raw(d).header;
        lwdata.data = dataset(subject_idx).raw(d).data;
    end

    % assign electrode coordinates
    fprintf('assigning electrode coordinates...')
    option = struct('filepath', sprintf('%s\\letswave 7\\res\\electrodes\\spherical_locations\\Standard-10-20-Cap81.locs', folder.toolbox), ...
        'suffix', '', 'is_save', 0);
    lwdata = FLW_electrode_location_assign.get_lwdata(lwdata, option);
    if d == 1
        RFSxLASER_info(subject_idx).preprocessing(3).process = sprintf('electrode coordinates assigned');
        RFSxLASER_info(subject_idx).preprocessing(3).params.layout = sprintf('standard 10-20-cap81');
        RFSxLASER_info(subject_idx).preprocessing(3).suffix = [];
        RFSxLASER_info(subject_idx).preprocessing(3).date = sprintf('%s', date);
    end
    
    % remove DC + linear detrend
    fprintf('removing DC and applying linear detrend...')
    option = struct('linear_detrend', 1, 'suffix', params.suffix{1}, 'is_save', 0);
    lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
    if d == 1
        RFSxLASER_info(subject_idx).preprocessing(4).process = sprintf('DC + linear detrend on all continuous data');
        RFSxLASER_info(subject_idx).preprocessing(4).suffix = params.suffix{1};
        RFSxLASER_info(subject_idx).preprocessing(4).date = sprintf('%s', date);
    end
    
    % bandpass
    fprintf('applying Butterworth bandpass filter...')
    option = struct('filter_type', 'bandpass', 'high_cutoff', params.bandpass(2),'low_cutoff', params.bandpass(1),...
        'filter_order', 4, 'suffix', params.suffix{2}, 'is_save', 0);
    lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);
    if d == 1
        RFSxLASER_info(subject_idx).preprocessing(5).process = sprintf('bandpass filtered');
        RFSxLASER_info(subject_idx).preprocessing(5).params.method = sprintf('Butterworth');
        RFSxLASER_info(subject_idx).preprocessing(5).params.order = 4;
        RFSxLASER_info(subject_idx).preprocessing(5).params.limits = params.bandpass;
        RFSxLASER_info(subject_idx).preprocessing(5).suffix = params.suffix{2};
        RFSxLASER_info(subject_idx).preprocessing(5).date = sprintf('%s', date);
    end
    
    % 50 Hz notch
    fprintf('applying FFT notch filter...')
    option = struct('filter_type', 'notch', 'notch_fre', 50, 'notch_width', 2, 'slope_width', 2,...
        'harmonic_num', 2,'suffix', params.suffix{3},'is_save', 0);
    lwdata = FLW_FFT_filter.get_lwdata(lwdata, option);
    if d == 1
        RFSxLASER_info(subject_idx).preprocessing(6).process = sprintf('notch filtered at 50 Hz');
        RFSxLASER_info(subject_idx).preprocessing(6).params.method = sprintf('FFT');
        RFSxLASER_info(subject_idx).preprocessing(6).params.width = 2;
        RFSxLASER_info(subject_idx).preprocessing(6).params.slope = 2;
        RFSxLASER_info(subject_idx).preprocessing(6).suffix = params.suffix{3};
        RFSxLASER_info(subject_idx).preprocessing(6).date = sprintf('%s', date);
    end
    
    % downsample and save
    fprintf('downsampling...\n')
    option = struct('x_dsratio', params.downsample, 'suffix', params.suffix{4}, 'is_save',1);
    lwdata = FLW_downsample.get_lwdata(lwdata, option);
    if d == 1
        RFSxLASER_info(subject_idx).preprocessing(7).process = sprintf('downsampled');
        RFSxLASER_info(subject_idx).preprocessing(7).params.ratio = params.downsample;
        RFSxLASER_info(subject_idx).preprocessing(7).params.fs_orig = 1/lwdata.header.xstep * params.downsample;
        RFSxLASER_info(subject_idx).preprocessing(7).params.fs_final = 1/lwdata.header.xstep;
        RFSxLASER_info(subject_idx).preprocessing(7).suffix = params.suffix{4};
        RFSxLASER_info(subject_idx).preprocessing(7).date = sprintf('%s', date);
    end

    % update dataset
    dataset(subject_idx).preprocessed(d).header = lwdata.header;
    dataset(subject_idx).preprocessed(d).data = lwdata.data;   
end
fprintf('done.\n')
fprintf('\n')

% check in letswave for bad channels
letswave
ok = 0; 
while ok == 0    
    key = input('Press ENTER to start', 's');
    if isempty(key)  
        ok = 1;  
    end
end    
clear ok key 

% interpolate channels if needed
prompt = {'interpolate channels:'};
dlgtitle = 'channel interpolation';
dims = [1 220];
definput = {strjoin({lwdata.header.chanlocs.labels}, ' ')};
chans2interpolate = inputdlg(prompt,dlgtitle,dims,definput);
chans2interpolate = split(chans2interpolate{1}, ' ');
clear prompt dlgtitle dims definput 
if ~isempty(chans2interpolate{1})
    % provide update
    fprintf('interpolating channel ')

    % loop through channels to interpolate
    for c = 1:length(chans2interpolate)
        fprintf('%s ...', chans2interpolate{c})

        % indentify the channel to interpolate
        chan_n = find(strcmp({lwdata.header.chanlocs.labels}, chans2interpolate{c}));

        % calculate distances with other electrodes
        chan_dist = -ones(length(lwdata.header.chanlocs), 1);
        for b = setdiff(1:length(lwdata.header.chanlocs), chan_n)
            if lwdata.header.chanlocs(b).topo_enabled == 1
                chan_dist(b) = sqrt((lwdata.header.chanlocs(b).X - lwdata.header.chanlocs(chan_n).X)^2 + ...
                    (lwdata.header.chanlocs(b).Y - lwdata.header.chanlocs(chan_n).Y)^2 + ...
                    (lwdata.header.chanlocs(b).Z - lwdata.header.chanlocs(chan_n).Z)^2);
            end
        end
        chan_dist((chan_dist==-1)) = max(chan_dist);
        [~,chan_dist] = sort(chan_dist);

        % identify neighbouring channels
        chan_dist = chan_dist(1:params.interp_chans);
        chans2use = {lwdata.header.chanlocs.labels};
        chans2use = chans2use(chan_dist);

        % cycle through all datasets
        for d = 1:length(dataset(subject_idx).raw)
            % select data
            lwdata.header = dataset(subject_idx).preprocessed(d).header;
            lwdata.data = dataset(subject_idx).preprocessed(d).data;

            % interpolate using the neighboring electrodes
            option = struct('channel_to_interpolate', chans2interpolate{c}, 'channels_for_interpolation_list', {chans2use}, ...
                'suffix', 'chan_interp', 'is_save', 0);
            lwdata = FLW_interpolate_channel.get_lwdata(lwdata, option);

            % update dataset
            dataset(subject_idx).preprocessed(d).header = lwdata.header;
            dataset(subject_idx).preprocessed(d).data = lwdata.data;  
        end

        % update info structure
        if c == 1
            RFSxLASER_info(subject_idx).preprocessing(8).process = sprintf('visual check - bad channels interpolated');
            RFSxLASER_info(subject_idx).preprocessing(8).suffix = [];
            RFSxLASER_info(subject_idx).preprocessing(8).date = sprintf('%s', date);
        end
        RFSxLASER_info(subject_idx).preprocessing(8).params.bad{c} = chans2interpolate{c};
        RFSxLASER_info(subject_idx).preprocessing(8).params.chans_used{c} = strjoin(chans2use, ' ');    
    end
else
    RFSxLASER_info(subject_idx).preprocessing(8).process = sprintf('visual check');
    RFSxLASER_info(subject_idx).preprocessing(8).suffix = [];
    RFSxLASER_info(subject_idx).preprocessing(8).date = sprintf('%s', date);
end

% segment ERP data
fprintf('*********** subject %d: segmenting ERP data ***********\n', subject_idx)
bad_counter = 1;
for d = 1:length(dataset(subject_idx).preprocessed)
    if ~contains(dataset(subject_idx).preprocessed(d).header.name, 'RS') 
        % provide update
        fprintf('%s:\n', dataset(subject_idx).preprocessed(d).header.name(27:end))

        % select data
        lwdata.header = dataset(subject_idx).preprocessed(d).header;
        lwdata.data = dataset(subject_idx).preprocessed(d).data;

        % re-label and filter events
        fprintf('checking events... ')
        event_idx = logical([]);
        for a = 1:length(lwdata.header.events)
            if strcmp(lwdata.header.events(a).code, params.eventcode{1})
                lwdata.header.events(a).code = params.eventcode_new{1};
                event_idx(a) = false; 
            elseif strcmp(lwdata.header.events(a).code, params.eventcode{2})
                lwdata.header.events(a).code = params.eventcode_new{2};
                event_idx(a) = false; 
            else
                event_idx(a) = true; 
            end
        end
        lwdata.header.events(event_idx) = [];

        % check event number
        fprintf('%d %s stimuli found.\n', length(lwdata.header.events), lwdata.header.events(1).code)
        event_idx = false(1, length(lwdata.header.events));
        if length(lwdata.header.events) > params.event_n
            % ask which events should be removed
            prompt = {sprintf('More than %d events were found\n. Which should be removed?', params.event_n)};
            dlgtitle = sprintf('%s', dataset(subject_idx).preprocessed(d).header.name(27:end));
            dims = [1 40];
            definput = {''};
            input = inputdlg(prompt,dlgtitle,dims,definput);

            % update index
            event_idx(str2num(input{1,1})) = true;

            % remove faulty events
            lwdata.header.events(event_idx) = [];

            % update info structure
            if bad_counter == 1 
                RFSxLASER_info(subject_idx).preprocessing(9).process = sprintf('faulty triggers removed');
                RFSxLASER_info(subject_idx).preprocessing(9).suffix = [];
                RFSxLASER_info(subject_idx).preprocessing(9).date = sprintf('%s', date);
            end
            RFSxLASER_info(subject_idx).preprocessing(9).params(bad_counter).dataset = extractAfter(dataset(subject_idx).preprocessed(d).header.name, sprintf('%s ', RFSxLASER_info(subject_idx).ID));    
            if ~isempty(str2num(input{1,1}))
                RFSxLASER_info(subject_idx).preprocessing(9).params(bad_counter).trig_removed = find(event_idx);
            end

            % update bad counter
            bad_counter = bad_counter + 1;
        end

        % re-reference to common average
        fprintf('re-referencing to common average...')
        option = struct('reference_list', {{lwdata.header.chanlocs(1:length(lwdata.header.chanlocs)).labels}}, ...
            'apply_list', {{lwdata.header.chanlocs(1:length(lwdata.header.chanlocs)).labels}}, 'suffix', params.suffix{5}, 'is_save', 0);
        lwdata = FLW_rereference.get_lwdata(lwdata, option);
        if d == 3
            RFSxLASER_info(subject_idx).preprocessing(10).process = sprintf('ERP data re-referenced to common average');
            RFSxLASER_info(subject_idx).preprocessing(10).suffix = params.suffix{5};
            RFSxLASER_info(subject_idx).preprocessing(10).date = sprintf('%s', date);
        end

        % segment
        fprintf('epoching ...')
        option = struct('event_labels', {lwdata.header.events(1).code}, 'x_start', params.epoch(1), 'x_end', params.epoch(2), ...
            'x_duration', params.epoch(2)-params.epoch(1), 'suffix', params.suffix{6}, 'is_save', 0);
        lwdata = FLW_segmentation.get_lwdata(lwdata, option);
        if d == 3
            RFSxLASER_info(subject_idx).preprocessing(11).process = sprintf('ERP data segmented');
            RFSxLASER_info(subject_idx).preprocessing(11).params.limits = params.epoch;
            RFSxLASER_info(subject_idx).preprocessing(11).suffix = params.suffix{6};
            RFSxLASER_info(subject_idx).preprocessing(11).date = sprintf('%s', date);
        end
    
        % remove DC + linear detrend
        fprintf('removing DC and applying linear detrend...\n')
        option = struct('linear_detrend', 1, 'suffix', params.suffix{7}, 'is_save', 1);
        lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
        if d == 3
            RFSxLASER_info(subject_idx).preprocessing(12).process = sprintf('DC + linear detrend on ERP epochs');
            RFSxLASER_info(subject_idx).preprocessing(12).suffix = params.suffix{7};
            RFSxLASER_info(subject_idx).preprocessing(12).date = sprintf('%s', date);
        end

        % update dataset
        dataset(subject_idx).preprocessed(d).header = lwdata.header;
        dataset(subject_idx).preprocessed(d).data = lwdata.data; 
    end
end
if isempty(RFSxLASER_info(subject_idx).preprocessing(9).process)        
    RFSxLASER_info(subject_idx).preprocessing(9).process = sprintf('no faulty triggers removed');
    RFSxLASER_info(subject_idx).preprocessing(9).suffix = [];
    RFSxLASER_info(subject_idx).preprocessing(9).date = sprintf('%s', date);
end
fprintf('done.\n')

% save to the output file
save(output_file, 'RFSxLASER_info', '-append');

% clear and continue
fprintf('\nsection 2 finished.\n\n')
clear params a b c d d_rfs i shift_samples trigger visual event_idx lwdata chans2interpolate chan_n chan_dist chans2use...
    fig option bad_counter bad_idx answer prompt dlgtitle dims definput input

%% 3) discard bad trials 
% ----- section input -----
params.prefix = 'dc ep reref ds notch bandpass dc';
params.n_files = 8;
params.suffix = 'ar';
% -------------------------
fprintf('section 3: bad trials\n')

% update info structure
load(output_file, 'RFSxLASER_info');
cd(folder.processed)

% define the names of checked files
file2process = dir(sprintf('%s*%s*.mat', params.prefix, RFSxLASER_info(subject_idx).ID));
if length(file2process) == params.n_files
    for i = 1:length(file2process)
        filenames{i} = sprintf('%s %s',params.suffix, file2process(i).name);
    end
else
    error('ERROR: Incorrect number of datasets found in the directory: %d\n', length(file2process));
end

% open letswave if not already open
addpath(genpath([folder.toolbox '\letswave 6']));
fig_all = findall(0, 'Type', 'figure');
open = true;
for f = 1:length(fig_all)
    if contains(get(fig_all(f), 'Name'), 'Letswave', 'IgnoreCase', true)
        open = false;
        break;
    end
end
if open
    fprintf('opening letswave:\n')
    letswave
end

% wait until all files are processed
wait4files(filenames);

% load dataset with bad trials removed
fprintf('updating dataset... ')
if exist('dataset') ~= 1
    % load checked data
    data2load = dir(sprintf('%s*%s*', params.suffix, RFSxLASER_info(subject_idx).ID));
    dataset = reload_dataset(data2load, subject_idx, 'checked');
else
    % append checked dataset
    dataset_old = dataset;
    data2load = dir(sprintf('%s*%s*', params.suffix, RFSxLASER_info(subject_idx).ID));
    dataset = reload_dataset(data2load, subject_idx, 'checked');
    dataset_new = dataset;
    [dataset_old.checked] = dataset_new.checked;
    dataset = dataset_old;
    clear dataset_old dataset_new
end
fprintf('done.\n')

% encode bad trials
fprintf('encoding bad trials... ')
RFSxLASER_info(subject_idx).preprocessing(13).process = sprintf('bad trials discarded');
RFSxLASER_info(subject_idx).preprocessing(13).params.GUI = 'letswave';
RFSxLASER_info(subject_idx).preprocessing(13).suffix = params.suffix;
RFSxLASER_info(subject_idx).preprocessing(13).date = sprintf('%s', date);
for a = 1:length(dataset(subject_idx).checked)
    % subset header
    header = dataset(subject_idx).checked(a).header;

    % extract discarded expochs
    if ~isempty(header.history(end).configuration)
        if ~isempty(header.history(end).configuration.parameters.rejected_epochs)
            discarded = header.history(end).configuration.parameters.rejected_epochs;
        else
            discarded = [];
        end
    end

    % encode 
    RFSxLASER_info(subject_idx).preprocessing(13).params.discarded(a).dataset = extractAfter(header.name, sprintf('%s ', RFSxLASER_info(subject_idx).ID));
    RFSxLASER_info(subject_idx).preprocessing(13).params.discarded(a).trials = discarded;
    RFSxLASER_info(subject_idx).preprocessing(13).params.kept(a).dataset = extractAfter(header.name, sprintf('%s ', RFSxLASER_info(subject_idx).ID));
    RFSxLASER_info(subject_idx).preprocessing(13).params.kept(a).trials = header.datasize(1);
end
fprintf('done.\n')

% save to the output file
save(output_file, 'RFSxLASER_info', '-append');

% clear and continue
fprintf('\nsection 3 finished.\n\n')
clear params a f i file2process filenames fig_all open header data2load discarded

%% 4) compute ICA
% ----- section input -----
params.prefix = 'ar dc ep reref ds notch bandpass dc';
params.suffix = 'ica';
params.ICA_comp = 30;
% -------------------------
fprintf('section 4: compute ICA\n')

% update info structure
load(output_file, 'RFSxLASER_info');
cd(folder.processed)

% load dataset if needed
if exist('dataset') ~= 1
    fprintf('re-loading dataset... ')
    data2load = dir(sprintf('%s*%s*', params.prefix, RFSxLASER_info(subject_idx).ID));
    dataset = reload_dataset(data2load, subject_idx, 'checked');
    fprintf('done.\n')
end

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% select dataset
lwdataset = dataset(subject_idx).checked;

% compute ICA and save  
fprintf('computing ICA matrix:\n')
option = struct('ICA_mode', 2, 'algorithm', 1, 'num_ICs', params.ICA_comp, 'suffix', params.suffix, 'is_save', 1);
lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
fprintf('done.\n')

% extract ICA parameters
matrix.mix = lwdataset(1).header.history(end).option.mix_matrix;
matrix.unmix = lwdataset(1).header.history(end).option.unmix_matrix;    
params.ICA_chanlocs = lwdataset(1).header.chanlocs;
for i = 1:size(matrix.mix, 2)
    params.ICA_labels{i} = ['IC',num2str(i)];
end
params.ICA_SR = 1/lwdataset(1).header.xstep;

% update dataset and adjust for letswave 6
dataset(subject_idx).ica = lwdataset;
for a = 1:length(dataset(subject_idx).ica)
    dataset(subject_idx).ica(a).header.history(11).configuration.gui_info.function_name = 'LW_ICA_compute_merged';  
    dataset(subject_idx).ica(a).header.history(11).configuration.parameters = dataset(subject_idx).ica(a).header.history(11).option;  
    [dataset(subject_idx).ica(a).header.history(11).configuration.parameters.ICA_um] = dataset(subject_idx).ica(a).header.history(11).configuration.parameters.unmix_matrix; 
    [dataset(subject_idx).ica(a).header.history(11).configuration.parameters.ICA_mm] = dataset(subject_idx).ica(a).header.history(11).configuration.parameters.mix_matrix; 
    dataset(subject_idx).ica(a).header.history(11).configuration.parameters = rmfield(dataset(subject_idx).ica(a).header.history(11).configuration.parameters, {'unmix_matrix' 'mix_matrix'});
    header = dataset(subject_idx).ica(a).header;
    save(sprintf('%s.lw6', dataset(subject_idx).ica(a).header.name), 'header');
end

% unmix data
for b = 1:length(dataset(subject_idx).ica)
    for e = 1:size(dataset(subject_idx).ica(b).data, 1)
        dataset(subject_idx).unmixed(b).header = dataset(subject_idx).ica(b).header;
        dataset(subject_idx).unmixed(b).data(e, :, 1, 1, 1, :) = matrix.unmix * squeeze(dataset(subject_idx).ica(b).data(e, :, 1, 1, 1, :));        
    end
end

% update info structure
RFSxLASER_info(subject_idx).preprocessing(14).process = 'ICA matrix computed';
RFSxLASER_info(subject_idx).preprocessing(14).params.method = 'runica';
RFSxLASER_info(subject_idx).preprocessing(14).params.components = params.ICA_comp;
RFSxLASER_info(subject_idx).preprocessing(14).params.chanlocs = params.ICA_chanlocs;
RFSxLASER_info(subject_idx).preprocessing(14).params.labels = params.ICA_labels;
RFSxLASER_info(subject_idx).preprocessing(14).params.SR = params.ICA_SR;
RFSxLASER_info(subject_idx).preprocessing(14).params.matrix = matrix;
RFSxLASER_info(subject_idx).preprocessing(14).suffix = params.suffix;
RFSxLASER_info(subject_idx).preprocessing(14).date = sprintf('%s', date);

% calculate PSD across all timepoints, components and trials 
fprintf('estimating spectral content...\n')
for c = 1:length(dataset(subject_idx).ica)
    for d = 1:params.ICA_comp
        for e = 1:size(dataset(subject_idx).unmixed(c).data, 1)
            [psd(c, d, e, :), freq] = pwelch(squeeze(dataset(subject_idx).unmixed(c).data(e, d, 1, 1, 1, :)), ...
                [], [], [], RFSxLASER_info(subject_idx).preprocessing(14).params.SR);  
        end
    end
end
psd = squeeze(mean(psd, [1, 3]));
RFSxLASER_info(subject_idx).preprocessing(14).params.PSD = psd;
RFSxLASER_info(subject_idx).preprocessing(14).params.freq = freq;
fprintf('done.\n')

% plot component topographies and spectral content
addpath(genpath([folder.toolbox '\letswave 6']));
figure('units','normalized','outerposition',[0 0 1 1]);
hold on
for f = 1:params.ICA_comp
    % plot the topography
    labels = {dataset(subject_idx).ica(1).header.chanlocs.labels};
    subplot(ceil(params.ICA_comp/3), 6, (f-1)*2 + 1);
    topoplot(double(matrix.mix(:, f)'), params.ICA_chanlocs, 'maplimits', [-4 4], 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
    set(gca,'color',[1 1 1]);
    title(params.ICA_labels{f})

    % plot the psd
    subplot(ceil(params.ICA_comp/3), 6, (f-1)*2 + 2);
    plot(freq(1:33), psd(f, 1:33));
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
end
sgtitle(sprintf('%s - ICA components', RFSxLASER_info(subject_idx).ID))
saveas(gcf, sprintf('%s\\figures\\%s_ICA.png', folder.output, RFSxLASER_info(subject_idx).ID))
figure_counter = figure_counter + 1;

% open letswave if not already open
fig_all = findall(0, 'Type', 'figure');
open = true;
for f = 1:length(fig_all)
    if contains(get(fig_all(f), 'Name'), 'Letswave', 'IgnoreCase', true)
        open = false;
        break;
    end
end
if open
    fprintf('opening letswave:\n')
    letswave
end

% save to the output file
save(output_file, 'RFSxLASER_info', '-append');

% clear and continue
fprintf('\nsection 4 finished.\nproceed to ICA\n\n')
clear params data2load lwdataset option a b c d e f i matrix header psd freq labels fig_all open

%% 5) encode ICA
% ----- section input -----
params.suffix = 'icfilt';
params.ICA_comp = 30;
params.plot_toi = [-0.1 0.5];
params.eoi = 'Cz';
% -------------------------
fprintf('section 5: encode ICA\n')
 
% update info structure
load(output_file, 'RFSxLASER_info');
cd(folder.processed)

% ask for subject number if necessary
if ~exist('subject_idx')
    prompt = {'subject number:'};
    dlgtitle = 'subject';
    dims = [1 40];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    subject_idx = str2num(input{1,1});
end
clear prompt dlgtitle dims definput input

% ask for input 
prompt = {'blinks:', 'horizontal:', 'muscles:', 'electrode:'};
dlgtitle = 'ICA';  
dims = [1 60];
definput = {'', '', '', ''};
input = inputdlg(prompt,dlgtitle,dims,definput);

% encode & save 
RFSxLASER_info(subject_idx).preprocessing(15).process = 'artifactual ICs discarded';
RFSxLASER_info(subject_idx).preprocessing(15).suffix = params.suffix;
RFSxLASER_info(subject_idx).preprocessing(15).date = sprintf('%s', date);
RFSxLASER_info(subject_idx).preprocessing(15).params.kept = params.ICA_comp - length([str2num(input{1}), str2num(input{2}), str2num(input{3}), str2num(input{4})]);
RFSxLASER_info(subject_idx).preprocessing(15).params.removed.blinks = str2num(input{1});
RFSxLASER_info(subject_idx).preprocessing(15).params.removed.horizontal = str2num(input{2});
RFSxLASER_info(subject_idx).preprocessing(15).params.removed.muscles = str2num(input{3});
RFSxLASER_info(subject_idx).preprocessing(15).params.removed.electrode = str2num(input{4});
save(output_file, 'RFSxLASER_info', '-append')

% clear and continue
fprintf('\nsection 5 finished.\n\n')
clear params data2load prompt dlgtitle dims definput input 

%% 6) import and process PsychoPy data
% ----- section input -----
params.stimulus = {'laser', 'RFS'};
params.side = {'right', 'left'};
params.intensity = {'high' 'low'};
params.folder = 'PsychoPy_blocks';
% -------------------------
fprintf('section 6: PsychoPy data import & processing\n')

% update info structures
output_vars = who('-file', output_file);
if ismember('RFSxLASER_info', output_vars)
    load(output_file, 'RFSxLASER_info')
else
    error('ERROR: the output file does not contain the subject & session information!')
end
if ismember('RFSxLASER_measures', output_vars)
    load(output_file, 'RFSxLASER_measures')
else
    RFSxLASER_measures = struct;
    save(output_file, 'RFSxLASER_measures', '-append')
    fprintf('WARNING: the output file does not contain the measures output\n --> creating new strucutre now\n')
end

% ask for subject number, if not defined
if ~exist('subject_idx')
    prompt = {'subject number:'};
    dlgtitle = 'subject';
    dims = [1 40];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    subject_idx = str2num(input{1,1});
end
clear prompt dlgtitle dims definput input

% load the .CSV
fprintf('importing PsychoPy data...\n')
file2import = dir(sprintf('%s\\%s\\%s\\*.csv', folder.raw, RFSxLASER_info(subject_idx).ID, params.folder));
if size(file2import, 1) ~= 1
    error('ERROR: %d files were found instead of expected 1!\n', size(file2import, 1))
end
option = detectImportOptions(sprintf('%s\\%s', file2import.folder, file2import.name));
option.SelectedVariableNames = {'block' 'condition', 'stim1_trials_thisRepN', 'slider_response', 'slider_rt'};
ratings_table = readtable(sprintf('%s\\%s', file2import.folder, file2import.name), option);

% formate and append to the output structure
idx = logical([]);
for a = 1:height(ratings_table)
    if isnan(ratings_table.block(a))
        idx(a) = false;
    else
        idx(a) = true;
        ratings_table.trial(a) = ratings_table.stim1_trials_thisRepN(a) + 1;
    end
end
ratings_table = ratings_table(idx, [1,2,6,4,5]);
RFSxLASER_measures.individual(subject_idx).ratings.raw = ratings_table;

% calculate basic statistics and append to the output structure
RFSxLASER_measures.individual(subject_idx).ratings.stats = struct([]);
for a = 1:length(params.stimulus)
    for b = 1:length(params.side)
        for c = 1:length(params.intensity)
            % identify condition label
            if strcmp(params.stimulus{a}, 'laser')
                label = 'LS_';
            elseif strcmp(params.stimulus{a}, 'RFS')
                label = 'RF_';
            end
            if strcmp(params.side{b}, 'right')
                label = [label 'RH_'];
            elseif strcmp(params.side{b}, 'left')
                label = [label 'LH_'];
            end
            if strcmp(params.intensity{c}, 'high')
                label = [label 'high'];
            elseif strcmp(params.intensity{c}, 'low')
                label = [label 'low'];
            end

            % subset the data
            idx = logical([]);
            for d = 1:height(ratings_table)
                if strcmp(ratings_table.condition(d), label)
                    idx(d) = true;
                else
                    idx(d) = false;
                end
            end
            data = [ratings_table.slider_response(idx)];

            if ~isempty(data)
                % remove missed or faulty trials
                for d = 1:length(RFSxLASER_info(subject_idx).preprocessing(13).params.discarded)  
                    if contains(RFSxLASER_info(subject_idx).preprocessing(13).params.discarded(d).dataset, params.stimulus{a}) && ...
                            contains(RFSxLASER_info(subject_idx).preprocessing(13).params.discarded(d).dataset, params.side{b}) && ...
                            contains(RFSxLASER_info(subject_idx).preprocessing(13).params.discarded(d).dataset, params.intensity{c}) 
                        bad_trials = RFSxLASER_info(subject_idx).preprocessing(13).params.discarded(d).trials;
                    end
                end
                data(bad_trials) = [];

                % encode conditions
                RFSxLASER_measures.individual(subject_idx).ratings.stats(end + 1).stimulus = params.stimulus{a};
                RFSxLASER_measures.individual(subject_idx).ratings.stats(end).side = params.side{b};
                RFSxLASER_measures.individual(subject_idx).ratings.stats(end).intensity = params.intensity{c};

                % calculate and append stats
                % t_value = tinv(0.975, length(data) - 1); 
                RFSxLASER_measures.individual(subject_idx).ratings.stats(end).mean = mean(data, 1);
                RFSxLASER_measures.individual(subject_idx).ratings.stats(end).SD = std(data, 0, 1);
                RFSxLASER_measures.individual(subject_idx).ratings.stats(end).SEM = RFSxLASER_measures.individual(subject_idx).ratings.stats(end).SD / sqrt(length(data)); 
                % RFSxLASER_measures(subject_idx).ratings.stats(end).CI_upper = RFSxLASER_measures(subject_idx).ratings.stats(end).mean + t_value * RFSxLASER_measures(subject_idx).ratings.stats(end).SEM; 
                % RFSxLASER_measures(subject_idx).ratings.stats(end).CI_lower = RFSxLASER_measures(subject_idx).ratings.stats(end).mean - t_value * RFSxLASER_measures(subject_idx).ratings.stats(end).SEM; 
            end
        end
    end
end
save(output_file, 'RFSxLASER_measures', '-append')

% clear and continue
fprintf('\nsection 6 finished.\n\n')
clear a b c d output_vars file2import option ratings_table idx label data bad_trials t_value 

%% 7) extraction of the N1 component from ERPs
% ----- section input -----
params.prefix = 'icfilt ica ar dc ep reref ds notch bandpass dc';
params.ref = 'AFz';
params.suffix = {'N1_reref' 'N1_ica' 'icfilt' 'N1_ready'};
params.stimulus = {'laser', 'RFS'};
% -------------------------
fprintf('section 7: N1 component extraction\n')

% update info structures
output_vars = who('-file', output_file);
if ismember('RFSxLASER_info', output_vars)
    load(output_file, 'RFSxLASER_info')
else
    error('ERROR: the output file does not contain the subject & session information!')
end
if ismember('RFSxLASER_measures', output_vars)
    load(output_file, 'RFSxLASER_measures')
else
    RFSxLASER_measures = struct;
    save(output_file, 'RFSxLASER_measures', '-append')
    fprintf('WARNING: the output file does not contain the measures output\n --> creating new strucutre now\n')
end

% ask for subject number, if not defined
if ~exist('subject_idx')
    prompt = {'subject number:'};
    dlgtitle = 'subject';
    dims = [1 40];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    subject_idx = str2num(input{1,1});
end
clear prompt dlgtitle dims definput input

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));
cd(folder.processed)

% load the data
fprintf('loading... ')
data2load = dir(sprintf('%s\\%s*%s*', folder.raw, params.prefix, RFSxLASER_info(subject_idx).ID));
dataset = reload_dataset(data2load, 1, 'icfilt');
fprintf('%d datasets found and loaded.\n', length(dataset.icfilt))

% re-reference to chosen frontal central electrode
fprintf('re-referencing to %s...\n', params.ref)
for d = 1:length(dataset.icfilt)
    % subset the dataset
    lwdata = dataset.icfilt(d);

    % re-reference
    option = struct('reference_list', {{params.ref}}, 'apply_list', {{lwdata.header.chanlocs.labels}},...
    'suffix', params.suffix{1}, 'is_save', 0);
    lwdata = FLW_rereference.get_lwdata(lwdata, option);  

    % update dataset
    dataset.N1(d) = lwdata;

    % encode to info structure
    if d == 1
        RFSxLASER_info(subject_idx).preprocessing(16).process = 're-referenced for N1 extraction';
        RFSxLASER_info(subject_idx).preprocessing(16).params.reference = params.ref;
        RFSxLASER_info(subject_idx).preprocessing(16).suffix = params.suffix{1};
        RFSxLASER_info(subject_idx).preprocessing(16).date = sprintf('%s', date);
    end
end 

% prepare flip dictionary
params.labels = {dataset.icfilt(1).header.chanlocs.labels};
params.chanlocs = dataset.icfilt(1).header.chanlocs;
labels_flipped = params.labels;
for i = 1:length(params.labels)
    electrode_n = str2num(params.labels{i}(end));
    if isempty(electrode_n)
    else
        if electrode_n == 0
            label_new = params.labels{i}(1:end-2);
            label_new = [label_new num2str(9)];
        elseif mod(electrode_n,2) == 1              % odd number --> left hemisphere                    
            label_new = params.labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n + 1)];
        else                                    % even number --> right hemisphere 
            label_new = params.labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n - 1)];
        end
        a = find(strcmpi(params.labels, label_new));
        if isempty(a)
        else
            labels_flipped{i} = label_new;
        end
    end
end
labels_dict = cat(1, params.labels, labels_flipped)';

% compute ICA separately for RFS and laser
IC_n = RFSxLASER_info(subject_idx).preprocessing(15).params.kept - 1;
for s = 1:length(params.stimulus)
    % update
    fprintf('%s stimulation:\n', params.stimulus{s})

    % select datasets for this stimulus
    data_idx = logical([]);
    for d = 1:length(dataset.icfilt)
        if contains(dataset.icfilt(d).header.name, params.stimulus{s})
            data_idx(d) = true;
        else
            data_idx(d) = false;
        end
    end
    lwdataset = dataset.N1(data_idx);

    % flip datasets that were stimulated on the left hand 
    addpath(genpath([folder.toolbox '\letswave 6']));
    for d = 1:length(lwdataset)
        % flip if leftsided
        if contains(dataset.icfilt(d).header.name, 'left')
            [lwdataset(d).header, lwdataset(d).data, ~] = RLW_flip_electrodes(lwdataset(d).header, lwdataset(d).data, labels_dict);
        end
    end

    % calculate ICA matrix
    addpath(genpath([folder.toolbox '\letswave 7']));
    fprintf('computing ICA matrix:\n')
    option = struct('ICA_mode', 2, 'algorithm', 1, 'num_ICs', IC_n, 'suffix', params.suffix{2}, 'is_save', 1);
    lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
    fprintf('done.\n')

    % extract ICA parameters
    matrix.mix = lwdataset(1).header.history(end).option.mix_matrix;
    matrix.unmix = lwdataset(1).header.history(end).option.unmix_matrix;    
    for i = 1:size(matrix.mix, 2)
        params.ICA_labels{i} = ['IC',num2str(i)];
    end
    params.ICA_SR = 1/lwdataset(1).header.xstep;

    % adjust for letswave 6 and update dataset  
    for a = 1:length(lwdataset)
        lwdataset(a).header.history(end).configuration.gui_info.function_name = 'LW_ICA_compute_merged';  
        lwdataset(a).header.history(end).configuration.parameters = lwdataset(a).header.history(end).option;  
        [lwdataset(a).header.history(end).configuration.parameters.ICA_um] = lwdataset(a).header.history(end).configuration.parameters.unmix_matrix; 
        [lwdataset(a).header.history(end).configuration.parameters.ICA_mm] = lwdataset(a).header.history(end).configuration.parameters.mix_matrix; 
        lwdataset(a).header.history(end).configuration.parameters = rmfield(lwdataset(a).header.history(end).configuration.parameters, {'unmix_matrix' 'mix_matrix'});
        header = lwdataset(a).header;
        save(sprintf('%s.lw6', lwdataset(a).header.name), 'header');
    end
    dataset.N1(data_idx) = lwdataset;

    % unmix data
    for a = 1:length(lwdataset)
        lwdataset_unmixed(a).header = lwdataset(a).header;
        for e = 1:size(lwdataset(a).data, 1)
            lwdataset_unmixed(a).data(e, :, 1, 1, 1, :) = matrix.unmix * squeeze(lwdataset(a).data(e, :, 1, 1, 1, :));        
        end
    end
    dataset.N1_unmixed(data_idx) = lwdataset_unmixed;

    % encode to info structure
    if s == 1
        RFSxLASER_info(subject_idx).preprocessing(17).process = 'ICA-based removal of N2P2 component';
        RFSxLASER_info(subject_idx).preprocessing(17).params.method = 'runica';
        RFSxLASER_info(subject_idx).preprocessing(17).params.components = IC_n;
        RFSxLASER_info(subject_idx).preprocessing(17).params.chanlocs = params.chanlocs;
        RFSxLASER_info(subject_idx).preprocessing(17).params.labels = params.ICA_labels;
        RFSxLASER_info(subject_idx).preprocessing(17).params.SR = params.ICA_SR;
        RFSxLASER_info(subject_idx).preprocessing(17).params.matrix.laser = matrix;
        RFSxLASER_info(subject_idx).preprocessing(17).suffix = params.suffix{2};
        RFSxLASER_info(subject_idx).preprocessing(17).date = sprintf('%s', date);
    elseif s == 2
        RFSxLASER_info(subject_idx).preprocessing(17).params.matrix.RFS = matrix;
    end

    % plot component topographies and mean signal at Cz
    addpath(genpath([folder.toolbox '\letswave 6']));
    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    for i = 1:IC_n
        % plot the topography
        subplot(ceil(IC_n/3), 6, (i-1)*2 + 1);
        topoplot(double(matrix.mix(:, i)'), params.chanlocs, 'maplimits', [-4 4], 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
        set(gca,'color',[1 1 1]);
        title(params.ICA_labels{i})

        % calculate mean timecourse and x
        for a = 1:length(lwdataset_unmixed)
            visual.data(a, :) = squeeze(mean(lwdataset_unmixed(a).data(:, i, 1, 1, 1, :), 1))'; 
        end
        visual.data = mean(visual.data, 1);
        visual.x = (0:lwdataset(a).header.datasize(6)-1)*lwdataset(a).header.xstep + lwdataset(a).header.xstart;
    
        % plot the psd
        subplot(ceil(IC_n/3), 6, (i-1)*2 + 2);
        plot_ERP(visual, 'xlim', [-0.1 0.8], 'colours', [0 0 0], 'legend', 'off', 'shading', 'off')
    end
    sgtitle(sprintf('%s - %s: N1 ICA components', RFSxLASER_info(subject_idx).ID, params.stimulus{s}))
    saveas(gcf, sprintf('%s\\figures\\%s_%s_N1-ICA.png', folder.output, RFSxLASER_info(subject_idx).ID, params.stimulus{s}))
    figure_counter = figure_counter + 1;
end

% open letswave if not already open
fig_all = findall(0, 'Type', 'figure');
open = true;
for f = 1:length(fig_all)
    if contains(get(fig_all(f), 'Name'), 'Letswave', 'IgnoreCase', true)
        open = false;
        break;
    end
end
if open
    fprintf('opening letswave:\n')
    addpath(genpath([folder.toolbox '\letswave 6']));
    letswave
end

% identify outcome filenames and wait for them to appear
for d = 1:length(data2load)
    filenames{d} = sprintf('%s %s',[params.suffix{3} ' ' params.suffix{2} ' ' params.suffix{1}], data2load(d).name);
end
wait4files(filenames);
close(gcf)

% load filtered data, flip back the datasets and save for letswave
fprintf('saving the final dataset... ')
for d = 1:length(dataset.icfilt)
    % load filtered data
    data_name = sprintf('%s %s',[params.suffix{3} ' ' params.suffix{2} ' ' params.suffix{1}], dataset.icfilt(d).header.name);
    load(sprintf('%s.lw6', data_name), '-mat');
    dataset.N1_filtered(d).header = header; 
    load(sprintf('%s.mat', data_name));
    dataset.N1_filtered(d).data = data; 

    % flip back if leftsided
    if contains(dataset.icfilt(d).header.name, 'left')
        [dataset.N1_filtered(d).header, dataset.N1_filtered(d).data, ~] = RLW_flip_electrodes(dataset.N1_filtered(d).header, dataset.N1_filtered(d).data, labels_dict);
    end

    % save with final suffix
    dataset.N1_filtered(d).header.name =sprintf('%s %s',params.suffix{4}, dataset.icfilt(d).header.name);
    header = dataset.N1_filtered(d).header;
    save(sprintf('%s.lw6', dataset.N1_filtered(d).header.name), 'header')
    data = dataset.N1_filtered(d).data;
    save(sprintf('%s.mat', dataset.N1_filtered(d).header.name), 'data')
end
fprintf('done.\n')

% encode ICA 
prompt = {'ICs removed for laser:' 'ICs removed for RFS:'};
dlgtitle = 'N1 ICA';
dims = [1 40];
definput = {'' ''};
input = inputdlg(prompt,dlgtitle,dims,definput);
RFSxLASER_info(subject_idx).preprocessing(17).params.N2P2_component.laser = split(input{1}, ',');
RFSxLASER_info(subject_idx).preprocessing(17).params.N2P2_component.RFS = split(input{2}, ',');
save(output_file, 'RFSxLASER_info', '-append')

% clear and continue
fprintf('\nsection 7 finished.\n\n')
clear a d e f i s output_vars data2load dataset lwdata option labels_flipped electrode_n label_new labels_dict data_idx lwdataset ...
        IC_n matrix lwdataset_unmixed fig_all open filenames data_name data header prompt dlgtitle dims definput input visual
 
% ask if the subject is done
answer = questdlg(sprintf('Do you want to continue to subject %d?', subject_idx + 1), 'Continue to next subject?', 'YES', 'NO', 'NO'); 
switch answer
    case 'NO'
    case 'YES'
    	subject_idx = subject_idx + 1;
        clear dataset answer
end

%% ===================== PART 3: group visualization & export for statistics ================
% directories
if ~exist('folder') 
    folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');    % MATLAB toolboxes
    folder.output = uigetdir(pwd, 'Choose the OneDrive folder');    % output folder --> figures, loutput file, exports 
end
folder.processed = uigetdir(pwd, 'Choose the data folder');         % processed data --> wherever you want to store the voluminous EEG data
cd(folder.output)

% output
study = 'RFSxLASER';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);
if ~exist('figure_counter')  
    figure_counter = 1;
end

% load the info structure
if exist(output_file) == 2
    output_vars = who('-file', output_file);
    % info
    if ismember('RFSxLASER_info', output_vars)
        load(output_file, 'RFSxLASER_info')
    else
        error('ERROR: the output file does not contain the subject & session information!')
    end

    % data
    if ismember('RFSxLASER_data', output_vars)
        load(output_file, 'RFSxLASER_data')
    else
        RFSxLASER_data = struct;
        save(output_file, 'RFSxLASER_data', '-append')
    end

    % measures
    if ismember('RFSxLASER_measures', output_vars)
        load(output_file, 'RFSxLASER_measures')
    else
        RFSxLASER_measures = struct;
        save(output_file, 'RFSxLASER_measures', '-append')
    end
else
    error('ERROR: output file not found!')
end
clear output_vars

%% 1) load data
% ----- section input -----
params.prefix = 'icfilt ica ar dc ep reref ds notch bandpass dc';
params.prefix_N1 = 'N1_ready';
params.side = {'right' 'left'};
params.stimulus = {'laser' 'RFS'};
params.intensity = {'high' 'low'};
params.baseline = [-0.3 -0.01];
% -------------------------
fprintf('section 1: load and prepare data\n')

% add letswave 6 to the top of search path
addpath(genpath([folder.toolbox '\letswave 6']));

% load data, normalize to baseline
% create subject-average original dataset
params.subjects = length(RFSxLASER_info);
subject_idx = logical([]);
dataset = [];
for s = 1:params.subjects
    fprintf('subject %d: ', s)

    % check available data
    data2load = dir(sprintf('%s\\%s*%s*', folder.processed, params.prefix, RFSxLASER_info(s).ID));
    
    % loop through datasets
    if ~isempty(data2load)
        fprintf('%d datasets found.\nloading... ', length(data2load))

        % update index
        subject_idx(s) = true;

        % loop through datasets
        for d = 1:length(data2load)
            if contains(data2load(d).name, 'lw6')
                % load the dataset
                [header, data] = CLW_load(sprintf('%s\\%s', data2load(d).folder, data2load(d).name));

                % normalize as z-score of the baseline 
                [header, data, ~] = RLW_baseline(header, data, 'operation', 'zscore', 'xstart', params.baseline(1), 'xend', params.baseline(2));

                % identify dataset
                if contains(data2load(d).name, params.stimulus{1})
                    thisdata.stimulus = params.stimulus{1};
                elseif contains(data2load(d).name, params.stimulus{2})
                    thisdata.stimulus = params.stimulus{2};
                end
                if contains(data2load(d).name, params.intensity{1})
                    thisdata.intensity = params.intensity{1};
                elseif contains(data2load(d).name, params.intensity{2})
                    thisdata.intensity = params.intensity{2};
                end
                if contains(data2load(d).name, params.side{1})
                    thisdata.side = params.side{1};
                elseif contains(data2load(d).name, params.side{2})
                    thisdata.side = params.side{2};
                end

                % save subject averages to the dataset structure
                statement = sprintf('dataset.%s.%s.%s(s, :, :) =  squeeze(mean(data, 1));', thisdata.stimulus, thisdata.intensity, thisdata.side);
                eval(statement)
            end
        end
        fprintf('done.\n')
    else
        % update index
        subject_idx(s) = false;
        fprintf('WARNING: no data were found!\nskipping this dataset...\n')
    end
end
RFSxLASER_data.original = dataset;

% prepare flip dictionary
params.labels = {header.chanlocs.labels};
params.chanlocs = header.chanlocs;
labels_flipped = params.labels;
for i = 1:length(params.labels)
    electrode_n = str2num(params.labels{i}(end));
    if isempty(electrode_n)
    else
        if electrode_n == 0
            label_new = params.labels{i}(1:end-2);
            label_new = [label_new num2str(9)];
        elseif mod(electrode_n,2) == 1              % odd number --> left hemisphere                    
            label_new = params.labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n + 1)];
        else                                    % even number --> right hemisphere 
            label_new = params.labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n - 1)];
        end
        a = find(strcmpi(params.labels, label_new));
        if isempty(a)
        else
            labels_flipped{i} = label_new;
        end
    end
end
labels_dict = cat(1, params.labels, labels_flipped)';

% flip normalized data to homogenize side of stimulation --> if right, flip
% create subject-average flipped dataset
% + check for missing datasets
fprintf('flipping data: ')
addpath(genpath([folder.toolbox '\letswave 6']));
header.datasize(1) = 1;
subject_idx = logical([]);
data_all = [];
for s = 1:params.subjects 
    fprintf('. ')
    for a = 1:length(params.stimulus)
        for b = 1:length(params.intensity)
            % check if data exist
            for c = 1:length(params.side)
                statement = sprintf('data = squeeze(dataset.%s.%s.%s(%d, :, :));', params.stimulus{a}, params.intensity{b}, params.side{c}, s);
                eval(statement)
                if sum(data, 'all') == 0
                    subject_idx(a, b, c, s) = false;
                else
                    subject_idx(a, b, c, s) = true;
                end
            end
            data = [];

            % load data from right-sided stimulation
            statement = sprintf('data_all(1, :, :) = squeeze(dataset.%s.%s.right(s, :, :));', params.stimulus{a}, params.intensity{b});
            eval(statement)
            
            % flip data from left-sided stimulation
            if subject_idx(a, b, 2, s)
                % load data
                statement = sprintf('data(1, :, 1, 1, 1, :) = squeeze(dataset.%s.%s.left(s, :, :));', params.stimulus{a}, params.intensity{b});
                eval(statement)
    
                % flip the data and append
                [header, data, ~] = RLW_flip_electrodes(header, data, labels_dict);
                data_all(2, :, :) = squeeze(data);     
            end

            % save to new dataset
            statement = sprintf('dataset.flipped.%s.%s(s, :, :) = squeeze(mean(data_all, 1));', params.stimulus{a}, params.intensity{b});
            eval(statement)
        end
    end
end
fprintf('done.\n')
RFSxLASER_data.original.subject_idx = subject_idx;
RFSxLASER_data.flipped = dataset.flipped;

% load N1-preprocessed data, normalize to baseline
% create subject-average original_N1 dataset
subject_idx = logical([]);
dataset = [];
for s = 19:params.subjects
    fprintf('subject %d: ', s)

    % check available data
    data2load = dir(sprintf('%s\\%s*%s*', folder.processed, params.prefix_N1, RFSxLASER_info(s).ID));
    
    % loop through datasets
    if ~isempty(data2load)
        fprintf('%d datasets found.\nloading... ', length(data2load))

        % update index
        subject_idx(s) = true;

        % loop through datasets
        for d = 1:length(data2load)
            if contains(data2load(d).name, 'lw6')
                % load the dataset
                [header, data] = CLW_load(sprintf('%s\\%s', data2load(d).folder, data2load(d).name));

                % normalize as z-score of the baseline 
                [header, data, ~] = RLW_baseline(header, data, 'operation', 'zscore', 'xstart', params.baseline(1), 'xend', params.baseline(2));

                % identify dataset
                if contains(data2load(d).name, params.stimulus{1})
                    thisdata.stimulus = params.stimulus{1};
                elseif contains(data2load(d).name, params.stimulus{2})
                    thisdata.stimulus = params.stimulus{2};
                end
                if contains(data2load(d).name, params.intensity{1})
                    thisdata.intensity = params.intensity{1};
                elseif contains(data2load(d).name, params.intensity{2})
                    thisdata.intensity = params.intensity{2};
                end
                if contains(data2load(d).name, params.side{1})
                    thisdata.side = params.side{1};
                elseif contains(data2load(d).name, params.side{2})
                    thisdata.side = params.side{2};
                end

                % save subject averages to the dataset structure
                statement = sprintf('dataset.%s.%s.%s(s, :, :) =  squeeze(mean(data, 1));', thisdata.stimulus, thisdata.intensity, thisdata.side);
                eval(statement)
            end
        end
        fprintf('done.\n')
    else
        % update index
        subject_idx(s) = false;
        fprintf('WARNING: no data were found!\nskipping this dataset...\n')
    end
end
RFSxLASER_data.original_N1 = dataset;

% flip normalized data to homogenize side of stimulation --> if right, flip
% create subject-average flipped dataset
% + check for missing datasets
fprintf('flipping data: ')
addpath(genpath([folder.toolbox '\letswave 6']));
header.datasize(1) = 1;
subject_idx = logical([]);
data_all = [];
for s = 1:params.subjects 
    fprintf('. ')
    for a = 1:length(params.stimulus)
        for b = 1:length(params.intensity)
            % check if data exist
            for c = 1:length(params.side)
                statement = sprintf('data = squeeze(dataset.%s.%s.%s(%d, :, :));', params.stimulus{a}, params.intensity{b}, params.side{c}, s);
                eval(statement)
                if sum(data, 'all') == 0
                    subject_idx(a, b, c, s) = false;
                else
                    subject_idx(a, b, c, s) = true;
                end
            end
            data = [];

            % load data from right-sided stimulation
            statement = sprintf('data_all(1, :, :) = squeeze(dataset.%s.%s.right(s, :, :));', params.stimulus{a}, params.intensity{b});
            eval(statement)
            
            % flip data from left-sided stimulation
            if subject_idx(a, b, 2, s)
                % load data
                statement = sprintf('data(1, :, 1, 1, 1, :) = squeeze(dataset.%s.%s.left(s, :, :));', params.stimulus{a}, params.intensity{b});
                eval(statement)
    
                % flip the data and append
                [header, data, ~] = RLW_flip_electrodes(header, data, labels_dict);
                data_all(2, :, :) = squeeze(data);     
            end

            % save to new dataset
            statement = sprintf('dataset.flipped.%s.%s(s, :, :) = squeeze(mean(data_all, 1));', params.stimulus{a}, params.intensity{b});
            eval(statement)
        end
    end
end
fprintf('done.\n')
RFSxLASER_data.original_N1.subject_idx = subject_idx;
RFSxLASER_data.flipped_N1 = dataset.flipped;

% save and continue
save(output_file, 'RFSxLASER_data', '-append')
clear a b c d i s data2load data header thisdata dataset statement subject_idx labels_flipped electrode_n labels_dict data_all label_new
fprintf('section 1 finished.\n\n')

%% 2) N2P2 ERP visualization and peak measure extraction
% ----- section input -----
params.eoi = 'Cz';
params.prefix = 'N1_ready';
params.stimulus = {'laser' 'RFS'};
params.intensity = {'high' 'low'};
params.ERP_measures = {'N2_latency' 'N2_amplitude' 'P2_latency' 'P2_amplitude'};
params.search_window = [0.05 0.75];
params.outliers = [16];
% -------------------------
fprintf('section 2: group statistics - ERP\n')

% select dataset
dataset = RFSxLASER_data.flipped;

% define subject filter
params.subjects = length(RFSxLASER_info);
subject_filter = logical(ones(1, params.subjects)); 
subject_filter(params.outliers) = false;

% define common visualization & stats parameters
data2load = dir(sprintf('%s\\%s*%s*', folder.processed, params.prefix));
load(sprintf('%s\\%s', data2load(1).folder, data2load(1).name), '-mat')
t_value = tinv(0.975, sum(subject_filter) - 1); 
visual.x = (0:header.datasize(6)-1)*header.xstep + header.xstart;
visual.labels = params.stimulus;
visual.chanlocs = header.chanlocs;
visual.chanlabels = {header.chanlocs.labels};
visual.eoi = find(strcmp(visual.chanlabels, params.eoi));
visual.colors = [0.3333    0.4471    0.9020;
    1.0000    0.0745    0.6510];
screen_size = get(0, 'ScreenSize');

% plot N2P2 ERPs for each intensity
for b = 1:length(params.intensity)
    % select laser data
    statement = sprintf('data = squeeze(dataset.laser.%s(subject_filter, visual.eoi, :));', params.intensity{b});
    eval(statement)
    visual.data(1, :) = mean(data, 1);
    visual.sem(1, :) = squeeze(std(data, 0, 1)) / sqrt(sum(subject_filter)); 
    visual.CI_upper(1, :) = visual.data(1, :) + t_value * visual.sem(1, :); 
    visual.CI_lower(1, :) = visual.data(1, :) - t_value * visual.sem(1, :);

    % select RFS data
    statement = sprintf('data = squeeze(dataset.RFS.%s(subject_filter, visual.eoi, :));', params.intensity{b});
    eval(statement)
    visual.data(2, :) = mean(data, 1);
    visual.sem(2, :) = squeeze(std(data, 0, 1)) / sqrt(sum(subject_filter)); 
    visual.CI_upper(2, :) = visual.data(2, :) + t_value * visual.sem(2, :); 
    visual.CI_lower(2, :) = visual.data(2, :) - t_value * visual.sem(2, :);

    % launch the figure
    fig = figure(figure_counter);    
    set(fig, 'Position', [screen_size(3)/4, screen_size(4)/4, 2*screen_size(3)/5, screen_size(4) / 2])

    % plot
    plot_ERP(visual, 'colours', visual.colors, 'labels', visual.labels, 'xlim', [-0.1 0.6], 'ylim', [-3.5 3.5], 'reverse', 'on')
    title(sprintf('%s stimulation intensity', params.intensity{b}))

    % save figure and update counter
    saveas(fig, sprintf('%s\\figures\\preliminary_%s_%s.png', folder.output, params.eoi, params.intensity{b}))
    figure_counter = figure_counter + 1;
end

% plot topographies at group average peak latencies
search_window = find(visual.x >= params.search_window(1) & visual.x <= params.search_window(2));
addpath(genpath([folder.toolbox '\letswave 6']));
for a = 1:length(params.stimulus)
    for b = 1:length(params.intensity)
        % subset dataset
        statement = sprintf('data = squeeze(mean(dataset.%s.%s(subject_filter, :, :), 1));', params.stimulus{a}, params.intensity{b});
        eval(statement)

        % indentify N2 idx
        [N2_amplitude, idx] = min(data(visual.eoi, search_window));
        N2_idx = search_window(1) + idx;

        % prepare data for visualization
        visual.data = data(:, N2_idx)';

        % plot & save N2 topography
        figure(figure_counter);
        topoplot(double(visual.data), visual.chanlocs, 'maplimits', [-2 2], 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
        set(gca,'color',[1 1 1]);
        title(sprintf('%s, %s intensity - N2', params.stimulus{a}, params.intensity{b}))
        saveas(gcf, sprintf('%s\\figures\\preliminary_topo_%s_%s_N2.png', folder.output, params.stimulus{a}, params.intensity{b}))
        figure_counter = figure_counter + 1;

        % indentify P2 idx
        [P2_amplitude, idx] = max(data(visual.eoi, search_window));
        P2_idx = search_window(1) + idx;

        % prepare data for visualization
        visual.data = data(:, P2_idx)';

        % plot & save P2 topography
        figure(figure_counter);
        topoplot(double(visual.data), visual.chanlocs, 'maplimits', [-2 2], 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
        set(gca,'color',[1 1 1]);
        title(sprintf('%s, %s intensity - P2', params.stimulus{a}, params.intensity{b}))
        saveas(gcf, sprintf('%s\\figures\\preliminary_topo_%s_%s_P2.png', folder.output, params.stimulus{a}, params.intensity{b}))
        figure_counter = figure_counter + 1;
    end
end

% extract N2P2 peak amplitudes and latencies
for a = 1:length(params.stimulus)
    for b = 1:length(params.intensity)
        % extract single-subject measures
        for s = 1:params.subjects
            % select data
            statement = sprintf('data = squeeze(dataset.%s.%s(%d, visual.eoi, :))'';', params.stimulus{a}, params.intensity{b}, s);
            eval(statement)

            % extract N2 measures
            [N2_amplitude, idx] = min(data(search_window));
            N2_latency = params.search_window(1) + idx*header.xstep;

            % extract P2 measures
            [P2_amplitude, idx] = max(data(search_window));
            P2_latency = params.search_window(1) + idx*header.xstep;

            % encode to the output structure
            statement = sprintf('RFSxLASER_measures.individual(s).ERP.%s.%s.N2_latency = N2_latency;', params.stimulus{a}, params.intensity{b});
            eval(statement)
            statement = sprintf('RFSxLASER_measures.individual(s).ERP.%s.%s.N2_amplitude = N2_amplitude;', params.stimulus{a}, params.intensity{b});
            eval(statement)
            statement = sprintf('RFSxLASER_measures.individual(s).ERP.%s.%s.P2_latency = P2_latency;', params.stimulus{a}, params.intensity{b});
            eval(statement)
            statement = sprintf('RFSxLASER_measures.individual(s).ERP.%s.%s.P2_amplitude = P2_amplitude;', params.stimulus{a}, params.intensity{b});
            eval(statement)
        end
    end
end

% calculate group statistics
row_counter = 1;
for a = 1:length(params.stimulus)
    for b = 1:length(params.intensity)
        for c = 1:length(params.ERP_measures)
            % extract the data
            data = [];
            for s = 1:params.subjects
                statement = sprintf('data(s) = [RFSxLASER_measures.individual(s).ERP.%s.%s.%s];', params.stimulus{a}, params.intensity{b}, params.ERP_measures{c});
                eval(statement)                
            end
            data = data(subject_filter);

            % calculate the stats
            RFSxLASER_measures.group_stats.ERP(row_counter).stimulus = params.stimulus{a}; 
            RFSxLASER_measures.group_stats.ERP(row_counter).intensity = params.intensity{b};
            RFSxLASER_measures.group_stats.ERP(row_counter).measure = params.ERP_measures{c};
            RFSxLASER_measures.group_stats.ERP(row_counter).mean = mean(data);
            RFSxLASER_measures.group_stats.ERP(row_counter).min = min(data);
            RFSxLASER_measures.group_stats.ERP(row_counter).max = max(data);
            RFSxLASER_measures.group_stats.ERP(row_counter).SD = std(data);
            RFSxLASER_measures.group_stats.ERP(row_counter).SEM = std(data) / sqrt(sum(subject_filter)); 
            RFSxLASER_measures.group_stats.ERP(row_counter).CI_upper = mean(data) + t_value * RFSxLASER_measures.group_stats.ERP(row_counter).SEM; 
            RFSxLASER_measures.group_stats.ERP(row_counter).CI_lower = mean(data) - t_value * RFSxLASER_measures.group_stats.ERP(row_counter).SEM;

            % update row counter
            row_counter = row_counter + 1;
        end
    end
end

% plot boxplots for group statistics
for b = 1:length(params.intensity)
    for c = 1:length(params.ERP_measures)
        % extract the data
        data = [];
        for s = 1:params.subjects
            statement = sprintf('data(s) = [RFSxLASER_measures.individual(s).ERP.laser.%s.%s];', params.intensity{b}, params.ERP_measures{c});
            eval(statement)
            statement = sprintf('data(params.subjects + s) = [RFSxLASER_measures.individual(s).ERP.RFS.%s.%s];', params.intensity{b}, params.ERP_measures{c});
            eval(statement)
        end
        visual.data = data([subject_filter, subject_filter])';

        % plot the boxplot
        visual.group = [ones(sum(subject_filter), 1); 2*ones(sum(subject_filter), 1)];
        figure(figure_counter)
        if contains(params.ERP_measures{c}, 'latency')    
            if contains(params.ERP_measures{c}, 'N2')  
                plot_box(visual, 'orientation', 'horizontal', 'xlim', [0.1 0.25], 'xlabel', 'time (s)')
            elseif contains(params.ERP_measures{c}, 'P2')
                plot_box(visual, 'orientation', 'horizontal', 'xlim', [0.18 0.48], 'xlabel', 'time (s)')
            end
        elseif contains(params.ERP_measures{c}, 'amplitude')
            plot_box(visual, 'ylabel', 'amplitude (\muV)')
        end
        title(sprintf('%s stimulation intensity - %s',  params.intensity{b}, strrep(params.ERP_measures{c}, '_', ' ')))

        % save and update counter
        saveas(gcf, sprintf('%s\\figures\\preliminary_box_%s_%s.png', folder.output, params.intensity{b}, params.ERP_measures{c}))
        figure_counter = figure_counter + 1;
    end
end

% paired t-test for high intensity stimulaiton
for c = 1:length(params.ERP_measures)
    % extract laser data
    data = [];
    for s = 1:params.subjects
        statement = sprintf('data(s) = [RFSxLASER_measures.individual(s).ERP.laser.high.%s];', params.ERP_measures{c});
        eval(statement)                
    end
    data_all(1, :) = data(subject_filter);

    % extract laser data
    data = [];
    for s = 1:params.subjects
        statement = sprintf('data(s) = [RFSxLASER_measures.individual(s).ERP.RFS.high.%s];', params.ERP_measures{c});
        eval(statement)                
    end
    data_all(2, :) = data(subject_filter);

    % run the paired t-test
    fprintf('%s:\n', strrep(params.ERP_measures{c}, '_', ' '))
    output = paired_t(data_all);

    % append the output
    statement = sprintf('ttest.%s = output;', params.ERP_measures{c});
    eval(statement)
    fprintf('\n')
end
RFSxLASER_measures.group_stats.ERP(1).ttest = ttest;

% save and continue
save(output_file, 'RFSxLASER_measures', '-append')
clear a b c s dataset data2load header data visual t_value fig statement screen_size search_window ttest subject_filter ...
    idx N2_amplitude N2_latency P2_amplitude P2_latency row_counter N2_idx P2_idx group fig data_all output
fprintf('section 2 finished.\n\n')

%% 3) ratings group statistics
% ----- section input -----
params.stimulus = {'laser' 'RFS'};
params.intensity = {'high' 'low'};
params.outliers = [3, 16];
% -------------------------
fprintf('section 3: group statistics - ratings\n')

% define subject filter
params.subjects = length(RFSxLASER_info);
subject_filter = logical(ones(1, params.subjects)); 
subject_filter(params.outliers) = false;

% calculate group statistics
row_counter = 1;
t_value = tinv(0.975, sum(subject_filter) - 1); 
for a = 1:length(params.stimulus)
    for b = 1:length(params.intensity)
        % identify the rows
        row = [];
        for c = 1:length(RFSxLASER_measures.individual(1).ratings.stats)
            if strcmp(RFSxLASER_measures.individual(1).ratings.stats(c).stimulus, params.stimulus{a})
                if strcmp(RFSxLASER_measures.individual(1).ratings.stats(c).intensity, params.intensity{b})
                    row(end + 1) = c;
                end
            end
        end

        % extract the data
        data = [];
        for s = find(subject_filter)
            for r = 1:length(row)
                data_subj(r) = RFSxLASER_measures.individual(s).ratings.stats(row(r)).mean;
            end
            data(end + 1) = mean(data_subj);
        end

        % save for t-test
        data_ttest(a, b, :) = data;

        % calculate the stats
        RFSxLASER_measures.group_stats.ratings(row_counter).stimulus = params.stimulus{a}; 
        RFSxLASER_measures.group_stats.ratings(row_counter).intensity = params.intensity{b};
        RFSxLASER_measures.group_stats.ratings(row_counter).mean = mean(data);
        RFSxLASER_measures.group_stats.ratings(row_counter).min = min(data);
        RFSxLASER_measures.group_stats.ratings(row_counter).max = max(data);
        RFSxLASER_measures.group_stats.ratings(row_counter).SD = std(data);
        RFSxLASER_measures.group_stats.ratings(row_counter).SEM = std(data) / sqrt(sum(subject_filter)); 
        RFSxLASER_measures.group_stats.ratings(row_counter).CI_upper = mean(data) + t_value * RFSxLASER_measures.group_stats.ratings(row_counter).SEM; 
        RFSxLASER_measures.group_stats.ratings(row_counter).CI_lower = mean(data) - t_value * RFSxLASER_measures.group_stats.ratings(row_counter).SEM;

        % update row counter
        row_counter = row_counter + 1;
    end
end

% paired t-test for high intensity stimulaiton
for b = 1:length(params.intensity)
    % extract data
    data(1, :) = squeeze(data_ttest(1, b, :))';
    data(2, :) = squeeze(data_ttest(2, b, :))';

    % run the paired t-test
    fprintf('%s intensity:\n', params.intensity{b})
    output = paired_t(data);

    % append the output
    statement = sprintf('ttest.%s = output;', params.intensity{b});
    eval(statement)
    fprintf('\n')
end
RFSxLASER_measures.group_stats.ratings(1).ttest = ttest;

% save and continue
save(output_file, 'RFSxLASER_measures', '-append')
clear a b c s r subject_filter row_counter row data data_subj data_ttest output statement ttest t_value
fprintf('section 3 finished.\n\n')

%% 4) N1 ERP visualization and peak measure extraction
% ----- section input -----
params.eoi = 'C3';
params.prefix = 'N1_ready';
params.stimulus = {'laser' 'RFS'};
params.ERP_measures = {'N1_latency' 'N1_amplitude'};
params.search_window = [0.05 0.3];
params.outliers = [16];
% -------------------------
fprintf('section 4: ERP - N1 evaluation\n')

% select dataset
dataset = RFSxLASER_data.flipped_N1;

% define subject filter
params.subjects = length(RFSxLASER_info);
subject_filter = logical(ones(1, params.subjects)); 
subject_filter(params.outliers) = false;

% define common visualization & stats parameters
data2load = dir(sprintf('%s\\%s*%s*', folder.processed, params.prefix));
load(sprintf('%s\\%s', data2load(1).folder, data2load(1).name), '-mat')
t_value = tinv(0.975, sum(subject_filter) - 1); 
visual.x = (0:header.datasize(6)-1)*header.xstep + header.xstart;
visual.labels = params.stimulus;
visual.chanlocs = header.chanlocs;
visual.chanlabels = {header.chanlocs.labels};
visual.eoi = find(strcmp(visual.chanlabels, params.eoi));
visual.colors = [0.3333    0.4471    0.9020;
    1.0000    0.0745    0.6510];
screen_size = get(0, 'ScreenSize');

% visualization: select laser data
data = squeeze(dataset.laser.high(subject_filter, visual.eoi, :));
visual.data(1, :) = mean(data, 1);
visual.sem(1, :) = squeeze(std(data, 0, 1)) / sqrt(sum(subject_filter)); 
visual.CI_upper(1, :) = visual.data(1, :) + t_value * visual.sem(1, :); 
visual.CI_lower(1, :) = visual.data(1, :) - t_value * visual.sem(1, :);

% visualization: select RFS data
data = squeeze(dataset.RFS.high(subject_filter, visual.eoi, :));
visual.data(2, :) = mean(data, 1);
visual.sem(2, :) = squeeze(std(data, 0, 1)) / sqrt(sum(subject_filter)); 
visual.CI_upper(2, :) = visual.data(2, :) + t_value * visual.sem(2, :); 
visual.CI_lower(2, :) = visual.data(2, :) - t_value * visual.sem(2, :);

% launch the figure
fig = figure(figure_counter);    
set(fig, 'Position', [screen_size(3)/4, screen_size(4)/4, 2*screen_size(3)/5, screen_size(4) / 2])

% plot
plot_ERP(visual, 'colours', visual.colors, 'labels', visual.labels, 'xlim', [-0.1 0.6], 'ylim', [-2 2], 'reverse', 'on')
title('high stimulation intensity - N1')

% save figure and update counter
saveas(fig, sprintf('%s\\figures\\preliminary_N1.png', folder.output))
figure_counter = figure_counter + 1;

% plot topographies at group average peak latencies
search_window = find(visual.x >= params.search_window(1) & visual.x <= params.search_window(2));
addpath(genpath([folder.toolbox '\letswave 6']));
for a = 1:length(params.stimulus)
    % subset dataset
    statement = sprintf('data = squeeze(mean(dataset.%s.high(subject_filter, :, :), 1));', params.stimulus{a});
    eval(statement)

    % indentify N1 idx
    [N1_amplitude, idx] = min(data(visual.eoi, search_window));
    N1_idx = search_window(1) + idx;

    % prepare data for visualization
    visual.data = data(:, N1_idx)';

    % plot & save N1 topography
    figure(figure_counter);
    topoplot(double(visual.data), visual.chanlocs, 'maplimits', [-0.8 0.8], 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
    set(gca,'color',[1 1 1]);
    title(sprintf('%s, high intensity - N1', params.stimulus{a}))
    saveas(gcf, sprintf('%s\\figures\\preliminary_topo_%s_N1.png', folder.output, params.stimulus{a}))
    figure_counter = figure_counter + 1;
end

% extract N1 peak amplitudes and latencies
for a = 1:length(params.stimulus)
    % extract single-subject measures
    for s = 1:params.subjects
        % select data
        statement = sprintf('data = squeeze(dataset.%s.high(%d, visual.eoi, :))'';', params.stimulus{a}, s);
        eval(statement)

        % extract N2 measures
        [N1_amplitude, idx] = min(data(search_window));
        N1_latency = params.search_window(1) + idx*header.xstep;

        % encode to the output structure
        statement = sprintf('RFSxLASER_measures.individual(s).ERP.%s.high.N1_latency = N1_latency;', params.stimulus{a});
        eval(statement)
        statement = sprintf('RFSxLASER_measures.individual(s).ERP.%s.high.N1_amplitude = N1_amplitude;', params.stimulus{a});
        eval(statement)
    end
end

% calculate group statistics
row_counter = length(RFSxLASER_measures.group_stats.ERP) + 1;
for a = 1:length(params.stimulus)
    for c = 1:length(params.ERP_measures)
        % extract the data
        data = [];
        for s = 1:params.subjects
            statement = sprintf('data(s) = [RFSxLASER_measures.individual(s).ERP.%s.high.%s];', params.stimulus{a}, params.ERP_measures{c});
            eval(statement)                
        end
        data = data(subject_filter);

        % calculate the stats
        RFSxLASER_measures.group_stats.ERP(row_counter).stimulus = params.stimulus{a}; 
        RFSxLASER_measures.group_stats.ERP(row_counter).intensity = 'high';
        RFSxLASER_measures.group_stats.ERP(row_counter).measure = params.ERP_measures{c};
        RFSxLASER_measures.group_stats.ERP(row_counter).mean = mean(data);
        RFSxLASER_measures.group_stats.ERP(row_counter).min = min(data);
        RFSxLASER_measures.group_stats.ERP(row_counter).max = max(data);
        RFSxLASER_measures.group_stats.ERP(row_counter).SD = std(data);
        RFSxLASER_measures.group_stats.ERP(row_counter).SEM = std(data) / sqrt(sum(subject_filter)); 
        RFSxLASER_measures.group_stats.ERP(row_counter).CI_upper = mean(data) + t_value * RFSxLASER_measures.group_stats.ERP(row_counter).SEM; 
        RFSxLASER_measures.group_stats.ERP(row_counter).CI_lower = mean(data) - t_value * RFSxLASER_measures.group_stats.ERP(row_counter).SEM;

        % update row counter
        row_counter = row_counter + 1;
    end
end

% plot boxplots for group statistics
for c = 1:length(params.ERP_measures)
    % extract the data
    data = [];
    for s = 1:params.subjects
        statement = sprintf('data(s) = [RFSxLASER_measures.individual(s).ERP.laser.high.%s];', params.ERP_measures{c});
        eval(statement)
        statement = sprintf('data(params.subjects + s) = [RFSxLASER_measures.individual(s).ERP.RFS.high.%s];', params.ERP_measures{c});
        eval(statement)
    end
    visual.data = data([subject_filter, subject_filter])';

    % plot the boxplot  
    figure(figure_counter)
    visual.group = [ones(sum(subject_filter), 1); 2*ones(sum(subject_filter), 1)];
    if contains(params.ERP_measures{c}, 'latency')        
        plot_box(visual, 'orientation', 'horizontal', 'xlim', [0.02 0.27], 'xlabel', 'time (s)')
    elseif contains(params.ERP_measures{c}, 'amplitude')
        plot_box(visual, 'ylabel', 'amplitude (\muV)')
    end
    title(sprintf('high stimulation intensity - %s',  strrep(params.ERP_measures{c}, '_', ' ')))

    % save and update counter
    saveas(gcf, sprintf('%s\\figures\\preliminary_box_%s.png', folder.output, params.ERP_measures{c}))
    figure_counter = figure_counter + 1;
end

% paired t-test 
for c = 1:length(params.ERP_measures)
    % extract laser data
    data = [];
    for s = 1:params.subjects
        statement = sprintf('data(s) = [RFSxLASER_measures.individual(s).ERP.laser.high.%s];', params.ERP_measures{c});
        eval(statement)                
    end
    data_all(1, :) = data(subject_filter);

    % extract laser data
    data = [];
    for s = 1:params.subjects
        statement = sprintf('data(s) = [RFSxLASER_measures.individual(s).ERP.RFS.high.%s];', params.ERP_measures{c});
        eval(statement)                
    end
    data_all(2, :) = data(subject_filter);

    % run the paired t-test
    fprintf('%s:\n', strrep(params.ERP_measures{c}, '_', ' '))
    output = paired_t(data_all);

    % append the output
    statement = sprintf('ttest.%s = output;', params.ERP_measures{c});
    eval(statement)
    fprintf('\n')
end
RFSxLASER_measures.group_stats.ERP(1).ttest.N1_amplitude = ttest.N1_amplitude;
RFSxLASER_measures.group_stats.ERP(1).ttest.N1_latency = ttest.N1_latency;

% save and continue
save(output_file, 'RFSxLASER_measures', '-append')
clear a b c o s dataset data2load header data visual t_value fig statement screen_size search_window ttest subject_filter ...
    N1_idx idx N1_amplitude N1_latency group h ax fig data_all output outliers
fprintf('section 4 finished.\n\n')

%% 5) conduction velocity calculation
% ----- section input -----
params.stimulus = {'laser' 'RFS'};
params.outliers = [16];
% -------------------------
fprintf('section 5: conduction velocity\n')

% define subject filter
params.subjects = length(RFSxLASER_info);
subject_filter = logical(ones(1, params.subjects)); 
subject_filter(params.outliers) = false;

% extract conduction distance per subject
subject = struct([]);
fprintf('extracting individual metrics:\n')
for s = 1:params.subjects
    % extract individual measures
    subject(s).arm = ((RFSxLASER_info(s).body.arm_left + RFSxLASER_info(s).body.arm_right)/2)/100;
    subject(s).head = (RFSxLASER_info(s).body.height / 8)/100;
    subject(s).distance = subject(s).arm + subject(s).head;
    subject(s).N1.laser = RFSxLASER_measures.individual(s).ERP.laser.high.N1_latency;
    subject(s).N1.RFS = RFSxLASER_measures.individual(s).ERP.RFS.high.N1_latency;

    % calculate individual conduction velocity
    subject(s).velocity.laser = subject(s).distance / subject(s).N1.laser;
    subject(s).velocity.RFS = subject(s).distance / subject(s).N1.RFS;

    % append to the output structure
    RFSxLASER_measures.individual(s).CV = subject(s);
end
fprintf('mean subject height: %d +- %d m\n', mean([subject(subject_filter).head] * 8), std([subject(subject_filter).head] * 8))
fprintf('mean subject arm length: %d +- %d m\n', mean([subject(subject_filter).arm]), std([subject(subject_filter).arm]))
fprintf('mean subject conduction distance: %d +- %d m\n\n', mean([subject(subject_filter).distance]), std([subject(subject_filter).distance]))

% calculate group statistics
t_value = tinv(0.975, sum(subject_filter) - 1); 
for a = 1:length(params.stimulus)
    data = [];
    for s = 1:params.subjects
        statement = sprintf('data(s) = [RFSxLASER_measures.individual(s).CV.velocity.%s];', params.stimulus{a});
        eval(statement)                
    end
    data = data(subject_filter);

    % calculate the stats
    RFSxLASER_measures.group_stats.CV(a).stimulus = params.stimulus{a}; 
    RFSxLASER_measures.group_stats.CV(a).intensity = 'high';
    RFSxLASER_measures.group_stats.CV(a).mean = mean(data);
    RFSxLASER_measures.group_stats.CV(a).min = min(data);
    RFSxLASER_measures.group_stats.CV(a).max = max(data);
    RFSxLASER_measures.group_stats.CV(a).SD = std(data);
    RFSxLASER_measures.group_stats.CV(a).SEM = std(data) / sqrt(sum(subject_filter)); 
    RFSxLASER_measures.group_stats.CV(a).CI_upper = mean(data) + t_value * RFSxLASER_measures.group_stats.CV(a).SEM; 
    RFSxLASER_measures.group_stats.CV(a).CI_lower = mean(data) - t_value * RFSxLASER_measures.group_stats.CV(a).SEM;
end

% extract the data for plotting
data = [];
for s = 1:params.subjects
    data(s) = RFSxLASER_measures.individual(s).CV.velocity.laser;
    data(params.subjects + s) = RFSxLASER_measures.individual(s).CV.velocity.RFS;
end
visual.data = data([subject_filter, subject_filter])';

% plot the boxplot
visual.labels = params.stimulus;
visual.colors = [0.3333    0.4471    0.9020;
    1.0000    0.0745    0.6510];
figure(figure_counter)
plot_box(visual, 'orientation', 'horizontal', 'xlabel', 'conduction velocity (m/s)')
saveas(gcf, sprintf('%s\\figures\\preliminary_box_CV.png', folder.output))
figure_counter = figure_counter + 1;

% save and continue
save(output_file, 'RFSxLASER_measures', '-append')
clear a s subject_filter subject row_counter data visual statement t_value
fprintf('section 5 finished.\n\n')

%% 6) ERP magnitude vs. pain intensity
% ----- section input -----
params.stimulus = {'laser' 'RFS'};
params.ERP_measures = {'N1_amplitude' 'N2_amplitude' 'P2_amplitude'};
params.outliers = [7, 16];
% -------------------------
fprintf('section 6: ERP magnitude vs. pain intensity\n')
filter
% define subject 
params.subjects = length(RFSxLASER_info);
subject_filter = logical(ones(1, params.subjects)); 
subject_filter(params.outliers) = false;

% select data
row_counter = 1;
data = struct([]);
for a = 1:length(params.stimulus)
    for s = 1:params.subjects
        if subject_filter(s)
            % extract subject & session info
            data(row_counter).ID = RFSxLASER_info(s).ID;
            data(row_counter).stimulus = params.stimulus{a};
    
            % extract ERP variables
            statement = sprintf('data(row_counter).N1_amplitude = RFSxLASER_measures.individual(s).ERP.%s.high.N1_amplitude;', params.stimulus{a});
            eval(statement)
            statement = sprintf('data(row_counter).N2_amplitude = RFSxLASER_measures.individual(s).ERP.%s.high.N2_amplitude;', params.stimulus{a});
            eval(statement)
            statement = sprintf('data(row_counter).P2_amplitude = RFSxLASER_measures.individual(s).ERP.%s.high.P2_amplitude;', params.stimulus{a});
            eval(statement)
    
            % extract pain rating
            idx = logical(zeros(1, length(RFSxLASER_measures.individual(s).ratings.stats)));
            for b = 1:length(RFSxLASER_measures.individual(s).ratings.stats)
                if strcmp(RFSxLASER_measures.individual(s).ratings.stats(b).stimulus, params.stimulus{a})
                    if strcmp(RFSxLASER_measures.individual(s).ratings.stats(b).intensity, 'high')
                        idx(b) = true;
                    end
                end
            end
            data(row_counter).rating = mean([RFSxLASER_measures.individual(s).ratings.stats(idx).mean]);

            % update the counter
            row_counter = row_counter + 1;
        end
    end
end

% plot linear regression
colors = [0.3333    0.4471    0.9020;
    1.0000    0.0745    0.6510];
for a = 1:length(params.stimulus)
    for b = 1:length(params.ERP_measures)
        % subset the data
        idx = logical(zeros(1, length(data)));
        for c = 1:length(data)
            if strcmp(data(c).stimulus, params.stimulus{a})
                idx(c) = true;
            end
        end
        statement = sprintf('data_reg(1, :) = [data(idx).%s];', params.ERP_measures{b});
        eval(statement)
        data_reg(2, :) = [data(idx).rating];

        % run simple linear regression
        model = fitlm(data_reg(1, :), data_reg(2, :));
        disp(model);
        x_range = linspace(min(data_reg(1, :)), max(data_reg(1, :)), 100)';
        [y_fit, y_ci] = predict(model, x_range);

        % plot scatter plot with regression line
        fig = figure(figure_counter);
        scatter(data_reg(1, :), data_reg(2, :), 60, 'MarkerFaceColor', colors(a, :), 'MarkerEdgeColor', 'none'); 
        hold on;
        plot(x_range, y_fit, 'k', 'LineWidth', 2);  
        fill([x_range; flipud(x_range)], [y_ci(:,1); flipud(y_ci(:,2))], colors(a, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');  
        hold off;
        
        % formatting
        box off;
        ax = gca;
        ax.XAxisLocation = 'bottom';
        ax.YAxisLocation = 'left';
        ax.TickDir = 'out'; 
        set(gca, 'FontSize', 12)
        xlabel(sprintf('%s (µV)', strrep(params.ERP_measures{b}, '_', ' ')));
        ylabel('pain intensity');
        title(sprintf('%s: effect of %s', params.stimulus{a}, strrep(params.ERP_measures{b}, '_', ' ')));
        
        % save and update the counter
        saveas(fig, sprintf('%s\\figures\\preliminary_regression_%s, %s.png', folder.output, params.stimulus{a}, params.ERP_measures{b}))
        figure_counter = figure_counter + 1;
    end
end

% save and continue
clear a b s subject_filter row_counter data idx data_reg statement data_reg model coeffs fig colors ...
    x_range y_fit y_ci
fprintf('section 6 finished.\n\n')

%% functions
function plot_thresholds(visual, threshold)   
% =========================================================================
% plots the output of threshold tracking
% ========================================================================= 
% plot the line
plot(visual.x, visual.intensity, ':', 'color', [0.6510    0.6510    0.6510], 'linewidth', 3)
hold on

% plot the trials
for e = 1:length(visual.x)
    plot(visual.x(e), visual.intensity(e), 'o', 'MarkerSize', 12, 'MarkerEdgeColor','none','MarkerFaceColor', visual.color(e, :))
    hold on
end

% plot three dummy points (one for each color) for the legend
h(1) = plot(NaN, NaN, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.6510    0.6510    0.6510], 'MarkerEdgeColor', 'none');  
h(2) = plot(NaN, NaN, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0    0.4471    0.7412], 'MarkerEdgeColor', 'none');  
h(3) = plot(NaN, NaN, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.8000    0.0157    0.0157], 'MarkerEdgeColor', 'none');  
h(4) = plot(NaN, NaN, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none'); 

% identify spans of the axes
y_lim = get(gca, "YLim");
y_span = y_lim(2) - y_lim(1);
x_lim = get(gca, "XLim");
x_span = x_lim(2) - x_lim(1);

% plot thresholds
line(0:visual.x(end)+1, repelem(round(mean(threshold.perception), 2), 1, length(visual.x)+2), ...
    'color', [0    0.4471    0.7412], 'linewidth', 1.5)
line(0:visual.x(end)+1, repelem(round(mean(threshold.pain), 2), 1, length(visual.x)+2), ...
    'color', [0.8000    0.0157    0.0157], 'linewidth', 1.5)
line(0:visual.x(end)+1, repelem(round(mean(threshold.subjective), 2), 1, length(visual.x)+2), ...
    'color', [0.8000    0.0157    0.0157], 'linewidth', 1.5, 'linestyle', '--')

% plot annotations
text(0.01*x_span, round(mean(threshold.perception), 2) + 0.02*y_span, ...
    sprintf('perception: %.2f %s', round(mean(threshold.perception), 2), visual.unit), 'fontsize', 14);
text(0.01*x_span, round(mean(threshold.pain), 2) + 0.02*y_span, ...
    sprintf('pain: %.2f %s', round(mean(threshold.pain), 2), visual.unit), 'fontsize', 14);
text(0.01*x_span, round(mean(threshold.subjective), 2) + 0.02*y_span, ...
    sprintf('subjective pain: %.2f %s', round(mean(threshold.subjective), 2), visual.unit), 'fontsize', 14);

% adjust visuals
xlim([0, visual.x(end)+1]); 
ylim([0.9*min(visual.intensity), 1.1*max(visual.intensity)])
xlabel('trials'); 
ylabel(sprintf('stimulation intensity (%s)', visual.unit));
set(gca, 'fontsize', 14)

% add the legend
legend(h, {'not perceived', 'not painful', 'painful', 'electric'}, 'Location', 'northwest', 'Box', 'off');
end
function dataset = reload_dataset(data2load, subject_idx, fieldname)
% =========================================================================
% Reloads pre-processed EEG data of a single subject for following 
% processing steps. 
% Input:    - list of datasets to loads
%           - subject index
%           - fieldname
% =========================================================================  
% initiate output
dataset = struct;

% subset header and data files
header_idx = logical([]);
data_idx = logical([]);
for d = 1:length(data2load)
    if contains(data2load(d).name, 'lw6') 
        header_idx(d) = true;
        data_idx(d) = false;
    elseif contains(data2load(d).name, 'mat') 
        header_idx(d) = false;
        data_idx(d) = true;
    end
end
headers = data2load(header_idx);
datas = data2load(data_idx);

% load all dataset for this condition
if length(datas) == length(headers) 
    for d = 1:length(datas)
        % load header
        load(sprintf('%s\\%s', headers(d).folder, headers(d).name), '-mat')
        statement = sprintf('dataset(%d).%s(d).header = header;', subject_idx, fieldname);
        eval(statement) 

        % load data
        load(sprintf('%s\\%s', datas(d).folder, datas(d).name))
        statement = sprintf('dataset(%d).%s(d).data = data;', subject_idx, fieldname);
        eval(statement) 
    end
else
    error('ERROR: Number of data files (.mat) does not match the number of header files (.lw6)! Check manually.')
end
end
function wait4files(filenames)
% =========================================================================
% waitForFiles pauses the script until all required files appear in the
% working working directory
% --> file names are specified in a cell array
% =========================================================================    
% loop to wait for files
while true
    allFilesExist = true;
    
    % check for each file
    for i = 1:length(filenames)
        if isempty(dir(filenames{i}))
            allFilesExist = false;
            break;
        end
    end
    
    % if all files exist, break the loop
    if allFilesExist
        break;
    end
    
    % pause for 2s
    pause(2);
end
end
function visual = create_visual(data, header, varargin)
% =========================================================================
% prepares a structure with data ready for ERP plotting
% --> data must be in matrix format: condition * subject * timepoint
% =========================================================================  
% initiate output
visual = struct;

% extract header parameters
visual.xstart = header.xstart;
visual.xstep = header.xstep;

% prepare original x 
visual.x = header.xstart : header.xstep : (header.datasize(6) * header.xstep)+header.xstart-header.xstep;

% loop through conditions
for c = 1:size(data, 1)
    % extract mean data and statistics
    visual.y(c, :) = squeeze(mean(data(c, :, :), 2));

    % calculate statistics
    % visual.t = tinv(0.975, size(data, 2) - 1); 
    visual.SD(c, :) = std(squeeze(data(c, :, :)));
    visual.SEM(c, :) = visual.SD(c, :) / sqrt(size(visual.y, 2)); 
    % visual.CI(c, :) = visual.t * visual.SEM(c, :); 
end
end
function plot_ERP(input, varargin)
% =========================================================================
% plots an event-related potential
% input = structure with fields:    
%           data --> condition/electrode * sample
%           x --> vector with time samples
%           CI_upper --> condition/electrode * sample
%           CI_lower --> condition/electrode * sample
% varargins = name-value pairs: 
%           xlim --> 2-element vector (min, max)     
%           ylim --> 2-element vector (min, max) 
%           colours --> n*3 matrix of RGB values
%           shading --> 'on'(default)/'off'
%           alpha --> a float (default 0.2)           
%           legend --> 'on'(default)/'off'
%           labels --> cell array with labels for the legend  
%           legend_loc --> legend location (default 'southeast')
%           eoi --> label of a channel to be highlighted
%           reverse --> 'on'/'off'(default) - flips y axis
%           interpolated --> time window that was interpolated
% =========================================================================  
% set defaults
x_limits = [0,0];
y_limits = [0,0];
col = prism(size(input.data, 1));
shading = true;
alpha = 0.2;
plot_legend = true;
for c = 1:size(input.data, 1)
    labels{c} = sprintf('condition %d', c);
end
legend_loc = 'southeast';
highlight = false;
reverse = false;
interpolate = false;

% check for varargins
if ~isempty(varargin)
    % x limits
    a = find(strcmpi(varargin, 'xlim'));
    if ~isempty(a)
        x_limits = varargin{a + 1};
    end

    % y limits
    b = find(strcmpi(varargin, 'ylim'));
    if ~isempty(b)
        y_limits = varargin{b + 1};
    end

    % colours
    c = find(strcmpi(varargin, 'colours'));
    if ~isempty(c)
        col = varargin{c + 1};
    end

    % shading - default on
    d = find(strcmpi(varargin, 'shading'));
    if ~isempty(d) && strcmp(varargin{d + 1}, 'off')
        shading = false;
    end

    % alpha
    e = find(strcmpi(varargin, 'alpha'));
    if ~isempty(e)
        alpha = varargin{e + 1};
    end

    % legend - default on
    f = find(strcmpi(varargin, 'legend'));
    if ~isempty(f) && strcmp(varargin{f + 1}, 'off')
        plot_legend = false;
    end    

    % labels
    g = find(strcmpi(varargin, 'labels'));
    if ~isempty(g)
        labels = varargin{g + 1};
    end

    % legend location
    h = find(strcmpi(varargin, 'legend_loc'));
    if ~isempty(h) 
        legend_loc = varargin{h + 1};
    end  

    % highlighted channel - default off
    i = find(strcmpi(varargin, 'eoi'));
    if ~isempty(i)
        eoi = varargin{i + 1};
        eoi_n = find(contains(input.chanlabels, eoi));
        highlight = true;
    end 

    % interpolated interval - default off
    j = find(strcmpi(varargin, 'interpolated'));
    if ~isempty(j)
        interpolate_toi = varargin{j + 1};
        interpolate = true;
    end 

    % reverse y axis - default off
    r = find(strcmpi(varargin, 'reverse'));
    if ~isempty(r) && strcmp(varargin{r + 1}, 'on')
        reverse = true;
    end
end

% loop through datasets to plot
for t = 1:size(input.data, 1) 
    P(t) = plot(input.x, input.data(t, :), 'Color', col(t, :), 'LineWidth', 2);
    hold on
    if shading
        F(t) = fill([input.x fliplr(input.x)],[input.CI_upper(t, :) fliplr(input.CI_lower(t, :))], ...
            col(t, :), 'FaceAlpha', alpha, 'linestyle', 'none');
        hold on
    end
end

% check y limits
if y_limits(1) == 0 && y_limits(2) == 0
    y_limits = ylim;
end

% highlight channel if required
if highlight
    P(end + 1) = plot(input.x, input.data(eoi_n, :), 'Color', [0.9216    0.1490    0.1490], 'LineWidth', 4);
    hold on
end

% shade interpolated window if required
if interpolate
    interpolate_x = [interpolate_toi(1), interpolate_toi(2), interpolate_toi(2), interpolate_toi(1)];
    interpolate_y = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
    fill(interpolate_x, interpolate_y, [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% plot stimulus
line([0, 0], y_limits, 'Color', 'black', 'LineWidth', 2.5, 'LineStyle', '--')

% plot legend if required
if plot_legend 
    legend(P, labels, 'Location', legend_loc, 'fontsize', 14)
    legend('boxoff');
else
    legend('off')
end

% axes
box off;
ax = gca;
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.TickDir = 'out'; 
ax.XColor = [0.5020    0.5020    0.5020]; 
ax.YColor = [0.5020    0.5020    0.5020]; 

% set x limits 
if x_limits(1) == 0 && x_limits(2) == 0
    xlim([input.x(1), input.x(end)]) 
else
    xlim(x_limits)
end

% referse y axis if required
if reverse
    set(gca, 'YDir', 'reverse');
end

% other parameters
xlabel('time (s)')
ylabel('amplitude (\muV)')
set(gca, 'FontSize', 14)
ylim(y_limits)
set(gca, 'Layer', 'Top')
end
function output = paired_t(data)
% =========================================================================
% calculates a simple paired t-test, outputs statistics and prints them
% data  =   MxN matrix
%           M --> n. of conditions (=2)
%           N --> n. of observations
% ========================================================================= 
% extract conditions
condition1 = data(1, :);
condition2 = data(2, :);

% run the t-test
[output.h, output.p, output.ci, output.stats] = ttest(condition1, condition2);

% display results
fprintf('t-statistic: %.3f\n', output.stats.tstat);
fprintf('p-value: %.5f\n', output.p);
fprintf('degrees of freedom: %d\n', output.stats.df);
fprintf('confidence interval: [%.3f, %.3f]\n', output.ci(1), output.ci(2));

% check significance
if output.h == 1
    disp('--> result: significant difference (p < 0.05)');
else
    disp('--> result: no significant difference (p >= 0.05)');
end
end
function plot_box(visual, varargin)
% =========================================================================
% plots an event-related potential
% visual = structure with fields:    
%           data --> n*2 matrix, n = n. of observations
%           group --> grouping vector
%           labels --> cell array with group latebls (optional field)
%           colors --> 2*3 matrix with RGB color codes (optional field)
% varargin = name-value pairs: 
%           orientation --> horizontal/vertical (default)  
%           xlim --> x axis limits 
%           xlabel --> x axis label 
%           ylim --> y axis limits 
%           ylabel --> y axis label 
% =========================================================================  
% set defaults
colors = prism(length(unique(visual.group)));
horizontal = false;

% check for optional plotting parameters
fields = fieldnames(visual);
if ~isempty(varargin)
    % oritentation
    a = find(strcmpi(varargin, 'orientation'));
    if ~isempty(a) && strcmp(varargin{a + 1}, 'horizontal')
        horizontal = true;
    end
end
if sum(ismember(fields, 'colors')) > 0
    colors = visual.colors;
end

% plot boxplots
if horizontal
    boxplot(visual.data, visual.group, 'Labels', visual.labels, 'colors', colors, 'orientation', 'horizontal','symbol', 'o');
else
    boxplot(visual.data, visual.group, 'Labels', visual.labels, 'colors', colors, 'symbol', 'o');
end

% recolor outliers
outliers = findobj(gca, 'Tag', 'Outliers');
outliers = outliers([2,1]);
for o = 1:length(outliers)
    set(outliers(o), 'Marker', 'o', ...  
                           'MarkerFaceColor', visual.colors(o, :), ... 
                           'MarkerEdgeColor', visual.colors(o, :), ...  
                           'MarkerSize', 6); 
end

% add lebels
if ~isempty(varargin)
    % x limits
    a = find(strcmpi(varargin, 'xlim'));
    if ~isempty(a) 
        xlim(varargin{a + 1});
    end

    % x label
    b = find(strcmpi(varargin, 'xlabel'));
    if ~isempty(b) 
        xlabel(varargin{b + 1});
    end
        
    % y limits
    c = find(strcmpi(varargin, 'ylim'));
    if ~isempty(c) 
        ylim(varargin{c + 1});
    end

    % y label
    d = find(strcmpi(varargin, 'ylabel'));
    if ~isempty(d) 
        ylabel(varargin{d + 1});
    end
end

% adjust visuals
h = findobj(gca, 'Type', 'Line');
set(h, 'LineWidth', 1.5); 
box off;
ax = gca;
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.TickDir = 'out'; 
set(gca, 'FontSize', 12)
end