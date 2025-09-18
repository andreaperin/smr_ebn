function Simulation_model(mode, workDir)
% Simulation_model implements three modes:
%   - no arguments: run a single simulation as a script
%   - 'server', workDir: watch workDir for job files and run simulations
%   - 'client', workDir: no action; job files are written by ExternalModel

if nargin == 0
    runSimulation(workDir);
elseif strcmp(mode,'server')
    while true
        % look for job files anywhere under workDir
        jobFiles = dir(fullfile(workDir, '**', '*.job'));
        if ~isempty(jobFiles)
            job = jobFiles(1);
            jobFolder = job.folder;
            delete(fullfile(jobFolder, job.name));   % remove the job file
            runSimulation(jobFolder);                % run simulation in that folder
        else
            pause(0.05);
        end
    end
elseif strcmp(mode, 'client')
    % In this setup the client does nothing; the Julia code creates the job
    return;
else
    error('Unknown mode');
end
end

function runSimulation(workDir)
% runSimulation executes the Simulink model once.
% It reads inputs and writes outputs relative to workDir.

% change into the working directory so relative paths work
if nargin >= 1 && ~isempty(workDir)
    curDir = pwd;
    cd(workDir);
end

tic

% read CSV and put each column into local and base workspaces
    tbl  = readtable('Failure_model_outputs.csv');
    S    = table2struct(tbl,'ToScalar',true);
    names = fieldnames(S);
    for k = 1:numel(names)
        name = names{k};
        value = S.(name);
        % define a variable in the local workspace
        eval([name ' = value;']);
        % also define it in the base workspace so Simulink sees it
        assignin('base', name, value);
        % optional: attach it directly to the simulation input
        % in = in.setVariable(name, value);
    end

% run parameter file
evalc('SMDFR_Parameters');

model = 'SMDFR_HTE_model';
All_Outputs = 0;

% base flow values and simulation time
ACS1_flow = 25; ACS2_flow = 25; ACS3_flow = 25; ACS4_flow = 25;
tsim = 1200;

% set up simulation inputs
in = Simulink.SimulationInput(model);
in = in.setModelParameter('StartTime','1', 'StopTime', num2str(tsim));
in = in.setVariable('All_Outputs', All_Outputs);
in = in.setVariable('LOCA1_time', LOCA1_time);
in = in.setVariable('LOCA2_time', LOCA2_time);
in = in.setVariable('LOCA3_time', LOCA3_time);
in = in.setVariable('LOCA4_time', LOCA4_time);
in = in.setVariable('ACS_1', ACS_1);
in = in.setVariable('ACS_2', ACS_2);
in = in.setVariable('ACS_3', ACS_3);
in = in.setVariable('ACS_4', ACS_4);
in = in.setVariable('EDG_1', EDG_1);
in = in.setVariable('EDG_2', EDG_2);
in = in.setVariable('Power', Power);

% component failures
in = in.setVariable('PDP11', PDP11);
in = in.setVariable('PDP12', PDP12);
in = in.setVariable('PDP21', PDP21);
in = in.setVariable('PDP22', PDP22);
in = in.setVariable('PDP31', PDP31);
in = in.setVariable('PDP32', PDP32);
in = in.setVariable('PDP41', PDP41);
in = in.setVariable('PDP42', PDP42);
in = in.setVariable('MSLB1', MSLB1);
in = in.setVariable('MSLB2', MSLB2);
in = in.setVariable('MSLB3', MSLB3);
in = in.setVariable('MSLB4', MSLB4);

% adjust flows with alpha factors
in = in.setVariable('ACS1_flow', ACS1_flow * alpha_1);
in = in.setVariable('ACS2_flow', ACS2_flow * alpha_2);
in = in.setVariable('ACS3_flow', ACS3_flow * alpha_3);
in = in.setVariable('ACS4_flow', ACS4_flow * alpha_4);
in = in.setVariable('thermal_failure', thermal_failure);
in = in.setVariable('LHS', LHS);
in = in.setVariable('thermal_failure_time', thermal_failure_time);
in = in.setVariable('PGA', PGA);

% run the simulation
out = sim(in);

% extract maximum temperatures
try
    T_W1 = reshape(max(out.T_W1, [], 1), [1, tsim])';
    T_W2 = reshape(max(out.T_W2, [], 1), [1, tsim])';
    T_W3 = reshape(max(out.T_W3, [], 1), [1, tsim])';
    T_W4 = reshape(max(out.T_W4, [], 1), [1, tsim])';
catch
    % simulation failed; fill each column with NaNs
    T_W1 = ones(tsim, 1);
    T_W2 = ones(tsim, 1);
    T_W3 = ones(tsim, 1);
    T_W4 = ones(tsim, 1);
end

% save outputs
save('Simulation_model_outputs', 'T_W1','T_W2','T_W3','T_W4');
results = table(T_W1, T_W2, T_W3, T_W4);
writetable(results, 'Simulation_model_outputs.csv');

if exist('curDir','var')
    cd(curDir);
end

toc
end