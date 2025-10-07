function [T_W1, T_W2, T_W3, T_W4] = Simulation_model(LOCA1_time_no_csv)

tic

tbl  = readtable('Failure_model_outputs.csv');
S    = table2struct(tbl,'ToScalar',true);
names = fieldnames(S);

for k = 1:numel(names)
    name = names{k};
    value = S.(name);

    if iscell(value), value = cell2mat(value); end
    if ischar(value), value = str2double(value); end

    assignin('base', name, value);
end

% Override LOCA1_time with user input
assignin('base', 'LOCA1_time', LOCA1_time_no_csv);

% Load parameter script
evalc('SMDFR_Parameters');

% === Model configuration ===
model = 'SMDFR_HTE_model';
load_system(model);
set_param(model, 'UnconnectedOutputMsg', 'none');

All_Outputs = 0;
ACS1_flow = 25; ACS2_flow = 25; ACS3_flow = 25; ACS4_flow = 25;
tsim = 1200;
alpha_1 = 1; alpha_2 = 1; alpha_3 = 1; alpha_4 = 1;

% === Build simulation input ===
in = Simulink.SimulationInput(model);
in = in.setModelParameter('StartTime', '1', 'StopTime', num2str(tsim));

in = in.setVariable('All_Outputs', All_Outputs);
in = in.setVariable('LOCA1_time', LOCA1_time_no_csv);

% Set remaining variables from base workspace
varlist = {'LOCA2_time','LOCA3_time','LOCA4_time','ACS_1','ACS_2','ACS_3','ACS_4', ...
           'EDG_1','EDG_2','Power','PDP11','PDP12','PDP21','PDP22','PDP31','PDP32', ...
           'PDP41','PDP42','MSLB1','MSLB2','MSLB3','MSLB4','thermal_failure','LHS', ...
           'thermal_failure_time','PGA'};
for k = 1:numel(varlist)
    if evalin('base', sprintf('exist(''%s'', ''var'')', varlist{k}))
        in = in.setVariable(varlist{k}, evalin('base', varlist{k}));
    end
end

in = in.setVariable('ACS1_flow', ACS1_flow * alpha_1);
in = in.setVariable('ACS2_flow', ACS2_flow * alpha_2);
in = in.setVariable('ACS3_flow', ACS3_flow * alpha_3);
in = in.setVariable('ACS4_flow', ACS4_flow * alpha_4);

out = sim(in);

% === Extract temperature results ===
try
    T_W1 = reshape(max(out.T_W1, [], 1), [1, tsim])';
    T_W2 = reshape(max(out.T_W2, [], 1), [1, tsim])';
    T_W3 = reshape(max(out.T_W3, [], 1), [1, tsim])';
    T_W4 = reshape(max(out.T_W4, [], 1), [1, tsim])';
catch
    warning('Temperature extraction failed; filling with NaNs.');
    T_W1 = NaN(tsim, 1);
    T_W2 = NaN(tsim, 1);
    T_W3 = NaN(tsim, 1);
    T_W4 = NaN(tsim, 1);
end

toc
end
