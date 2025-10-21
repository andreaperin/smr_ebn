function [T_W1, T_W2, T_W3, T_W4] = Simulation_model_hydrogen( ...
    LOCA_time, ...
    LOOP_time, ...
    LHS_time, ...
    MSLB_time, ...
    MSLBH2_time, ...
    LOOPH2_time, ...
    ACS_time, ...
    ACS_rtime, ...
    EDG_time, ...
    EDG_rtime, ...
    pdp_time, ...
    pdp_rtime)

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
assignin('base', 'LOCA1_time', LOCA_time);
assignin('base', 'ACS_1', ACS_time);
assignin('base', 'Power', min(LOOPH2_time,LOOP_time));
assignin('base', 'MSLB1', min(MSLB_time, MSLBH2_time));
assignin('base', 'ACS1_response_time', ACS_rtime);
assignin('base', 'EDG_1', EDG_time);
assignin('base', 'EDG1_response_time', EDG_rtime);
assignin('base', 'PDP11', pdp_time);
assignin('base', 'PDP11_response_time', pdp_rtime);
assignin('base', 'LHS', LHS_time);
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
in = in.setVariable('LOCA1_time', LOCA_time);
in = in.setVariable('ACS_1', ACS_time);
in = in.setVariable('Power', min(LOOPH2_time,LOOP_time));
in = in.setVariable('MSLB1', min(MSLB_time, MSLBH2_time));
in = in.setVariable('ACS1_response_time', ACS_rtime);
in = in.setVariable('EDG_1', EDG_time);
in = in.setVariable('EDG1_response_time', EDG_rtime);
in = in.setVariable('PDP11', pdp_time);
in = in.setVariable('PDP11_response_time', pdp_rtime);
in = in.setVariable('LHS', LHS_time);
% Set remaining variables from base workspace
varlist = {'LOCA2_time','LOCA3_time','LOCA4_time','ACS_2','ACS_3','ACS_4', ...
           'EDG_2','PDP12','PDP21','PDP22','PDP31','PDP32', ...
           'PDP41','PDP42','MSLB2','MSLB3','MSLB4','thermal_failure', ...
           'thermal_failure_time','PGA',"ACS2_response_time", ...
           "ACS3_response_time","ACS4_response_time","EDG2_response_time",...
           "PDP12_response_time","PDP21_response_time","PDP22_response_time",...
           "PDP31_response_time","PDP32_response_time","PDP41_response_time","PDP42_response_time" };
for k = 1:numel(varlist)
    if evalin('base', sprintf('exist(''%s'', ''var'')', varlist{k}))
        in = in.setVariable(varlist{k}, evalin('base', varlist{k}));
    end
end

in = in.setVariable('ACS1_flow', ACS1_flow * alpha_1);
in = in.setVariable('ACS2_flow', ACS2_flow * alpha_2);
in = in.setVariable('ACS3_flow', ACS3_flow * alpha_3);
in = in.setVariable('ACS4_flow', ACS4_flow * alpha_4);

vars = in.Variables;
fprintf('Simulation inputs:\n');
for idx = 1:numel(vars)
    varName = vars(idx).Name;
    varValue = vars(idx).Value;
    if isnumeric(varValue) || islogical(varValue)
        valueStr = mat2str(varValue);
    elseif isstring(varValue)
        valueStr = strjoin(cellstr(varValue), ', ');
    elseif ischar(varValue)
        valueStr = varValue;
    else
        valueStr = strtrim(evalc('disp(varValue)'));
    end
    fprintf('  %s: %s\n', varName, valueStr);
end


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
