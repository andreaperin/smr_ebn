function [T_W1, T_W2, T_W3, T_W4] = Simulation_model_hydrogen_parallel( ...
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

% Load parameter script
evalc('SMDFR_Parameters');

% === Model configuration ===
model = 'SMDFR';
load_system(model);
set_param(model, 'UnconnectedOutputMsg', 'none');

% Defaults (used if not found in base workspace or not provided as vectors)
All_Outputs = 0;
ACS1_flow = 25; ACS2_flow = 25; ACS3_flow = 25; ACS4_flow = 25;
tsim = 1200;
alpha_1 = 1; alpha_2 = 1; alpha_3 = 1; alpha_4 = 1;

% Additional params (can be overridden by base workspace values below)
ACS_2 = 1; ACS_3 = 1; ACS_4 = 1;
LOCA2_time = 1200; LOCA3_time = 1200; LOCA4_time = 1200;
EDG_2 = 1;
PDP12 = 0; PDP21 = 0; PDP22 = 0; PDP31 = 0; PDP32 = 0; PDP41 = 0; PDP42 = 0;
MSLB2 = 0; MSLB3 = 0; MSLB4 = 0;
thermal_failure = 0; thermal_failure_time = 1200; PGA = 0;
ACS2_response_time = 0; ACS3_response_time = 0; ACS4_response_time = 0;
EDG2_response_time = 0;
PDP12_response_time = 0; PDP21_response_time = 0; PDP22_response_time = 0;
PDP31_response_time = 0; PDP32_response_time = 0; PDP41_response_time = 0; PDP42_response_time = 0;

N = numel(LOCA_time);

% List of the "missing inputs" you wanted to also push into each in(i)
extraVarList = { ...
    'LOCA2_time','LOCA3_time','LOCA4_time', ...
    'ACS_2','ACS_3','ACS_4', ...
    'EDG_2', ...
    'PDP12','PDP21','PDP22','PDP31','PDP32','PDP41','PDP42', ...
    'MSLB2','MSLB3','MSLB4', ...
    'thermal_failure','thermal_failure_time','PGA', ...
    'ACS2_response_time','ACS3_response_time','ACS4_response_time', ...
    'EDG2_response_time', ...
    'PDP12_response_time','PDP21_response_time','PDP22_response_time', ...
    'PDP31_response_time','PDP32_response_time','PDP41_response_time','PDP42_response_time' ...
};

% Build SimulationInput array
in(1,N) = Simulink.SimulationInput(model);
for i = 1:N
    in(i) = Simulink.SimulationInput(model);
    in(i) = in(i).setModelParameter('StartTime', '1', 'StopTime', num2str(tsim));

    % Function-argument inputs (index if vectors)
    in(i) = in(i).setVariable('All_Outputs', All_Outputs);
    in(i) = in(i).setVariable('LOCA1_time',    idxOrScalar(LOCA_time,    i));
    in(i) = in(i).setVariable('ACS_1',         idxOrScalar(ACS_time,     i));
    in(i) = in(i).setVariable('Power',         min(idxOrScalar(LOOPH2_time,i), idxOrScalar(LOOP_time,i)));
    in(i) = in(i).setVariable('MSLB1',         min(idxOrScalar(MSLB_time, i), idxOrScalar(MSLBH2_time,i)));
    in(i) = in(i).setVariable('ACS1_response_time', idxOrScalar(ACS_rtime, i));
    in(i) = in(i).setVariable('EDG_1',             idxOrScalar(EDG_time,  i));
    in(i) = in(i).setVariable('EDG1_response_time', idxOrScalar(EDG_rtime, i));
    in(i) = in(i).setVariable('PDP11',              idxOrScalar(pdp_time,  i));
    in(i) = in(i).setVariable('PDP11_response_time',idxOrScalar(pdp_rtime, i));
    in(i) = in(i).setVariable('LHS',                idxOrScalar(LHS_time,  i));

    % Flows with alphas (index alphas if you later make them vectors)
    in(i) = in(i).setVariable('ACS1_flow', ACS1_flow * idxOrScalar(alpha_1,i));
    in(i) = in(i).setVariable('ACS2_flow', ACS2_flow * idxOrScalar(alpha_2,i));
    in(i) = in(i).setVariable('ACS3_flow', ACS3_flow * idxOrScalar(alpha_3,i));
    in(i) = in(i).setVariable('ACS4_flow', ACS4_flow * idxOrScalar(alpha_4,i));

    % Push all the "missing" variables too:
    for k = 1:numel(extraVarList)
        vn = extraVarList{k};
        val = preferLocalElseBase(vn);       % get scalar or vector
        in(i) = in(i).setVariable(vn, idxOrScalar(val, i));  % index if vector
    end
end

% Decide which cases need an actual simulation (rest get T_W1 fixed at 1100)
conditionThreshold = 1200;
simulateMask = false(1, N);
for i = 1:N
    simulateMask(i) = (LOCA_time(i) < conditionThreshold) || ...
                      (MSLB_time(i) < conditionThreshold) || ...
                      (LHS_time(i) < conditionThreshold);
end

% Run only the selected cases in parallel
save_system(model);
out = [];
if any(simulateMask)
    out = parsim(in(simulateMask), 'ShowProgress', 'on');
end

% === Extract temperature results (index outputs properly) ===
T_W1 = 1100 * ones(N, tsim);  % default constant for non-simulated cases
T_W2 = NaN(N, tsim); T_W3 = NaN(N, tsim); T_W4 = NaN(N, tsim);
outIdx = 1;
for i = 1:N
    if ~simulateMask(i)
        continue
    end
    try
        T_W1(i,:) = reshape(max(out(outIdx).T_W1, [], 1), [1, tsim])';
        T_W2(i,:) = reshape(max(out(outIdx).T_W2, [], 1), [1, tsim])';
        T_W3(i,:) = reshape(max(out(outIdx).T_W3, [], 1), [1, tsim])';
        T_W4(i,:) = reshape(max(out(outIdx).T_W4, [], 1), [1, tsim])';
    catch
        warning('Temperature extraction failed for run %d; filling with NaNs.', i);
        % leave NaNs as initialized
    end
    outIdx = outIdx + 1;
end

toc
end

% ----- helpers -----
function v = idxOrScalar(x, i)
% If x is a vector with length >= i, return x(i); otherwise return x unchanged.
    if isnumeric(x) || islogical(x)
        if ~isscalar(x) && numel(x) >= i
            v = x(i);
        else
            v = x;
        end
    else
        % pass-through for non-numeric types (strings, etc.)
        try
            if ~isscalar(x) && numel(x) >= i
                v = x(i);
            else
                v = x;
            end
        catch
            v = x;
        end
    end
end

function val = preferLocalElseBase(varName)
% Return local variable value if it exists in caller; else take from base workspace;
% if neither exists, create an empty [] so setVariable still succeeds.
    if evalin('caller', sprintf('exist(''%s'',''var'')', varName))
        val = evalin('caller', varName);
    elseif evalin('base', sprintf('exist(''%s'',''var'')', varName))
        val = evalin('base', varName);
    else
        val = [];
    end
end
 
