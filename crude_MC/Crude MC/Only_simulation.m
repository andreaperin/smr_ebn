clear; clc; close all;

load("Input_single_reactor_0.mat")

% === Define fixed parameters (same for all simulations) ===
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

% === Assign fixed variables to all SimulationInput objects ===
vars_to_set = who;  % all variables in workspace
% Remove ones you don’t want to set (e.g., 'in_filtered', 'out', etc.)
vars_to_set = setdiff(vars_to_set, {'in_filtered','out','vars_to_set','ans'});

for i = 1:length(in_filtered)
    for v = 1:length(vars_to_set)
        val = eval(vars_to_set{v});
        if length(val)==length(in_filtered)
            in_filtered(i) = setVariable(in_filtered(i), vars_to_set{v}, val(i));
        else
            in_filtered(i) = setVariable(in_filtered(i), vars_to_set{v}, val);
        end
    end
end

% === Run simulations in parallel ===
out = parsim(in_filtered);

% === Post-processing: count failures ===
fail = 0;
for i = 1:length(out)
    if max(max(out(i).T_W1)) > 1244
        fail = fail + 1;
    end
end

fprintf('Total failed simulations: %d out of %d\n', fail, length(out));
FP = fail/N;
fprintf('Failure probability: %d\n', FP);

save(strcat('Results_single_reactor_hydrogen_',num2str(round(age_considered))))
