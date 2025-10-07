function [T_W1, T_W2, T_W3, T_W4] = Simulation_model( ...
    t_LOCA1, t_LOCA2, t_LOCA3, t_LOCA4, ...
    t_ACS1,  t_ACS2,  t_ACS3,  t_ACS4,  ...
    tr_ACS1, tr_ACS2, tr_ACS3, tr_ACS4, ...
    t_EDG1,  t_EDG2,  ...
    tr_EDG1, tr_EDG2, ...
    t_LOOP)

    tic

    % === Load defaults from CSV to base workspace (kept as in original) ===
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

    % === Load parameter script ===
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

    % Generic controls
    in = in.setVariable('All_Outputs', All_Outputs);

    % === Map EBN time nodes to Simulink workspace variables ===
    % LOCA times
    in = in.setVariable('LOCA1_time', t_LOCA1);
    in = in.setVariable('LOCA2_time', t_LOCA2);
    in = in.setVariable('LOCA3_time', t_LOCA3);
    in = in.setVariable('LOCA4_time', t_LOCA4);

    % ACS actuation times (note: model variables are named ACS_1..ACS_4)
    in = in.setVariable('ACS_1', t_ACS1);
    in = in.setVariable('ACS_2', t_ACS2);
    in = in.setVariable('ACS_3', t_ACS3);
    in = in.setVariable('ACS_4', t_ACS4);

    % ACS recovery/response times
    in = in.setVariable('tr_ACS1', tr_ACS1);
    in = in.setVariable('tr_ACS2', tr_ACS2);
    in = in.setVariable('tr_ACS3', tr_ACS3);
    in = in.setVariable('tr_ACS4', tr_ACS4);

    % EDG start/availability and transients
    in = in.setVariable('t_EDG1', t_EDG1);
    in = in.setVariable('t_EDG2', t_EDG2);
    in = in.setVariable('tr_EDG1', tr_EDG1);
    in = in.setVariable('tr_EDG2', tr_EDG2);

    % LOOP time
    in = in.setVariable('t_LOOP', t_LOOP);

    % === Keep remaining variables from base workspace if present ===
    varlist = {'EDG_1','EDG_2','Power','PDP11','PDP12','PDP21','PDP22','PDP31','PDP32', ...
               'PDP41','PDP42','MSLB1','MSLB2','MSLB3','MSLB4','thermal_failure','LHS', ...
               'thermal_failure_time','PGA'};
    for k = 1:numel(varlist)
        if evalin('base', sprintf('exist(''%s'', ''var'')', varlist{k}))
            in = in.setVariable(varlist{k}, evalin('base', varlist{k}));
        end
    end

    % ACS flows (affected by alphas)
    in = in.setVariable('ACS1_flow', ACS1_flow * alpha_1);
    in = in.setVariable('ACS2_flow', ACS2_flow * alpha_2);
    in = in.setVariable('ACS3_flow', ACS3_flow * alpha_3);
    in = in.setVariable('ACS4_flow', ACS4_flow * alpha_4);

    % === Run simulation ===
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
