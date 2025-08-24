function msfcn_indirectps_v1(block)
% Level-2 MATLAB file S-Function for limited integrator demo.

%   Copyright 1990-2009 The MathWorks, Inc.

  setup(block);
  
  
%endfunction

function setup(block)
  
  %% Register number of dialog parameters   
  % block.NumDialogPrms = 3; -> check Derivative(block) function

  %% Register number of input and output ports
  % InputPorts
  % - Current simulation time (t)
  % - Number of components (Ncomps)
  % - Number of states (Nstates)
  % - Matrix of states (mat_states_index)
  % - Transition cell (mat_trans)
  % - failed states, from MCS (failed_states)
  % - initial state (initial_state)
  block.NumInputPorts = 7;

  % Internal states (and eventually output ports):
  % 1 Current system state (discrete)
  % 2 Current transition component (discrete)
  % 3 Future system state (discrete)
  % 4 Future transition component (discrete)
  % 5 Future transition time (continuous)
  % 6 Previous transition time (continuous, just Dwork, and not and outputport)
  % 7 Failure flag (discrete)
  % obs: the state variables can become an output variable
  block.NumOutputPorts = 6;

  %% Setup functional port properties to dynamically inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;

  %% Allowing signals to be arrays
  block.AllowSignalsWithMoreThan2D = 1;

  %% Setting inputs
  % The field '.Dimensions' can be ommited if this is either dinamically or inherited
  block.InputPort(1).Dimensions = 1;
  % Direct feedthrough FOR EACH INPUT must be set to true if the input is used in the mldOutput, or mdlGetTimeOfNextVarHit routine.
  block.InputPort(1).DirectFeedthrough = true;

  block.InputPort(2).Dimensions = 1;
  block.InputPort(3).Dimensions = 3;
  block.InputPort(4).Dimensions = [3 2];
  block.InputPort(5).Dimensions = [2 2 3];
  block.InputPort(6).Dimensions = [3 3];
  block.InputPort(7).Dimensions = 3;
  
  %% Setting outputs
  block.OutputPort(1).Dimensions = 3;
  block.OutputPort(2).Dimensions = 1;
  block.OutputPort(3).Dimensions = 3;
  block.OutputPort(4).Dimensions = 1;
  block.OutputPort(5).Dimensions = 1;
  % For failure flag 
  block.OutputPort(6).Dimensions = 1;
 
  %% Set block sample time to continuous
  % block.SampleTimes = [0 0];

  %% Set block sample time to a fixed step ([step, offset])
  %block.SampleTimes = [0.1 0];

  %% Set block sample time to inherited
  block.SampleTimes = [-1 0];
  
  %% Setup Dwork
  % obs: one state variable for each output
  block.NumContStates = 2;
  % Do not initialize discrete states in the setup method (do it in
  % DoPostPropSetup method)
  % block.NumDiscStates = 5;

  %% Set the block simStateCompliance to default (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Register methods
  block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Update',                  @Update);
  block.RegBlockMethod('Derivatives',             @Derivative);  
  block.RegBlockMethod('Terminate',               @Terminate); % Required
  block.RegBlockMethod('SetInputPortSamplingMode',@SetInpPortFrameData); 
%endfunction

%% Sampling mode
function SetInpPortFrameData(block, idx, fd)

  block.InputPort(idx).SamplingMode = fd;
  block.OutputPort(1).SamplingMode  = fd;
  block.OutputPort(2).SamplingMode  = fd;
  block.OutputPort(3).SamplingMode  = fd;
  block.OutputPort(4).SamplingMode  = fd;
  block.OutputPort(5).SamplingMode  = fd;
  block.OutputPort(6).SamplingMode  = fd;

%endfunction

%% Post Propagation (Dwork)
function DoPostPropSetup(block)

  block.NumDworks = 7;

  block.Dwork(1).Name            = 'xS_1'; 
  block.Dwork(1).Dimensions      = 3;
  block.Dwork(1).DatatypeID      = 0;
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = true;

  block.Dwork(2).Name            = 'xS_2'; 
  block.Dwork(2).Dimensions      = 1;
  block.Dwork(2).DatatypeID      = 0;
  block.Dwork(2).Complexity      = 'Real';
  block.Dwork(2).UsedAsDiscState = true;

  block.Dwork(3).Name            = 'xS_3';
  block.Dwork(3).Dimensions      = 3;
  block.Dwork(3).DatatypeID      = 0;
  block.Dwork(3).Complexity      = 'Real';
  block.Dwork(3).UsedAsDiscState = true;

  block.Dwork(4).Name            = 'xS_4';
  block.Dwork(4).Dimensions      = 1;
  block.Dwork(4).DatatypeID      = 0;
  block.Dwork(4).Complexity      = 'Real';
  block.Dwork(4).UsedAsDiscState = true;

  block.Dwork(5).Name            = 'xS_5';
  block.Dwork(5).Dimensions      = 1;
  block.Dwork(5).DatatypeID      = 0;
  block.Dwork(5).Complexity      = 'Real';
  block.Dwork(5).UsedAsDiscState = false;

  block.Dwork(6).Name            = 'xS_6';
  block.Dwork(6).Dimensions      = 1;
  block.Dwork(6).DatatypeID      = 0;
  block.Dwork(6).Complexity      = 'Real';
  block.Dwork(6).UsedAsDiscState = false;

  block.Dwork(7).Name            = 'xS_7';
  block.Dwork(7).Dimensions      = 1;
  block.Dwork(7).DatatypeID      = 0;
  block.Dwork(7).Complexity      = 'Real';
  block.Dwork(7).UsedAsDiscState = true;

%endfunction

%% Initial Conditions
function InitConditions(block)
  % Internal states (and output ports):
  % 1 Current system state
  block.Dwork(1).Data = [1 1 1];
  % 2 Current state of the component under transition
  block.Dwork(2).Data = 1;
  % 3 Future system state 
  block.Dwork(3).Data = [0 0 0];
  % 4 Future state of the component under transition
  block.Dwork(4).Data = 1;
  % 5 Future transition time
  block.Dwork(5).Data = 0;
  % 6 Previous transition time
  block.Dwork(6).Data = 0;
  % 7 Failure flag
  block.Dwork(7).Data = 0;

  % Initialize Dwork
  %block.ContStates.Data(1) = block.InputPort(1).Data;
  %block.DiscStates.Data(1) = 0;

%endfunction

%% Output
function Output(block)
  % Declare output ports as a function of Dworks. Define Dworks at next step
  % Define explicitly output ports as function of input ports, if any

  block.OutputPort(1).Data = block.Dwork(1).Data;
  block.OutputPort(2).Data = block.Dwork(2).Data;
  block.OutputPort(3).Data = block.Dwork(3).Data;
  block.OutputPort(4).Data = block.Dwork(4).Data;
  block.OutputPort(5).Data = block.Dwork(5).Data;
  % For failure flag
  block.OutputPort(6).Data = block.Dwork(7).Data;

%endfunction

%% Update: for major integration steps. Discrete states are updated here!
function Update(block)
    % InputPorts
    % - Current simulation time (t)
    % - Number of components (Ncomps)
    % - Number of states (Nstates)
    % - Matrix of states (mat_states_index)
    % - Transition cell (cell_trans)
    % - failed states, from MCS (failed_states)
    % - initial state (initial_state)
    t = block.InputPort(1).Data;
    Ncomps = block.InputPort(2).Data;
    Nstates = block.InputPort(3).Data;
    mat_states_index = block.InputPort(4).Data;
    mat_trans = block.InputPort(5).Data;
    failed_states = block.InputPort(6).Data;
    initial_sys_state = block.InputPort(7).Data;
    
    % Internal states (and output ports):
    % 1 Current system state
    current_sys_state = block.Dwork(1).Data;
    % 2 Current state of the component under transition
    current_state_comp = block.Dwork(2).Data;
    % 3 Future system transition
    trans_sys = block.Dwork(3).Data;
    % 4 Future state of the component under transition
    new_trans = block.Dwork(4).Data;
    % 5 Future transition time
    new_ttf = block.Dwork(5).Data;
    % 6 Previous ttf (not an output port, but rather a state variable
    prev_ttf = block.Dwork(6).Data;
    % 7 Failure flag
    failure_flag = block.Dwork(7).Data;

    if new_ttf == 0
    % sample the transition times, based on current component state
    % [mat_ttf, min_ttf, comp_failure, new_trans] = indirect_method(Ncomps, Nstates, cell_trans, current_state, mat_states_index);
        mat_ttf = zeros(Ncomps, max(Nstates));
        % iterate over components
        for c = 1:Ncomps
            current_state_comp = current_sys_state(c);
            states_c = mat_states_index(c,:);
            trans_c_vec = mat_trans(current_state_comp,:,c);
            % iterate and sample over states
            for j = 1:length(states_c)
                if (j ~= current_state_comp) && (states_c(j) ~= 0)
                    lambda_trans_j = trans_c_vec(j);
                    ttf_j = (-1/lambda_trans_j)*log(rand);
                    mat_ttf(c,j) = ttf_j;
                end
            end
        end
        % indirect method: find the minimum sampled time to transition
        [r, c] = find(mat_ttf > 0);
        % comp, trans_state, non-negative sampled ttf
        mat_ttf_index = [r, c, mat_ttf(mat_ttf > 0)];
        [a, b] = find(mat_ttf_index(:,3) == min(mat_ttf_index(:,3)));
        trans_c = mat_ttf_index(a,1);
        new_trans = mat_ttf_index(a,2);
        min_ttf = mat_ttf_index(a,3);
        new_ttf = prev_ttf + min_ttf;
        %display(mat_ttf);   % !!!
        %display(new_trans); % !!!
        trans_sys = repmat(current_sys_state, 1);
        % updating new system state (NEXT TIME TRANSITION)
        trans_sys(trans_c) = new_trans;
        block.Dwork(3).Data = trans_sys;
        block.Dwork(4).Data = new_trans;
        block.Dwork(5).Data = new_ttf;
    end

    % if cst is less than the next transition time:
    if t < new_ttf
        % 1 Current system state
        %block.Dwork(1).Data = current_sys_state;
        % 2 Current state of component
        %block.Dwork(2).Data = current_state_comp;
        % 3 Future system transition
        %block.Dwork(3).Data = trans_sys;
        % 4 Future transition component
        %block.Dwork(4).Data = new_trans;
        % 5 Future transition time
        %block.Dwork(5).Data = new_ttf;
        % 6 Previous transition time
        %block.Dwork(6).Data = prev_ttf;
        % 7 Previous transition time
        %block.Dwork(7).Data = failure_flag;
    % else (i.e., reaching the transition time)
    else
    % 1) assign new states in the current states
        % 1 Current system state
        current_sys_state = repmat(trans_sys,1);
        block.Dwork(1).Data = current_sys_state;
        % 2 Current state of component
        current_state_comp = new_trans;
        block.Dwork(2).Data = current_state_comp;
    
        % detecting failure state
        % bool_failure_state = all(any(failed_states == current_sys_state));
        bool_array = failed_states == current_sys_state';
        if any(all(bool_array'))
            failure_flag = 1;
        else
            failure_flag = 0;
        end
        block.Dwork(7).Data = failure_flag;
        % detecting first failure state (for reliability analysis)
        % if failure_flag == 1 && first_failure_flag == 0
        %     first_failure_flag = 1;
        % end
    % and 2) store the next ones 
        mat_ttf = zeros(Ncomps, max(Nstates));
        % iterate over components
        for c = 1:Ncomps
            current_state_comp = current_sys_state(c);
            states_c = mat_states_index(c,:);
            trans_c_vec = mat_trans(current_state_comp,:,c);
            % iterate and sample over states
            for j = 1:length(states_c)
                if (j ~= current_state_comp) && (states_c(j) ~= 0)
                    lambda_trans_j = trans_c_vec(j);
                    ttf_j = (-1/lambda_trans_j)*log(rand);
                    mat_ttf(c,j) = ttf_j;
                end
            end
        end
        % indirect method: find the minimum sampled time to transition
        [r, c] = find(mat_ttf > 0);
        % comp, trans_state, non-negative sampled ttf
        mat_ttf_index = [r, c, mat_ttf(mat_ttf > 0)];
        [a, b] = find(mat_ttf_index(:,3) == min(mat_ttf_index(:,3)));
        trans_c = mat_ttf_index(a, 1);
        new_trans = mat_ttf_index(a,2);
        min_ttf = mat_ttf_index(a,3);
        prev_ttf = t;
        block.Dwork(6).Data = prev_ttf;
        new_ttf = prev_ttf + min_ttf;
        %display(mat_ttf);   % !!!
        %display(new_trans); % !!!
        trans_sys = repmat(current_sys_state, 1);
        % updating new system state (NEXT TIME TRANSITION)
        trans_sys(trans_c) = new_trans;
        block.Dwork(3).Data = trans_sys;
        block.Dwork(4).Data = new_trans;
        block.Dwork(5).Data = new_ttf;
    end
%endfunction

%% Derivative. Update of continuous states
function Derivative(block)
  block.Derivatives.Data = [0 0];
     %block.Output
%      lb = block.DialogPrm(1).Data;
%      ub = block.DialogPrm(2).Data;
%      u =  block.InputPort(1).Data;
%  
%      if (block.ContStates.Data <= lb && u < 0) || (block.ContStates.Data >= ub && u > 0)
%        block.Derivatives.Data = 0;
%      else
%        block.Derivatives.Data = u;
%      end

%endfunction

%% Terminate
function Terminate(block)

%endfunction
