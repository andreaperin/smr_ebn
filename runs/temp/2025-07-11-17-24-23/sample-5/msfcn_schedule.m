function msfcn_schedule(block)
% Level-2 MATLAB file S-Function for limited integrator demo.

%   Copyright 1990-2009 The MathWorks, Inc.

  setup(block);
  
%endfunction

function setup(block)
  
  %% Register number of dialog parameters   
  % block.NumDialogPrms = 3; -> check Derivative(block) function

  %% Register number of input and output ports
  block.NumInputPorts  = 4;
  block.NumOutputPorts = 2;

  %% Setup functional port properties to dynamically inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;

  %% Allowing signals to be arrays
  block.AllowSignalsWithMoreThan2D = 1;
 
  % The field '.Dimensions' can be ommited if this is either dinamically or inherited 
  %% Setting inputs
  block.InputPort(1).Dimensions = 1;
  % Direct feedthrough FOR EACH INPUT must be set to true if the input is used in the mldOutput, or mdlGetTimeOfNextVarHit routine.
  block.InputPort(1).DirectFeedthrough = true;

  block.InputPort(2).Dimensions = 1;
  block.InputPort(2).DirectFeedthrough = true;

  % For threshold of the activation function
  block.InputPort(3).Dimensions = 1;
  %block.InputPort(3).DirectFeedthrough = true;

  % FOR past DISCRETE STATE VECTOR
  block.InputPort(4).Dimensions = 2;
  %block.InputPort(4).DirectFeedthrough = true;
  
  %% Setting outputs
  % Activation of PS
  block.OutputPort(1).Dimensions = 1;
  % FOR current DISCRETE STATE VECTOR
  block.OutputPort(2).Dimensions = block.InputPort(4).Dimensions;
  
 
  %% Set block sample time to continuous
  % block.SampleTimes = [0 0];
  %% Set block sample time to a fixed step ([step, offset])
  %block.SampleTimes = [0.1 0];
  %% Set block sample time to inherited
  block.SampleTimes = [-1 0];
  
  %% Setup Dwork
  %block.NumContStates = 2;
  %block.NumDiscStates = 2;

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

%endfunction

%% Post Propagation (Dwork)
function DoPostPropSetup(block)

  block.NumDworks = 2;
  block.Dwork(1).Name = 'act_PS1'; 
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = true;

  %block.NumDworks = 1;
  block.Dwork(2).Name = 'state_PS1'; 
  block.Dwork(2).Dimensions      = block.InputPort(4).Dimensions;
  block.Dwork(2).DatatypeID      = 0;
  block.Dwork(2).Complexity      = 'Real';
  block.Dwork(2).UsedAsDiscState = true;

%endfunction

%% Initial Conditions
function InitConditions(block)

  % Initialize Dwork
  %block.ContStates.Data(1) = block.InputPort(1).Data;
  %block.DiscStates.Data(1) = 0;
  block.Dwork(1).Data = 0;
  block.Dwork(2).Data = [0 1];

%endfunction

%% Output
function Output(block)

  block.OutputPort(1).Data = block.Dwork(1).Data;
  block.OutputPort(2).Data = block.Dwork(2).Data;

%endfunction

%% Update: for major integration steps. Discrete states are updated here!
function Update(block)
    
     %state_var_ps: which protection system is being demanded
     state_var_ps = block.InputPort(1).Data;
     %failure_flag: if the system is not faulty then it's available
     failure_flag = block.InputPort(2).Data;
     if failure_flag == 0
         ava_ps = 1;
     else
         ava_ps = 0;
     end

     %treshold for demanding of protection system
     tres_demand_ps = block.InputPort(3).Data;
     if state_var_ps  == tres_demand_ps
         demand_ps = 1;
     else
         demand_ps = 0;
     end
     
     current_state = [demand_ps, ava_ps];
     prev_state = block.InputPort(4).Data';

     if demand_ps == 1
         if ava_ps == 1
             act_ps = 1;
             % cases: event list criteria
             if ~ all(current_state == prev_state)
                 nombre = get_param(bdroot, 'Name');
                 %set_param(nombre,"SimulationCommand","stop");
             end
         else
             act_ps = 0;
             if ~ all(current_state == prev_state)
                 nombre = get_param(bdroot, 'Name');
                 %set_param(nombre,"SimulationCommand","stop");
             end
         end
     else
         act_ps = 0;
         % saving event
         if prev_state(1) == 1
              nombre = get_param(bdroot, 'Name');
              %set_param(nombre,"SimulationCommand","stop");
         end
     end 

     block.Dwork(1).Data = act_ps;
     block.Dwork(2).Data = current_state;
     
%endfunction

%% Derivative
function Derivative(block)
%  block.Derivatives.Data = 0;
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

