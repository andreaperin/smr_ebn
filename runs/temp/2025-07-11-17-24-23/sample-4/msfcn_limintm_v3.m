function msfcn_limintm_v3(block)
% Level-2 MATLAB file S-Function for limited integrator demo.

%   Copyright 1990-2009 The MathWorks, Inc.

  setup(block);
  
%endfunction

function setup(block)
  
  %% Register number of dialog parameters   
  % block.NumDialogPrms = 3; -> check Derivative(block) function

  %% Register number of input and output ports
  block.NumInputPorts  = 4;
  block.NumOutputPorts = 5;

  %% Setup functional port properties to dynamically inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  % The field '.Dimensions' can be ommited if this is either dinamically or inherited 
  %% Setting inputs
  block.InputPort(1).Dimensions = 1;
  % Direct feedthrough FOR EACH INPUT must be set to true if the input is used in the mldOutput, or mdlGetTimeOfNextVarHit routine.
  block.InputPort(1).DirectFeedthrough = true;

  % FOR CLOCK!
  block.InputPort(3).Dimensions = 1;
  block.InputPort(3).DirectFeedthrough = true;

  % LOWER BOUND - DEACTIVATION OF PROTECTION SYSTEM
  block.InputPort(4).Dimensions = 1;
  block.InputPort(4).DirectFeedthrough = true;
  
  %% Setting outputs
  block.OutputPort(1).Dimensions = block.InputPort(1).Dimensions;
  block.OutputPort(2).Dimensions = 1;
  block.OutputPort(3).Dimensions = 1;
  %TEST
  block.OutputPort(4).Dimensions = 1;
  %block.OutputPort(5).Dimensions = 2;
 
  %% Set block sample time to continuous
  % block.SampleTimes = [0 0];
  %% Set block sample time to a fixed step ([step, offset])
  %block.SampleTimes = [0.1 0];
  %% Set block sample time to inherited
  block.SampleTimes = [-1 0];
  
  %% Setup Dwork
  block.NumContStates = 1;
  %block.NumDiscStates = 1;

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
  %TEST
  block.OutputPort(4).SamplingMode  = fd;
  block.OutputPort(5).SamplingMode = fd;

%endfunction

%% Post Propagation (Dwork)
function DoPostPropSetup(block)

  block.NumDworks = 3;
  block.Dwork(1).Name = 'xD'; 
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = true;

  %block.NumDworks = 1;
  block.Dwork(2).Name = 'xD_2'; 
  block.Dwork(2).Dimensions      = 1;
  block.Dwork(2).DatatypeID      = 0;
  block.Dwork(2).Complexity      = 'Real';
  block.Dwork(2).UsedAsDiscState = true;

  % for the first transition flag
  block.Dwork(3).Name            = 'xD_3';
  block.Dwork(3).Dimensions      = 1;
  block.Dwork(3).DatatypeID      = 0;
  block.Dwork(3).Complexity      = 'Real';
  block.Dwork(3).UsedAsDiscState = true;

%endfunction

%% Initial Conditions
function InitConditions(block)

  % Initialize Dwork
  block.ContStates.Data(1) = block.InputPort(1).Data;
  %block.DiscStates.Data(1) = 0;
  block.Dwork(1).Data = 0;
  block.Dwork(2).Data = 0;
  block.Dwork(3).Data = 0;
%endfunction

%% Output
function Output(block)

  block.OutputPort(2).Data = block.Dwork(1).Data;
  block.OutputPort(1).Data = block.InputPort(1).Data;
  block.OutputPort(3).Data = block.Dwork(2).Data;
  %TEST
  reloj = block.InputPort(3).Data;
  if reloj <= 5
    block.OutputPort(4).Data = 0;
  else
    block.OutputPort(4).Data = length(block.InputPort(1).Data);
  end 
  %block.OutputPort(5).Data = ?

%endfunction

%% Update: for major integration steps. Discrete states are updated here!
function Update(block)
    
     a = block.InputPort(2).Data;
     u = block.InputPort(1).Data;
     lb_ps = block.InputPort(4).Data;
     % !!!
     %arreglo_var = block.InputPort(4).Data;
     ub = min(a);
     lb = 0.9995;
     states = 0:1:length(a)+1;
     states_max = block.Dwork(2).Data;
     for i = 1:length(a)
         if i == 1
             if u <= ub && u >= lb
                 block.Dwork(1).Data = states(i);
                 states_max = max(states_max, states(i)); 
             end
         elseif i == length(a)
             if u <= a(i) && u >= a(i-1)
                 block.Dwork(1).Data = states(i);
                 states_max = max(states_max, states(i));
             elseif u > a(i)
                 %first_transition_flag = block.Dwork(3).Data;
                 %if first_transition_flag  == 0
                 block.Dwork(1).Data = states(i+1);
                 states_max = max(states_max, states(i+1));
                 %first_transition_flag = first_transition_flag + 1;
                 %block.Dwork(3).Data = first_transition_flag;
                 %end
             end
         else
             if u <= a(i) && u >= a(i-1)
                 block.Dwork(1).Data = states(i);
                 states_max = max(states_max, states(i));
             end
         end
     end

     % CONDITION: DEACTIVATION OF PROTECTION SYSTEM
     if u < lb_ps && states_max > min(states)
         states_max = 0;
         % stopping criteria?
     end

     block.Dwork(2).Data = states_max;
%endfunction

%% Derivative
function Derivative(block)
  block.Derivatives.Data = 0;
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

