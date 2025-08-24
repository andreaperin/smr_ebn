clear
close all

tic

load('Failure_model_outputs')
evalc('SMDFR_Parameters');

model = 'SMDFR_HTE_model';
All_Outputs = 0;

%ADD UNCERTAINTY???
ACS1_flow = 25;
ACS2_flow = 25;
ACS3_flow = 25;
ACS4_flow = 25;
tsim=1200;

in = Simulink.SimulationInput(model);
in = in.setModelParameter('StartTime','1','StopTime',num2str(tsim));
in = in.setVariable('All_Outputs',All_Outputs);
in = in.setVariable('LOCA1_time',LOCA1_time);
in = in.setVariable('LOCA2_time',LOCA2_time);
in = in.setVariable('LOCA3_time',LOCA3_time);
in = in.setVariable('LOCA4_time',LOCA4_time);
in = in.setVariable('ACS_1',ACS_1);
in = in.setVariable('ACS_2',ACS_2);
in = in.setVariable('ACS_3',ACS_3);
in = in.setVariable('ACS_4',ACS_4);
in = in.setVariable('EDG_1',EDG_1);
in = in.setVariable('EDG_2',EDG_2);
in = in.setVariable('Power',Power);
% in = in.setVariable('ACS1_response_time',ACS1_response_time);
% in = in.setVariable('ACS2_response_time',ACS2_response_time);
% in = in.setVariable('ACS3_response_time',ACS3_response_time);
% in = in.setVariable('ACS4_response_time',ACS4_response_time);
% in = in.setVariable('EDG1_response_time',EDG1_response_time);
% in = in.setVariable('EDG2_response_time',EDG2_response_time);
% in = in.setVariable('PDP11_response_time',PDP11_response_time);
% in = in.setVariable('PDP12_response_time',PDP12_response_time);
% in = in.setVariable('PDP21_response_time',PDP21_response_time);
% in = in.setVariable('PDP22_response_time',PDP22_response_time);
% in = in.setVariable('PDP31_response_time',PDP31_response_time);
% in = in.setVariable('PDP32_response_time',PDP32_response_time);
% in = in.setVariable('PDP41_response_time',PDP41_response_time);
% in = in.setVariable('PDP42_response_time',PDP42_response_time);
in = in.setVariable('PDP11',PDP11);
in = in.setVariable('PDP12',PDP12);
in = in.setVariable('PDP21',PDP21);
in = in.setVariable('PDP22',PDP22);
in = in.setVariable('PDP31',PDP31);
in = in.setVariable('PDP32',PDP32);
in = in.setVariable('PDP41',PDP41);
in = in.setVariable('PDP42',PDP42);
in = in.setVariable('MSLB1',MSLB1);
in = in.setVariable('MSLB2',MSLB2);
in = in.setVariable('MSLB3',MSLB3);
in = in.setVariable('MSLB4',MSLB4);
in = in.setVariable('ACS1_flow',ACS1_flow*alpha_1);
in = in.setVariable('ACS2_flow',ACS2_flow*alpha_2);
in = in.setVariable('ACS3_flow',ACS3_flow*alpha_3);
in = in.setVariable('ACS4_flow',ACS4_flow*alpha_4);
in = in.setVariable('thermal_failure',thermal_failure);
in = in.setVariable('LHS',LHS);
in = in.setVariable('thermal_failure_time',thermal_failure_time);
in = in.setVariable('PGA',PGA);

out = sim(in);

T_W1 = reshape(max(out.T_W1,[],1),[1,tsim]);
T_W2 = reshape(max(out.T_W2,[],1),[1,tsim]);
T_W3 = reshape(max(out.T_W3,[],1),[1,tsim]);
T_W4 = reshape(max(out.T_W4,[],1),[1,tsim]);

save('Simulation_model_outputs','T_W1','T_W2','T_W3','T_W4')

toc