%Function to calculate the overpressure from a delayed ignition of H2
function [P,m_dot_release] = Overpressure (r,d,t,m_dot)

section = randi(5); %sampling where the leakage occurs to get op. conditions
v = [3306,1899,2197,2230,2537]; %from python code CoolProp
rho = [0.0385,0.111,0.139,0.217,0.3]; %from python code CoolProp
rho = rho(section);
v = v(section);
m_dot_release = min(m_dot,v*rho*pi*d^2/4);
m_H = m_dot_release*t;
P_atm = 101325;
E = 119*10^6*m_H*2;
R = r*P_atm^(1/3)/E^(1/3);
P = 0.34.*R.^(-4/3)+0.062.*R.^(-2)+0.0033.*R.^(-3); %from "SANDIA - Final Report on Hydrogen Plant Hazards and Risk Analysis Supporting Hydrogen Plant Siting near Nuclear Power Plants"
P = P*P_atm;