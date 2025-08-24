function [LOP,m_dot_release,sw_f,tr_f,P] = Loss_of_Power(d,t,sw_fun,tr_fun,demand,distance)

%Overpressure
r=distance;
m_dot = 3.3; %from "US Department of Energy, Preconceptual Designs of Coupled Power Delivery between a 4-Loop PWR and 100-500 MWe HTSE Plants. April 2023."
m_dot=m_dot*demand;
[P,m_dot_release] = Overpressure (r,d,t,m_dot);

%Power network fragility
sw_f = min(1,max(0,sw_fun(P)));
tr_f = min(1,max(0,tr_fun(P)));
r1 = rand();
r2 = rand();
r3 = rand();

if r1<sw_f || (r2<tr_f && r3<tr_f)
    LOP = 1;
else
    LOP=0;
end