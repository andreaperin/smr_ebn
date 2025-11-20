function p_fail = loop_prob(PGA)

p_fail = normcdf((log(PGA/(9.81))+0.25*normcdf(0.95))/0.2);

end