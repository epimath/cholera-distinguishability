% Cost function calculator
% 1/14/2016 Marisa Eisenberg (marisae@umich.edu)

function negativeLL = cost_ML(tspan, data, xfcn, x0fcn, param, yfcn, costfcn)
param = abs(param); %forces all parameters to effectively be positive (so that if you try a negative value it just uses the postive one anyhow)

x0 = x0fcn(data,param(end));

[t,x] = ode45(xfcn,tspan,x0,[],param);

y=yfcn(param(end),x);

negativeLL = costfcn(data, y);
