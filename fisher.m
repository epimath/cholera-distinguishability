% Simplified FIM assuming unweighted least squares
% 9/8/2015 Marisa Eisenberg (marisae@umich.edu)

function FIM = fisher(tspan,x0,params,xfcn,yfcn)
% using 1% perturbations in the parameter to calculate the parameter
% sensitivities, with the simplified FIM = X'X. (note we've dropped the
% weighting term for the structural identifiability analysis. For more
% information, please see:
% Eisenberg MC, Hayashi MAL. 2014. Determining Identifiable Parameter 
% Combinations Using Subset Profiling. Mathematical Biosciences 256: 116?126.
 
X = [];

for j=1:length(params)
    params1 = params;
    params2 = params;
    params1(j) = 1.01*params(j);
    params2(j) = 0.99*params(j);
    [t x1] = ode45(xfcn,tspan,x0,[],params1);
    [t x2] = ode45(xfcn,tspan,x0,[],params2);
    X = [X, (yfcn(params1(end),x1) - yfcn(params2(end),x2))./(0.02*params(j))];
    %this fills in the jth column of the design matrix with the sensitivities to parameter j 
    %at each time point.  Note x1(:,2) is the infected compartment.
end

FIM = X'*X;




