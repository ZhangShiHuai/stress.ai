function [w, cost, MLE] = Logistic_Regression(Slip_Tendency,criticality)
%
% Input parameter
x = Slip_Tendency;
y = criticality;
% Input of sigmiod function
[m, n] = size(x);
z = [ones(m, 1) x];
%--------------------------Minimize the cost function----------------------
% Initialize fitting parameters
initial_w = zeros(n + 1, 1);
%  Set options for fminunc
options = optimset('GradObj', 'on', 'MaxIter', 400);
%  Run fminunc to obtain the optimal w
%  This function will return w and the cost 
[w, cost] = fminunc(@(t)(costFunction(t, z, y)), initial_theta, options);
MLE = exp(-cost);
end

