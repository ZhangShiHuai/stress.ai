function [theta, cost, MLE] = Logistic_Regression(Slip_Tendency,criticality)
%
% Input parameter
X = Slip_Tendency;
y = criticality;
% Add intercept term to X
[m, n] = size(X);
X = [ones(m, 1) X];
%--------------------------Minimize the cost function----------------------
% Initialize fitting parameters
initial_theta = zeros(n + 1, 1);
%  Set options for fminunc
options = optimset('GradObj', 'on', 'MaxIter', 400);
%  Run fminunc to obtain the optimal theta
%  This function will return theta and the cost 
[theta, cost] = fminunc(@(t)(costFunction(t, X, y)), initial_theta, options);
MLE = exp(-cost);
end

