function [J, grad] = costFunction(w, x, y)
%COSTFUNCTION Compute cost and gradient for logistic regression
%
%
% Initialize some useful values
m = length(y); % number of fractures

%
%=====================Compute J (cost) ====================

H = sigmoid(x*w);
T = y.*log(H) + (1 - y).*log(1 - H);
J = -1/m*sum(T);

% ====================Compute grad (gradient)==================

for i = 1 : m
	grad = grad + (H(i) - y(i)) * x(i,:)';
end

grad = 1/m*grad;

% =============================================================

end
