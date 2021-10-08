function g = sigmoid(z)
% Sigmoid function computes the sigmoid of z.
% 
% g = zeros(size(z));
g = 1./(1 + exp(-z));
end
