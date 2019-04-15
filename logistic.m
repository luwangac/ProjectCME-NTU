function [f, df, y] = logistic(weights, data, targets, hyperparams)
% Calculate log likelihood and derivatives with respect to weights.
%
% Note: N is the number of examples and 
%       M is the number of features per example.
% Inputs:
% 	weights:    (M+1) x 1 vector of weights, where the last element
%               corresponds to bias (intercepts).
% 	data:       N x M data matrix where each row corresponds 
%               to one data point.
%	  targets:    N x 1 vector of targets class probabilities.
%   hyperparameters: The hyperparameter structure
%
% Outputs:
%	f:             The scalar error value.
%	df:            (M+1) x 1 vector of derivatives of error w.r.t. weights.
%   y:             N x 1 vector of probabilities. This is the output of the classifier.
%

%TODO: finish this function
N = size(data,1);
lambda  = hyperparams.weight_regularization; %get the hyperparameters
y = logistic_predict(weights,data); %prediction
error = y - targets;

%leave theta(final) alone; do not regularize bias. Convention.
temp = weights;
temp(end) = 0;
data = [data ones(size(data,1),1)];
grad = (1/N)*data'*error +(lambda/N)*temp; %when lambda = 0, there is no penalty
[f, frac_correct] = evaluate(targets,y);

% =============================================================
df = grad(:); %return grad as a column vector. For an n by m matrix
%the first n items are from first column, the next n items from 2nd column
%etc

end
