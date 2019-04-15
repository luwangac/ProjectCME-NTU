%% Clear workspace.
clear
close all;

%% Load data and visualize data.
load mnist_train;
load mnist_valid;
load mnist_train_small;
load mnist_test.mat;
m = size( train_inputs_small, 1);
rand_indices = randperm(m,5); %generate 5 random samples for visualizing (<=5 will be displayed, depending on 
%plot_digits)
sel =  train_inputs_small(rand_indices, :); 
plot_digits(sel);
fprintf('Program paused. Press enter to continue.\n');
pause;
%% TODO: Initialize hyperparameters.
% Learning rate
hyperparameters.learning_rate = 0.85;
% Weight regularization parameter

hyperparameters.weight_regularization = .1; %zero implies logistic regression without penalty
% Number of iterations
hyperparameters.num_iterations = 200;
% Logistics regression weights
% TODO: Set random weights.
num_weights = size( train_inputs_small, 2);
weights = rand(num_weights+1,1)*0.01;


%% Verify that your logistic function produces the right gradient, diff should be very close to 0
% this creates small random data with 80 examples and 10 dimensions and checks the gradient on
% that data.
nexamples = 80;
ndimensions = 10;
diff = checkgrad('logistic', ...
	             randn((ndimensions + 1), 1), ...   % weights
                 0.001,...                          % perturbation
                 randn(nexamples, ndimensions), ... % data        
                 rand(nexamples, 1), ...            % targets
                 hyperparameters);                        % other hyperparameters

N = size( train_inputs_small, 1);
%% Begin learning with gradient descent.
cross_entropy_trains = zeros(1, hyperparameters.num_iterations);
cross_entropy_valids = zeros(1, hyperparameters.num_iterations);
for t = 1:hyperparameters.num_iterations

	%% TODO: You will need to modify this loop to create plots etc.

	% Find the negative log likelihood and derivative w.r.t. weights.
	[f, df, predictions] = logistic(weights, ...
                                            train_inputs_small, ...
                                            train_targets_small, ...
                                           hyperparameters);

  [cross_entropy_train, frac_correct_train] = evaluate( train_targets_small, predictions);
                    cross_entropy_trains(t) = cross_entropy_train;

	% Find the fraction of correctly classified validation examples.
	[temp, temp2, frac_correct_valid] = logistic(weights, ...
                                                 valid_inputs, ...
                                                 valid_targets, ...
                                                 hyperparameters);

    if isnan(f) || isinf(f)
		error('nan/inf error');
    end

	%% Update parameters.
	weights = weights - hyperparameters.learning_rate .* df / N;

    predictions_valid = logistic_predict(weights, valid_inputs);
    [cross_entropy_valid, frac_correct_valid] = evaluate(valid_targets, predictions_valid);
                    cross_entropy_valids(t) = cross_entropy_valid;

	%% Print some stats.
	fprintf(1, 'ITERATION:%4i   NLOGL:%4.2f TRAIN CE %.6f TRAIN FRAC:%2.2f VALIC_CE %.6f VALID FRAC:%2.2f\n',...
			t, f/N, cross_entropy_train, frac_correct_train*100, cross_entropy_valid, frac_correct_valid*100);

end
plot(1:hyperparameters.num_iterations, cross_entropy_valids, '-r', 1:hyperparameters.num_iterations, cross_entropy_trains, '-b'); 
legend({'validation', 'training'}); xlabel('Number of Iterations'); ylabel('Cross-entropy');
%draw table 
%%
t = table();
t.CROSS_ENTROPY_TRAIN = cross_entropy_train;
t.CROSS_ENTROPY_VALID = cross_entropy_valid;
t.FRACTION_CORRECT_TRAIN = frac_correct_train;
t.FRACTION_CORRECT_VALID = frac_correct_valid;
t.HYPERPARAMETERS  = struct2table(hyperparameters);
