%% The test script to run knn

k = 1;

resource_train = load('mnist_train.mat');
train_data = resource_train.train_inputs;
train_labels = resource_train.train_targets;
 
resource_valid = load('mnist_valid.mat');
valid_data = resource_valid.valid_inputs;
valid_labels = resource_valid.valid_targets;

resource_valid = load('mnist_test.mat');
test_data = resource_valid.test_inputs;
test_labels = resource_valid.test_targets;


classification_rate = zeros(1, 5);
%Testing with validation data
for i =1:5 
	caculate_labels = run_knn(2*i-1, train_data, train_labels, valid_data);

	accurate = caculate_labels-valid_labels;

	classification_rate(i) = (size(caculate_labels,1)-sum(abs(accurate)))/size(caculate_labels,1);

end
%subplot(2,1,1)
plot(1:2:9, classification_rate,'ro'); ylabel('rate'); xlabel('k');
title('Classification rate with validation data');

%Testing with test data
%%
for i =1:5 
	caculate_labels = run_knn(2*i-1, train_data, train_labels, test_data);

	accurate = caculate_labels-test_labels;

	classification_rate(i) = (size(caculate_labels,1)-sum(abs(accurate)))/size(caculate_labels,1);
end
%subplot(2,1,2)
plot(1:2:9, classification_rate,'ro'); ylabel('rate'); xlabel('k');
title('Classification rate with test data');

