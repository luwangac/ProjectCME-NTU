%% Clear workspace.
clear
close all;

%% Data Import
Data = importdata('Q4_abalone_data.txt');
%separate categorical data
Sex = zeros(length(Data.textdata),1);
text_sex = cell2mat(Data.textdata);
for i = 1 : length(Data.textdata)
    if text_sex(i) == 'F'
        Sex(i) = 1;
    elseif text_sex(i) == 'I'
        Sex(i) = 0.75;
    end
    
end
abaloneData = [Sex Data.data];
%% Preprocessing data.
m = size(abaloneData,1);
rand_indices = randperm(m,1000); 
aba_train_inputs = abaloneData(rand_indices,1:end-1); 
aba_train_targets = abaloneData(rand_indices,end); %final column is target
%save('aba_train.mat','aba_train_inputs','aba_train_targets');
%
new_set = setdiff(abaloneData,abaloneData(rand_indices,:),'rows'); %data set excluding training set
rand_valid = randperm(size(new_set,1),1000);
aba_valid_inputs = new_set(rand_valid,1:end-1); 
aba_valid_targets = new_set(rand_valid,end);%final column is target
%save('aba_valid','aba_valid_inputs','aba_valid_targets');
%
newer_set = setdiff(new_set,new_set(rand_valid,:),'rows'); %data set excluding training and validation set
aba_test_inputs = newer_set(:,1:end-1); 
aba_test_targets = newer_set(:,end);%final column is target
%save('aba_test','aba_test_inputs','aba_test_targets');
hyperparameters.learning_rate = 0.85;
% Weight regularization parameter
hyperparameters.weight_regularization = 0.1; %zero implies  regression without penalty
% Number of iterations
hyperparameters.num_iterations = 50;
%% Linear Regression
X = [ones(size(aba_train_inputs,1),1) aba_train_inputs];
modelOLS = pinv(X)*aba_train_targets; %or
modelOLS2 = (X'*X)\X'*aba_train_targets;

%% Ridge Regression
lambda = 0:.01:20;
b = ridge(aba_train_targets, aba_train_inputs, lambda);
plot(lambda, b');
xlabel('Ridge parameter'); ylabel('Standardized coef.');
title('Ridge Trace for Abalone Data');
lambda_opt = 9;
X = [ones(size(aba_train_inputs,1),1) aba_train_inputs];
modelRidge = (X'*X+lambda_opt*eye(size(X,2)))\X'*aba_train_targets;

%% Lasso model
[BlassoAll, stats] = lasso(aba_train_inputs, aba_train_targets, 'CV', 20);
lassoPlot(BlassoAll,stats,'PlotType','CV'); 
modelLasso = [stats.Intercept(stats.Index1SE);...
    BlassoAll(:,stats.Index1SE)];
%% Comparisons
blinear = [ones(size(aba_valid_inputs,1),1) aba_valid_inputs]...
    *modelOLS;
cross_entropy_linear = evaluate(blinear,aba_valid_targets);
bridge = [ones(size(aba_valid_inputs,1),1) aba_valid_inputs]...
    *modelRidge;
cross_entropy_ridge = evaluate(bridge,aba_valid_targets);
blasso = [ones(size(aba_valid_inputs,1),1) aba_valid_inputs]...
    *modelLasso;
cross_entropy_lasso = evaluate(blasso,aba_valid_targets);
xx = linspace(min([aba_valid_targets;blinear;blasso;bridge]),max([aba_valid_targets;blinear;blasso;bridge]),100);
plot(aba_valid_targets,blinear,'o',aba_valid_targets,blasso,'x',aba_valid_targets,bridge,'.',xx,xx,'+');
xlabel('Actual target value'); ylabel('Estimates by model');
legend('OLS','LASSO','Ridge','Truth');

%% Comparison on Test Data
test_blinear = [ones(size(aba_test_inputs,1),1) aba_test_inputs]...
    *modelOLS;
test_bridge = [ones(size(aba_test_inputs,1) ,1) aba_test_inputs]...
    *modelRidge;
test_cross_entropy_linear = evaluate(test_blinear,aba_test_targets);
test_cross_entropy_ridge = evaluate(test_bridge,aba_test_targets);
xx = linspace(min([aba_test_targets;test_blinear;test_bridge]),max([aba_test_targets;test_blinear;test_bridge]),100);
plot(aba_test_targets,test_blinear,'ro',aba_test_targets,test_bridge,'m*',xx,xx,'+');
xlabel('Actual target value'); ylabel('Estimates by model');
legend('Linear','Ridge','Truth');