clear all
close all
clc

load('fashion_mnist.mat')

X_train = im2double(X_train);
X_test = im2double(X_test);

X_train = reshape(X_train,[60000 28 28 1]);
X_train = permute(X_train,[2 3 4 1]);

X_test = reshape(X_test,[10000 28 28 1]);
X_test = permute(X_test,[2 3 4 1]);

X_valid = X_train(:,:,:,1:5000);
X_train = X_train(:,:,:,5001:end);

y_valid = categorical(y_train(1:5000))';
y_train = categorical(y_train(5001:end))';
y_test = categorical(y_test)';
layer = reluLayer;
layers = [imageInputLayer([28 28 1])
        fullyConnectedLayer(700)
        layer
        fullyConnectedLayer(600)
        layer
        fullyConnectedLayer(500)
        layer
        fullyConnectedLayer(400)
        layer
        fullyConnectedLayer(300)
        layer
        fullyConnectedLayer(100)
        layer
        fullyConnectedLayer(10)
        softmaxLayer
        classificationLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',7,...
    'InitialLearnRate',1e-4, ...
    'L2Regularization',1e-5, ...
    'ValidationData',{X_valid,y_valid}, ...
    'GradientDecayFactor',0.9,...
    'ValidationFrequency',40,...
    'Verbose',false, ...
    'Plots','training-progress');

net = trainNetwork(X_train,y_train,layers,options);

y_pred = classify(net,X_test);

plotconfusion(y_test,y_pred)

%%
clear all
close all
clc

load('fashion_mnist.mat')

X_train = im2double(X_train);
X_test = im2double(X_test);

X_train = reshape(X_train,[60000 28 28 1]);
X_train = permute(X_train,[2 3 4 1]);

X_test = reshape(X_test,[10000 28 28 1]);
X_test = permute(X_test,[2 3 4 1]);

X_valid = X_train(:,:,:,1:5000);
X_train = X_train(:,:,:,5001:end);

y_valid = categorical(y_train(1:5000))';
y_train = categorical(y_train(5001:end))';
y_test = categorical(y_test)';

layers = [
    imageInputLayer([28 28 1],"Name","imageinput")
    convolution2dLayer([2 2],10,"Name","conv_1","Padding","same")
    tanhLayer("Name","tanh_1")
    maxPooling2dLayer([2 2],"Name","maxpool_1","Padding","same")
    convolution2dLayer([2 2],22,"Name","conv_2")
    tanhLayer("Name","tanh_2")
    maxPooling2dLayer([5 5],"Name","maxpool_2","Padding","same")
    convolution2dLayer([2 2],130,"Name","conv_3")
    tanhLayer("Name","tanh_3")
    fullyConnectedLayer(90,"Name","fc_1")
    tanhLayer("Name","tanh_4")
    fullyConnectedLayer(10,"Name","fc_2")
    softmaxLayer("Name","softmax")
    classificationLayer("Name","classoutput")];

options = trainingOptions('adam', ...
    'MaxEpochs',5,...
    'InitialLearnRate',1e-3, ...
    'L2Regularization',1e-4, ...
    'ValidationData',{X_valid,y_valid}, ...
    'Verbose',false, ...
    'Plots','training-progress');

net = trainNetwork(X_train,y_train,layers,options);

y_pred = classify(net,X_test);

plotconfusion(y_test,y_pred)
%%
net = resnet50;

y_pred = classify(net,X_test);