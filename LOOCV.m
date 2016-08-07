function [test_result, CSS, CrossTalkMtx, nonlinearCoef, fval, lambda, DeltaE_train, DeltaE_test] = LOOCV(RGB, SPD)
totalPatchesNum = size(RGB,1);
wlNum = size(SPD,2);
wlInterval = (780-380)/(wlNum-1);
patchesNum = 1;
const_geometry = (pi/4)*((1/4)^2); % #F = 4
ISO = 1*ones(patchesNum,1);
ExposureTime = 1/15*ones(patchesNum,1);
Const = ISO.*ExposureTime * const_geometry * wlInterval;
test_result = zeros(totalPatchesNum,1);
% prelocating
CSS = cell(totalPatchesNum,1);
CrossTalkMtx = cell(totalPatchesNum,1);
nonlinearCoef = zeros(totalPatchesNum,3);
DeltaE_train = zeros(totalPatchesNum,totalPatchesNum-1);
fval = zeros(totalPatchesNum,1);
lambda = zeros(totalPatchesNum,1);
DeltaE_test = zeros(totalPatchesNum,1);

for i = 1:totalPatchesNum
    RGB_train = RGB;
    RGB_train(i,:) = [];
    RGB_test = RGB(i,:);
    SPD_train = SPD;
    SPD_train(i,:) = [];
    SPD_test = SPD(i,:);
    disp(['Processing the testing of No.',num2str(i),' patch.']);
    [CSS{i}, ~, CrossTalkMtx{i}, nonlinearCoef(i,:), ~, DeltaE_train(i,:), ~,fval(i,:),~,~,lambda(i,:)] = CameraResponsePrediction_training(RGB_train, SPD_train, 0.1, 12000);
    RGB_reconst_testing = (real( ( diag(Const)*SPD_test*CSS{i}*CrossTalkMtx{i} + nonlinearCoef(i,1) ).^nonlinearCoef(i,3) )) + nonlinearCoef(i,2);
    DeltaE_test(i,:) = sRGB2CIEDeltaE(RGB_test.^(1/2.2),RGB_reconst_testing.^(1/2.2),'cie00');
    test_result(i,:) = mean(DeltaE_test(i)) + 0.1*max(DeltaE_test(i));
    clear RGB_train RGB_test RGB_reconst_testing DeltaE
    close all;clc;
end