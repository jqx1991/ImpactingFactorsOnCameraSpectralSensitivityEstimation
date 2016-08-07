load('M:\D3x\Central\data\RGB_mean_ranked.mat')
RGB_D65 = RGB_mean_ranked;
RGB_D65 = RGB_D65(1:96,:);
load('E:\Dropbox\Works\papers\ResponsePrediction\SPD_Central.mat')
for i = 1:96
    SPD_D65(i,:) = (SPD_Central(2*i-1,:)+SPD_Central(2*i,:))/2;
end
SPD_D65 = SPD_D65(:,1:10:end);
clear RGB_mean_ranked SPD_Central

load('M:\D3x\A\data\RGB_mean_ranked.mat')
RGB_A = RGB_mean_ranked;
RGB_A = RGB_A(1:24,:);
load('E:\Dropbox\Works\papers\ResponsePrediction\SPD_Central_A.mat')
for i = 1:24
    SPD_A(i,:) = (SPD_Central_A(2*i-1,:)+SPD_Central_A(2*i,:))/2;
end
SPD_A = SPD_A(:,1:10:end);
clear RGB_mean_ranked SPD_Central_A

const_geometry = (pi/4)*((1/4)^2); % #F = 4
% reestimation
[~,order] = sort(DeltaE_test);
trainNum = 4:4:48;
for i = 1:length(trainNum)
    trainPatchesOrder = order(1:trainNum(i));
    testPatchesOrder = setdiff([1:96],trainPatchesOrder);
    testPatchesNum = length(testPatchesOrder);
    RGB_train = RGB_D65(trainPatchesOrder,:);
    RGB_test = RGB_D65(testPatchesOrder,:);
    SPD_train = SPD_D65(trainPatchesOrder,:);
    SPD_test = SPD_D65(testPatchesOrder,:);
    [CSS{i},~, CrossTalkMtx{i}, nonlinearCoef{i},~, DeltaE{i},~,~,~,~,lambda0{i}] = CameraResponsePrediction_training(RGB_train, SPD_train, 0.1, 15000);
    ISO = 1*ones(testPatchesNum,1);
    ExposureTime = 1/15*ones(testPatchesNum,1);
    Const = ISO.*ExposureTime * const_geometry * 10;
    ISO_A = 1*ones(24,1);
    ExposureTime_A = 1/15*ones(24,1);
    Const_A = ISO_A.*ExposureTime_A * const_geometry * 10;
    RGB_reconst_test_D65 = (real( ( diag(Const)*SPD_test*CSS{i}*CrossTalkMtx{i} + nonlinearCoef{i}(1) ).^nonlinearCoef{i}(3) )) + nonlinearCoef{i}(2);
    RGB_reconst_test_A = (real( ( diag(Const_A)*SPD_A*CSS{i}*CrossTalkMtx{i} + nonlinearCoef{i}(1) ).^nonlinearCoef{i}(3) )) + nonlinearCoef{i}(2);
    DeltaE_final_D65{i} = sRGB2CIEDeltaE(RGB_test.^(1/2.2),RGB_reconst_test_D65.^(1/2.2),'cie00');
    DeltaE_final_A{i} = sRGB2CIEDeltaE(RGB_A.^(1/2.2),RGB_reconst_test_A.^(1/2.2),'cie00');
    clear trainPatchesOrder testPatchesOrder testPatchesNum RGB_train RGB_test SPD_train SPD_test ISO ISO_A ExposureTime ExposureTime_A Const Const_A RGB_reconst_test_D65 RGB_reconst_test_A
    close all
end