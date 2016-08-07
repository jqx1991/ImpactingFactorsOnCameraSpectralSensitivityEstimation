function [CSS, CSS0, CrossTalkMtx, nonlinearCoef, nonlinearCoef0, DeltaE, RelativeError,fval,exitflag,output,lambda0] = CameraResponsePrediction_training(RGB, SPD, maxWeight, maxEvals)
% function CameraSpectralSensitivityEstimate_weighted() calculate camera
% spectral sensitivity using nonnegtive-smoothness algorithm, given
% captured camera responses CapturedRGB and corresponding SPD captured by
% CS2000, each row of which denoting power distribution of a specified
% patch over 400-780nm.
%
% Instead of choosing all color patches equally as function
% CameraSpectralSensitivityEstimate() does,
% CameraSpectralSensitivityEstimate_weighted() assignes different weights
% to different patches(weight = 1 denotes picking a patch normally and
% weight = 0 denotes totally preventing a patch from calculation) so as to
% get optimal CSS.
%
% This function might be VERY slow, owing to time-consuming non-linear
% least squares method, as the number of patches is large.

if ~exist('max2thDerivative','var')
    max2thDerivative = 1000;
end
if ~exist('maxEvals','var')
    maxEvals = 30000;
end

% Colorimetric Characterazation Using Root-Polynomial Method
% XYZ = RGB2XYZ_RPCC(RGB)/100;

wlNum = size(SPD,2);
wlInterval = (780-380)/(wlNum-1);
patchesNum = size(SPD,1);

% Initial Guess for Camera Spectral Sensitivity Estimation Begin
S = zeros(wlNum-2,wlNum);
for i = 1:wlNum-2
    S(i,i) = -1;
    S(i,i+1) = 2;
    S(i,i+2) = -1;
end
Z = [ones(patchesNum,1);zeros(wlNum-2,1)];

const_geometry = (pi/4)*((1/4)^2); % #F = 4
ISO = 1*ones(patchesNum,1);
ExposureTime = 1/15*ones(patchesNum,1);
Const = ISO.*ExposureTime * const_geometry * wlInterval;
L_relative_R = ((diag(RGB(:,1)))^-1)*SPD;
L_relative_G = ((diag(RGB(:,2)))^-1)*SPD;
L_relative_B = ((diag(RGB(:,3)))^-1)*SPD;
% use RegTools
[U,sm,~,~] = cgsvd(blkdiag(diag(Const)*L_relative_R,diag(Const)*L_relative_G,diag(Const)*L_relative_B),blkdiag(S,S,S));
[lambda0,~,~,~] = l_curve(U,sm,ones(patchesNum*3,1));
clear U sm

lambda0 = lambda0/120;
A_Red_Ch  = [diag(Const)*L_relative_R;lambda0*S];
A_Green_Ch  = [diag(Const)*L_relative_G;lambda0*S];
A_Blue_Ch  = [diag(Const)*L_relative_B;lambda0*S];
% the CSS estimation formula
CSS_Red_Ch = lsqnonneg(A_Red_Ch, Z);
CSS_Green_Ch = lsqnonneg(A_Green_Ch, Z);
CSS_Blue_Ch = lsqnonneg(A_Blue_Ch, Z);
CSS0 = [CSS_Red_Ch, CSS_Green_Ch, CSS_Blue_Ch];
% Camera Spectral Sensitivity Estimation End
% Goodness of Estimation
RGB_reconst = diag(Const)*SPD*[CSS_Red_Ch, CSS_Green_Ch, CSS_Blue_Ch];

fitModel = @(x,t) real((t+x(1)).^x(3)) + x(2);
options = optimoptions(@lsqcurvefit ,'MaxFunEvals',20000);
nonlinearCoef0 = lsqcurvefit(fitModel,[0;0;1],RGB_reconst(:), RGB(:), [-1E-3;-100;0.9],[100;100;1],options);

% Initial Guess for Camera Spectral Sensitivity Estimation End
% Optimization Begin
CSS_Red_Ch_opt_temp = '';CSS_Green_Ch_opt_temp = '';CSS_Blue_Ch_opt_temp = '';
for i = 1:wlNum
    CSS_Red_Ch_opt_t = ['x(',num2str(i),') '];
    CSS_Red_Ch_opt_temp = [CSS_Red_Ch_opt_temp,CSS_Red_Ch_opt_t];
    CSS_Green_Ch_opt_t = ['x(',num2str(i+wlNum),') '];
    CSS_Green_Ch_opt_temp = [CSS_Green_Ch_opt_temp,CSS_Green_Ch_opt_t];
    CSS_Blue_Ch_opt_t = ['x(',num2str(i+2*wlNum),') '];
    CSS_Blue_Ch_opt_temp = [CSS_Blue_Ch_opt_temp,CSS_Blue_Ch_opt_t];
end
clear CSS_Red_Ch_opt_t CSS_Green_Ch_opt_t CSS_Blue_Ch_opt_t

eval(['CSS_opt = @(x) [[',CSS_Red_Ch_opt_temp,']'', [',CSS_Green_Ch_opt_temp,']'', [',CSS_Blue_Ch_opt_temp,']''];']);

CrossTalkMtx = @(x) [x(3*wlNum+1) x(3*wlNum+2) 0; x(3*wlNum+3) x(3*wlNum+4) x(3*wlNum+5); 0 x(3*wlNum+6) x(3*wlNum+7)];
RGB_reconst = @(x) (real( ( diag(Const)*SPD*CSS_opt(x)*CrossTalkMtx(x) + x(3*wlNum+8) ).^x(3*wlNum+10) )) + x(3*wlNum+9);
RGB_reconst = @(x) max(RGB_reconst(x),0);
RGB_reconst = @(x) min(RGB_reconst(x),1);

% XYZ_reconst = @(x) RGB2XYZ_RPCC(RGB_reconst(x))/100;

% DeltaE = @(x) sRGB2CIEDeltaE(XYZ,XYZ_reconst(x),'cie00','XYZ');
DeltaE = @(x) sRGB2CIEDeltaE(RGB.^(1/2.2),RGB_reconst(x).^(1/2.2),'cie00');
costfun = @(x) mean(DeltaE(x)) + maxWeight*max(DeltaE(x));

CSS_opt0 = [CSS_Red_Ch; CSS_Green_Ch; CSS_Blue_Ch; [1;0;0;1;0;0;1]; nonlinearCoef0];
% force the maximum of second-derivative to be less than max2thDerivative
% for 5nm wavelength interval, max2thDerivative = 0.15 is recommended
% A = [[blkdiag(S,S,S),zeros(3*wlNum-6,7)];[blkdiag(-S,-S,-S),zeros(3*wlNum-6,7)]];
% b = max2thDerivative*ones(6*wlNum-12,1);
lb = [0.6*CSS_Red_Ch; 0.6*CSS_Green_Ch; 0.6*CSS_Blue_Ch; [0.95;0;0;0.95;0;0;0.95];   -1E-1;-100;0.9];
ub = [1.5*CSS_Red_Ch; 1.5*CSS_Green_Ch ;1.5*CSS_Blue_Ch; [1;0.05;0.05;1;0.05;0.05;1]; 100;100;1];

options = optimoptions(@fmincon ,'MaxFunEvals',maxEvals,'Display','iter','PlotFcns',@optimplotfval);
[CSS_optimum,fval,exitflag,output] = fmincon(costfun,CSS_opt0,[],[],[],[],lb,ub,@myConfun,options);
clear CSS_opt CrossTalkMtx RGB_reconst DeltaE
% Optimization End
CSS_Red_Ch_opt = CSS_optimum(1:wlNum);
CSS_Green_Ch_opt = CSS_optimum(wlNum+1:2*wlNum);
CSS_Blue_Ch_opt = CSS_optimum(2*wlNum+1:3*wlNum);
CSS = [CSS_Red_Ch_opt, CSS_Green_Ch_opt, CSS_Blue_Ch_opt];
CrossTalkMtx = [CSS_optimum(3*wlNum+1) CSS_optimum(3*wlNum+2) 0; CSS_optimum(3*wlNum+3) CSS_optimum(3*wlNum+4) CSS_optimum(3*wlNum+5); 0 CSS_optimum(3*wlNum+6) CSS_optimum(3*wlNum+7)];

nonlinearCoef = CSS_optimum(3*wlNum+8:3*wlNum+10)';
C1 = nonlinearCoef(1); C2 = nonlinearCoef(2); alpha = nonlinearCoef(3);
RGB_reconst = (real( ( diag(Const)*SPD*CSS*CrossTalkMtx + C1 ).^alpha )) + C2;
% XYZ_reconst = RGB2XYZ_RPCC(RGB_reconst)/100;

% DeltaE = sRGB2CIEDeltaE(XYZ,XYZ_reconst,'cie00','XYZ');
DeltaE = sRGB2CIEDeltaE(RGB.^(1/2.2),RGB_reconst.^(1/2.2),'cie00')';
RelativeError = mean(abs(RGB-RGB_reconst)./RGB,2);
% set graphFlag = 1 to display color diffs of all patches, 0 to display max
% and mean values only, -1 to display nothing.

end
function [c,ceq]=myConfun(x)
c = -real(x(3*41+8)^x(3*41+10))-x(3*41+9);
ceq=[];
end
