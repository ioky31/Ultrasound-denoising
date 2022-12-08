clc;clf;clear;close all;
clean = imread('kidney_cut.bmp');
clean = imresize(clean, [128, 128]);
noise = imread('kidney_cut_sigma0.1.bmp');
noise = imresize(noise, [128, 128]);

%% denoise
% FROST
% denoised=fcnFrostFilter(noise);
% LEE
% denoised = myLee(noise);
% SRAD
denoised1 = specklefilt(noise,DegreeOfSmoothing=0.01,NumIterations=90);
% SARBM3D
denoised2 = SARBM3D_v10(noise,40);
% NPD
denoised3 = NPD(noise);
% OBNLM
% denoised = OBNLM(noise);
% NLM
% denoised = imnlmfilt(noise);
% FANS
denoised4 = FANS(noise,72);
% LF
% denoised = LF(noise,5, 5, 5);
% res = im - noise;
% var(res(:))
% WNNM
denoised5 = my_WNNM(noise,2.7);
% no_log
denoised6 = imread("D:\swimir_result\no_log\sigma0.1\test_images\kidney_cut.png");
denoised6 = imresize(denoised6, [128,128]);
% log
denoised7 = imread("D:\swimir_result\log\sigma0.1\test_images\kidney_cut.bmp");
denoised7 = imresize(denoised7, [128,128]);
%% imshow
figure('Name',"Denoise",'NumberTitle','off')
subplot(3,3,1)
imshow(uint8(clean))
title("Clean")
subplot(3,3,2)
imshow(uint8(noise))
title("Noise")
subplot(3,3,3)
imshow(uint8(denoised1))
title("SRAD")
subplot(3,3,4)
imshow(uint8(denoised2))
title("SARBM3D")
subplot(3,3,5)
imshow(uint8(denoised3))
title("NPD")
subplot(3,3,6)
imshow(uint8(denoised4))
title("FANS")
subplot(3,3,7)
imshow(uint8(denoised5))
title("WNNM")
subplot(3,3,8)
imshow(uint8(denoised6))
title("Non-Log")
subplot(3,3,9)
imshow(uint8(denoised7))
title("Log")
% [ssimval,peaksnr,snr,MSE,epi] = indicator(clean, denoised);