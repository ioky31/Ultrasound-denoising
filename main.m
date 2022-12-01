clc;clf;clear;close all;
im = imread('C:\Users\86153\Desktop\组会\数据\invivo_normal\test\clean\637_A.png');
im = imresize(im, [128, 128]);


% denoised = imread("D:\result\2\images\638_A\638_A_10000.png");
% denoised = imresize(denoised, [128,128]);
%% add noise
sigma = sqrt(0.5);
im = double(im);
noise = exp(log(im) + normrnd(0,sigma,size(im,1),size(im,2)));
imwrite(uint8(noise),'noise0.5.png')
%% denoise
% FROST
% denoised=fcnFrostFilter(noise);
% LEE
% denoised = myLee(noise);
% SRAD
% denoised = specklefilt(noise,DegreeOfSmoothing=0.3,NumIterations=50);
% SARBM3D
% denoised = SARBM3D_v10(noise,6);
% NPD
% denoised = NPD(noise);
% OBNLM
% denoised = OBNLM(noise);
% NLM
% denoised = imnlmfilt(noise);
% FANS
% denoised = FANS(noise,6);
% LF
% denoised = LF(noise,370);
% res = im - noise;
% var(res(:))
%% imshow
% my_imshow(im, noise, denoised)
% [ssimval,peaksnr,snr,MSE,epi] = indicator(im, denoised);