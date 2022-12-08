% Edge preservation index created on 28-11-2017, R.Sudhakar
% equation through paper from Farook Sataar, IEEE transactions on IP vol 6 no 6 June 1997
clc;
clear all;
close all;
% get the original image of your choice
s=imread('cameraman.tif');

% add with the type of noise along with density
ns=imnoise(s,'salt & pepper',0.2);

% a high pass laplacian filter is created
H = fspecial('laplacian',0.5) ;

% noisy image is restored  with median filter of size 3 X 3
y=medfilt2(ns);

% preparing the components of EPI

% input image is highpass filtered with laplacian filter
deltas=imfilter(s,H);
meandeltas=mean2(deltas);

% Restored Image is highpass filtered with laplacian filter
deltascap=imfilter(y,H);
meandeltascap=mean2(deltascap);

% computation of EPI
p1=deltas-meandeltas;
p2=deltascap-meandeltascap;

num=sum(sum(p1.*p2));
den=(sum(sum(p1.^2)))*(sum(sum(p2.^2)));
epi=num/sqrt(den)

% figure(1);
% imshow(s)
% figure(2);
% imshow(ns);
% figure(3);
% imshow(y);
% figure(4);
% imshow(deltas);
% figure(5);
% imshow(deltascap);


