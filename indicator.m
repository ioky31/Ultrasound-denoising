function [ssimval,peaksnr,snr,MSE,epi] = indicator(A, ref)
if class(A)~="double"
    A = double(A);
end
if class(ref)~="double"
    ref = double(ref);
end
[ssimval,~] = ssim(A,ref);
[peaksnr, snr] = psnr(uint8(A), uint8(ref));
MSE = mean((A - ref).^2,'all');
epi = EPI(A, ref);
end