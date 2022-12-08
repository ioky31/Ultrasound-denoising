function fimgd = OBNLM(img)
    addpath("OBNLMpackage")
    %params
    M = 7;      % search area size (2*M + 1)^2
    alpha = 3;  % patch size (2*alpha + 1)^2
    h = 0.7;    % smoothing parameter [0-infinite].
    % If you can see structures in the "Residual image" decrease this parameter
    % If you can see speckle in the "denoised image" increase this parameter
    
    offset = 100; % to avoid Nan in Pearson divergence computation
    [dimxy dimt] = size(size(img));
    if ( dimt > 2)
        img = rgb2gray(img);
    end
    
    % Intensity normalization
    imgd = double(img);
    mini = (min(imgd(:)));
    imgd = (imgd - mini);
    maxi = max(imgd(:));
    imgd = (imgd / maxi) * 255;
    imgd = imgd + offset; % add offset to enable the pearson divergence computation (i.e. avoid division by zero).
    s = size(imgd);
    
    % Padding
    imgd = padarray(imgd,[alpha alpha],'symmetric');
    fimgd=bnlm2D(imgd,M,alpha,h);
    fimgd = fimgd - offset;
    imgd = imgd - offset;
    imgd = imgd(alpha+1: s(1)+alpha, alpha+1: s(2)+alpha);
    fimgd = fimgd(alpha+1: s(1)+alpha, alpha+1: s(2)+alpha);
    fimgd = (fimgd - min(fimgd(:)))/max(fimgd(:))*255;
end