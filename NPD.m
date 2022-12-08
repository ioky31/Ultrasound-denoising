function img_denoised = NPD(img_noisy)

%
% This is a demo program of the paper J. Tian and L. Chen, "Image
% despeckling using a non-parametric statistical model of wavelet coefficients
% Biomedical Signal Processing and Control, Vol. 6, No. 4, Oct. 2011, pp. 432-437.
%
% Contact: eejtian@gmail.com
%
% Note that the PSNR values could be slightly different with that reported 
% in paper due to the random number generator used to generate noisy image.
%
% The input image should have a square size.
%
% Acknowledgement: This program needs Wavelet transform toolbox, which is 
% downloaded from http://www-stat.stanford.edu/~wavelab/.
addpath('Wave802Tool');   
img_noisy = double(img_noisy);
% wavelet parameters
wbase = 'Daubechies';
mom = 8;
dwt_level = 5;
% [~,J] = func_quadlength(img_truth);
% L = J-dwt_level;

% Perform image classification using ant colony optimization
img_mask = func_aco_classification(img_noisy, wbase, mom, dwt_level);

% Perform image denoising
win_size=3;
img_denoised = func_denoise_kde(img_noisy, img_mask, wbase, mom, dwt_level, win_size);

%--------------------------------------------------------------------------
%-------------------------- Inner Function --------------------------------
%--------------------------------------------------------------------------

function x_mask = func_aco_classification(x_in, wbase, mom, dwt_level)

% [nrow, ncol] = size(x_in);
L = log2(size(x_in,2))-dwt_level;
qmf = func_MakeONFilter(wbase,mom);
wx  = func_FWT2_PO(x_in, L, qmf);
[~,J] = func_dyadlength(wx);

x_mask = zeros(size(wx));
wx = abs(wx);
for j=(J-1):-1:L

    [t1,t2] = func_dyad2HH(j);
    ant_search_range_row_min = ones(size(wx(t1,t2)));
    ant_search_range_row_max = ones(size(wx(t1,t2))).*size(wx(t1,t2),1);
    ant_search_range_col_min = ones(size(wx(t1,t2)));
    ant_search_range_col_max = ones(size(wx(t1,t2))).*size(wx(t1,t2),2);
    x_mask(t1,t2) = func_aco_thresholding_estimate_mask(wx(t1,t2), ant_search_range_row_min,ant_search_range_row_max,ant_search_range_col_min,ant_search_range_col_max);
    [t1,t2] = func_dyad2HL(j);
    ant_search_range_row_min = ones(size(wx(t1,t2)));
    ant_search_range_row_max = ones(size(wx(t1,t2))).*size(wx(t1,t2),1);
    ant_search_range_col_min = ones(size(wx(t1,t2)));
    ant_search_range_col_max = ones(size(wx(t1,t2))).*size(wx(t1,t2),2);
    x_mask(t1,t2) = func_aco_thresholding_estimate_mask(wx(t1,t2), ant_search_range_row_min,ant_search_range_row_max,ant_search_range_col_min,ant_search_range_col_max);
    [t1,t2] = func_dyad2LH(j);
    ant_search_range_row_min = ones(size(wx(t1,t2)));
    ant_search_range_row_max = ones(size(wx(t1,t2))).*size(wx(t1,t2),1);
    ant_search_range_col_min = ones(size(wx(t1,t2)));
    ant_search_range_col_max = ones(size(wx(t1,t2))).*size(wx(t1,t2),2);
    x_mask(t1,t2) = func_aco_thresholding_estimate_mask(wx(t1,t2), ant_search_range_row_min,ant_search_range_row_max,ant_search_range_col_min,ant_search_range_col_max);
end   

%--------------------------------------------------------------------------
%-------------------------- Inner Function --------------------------------
%--------------------------------------------------------------------------

function result = func_aco_thresholding_estimate_mask(wx, ant_search_range_row_min,ant_search_range_row_max,ant_search_range_col_min,ant_search_range_col_max)

% System setup
ant_move_step_within_iteration = 8; % the numbe of iterations?
total_iteration_num = 4;
search_clique_mode = 8;  

wx = abs(wx);
[nrow, ncol] = size(wx);
wx_1D = func_2D_LexicoOrder(wx);

% initialization
h = zeros(size(wx));
h_1D = func_2D_LexicoOrder(h);
p = 1./(wx.*(wx~=0)+1.*(wx==0)).*(wx~=0)+(wx==0);
p_1D = func_2D_LexicoOrder(p);

%paramete setting
alpha = 1;      
beta = 2;       
rho = 0.1;      
w = 0.6;        
A = 5000;       
B = 10;         


%use one ant for each pixel position
ant_total_num = nrow*ncol;
ant_current_row = zeros(ant_total_num, 1); % record the location of ant
ant_current_col = zeros(ant_total_num, 1); % record the location of ant
ant_current_val = zeros(ant_total_num, 1); % record the location of ant
rr = (meshgrid(1:nrow))'; cc = meshgrid(1:ncol);
ant_current_row = rr(:);
ant_current_col = cc(:);
ant_current_val = wx_1D((ant_current_row-1).*ncol+ant_current_col);

ant_search_range_row_min = ant_search_range_row_min(:);
ant_search_range_row_min = padarray(ant_search_range_row_min, [0 search_clique_mode-1],'replicate','post');
ant_search_range_row_max = ant_search_range_row_max(:);
ant_search_range_row_max = padarray(ant_search_range_row_max, [0 search_clique_mode-1],'replicate','post');
ant_search_range_col_min = ant_search_range_col_min(:);
ant_search_range_col_min = padarray(ant_search_range_col_min, [0 search_clique_mode-1],'replicate','post');
ant_search_range_col_max = ant_search_range_col_max(:);
ant_search_range_col_max = padarray(ant_search_range_col_max, [0 search_clique_mode-1],'replicate','post');                

            
for nIteration = 1: total_iteration_num               
%     ant_current_path_val_mean = zeros(ant_total_num,1);
    ant_current_path_val_mean = ant_current_val;            

    for nMoveStep = 1: ant_move_step_within_iteration-1                

        if search_clique_mode == 4
            ant_search_range_row = [ant_current_row-1, ant_current_row, ant_current_row+1, ant_current_row];
            ant_search_range_col = [ant_current_col, ant_current_col+1, ant_current_col, ant_current_col-1];                    
        elseif search_clique_mode == 8
            ant_search_range_row = [ant_current_row-1, ant_current_row-1, ant_current_row-1, ant_current_row, ant_current_row,ant_current_row+1, ant_current_row+1, ant_current_row+1];
            ant_search_range_col = [ant_current_col-1, ant_current_col, ant_current_col+1, ant_current_col-1, ant_current_col+1, ant_current_col-1, ant_current_col, ant_current_col+1];
        end

        ant_current_row_extend = padarray(ant_current_row, [0 search_clique_mode-1],'replicate','post');
        ant_current_col_extend = padarray(ant_current_col, [0 search_clique_mode-1],'replicate','post');
        ant_search_range_val = zeros(ant_total_num,search_clique_mode);

        %replace the positions our of the image's range

        temp = (ant_search_range_row>=ant_search_range_row_min) & (ant_search_range_row<=ant_search_range_row_max) & (ant_search_range_col>=ant_search_range_col_min) & (ant_search_range_col<=ant_search_range_col_max);
        ant_search_range_row = temp.*ant_search_range_row + (~temp).*ant_current_row_extend;
        ant_search_range_col = temp.*ant_search_range_col + (~temp).*ant_current_col_extend;

        ant_search_range_transit_prob_h = zeros(size(ant_search_range_val));
        ant_search_range_transit_prob_p = zeros(size(ant_search_range_val));                                

        for ii=1:search_clique_mode
            ant_search_range_val(:,ii) = wx_1D((ant_search_range_row(:,ii)-1).*ncol+ant_search_range_col(:,ii));

            temp = ant_search_range_val(:,ii);
            temp = abs(temp - ant_current_path_val_mean);    
            temp = 1./(temp.*(temp~=0)+1.*(temp==0)) .* (temp~=0)+(temp==0);

            ant_search_range_transit_prob_h(:,ii) = temp;
            ant_search_range_transit_prob_p(:,ii) = p_1D((ant_search_range_row(:,ii)-1).*ncol+ant_search_range_col(:,ii));

        end
        temp = (ant_search_range_transit_prob_h.^alpha) .* (ant_search_range_transit_prob_p.^beta);

        temp_sum = sum(temp,2);
        temp_sum = padarray(temp_sum, [0 search_clique_mode-1],'replicate','post');

        ant_search_range_transit_prob = temp ./ temp_sum;

        % generate a random number to determine the next position.
        rand('state', sum(100*clock));
        temp = rand(ant_total_num,1);
        temp = padarray(temp, [0 search_clique_mode-1],'replicate','post');
        temp = cumsum(ant_search_range_transit_prob,2)>=temp;
        temp = padarray(temp, [0 1],'pre');
        temp = (diff(temp,1,2)==1);

        temp_row = (ant_search_range_row .* temp)';
        [ii, jj, vv] = find(temp_row);
        ant_next_row = vv;
        temp_col = (ant_search_range_col .* temp)';
        [ii, jj, vv] = find(temp_col);
        ant_next_col = vv;

        ant_current_row = ant_next_row;
        ant_current_col = ant_next_col;                
        ant_current_val = wx_1D((ant_current_row-1).*ncol+ant_current_col);
        ant_current_path_val_mean = (ant_current_path_val_mean.*nMoveStep + ant_current_val)/(nMoveStep+1);

        rr = ant_current_row;
        cc = ant_current_col;
        p_1D((rr-1).*ncol+cc,1) = w*p_1D((rr-1).*ncol+cc,1) + 1./(A+B.*ant_current_path_val_mean);
        p = func_LexicoOrder_2D(p_1D, nrow, ncol);

    end % end of nMoveStep

    p = (1-rho).*p;
    p_1D = func_2D_LexicoOrder(p);

end % end of nIteration

clear h h_1D temp
clear ant_current_row ant_current_col ant_current_val
clear ant_current_path_val_mean
clear ant_search_range_row search_range_path_col ant_search_range_val
clear ant_search_range_transit_prob_h ant_search_range_transit_prob_v ant_search_range_transit_prob
clear ant_search_range_row_min ant_search_range_row_max ant_search_range_col_min ant_search_range_col_max

result = func_determine_class_fcm(wx, p);

% ******************************************************************************
% **************************Inner Function *************************************
% ******************************************************************************
% Bi-class classification algorithm using Fuzzy C-mean
function idx = func_determine_class_fcm(wx, p)

p_1D = func_2D_LexicoOrder(p);
wx_1D = func_2D_LexicoOrder(wx);
[nrow, ncol] = size(wx);

nFeature(:,1) = p_1D(:)./(max(max(p_1D)));
nFeature(:,2) = wx_1D(:)./(max(max(wx_1D)));
[center,U,obj_fcn] = fcm(nFeature, 2,[2.0 100 1e-5 0]);

idx = zeros(nrow:ncol,1);    
if sum(sum(U(1,:) >= U(2,:))) >= sum(sum(U(1,:) < U(2,:)))
    idx = (U(1,:) < U(2,:));
else
    idx = (U(1,:) >= U(2,:));
end
idx = func_LexicoOrder_2D(idx, nrow, ncol);    


%--------------------------------------------------------------------------
%-------------------------- Inner Function --------------------------------
%--------------------------------------------------------------------------
function x_out = func_denoise_kde(x_in, x_mask, wbase, mom, dwt_level, win_size)
[nRow, nColumn] = size(x_in);
L = log2(size(x_in,2))-dwt_level;

%%%%%%%%%%%%%%%%% estimate the noise_sigma from the noisy signal
qmf = func_MakeONFilter(wbase,mom);
[temp, coef] = func_NormNoise_2d(x_in, qmf);
noise_sigma = 1/coef;
wx  = func_FWT2_PO(x_in, L, qmf);
[n,J] = func_dyadlength(wx);
ws = wx;
for j=(J-1):-1:L
    [t1,t2] = func_dyad2HH(j);    
    ws(t1,t2) = func_subband_denoise(wx(t1,t2), x_mask(t1,t2),noise_sigma, win_size);
    [t1,t2] = func_dyad2HL(j);
    ws(t1,t2) = func_subband_denoise(wx(t1,t2), x_mask(t1,t2),noise_sigma, win_size);
    [t1,t2] = func_dyad2LH(j);
    ws(t1,t2) = func_subband_denoise(wx(t1,t2), x_mask(t1,t2),noise_sigma, win_size);
end
x_out  = func_IWT2_PO(ws, L, qmf);

%-------------------------------------------------------------------------
%------------------------------Inner Function ----------------------------
%-------------------------------------------------------------------------
% Denoise for each subband
function result = func_subband_denoise(x_in, x_mask, noise_sigma, win_size)

AA = padarray(x_in, [win_size win_size], 'replicate', 'bot');
AA = im2col(AA, [win_size*2+1 win_size*2+1],'sliding');

x_in_1d = (x_in(:))';
result_1d = zeros(size(AA));

var_gaussian = mean(AA.^2)-noise_sigma^2;
var_gaussian(var_gaussian<0)=0;
sigma_gaussian = sqrt(var_gaussian);
judge = (sigma_gaussian~=0);

BB = mean(AA);
y = x_in_1d(judge);

MM = padarray(x_mask, [win_size win_size], -1, 'bot');
MM = im2col(MM, [win_size*2+1 win_size*2+1],'sliding');
center_class = (x_mask(:))';

result_1d_num = zeros(size(center_class));

for i=1:size(AA,1)
    temp_x = AA(i,:);
    temp_A = temp_x(judge);
    temp_BB = BB(judge);
    temp_sigma = sigma_gaussian(judge);
    temp_result = (y.*temp_sigma.^2+temp_A.*noise_sigma.^2-temp_BB.*noise_sigma.^2) ./ (temp_sigma.^2+noise_sigma.^2);            
    result_1d_temp = zeros(size(result_1d(i,:)));
    result_1d_temp(judge) = temp_result;            
    judge_mask = (center_class==MM(i,:));
    result_1d_temp(~judge_mask) = 0;
    result_1d_num = result_1d_num + double(judge_mask);
    result_1d(i,:) = result_1d_temp;
end
result_1d = sum(result_1d,1)./result_1d_num;

result = reshape(result_1d,size(x_in,1),size(x_in,2));


%-------------------------------------------------------------------------
%------------------------------Inner Function ----------------------------
%-------------------------------------------------------------------------
% Calculate the PSNR performance to two images
function result = func_psnr_gray(f, g)

f = double(f);
g = double(g);
Q=255; MSE=0;
[M,N]=size(f);
h = f - g;
MSE = sum(sum(h.*h));
MSE=MSE/M/N;
result=10*log10(Q*Q/MSE);

%-------------------------------------------------------------------------
%------------------------------Inner Function ----------------------------
%-------------------------------------------------------------------------
% Convert data from 2D format to 1D format
function result = func_2D_LexicoOrder(x)
temp = x';
result = temp(:);

%-------------------------------------------------------------------------
%------------------------------Inner Function ----------------------------
%-------------------------------------------------------------------------
% Convert data from 1D format to 2D format
function result = func_LexicoOrder_2D(x, nRow, nColumn)

result = reshape(x, nColumn, nRow)';