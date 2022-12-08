function denoised = my_WNNM(noise,nSig)
addpath('WNNM_CVPR2014-master');   
noise = double(noise);
Par   = ParSet(nSig); 
denoised = WNNM_DeNoising(noise, Par);   
end