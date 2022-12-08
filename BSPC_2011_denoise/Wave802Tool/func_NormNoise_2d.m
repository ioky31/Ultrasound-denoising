function [y,coef] = func_NormNoise_2d(x,qmf)
% NormNoise2 -- Estimates noise level, Normalize signal to noise level 1
%  Usage
%    [y,coef] = NormNoise2(x,qmf)
%  Inputs
%    x     2-d signal
%    qmf   quadrature mirror filter
%  Outputs
%    y     2-d signal, scaled so wavelet coefficients
%          at finest level have median absolute deviation 1.
%    coef  estimation of 1/sigma
%
%  Description
%    This is required pre-processing to use any of the DeNoising
%    tools on naturally-occurring data.
%
%  See Also
%    WaveShrink, CPDeNoise,WPDeNoise
%

	[n,J] = func_quadlength(x);
	wc = x; 
	nc = n;
	for jscal=J-1,
		top = (nc/2+1):nc; bot = 1:(nc/2);
		for ix=1:nc,
			row = wc(ix,1:nc);
			wc(ix,bot) = func_DownDyadLo(row,qmf);
			wc(ix,top) = func_DownDyadHi(row,qmf);
		end
		for iy=1:nc,
			row = wc(1:nc,iy)';
			wc(top,iy) = func_DownDyadHi(row,qmf)';
			wc(bot,iy) = func_DownDyadLo(row,qmf)'; 
		 end
		nc = nc/2;
	end   
    
    [t1,t2] = func_dyad2HH(jscal);
    
    temp = wc(t1,t2);
        
        
	s = median(abs(temp(:)));
	if s ~= 0
		y =  .6745 .* x ./s;
                coef = .6745 /s;
	else
		y = x;
                coef = 1;
	end
    