function epi = EPI(A,ref)
% a high pass laplacian filter is created
H = fspecial('laplacian',0.5) ;
% input image is highpass filtered with laplacian filter
deltas=imfilter(A,H);
meandeltas=mean2(deltas);

% Restored Image is highpass filtered with laplacian filter
deltascap=imfilter(ref,H);
meandeltascap=mean2(deltascap);

% computation of EPI
p1=deltas-meandeltas;
p2=deltascap-meandeltascap;

num=sum(sum(p1.*p2));
den=(sum(sum(p1.^2)))*(sum(sum(p2.^2)));
epi=num/sqrt(den)

end