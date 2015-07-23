function negdata = tirbm_reconstruct_LB_fixconv_1d(S, W, pars)

ws = size(W,1);
patch_M = size(S,1);
numchannels = size(W,2);
numbases = size(W,3);

% Note: Reconstruction was off by a few pixels in the original code (above
% versions).. I fixed this as below:
S2 = S;
negdata2 = zeros(patch_M+ws-1, 1, numchannels);

for b = 1:numbases,
    H = reshape(W(:,:,b),[ws,1,numchannels]);
    negdata2 = negdata2 + conv2_mult(S2(:,:,b), H, 'full');
end

negdata = pars.C_sigm*negdata2;
% imagesc(negdata); colormap gray

return