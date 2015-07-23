function [poshidexp2 poshidprobs2] = tirbm_inference_fixconv_1d(imdata, W, hbias_vec, pars)

ws = size(W,1);
numbases = size(W,3);
numchannel = size(W,2);

poshidexp2 = zeros(size(imdata,1)-ws+1, 1, numbases);
if nargout>1
    poshidprobs2 = zeros(size(poshidprobs2));
end

% poshidexp2 = zeros(size(imdata,1)-ws+1, size(imdata,2)-ws+1, numbases);
for c=1:numchannel
    H = reshape(W(end:-1:1, c, :),[ws,1,numbases]);
    poshidexp2 = poshidexp2 + reshape(conv2_mult(imdata(:,:,c), H, 'valid'), size(poshidexp2));
end

for b=1:numbases
    poshidexp2(:,:,b) = pars.C_sigm/(pars.std_gaussian^2).*(poshidexp2(:,:,b) + hbias_vec(b));
    if nargout>1
        poshidprobs2(:,:,b) = 1./(1 + exp(-poshidexp2(:,:,b)));
    end
end

return