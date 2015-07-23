function [ferr dW_total dh_total dv_total poshidprobs poshidstates negdata stat] = ...
    fobj_tirbm_CD_LB_sparse_audio(imdata, W, hbias_vec, vbias_vec, pars, CD_mode, bias_mode, spacing, l2reg)

ws = size(W,1);

poshidexp = tirbm_inference_fixconv_1d(imdata, W, hbias_vec, pars);
[poshidstates poshidprobs] = tirbm_sample_multrand_1d(poshidexp, spacing);
if strcmp(CD_mode, 'mf'), poshidstates = poshidprobs; end
    
posprods = tirbm_vishidprod_fixconv_1d(imdata, poshidprobs, ws);
poshidact = squeeze(sum(sum(poshidprobs,1),2));
posvisact = squeeze(sum(sum(imdata,1),2));


neghidstates = poshidstates;
for j=1:pars.K_CD  
    negdata = tirbm_reconstruct_LB_fixconv_1d(neghidstates, W, pars);
    neghidexp = tirbm_inference_fixconv_1d(negdata, W, hbias_vec, pars);
    [neghidstates neghidprobs] = tirbm_sample_multrand_1d(neghidexp, spacing);
    if strcmp(CD_mode, 'mf'), neghidstates = neghidprobs; end
end
negprods = tirbm_vishidprod_fixconv_1d(negdata, neghidprobs, ws);
neghidact = squeeze(sum(sum(neghidprobs,1),2));
negvisact = squeeze(sum(sum(negdata,1),2));

ferr = mean( (imdata(:)-negdata(:)).^2 );


if strcmp(bias_mode, 'none')
    dhbias = 0;
    dvbias = 0;
    dW = 0;
elseif strcmp(bias_mode, 'simple')
    dhbias = squeeze(mean(mean(poshidprobs,1),2)) - pars.pbias;
    dvbias = 0;
    dW = 0;
else 
    error('wrong adjust_bias mode!');
end

numcases1 = size(poshidprobs,1)*size(poshidprobs,2);
numcases2 = size(imdata,1)*size(imdata,2);

dW_total1 = (posprods-negprods)/numcases1;
dW_total2 = - l2reg*W;
dW_total3 = - pars.pbias_lambda*dW;
dW_total = dW_total1 + dW_total2 + dW_total3;

stat = [];
stat.dWnorm_CD = norm(vec(dW_total1));
stat.dWnorm_l2 = norm(vec(dW_total2));


dh_total = (poshidact-neghidact)/numcases1 - pars.pbias_lambda*dhbias;
dv_total = (posvisact-negvisact)/numcases2;

% fprintf('||W||=%g, ||dWprod|| = %g, ||dWl2|| = %g, ||dWsparse|| = %g\n', sqrt(sum(W(:).^2)), sqrt(sum(dW_total1(:).^2)), sqrt(sum(dW_total2(:).^2)), sqrt(sum(dW_total3(:).^2)));

return
