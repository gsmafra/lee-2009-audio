function vishidprod2 = tirbm_vishidprod_fixconv_1d(imdata, H, ws)

numchannels = size(imdata,3);
numbases = size(H,3);

selidx1 = size(H,1):-1:1;
vishidprod2 = zeros(ws,1,numchannels,numbases);

for b=1:numbases
    vishidprod2(:,:,:,b) = conv2_mult2(imdata, H(selidx1, :, b), 'valid');
end

vishidprod2 = reshape(vishidprod2, [ws, numchannels, numbases]);

return
