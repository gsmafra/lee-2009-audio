function y = conv2_mult(a, B, convopt)
y = [];
for i=1:size(B,3)
    y(:,:,i) = conv2(a, B(:,:,i), convopt);
end
return
