function y = conv2_mult_pairwise(A, B, convopt)
y = [];
assert(size(A,3)==size(B,3));
for i=1:size(A,3)
    y(:,:,i) = conv2(A(:,:,i), B(:,:,i), convopt);
end
return
