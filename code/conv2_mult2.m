function y = conv2_mult(A, b, convopt)
y = [];
for i=1:size(A,3)
    y(:,:,i) = conv2(A(:,:,i), b, convopt);
end
return
