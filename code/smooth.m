function X = smooth(X, basis_len, sharpness);

% SMOOTH smooths a signal so that it decays smoothly to
% zero along the first dimension. 

num_channels = size(X,3);

for k = 1:num_channels,
  patch_len = size(X,1);
  num_freqs = size(X,2);
  [C,R] = meshgrid(1:num_freqs, 1:patch_len);
  X(:,:,k) = X(:,:,k) ./ (1+exp(-(sharpness*2)*(R/basis_len)+sharpness));
  X(:,:,k) = X(:,:,k) ./ (1+exp(-(sharpness*2)*((patch_len+1-R)/basis_len)+sharpness));
end

