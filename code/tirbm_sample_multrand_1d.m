function [H HP Hc HPc] = tirbm_sample_multrand_1d(poshidexp, spacing)
% poshidexp is 3d array
poshidexp = max(min(poshidexp,20),-20); % DEBUG: This is for preventing NAN values
poshidprobs = exp(poshidexp);
poshidprobs_mult = zeros(spacing+1, size(poshidprobs,1)*size(poshidprobs,2)*size(poshidprobs,3)/spacing);
poshidprobs_mult(end,:) = 1;
% TODO: replace this with more realistic activation, bases..
for r=1:spacing
    temp = poshidprobs(r:spacing:end, :, :);
    poshidprobs_mult(r,:) = temp(:);
end

[S1 P1] = multrand2(poshidprobs_mult');
S = S1';
P = P1';
clear S1 P1

% convert back to original sized matrix
H = zeros(size(poshidexp));
HP = zeros(size(poshidexp));
for r=1:spacing
    H(r:spacing:end, :, :) = reshape(S(r,:), [size(H,1)/spacing, size(H,2), size(H,3)]);
    HP(r:spacing:end, :, :) = reshape(P(r,:), [size(H,1)/spacing, size(H,2), size(H,3)]);
end

if nargout >2
    Sc = sum(S(1:end-1,:));
    Pc = sum(P(1:end-1,:));
    Hc = reshape(Sc, [size(poshidexp,1)/spacing,size(poshidexp,2),size(poshidexp,3)]);
    HPc = reshape(Pc, [size(poshidexp,1)/spacing,size(poshidexp,2),size(poshidexp,3)]);
end

return
