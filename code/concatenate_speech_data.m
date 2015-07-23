% octave

function [Pconc startframe_list] = concatenate_speech_data(Pall, idx)

count = 0;
for i=1:length(idx)
    count = count + size(Pall{idx(i)},2);
end

numfeat = size(Pall{1},1);
Pconc = zeros(numfeat, count); 
startframe_list = [];

count = 0;
for i=1:length(idx)
    startframe_list(i) = count+1;
    Pconc(:, count+1:count+size(Pall{idx(i)},2)) = Pall{idx(i)};
    count = count + size(Pall{idx(i)},2);
end

return

