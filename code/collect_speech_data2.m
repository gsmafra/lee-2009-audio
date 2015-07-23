function Pall = collect_speech_data2()

Pall = {};
flist= get_speech_filenames2();
parfor i= 1:length(flist)
    if mod(i,10)==0, fprintf('.'); end
    if mod(i,1000)==0, fprintf('\n'); end
    [P, P2, y] = load_spectrogram(flist{i});
    Pall{i} = P2;
end

return

