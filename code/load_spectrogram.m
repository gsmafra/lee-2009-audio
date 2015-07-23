function [P, P2, y] = load_spectrogram(fname_prefix, opt_figure)

if ~exist('opt_figure', 'var')
    opt_figure = false;
end

fname_wav = [fname_prefix '.wav'];

% read wav file
% [y,fs,ffx]=readsph(fname_wav);
[y, fs, nbits] = wavread(fname_wav);

% convert to spectrogram
P = []; 
P2 = get_spectrogram_orig(y, 0, fs);

return
