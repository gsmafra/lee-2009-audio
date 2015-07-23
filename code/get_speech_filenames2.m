function flist = get_speech_filenames2()

flist = [];
fid = fopen('./flist_train_alert.txt', 'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    % disp(tline);
    flist{length(flist)+1} = tline;
end
fclose(fid);

return

