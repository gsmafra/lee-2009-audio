function val = get_constant(id);

fid = fopen('constants.txt');
while true,
  thisid = fscanf(fid, '%s', 1);
  if strcmp(thisid, ''), break; end;
  thisval = fscanf(fid, '%s', 1);

  if strcmp(thisid, id),
    if thisval(1) >= '0' && thisval(1) <= '9',
      thisval = str2num(thisval);
    end
    val = thisval;
    fclose(fid);
    return;
  end
  

end

fclose(fid);
error('Constant not found: %s\n', id);

