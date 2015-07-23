function fidx=writewav(d,fs,filename,mode,nskip)
%WRITEWAV Creates .WAV format sound files FIDX=(D,FS,FILENAME,MODE,NSKIP)
%
%   The input arguments for WRITEWAV are as follows:
%
%       D           The sampled data to save
%       FS          The rate at which the data was sampled
%       FILENAME    A string containing the name of the .WAV file to create or
%                        alternatively the FIDX output from a previous writewav call
%       MODE        String containing any reasonable mixture of the following (*=default):
%
%  Precision: 'a'    for 8-bit A-law PCM
%             'u'    for 8-bit mu-law PCM
%            '16' *	for 16 bit PCM data
%             '8'    for 8 bit PCM data
%             ...    any number in the range 2 to 32 for PCM
%	  Dither: 'w'    White triangular dither of amplitude +-1 LSB (PCM modes only)
%             'h'    High pass dither (filtered by 1-1/z) (PCM modes only)
%             'l'    Low pass dither (filtered by 1+1/z) (PCM modes only)
%    Scaling: 's' *  Auto scale to make data peak = +-1
%             'r'    Raw unscaled data (integer values)
%             'q'    Scaled to make 0dBm0 be unity mean square
%             'p'  	Scaled to make +-1 equal full scale
%             'o'    Scale to bin centre rather than bin edge (e.g. 127 rather than 127.5 for 8 bit values)
%                     (can be combined with n+p,r,s modes)
%             'n'    Scale to negative peak rather than positive peak (e.g. 128.5 rather than 127.5 for 8 bit values)
%                     (can be combined with o+p,r,s modes)
%     Offset: 'y' *	Correct for offset in <=8 bit PCM data
%             'z'    No offset correction
%   File I/O: 'f'    Do not close file on exit
%             'd'    Look in data directory: voicebox('dir_data')
%
%        NSKIP      Number of samples to skip before writing or -1[default] to continue from previous write
%                   Only valid if FIDX is specified for FILENAME 
%               
% Output Parameter:
%
%	FIDX     Information row vector containing the element listed below.
%
%           (1)  file id
%			(2)  current position in file (in samples, 0=start of file)
%           (3)  dataoff	length of file header in bytes
%           (4)  nsamp	number of samples
%           (5)  nchan	number of channels
%           (6)  nbyte	bytes per data value
%           (7)  bits	number of bits of precision
%           (8)  code	Data format: 1=PCM, 2=ADPCM, 6=A-law, 7=Mu-law
%           (9)  fs	sample frequency
%           (10) dither state variable
%
%   Note: WRITEWAV will create an 16-bit PCM, auto-scaled wave file by default.
%   For stereo data, d(:,1) is the left channel and d(:,2) the right
%
%   See also READWAV

%   *** Note on scaling ***
%   If we want to scale signal values in the range +-1 to an integer in the
%   range [-128,127] then we have four plausible choices corresponding to
%   scale factors of (a) 127, (b) 127.5, (c) 128 or (d) 128.5 but each choice
%   has disadvantages. 
%   For forward scaling: (c) and (d) cause clipping on inputs of +1.
%   For reverse scaling: (a) and (b) can generate output values < -1.
%   Any of these scalings can be selected via the mode input: (a) 'o', (b) default, (c) 'on', (d) 'n'

%	   Copyright (C) Mike Brookes 1998
%      Version: $Id: writewav.m,v 1.4 2007/05/04 07:01:39 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Acknowledgements
%   Thanks to Hugh Barnes for sorting out seek problems with MATLAB 6.5

if nargin<3 error('Usage: WRITEWAV(data,fs,filename,mode,nskip)'); end

info=zeros(1,9);
info(9)=fs;
if nargin<4 mode='s';
else mode = [mode(:).' 's'];
end
info(8)=1;     % default mode is PCM
mno=all(mode~='o');                      % scale to input limits not output limits
k=find((mode>='0') & (mode<='9'));
if k, info(7)=str2num(mode(k));
else info(7)=16;
end
info(6)=ceil(info(7)/8);
lo=-pow2(0.5,info(7));
hi=-1-lo;
pk=pow2(0.5,8*info(6))*(1-(mno/2-all(mode~='n'))/lo);  % use modes o and n to determine effective peak
if any(mode=='a') info(8)=6; pk=4032+mno*64; info(7)=8; end
if any(mode=='u') info(8)=7; pk=8031+mno*128; info(7)=8; end			% is this pk value correct ?
k=find((mode>='p') & (mode<='s'));
sc=mode(k(1)); 
z=128*all(mode~='z');
if any(mode=='w') di='w';                       % select dither mode
elseif any(mode=='h') di='h';
elseif any(mode=='l') di='l';
else di='n';
end


[n,nc]=size(d);
if n==1 n=nc; nc=1;
else d = d.';
end;
if nc>10 error('WRITEWAV: attempt to write a sound file with >10 channels'); end
nc=max(nc,1);
ncy=nc*info(6);                     % bytes per sample time
nyd=n*ncy;                          % bytes to write

if ischar(filename)
    if any(mode=='d')
        filename=fullfile(voicebox('dir_data'),filename);
    end
    ny=nyd;
    if isempty(findstr(filename,'.')) filename=[filename,'.wav']; end
    fid=fopen(filename,'wb+','l');
    if fid == -1 error(sprintf('Can''t open %s for output',filename)); end
    info(1)=fid;
    fwrite(fid,'RIFF','uchar');
    fwrite(fid,36+ny,'ulong');
    fwrite(fid,'WAVEfmt ','uchar');
    fwrite(fid,[16 0 info(8) nc],'ushort');
    fwrite(fid,[fs fs*ncy],'ulong');
    fwrite(fid,[ncy info(7)],'ushort');
    fwrite(fid,'data','uchar');
    fwrite(fid,ny,'ulong');
    nskip=0;
    info(3)=44;
    info(4)=n;
    info(2)=n;
    info(10)=rand(1);                       % seed for dither generation
else
    info=filename;
    fid=info(1);
    fseek(fid,0,1); % go to end of file
    if nargin<5 nskip=info(2); 
    elseif nskip<0 nskip=info(2);
    end
    info(2)=n+nskip;                         % file position following this write operation (in samples)
    ny=nyd+nskip*ncy;                        % file position following this write operation (in bytes following header)
    if n & (info(2)>info(4))                 % update high water mark
        if ~info(4)                           % if no data written previously
            fseek(fid,22,-1); fwrite(fid,nc,'ushort'); 
            fseek(fid,28,-1); fwrite(fid,fs*ncy,'ulong');
            fwrite(fid,ncy,'ushort');
        end
        fseek(fid,4,-1); fwrite(fid,36+ny,'ulong');
        fseek(fid,40,-1); fwrite(fid,ny,'ulong');
        info(4)=info(2);
    end
end
info(5)=nc;

if n
    if fseek(fid,0,-1) error(sprintf('Cannot rewind file')); end % MATLAB V6.5 fails if this is omitted
    if fseek(fid,info(3)+nskip*nc*info(6),-1) error(sprintf('Cannot seek to byte %d in output file',info(3)+nskip*nc*info(6))); end
    if sc~='r'
        if sc=='s' pd=max(abs(d(:))); pd=pd+(pd==0);
        elseif sc=='p' pd=1;
        else 
            if info(8)==7
                pd=2.03761563;
            else
                pd=2.03033976;
            end
        end
        d=pk/pd*d;
    end
    if info(8)<6
        if di=='n'
            d=round(d);
        else
            [d,info(10)]=ditherq(d,di,info(10));
        end
        d=min(max(d,lo),hi)*pow2(1,8*info(6)-info(7));       % clip data and shift to most significant bits
    else
        z=0;
        if info(8) < 7
            d=lin2pcma(d,213,1);
        else
            d=lin2pcmu(d,1);
        end
    end
    
    
    if info(6)<3
        if info(6)<2
            fwrite(fid,d+z,'uchar');
        else
            fwrite(fid,d,'short');
        end
    else
        if info(6)<4
            d=d(:)';
            d2=floor(d/65536);
            d=d-65536*d2;
            fwrite(fid,[rem(d,256); floor(d/256); d2],'uchar');
        else
            fwrite(fid,d,'long');
        end
    end
    if rem(ny,2) fwrite(fid,0,'uchar'); end
end
if all(mode~='f') fclose(fid); end
if nargout fidx=info; end
