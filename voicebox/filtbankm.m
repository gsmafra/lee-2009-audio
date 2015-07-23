function [x,cf,il,ih]=filtbankm(p,n,fs,fl,fh,w)
%FILTBANKM determine matrix for a linear/mel/erb/bark-spaced filterbank [X,MN,MX]=(P,N,FS,FL,FH,W)
%
% Inputs:
%       p   number of filters in filterbank or the filter spacing in k-mel/bark/erb [ceil(4.6*log10(fs))]
%		n   length of fft
%		fs  sample rate in Hz
%		fl  low end of the lowest filter in Hz [default = 0]
%		fh  high end of highest filter in Hz [default = fs/2]
%		w   any sensible combination of the following:
%
%             'b' = bark scale instead of mel
%             'e' = erb-rate scale
%             'l' = log10 Hz frequency scale
%             'f' = linear frequency scale
%             'm' = mel frequency scale
%
%             'c' = fl & fh specify centre of low and high filters instead of edges
%             'h' = fl & fh are in mel/erb/bark/log10 instead of Hz
%             'H' = cf outputs are in mel/erb/bark/log10 instead of Hz
%
%		      'y' = lowest filter remains at 1 down to 0 frequency and
%			        highest filter remains at 1 up to nyquist freqency
%		            The total power in the fft is preserved (unless 'u' is specified).
%
%             'u' = scale filters to sum to unity instead of to preserve power.
%
%             's' = single-sided input: do not include symmetric negative frequencies
%             'S' = single-sided output: do not mirror the non-DC filter characteristics
%
%             'g' = plot filter coefficients as graph
%             'G' = plot filter coefficients as image [default if no output arguments present]
%

%
% Outputs:	x     a sparse matrix containing the filterbank amplitudes
%		          If the mn and mx outputs are given then size(x)=[p,mx-mn+1]
%                 otherwise size(x)=[p,1+floor(n/2)]
%                 Note that the peak filter values equal 2 to account for the power
%                 in the negative FFT frequencies.
%           cf    the filterbank centre frequencies in Hz (see 'H' option)
%		    il    the lowest fft bin with a non-zero coefficient
%		    ih    the highest fft bin with a non-zero coefficient
%
% The routine perorms interpolation of the input spectrum by convolving the power spectrum
% with a triangular filter and then simulates a filterbank with asymetric triangular filters.
%
% Examples of use:
%
% (a) Calcuate the Mel-frequency Cepstral Coefficients
%
%       f=rfft(s);			        % rfft() returns only 1+floor(n/2) coefficients
%		x=filtbankm(p,n,fs,0,fs/2,'m');	        % n is the fft length, p is the number of filters wanted
%		z=log(x*abs(f).^2);         % multiply x by the power spectrum
%		c=dct(z);                   % take the DCT
%
% (b) Calcuate the Mel-frequency Cepstral Coefficients efficiently
%
%       f=fft(s);                        % n is the fft length, p is the number of filters wanted
%       [x,cf,na,nb]=filtbankm(p,n,fs,0,fs/2,'m');   % na:nb gives the fft bins that are needed
%       z=log(x*(f(na:nb)).*conj(f(na:nb)));
%
% (c) Plot the calculated filterbanks
%
%      plot((0:floor(n/2))*fs/n,filtbankm(p,n,fs,0,fs/2,'m')')   % fs=sample frequency
%
% (d) Plot the filterbanks
%
%      filtbankm(p,n,fs,0,fs/2,'m');
%
% References:
%
% [1] S. S. Stevens, J. Volkman, and E. B. Newman. A scale for the measurement
%     of the psychological magnitude of pitch. J. Acoust Soc Amer, 8: 185–19, 1937.
% [2] S. Davis and P. Mermelstein. Comparison of parametric representations for
%     monosyllabic word recognition in continuously spoken sentences.
%     IEEE Trans Acoustics Speech and Signal Processing, 28 (4): 357–366, Aug. 1980.


%      Copyright (C) Mike Brookes 1997-2009
%      Version: $Id: filtbankm.m,v 1.1 2010/01/02 16:28:38 dmb Exp $
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

% Note "FFT bin_0" assumes DC = bin 0 whereas "FFT bin_1" means DC = bin 1

if nargin < 6
    w='tz'; % default options
    if nargin < 5
        fh=0.5; % max freq is the nyquist
        if nargin < 4
            fl=0; % min freq is DC
        end
    end
end
wr=' ';   % default warping is mel
for i=1:length(w)
    if any(w(i)=='lebm');
        wr=w(i);
    end
end
mflh=[fl fh];
if ~any(w=='h')
    switch wr
        case 'm'
            mflh=frq2mel(mflh);       % convert frequency limits into mel
        case 'l'
            if fl<=0
                error('Low frequency limit must be >0 for l option');
            end
            mflh=log10(mflh);       % convert frequency limits into log10 Hz
        case 'e'
            mflh=frq2erb(mflh);       % convert frequency limits into erb-rate
        case 'b'
            mflh=frq2bark(mflh);       % convert frequency limits into bark
    end
end
melrng=mflh*(-1:2:1)';          % mel range
fn2=floor(n/2);     % bin index of highest positive frequency (Nyquist if n is even)
if isempty(p)
    p=ceil(4.6*log10(fs));         % default number of filters
end
if any(w=='c')              % c option: specify fiter centres not edges
    if p<1
        p=round(melrng/(p*1000))+1;
    end
    melinc=melrng/(p-1);
    mflh=mflh+(-1:2:1)*melinc;
else
    if p<1
        p=round(melrng/(p*1000))-1;
    end
    melinc=melrng/(p+1);
end
%
% Calculate the FFT bins0 corresponding to the filters
%
cf=mflh(1)+(0:p+1)*melinc;
switch wr
    case 'l'
        mb=10.^(cf)*n/fs;
    case 'e'
        mb=erb2frq(cf)*n/fs;
    case 'b'
        mb=bark2frq(cf)*n/fs;
    case 'm'
        mb=mel2frq(cf)*n/fs;
    otherwise
        mb=cf*n/fs;
end

% define mel filters (including both ends) in FFT bins0

nyqbin=n-fn2;
dm=mb(2:end)-mb(1:end-1); % lower triangle length of filter0 i

mf=floor(mb);  % low bin0 number for each filter
kl=mf(1);       % lowest bin0 needed to include all of mb
kh=mf(end)+1;   % highest bin0 needed to include all of mb
nk=kh-kl+1;     % number of bins needed to include all of mb
% kln=min(mf(2),kl+1);   % lowest bin0 with non-zero coefficient
% khn=max(mf(end-1)+1,kh-1);  % highest bin0 with non-zero coefficient
kf=cumsum(accumarray((mf-kl+2)',1,[nk 1]))'; % gives next highest filter0 for the bins0 (1:nk)+kl-1
mx=mb-mf;       % fractional part of each filter position


% integrate the product of triangles; each triangle pair is attached to the one
% whose rightmost edge has the lower value

% first take those attached to the filters: for j=1:p ...
% take lower triangle of filter0 j1 with lower triangle of bin0 mf(j1+1)+1
j1=1:p;
i1=mf(j1+1)+1;
dma=dm(j1); % length of lower filter triangle
mxa=mx(j1+1); % distance to next FFT bin downwards
l1=min(dma,mxa);  % integration length
s1=(mxa +(l1/3-(mxa+dma)/2).*l1./dma).*l1;
% now do lower triangle of filter0 i  with upper triangle of previous bin0: mf(i+1)
i5=i1-1;
s5=(1-mxa +((dma-1+mxa)/2-l1/3).*l1./dma).*l1;
% now do upper triangular of filter0 j1 with with the lower triangle of bin0 mf(j1+2)+1
i2=mf(j1+2)+1;
dmb=dm(j1+1); % length of upper filter triangle
mxb=mx(j1+2); % tail length of bin0 i2 lower triangle
l2=min(dmb,mxb); % integration length
s2=(mxb/2-l2/3).*l2.^2./dmb;
% now do upper triangle of filter0 i  with same upper triangle of bin0: mf(i+1)
i6=i2-1;
s6=((1-mxb)/2+l2/3).*l2.^2./dmb;

% now do those attached to the bins: for i= bin0 number
% take lower triangle of bin0 kl+i with lower triangle of filter0 kf(i-kl+1)
i3=kl+1:mf(p+1);           % range of bin0 numbers needed
j3=kf(i3-kl+1);
dmc=dm(j3);            % total length of filter triangle
l3f=i3-mb(kf(i3-kl+1)); % available length of filter triangle
l3=min(1,l3f);   % integration length
s3=(l3f+(l3/3-(1+l3f)/2).*l3).*l3./dmc;
% take upper triangle of previous bin0 kl+i-1 with same lower triangle of filter0 kf(i-kl+1)
i4=i3-1;
s4=(l3f/2-l3/3).*l3.^2./dmc;
% Note s3+s4 = (l3f-l3/2).*l3./dmc
% take lower triangle of bin0 kl+i with upper triangle of previous filter0 kf(i-kl+1)-1
i7=mf(2)+1:kh-1;   % range of bin0 numbers needed
j7=kf(i7-kl+1)-1;   % corresponding filter0 numbers
dmd=dm(j7+1);            % total length of filter upper triangle
l7f=i7-mb(j7+1); % available length of filter triangle below bin0 i
l7=min(l7f,1);
s7= (dmd-l7f+((1-dmd+l7f)/2-l7/3).*l7).*l7./dmd;
% take upper triangle of previous bin0 kl+i with same upper triangle of filter0 kf(i-kl+1)-1
i8=i7-1;            % previous bin0
s8=((dmd-l7f)/2+l7/3).*l7.^2./dmd;
% Note s7+s8 = (dmd-l7f+l7/2).*l7./dmd


r=[j1 j1 j3 j3 j1 j1 j7 j7];
c=[i1 i5 i3 i4 i2 i6 i7 i8];
v=[s1 s5 s3 s4 s2 s6 s7 s8];
v((abs(c)>fn2) | (c+nyqbin<=0))=0; % kill entries outside the Nyquist range

% ss=[repmat(1,1,p)  repmat(5,1,p) repmat(3,1,numel(j3))  repmat(4,1,numel(j3)) ...
%     repmat(2,1,p) repmat(6,1,p) repmat(7,1,numel(j7)) repmat(8,1,numel(j7))];
if any(w=='y')
    % add in FFT bins outside the range we have considered
    nbl=2*(p+numel(j3));
    v(find(r(1:nbl)==1))=0;  % delete entries using left triangle of filter 1
    r=[r repmat(1,1,mf(2)+2)];
    c=[c mf(2)+1 mf(2) 0:mf(2)-1];
    v=[v mx(2)^2/2 1-(1-mx(2))^2/2 ones(1,mf(2))];
    %     ss=[ss repmat(10,1,mf(2)+2)];
    v(nbl+find(r(nbl+1:end)==p))=0;  % delete entries using right triangle of filter p
    if fn2>=mf(p+1)
        r=[r p];
        c=[c mf(p+1)];
        v=[v (1-mx(p+1))^2/2];
        %         ss=[ss 11];
        if fn2>mf(p+1)
            r=[r repmat(p,1,fn2-mf(p+1))];
            c=[c mf(p+1)+1:fn2];
            v=[v 1-mx(p+1)^2/2 ones(1,fn2-mf(p+1)-1)];
            %             ss=[ss repmat(11,1,fn2-mf(p+1))];
        end
    end
end
ss=2*any(w=='s')+any(w=='S');
dcfilt=find(4*abs(mb(3:p+2)+mb(1:p))-mb(3:p+2)+mb(1:p)<0,1);     % find the DC filter bin, if any
if numel(dcfilt)
    switch ss
        case 0
            v(r~=dcfilt)=2*v(r~=dcfilt);      % double all outputs except the DC filter
        case 1
            v((r==dcfilt) & (c<0))=0;           % zero negative frequency contributions to DC filter
        case 3
            v(c<0)=0;           % zero negative frequency components
    end
else
    switch ss
        case 0
            v=2*v;                          % double all outputs since no DC filter exists
        case 3
            v(c<0)=0;           % zero negative frequency components
    end
end
c=abs(c);                           % reflect negative frequencies
msk=(v==0);     % kill all entries that are null
r(msk)=[];
c(msk)=[];
v(msk)=[];

%
% sort out the output argument options
%
if ~any(w=='H')
    cf=mb*fs/n;         % output Hz instead of mel/erb/...
end
cf=cf(2:p+1);
if nargout > 2
    il=min(c)+1;        % find range of FFT bins used (+1 because DC term is at 1)
    ih=max(c)+1;
else
    il=1;
    ih=1+fn2;
end
x=sparse(r,c-il+2,v,p,ih-il+1);
if any(w=='u')
    sx=sum(x,2);
    x=x./repmat(sx+(sx==0),1,size(x,2));
end
%
% plot results if no output arguments or g option given
%
if ~nargout || any(w=='g') || any(w=='G') % plot idealized filters
    if ~any(w=='g') && ~any(w=='G')
        w=[w 'G'];
    end
    newfig=0;
    if  any(w=='g')
        plot((il:ih)*fs/n,x');
        title(['filtbankm: mode = ' w]);
        xlabel(['Frequency (' xticksi 'Hz)']);
        ylabel('Weight');
        newfig=1;
    end

    if  any(w=='G')
        if newfig
            figure;
        end
        imagesc((il:ih)*fs/n,1:p,x);
        axis 'xy'
        colorbar;
        cblabel('Weight');
        switch wr
    case 'l'
        type='Log-spaced';
    case 'e'
        type='Erb-spaced';
    case 'b'
        type='Bark-spaced';
    case 'm'
        type='Mel-spaced';
    otherwise
        type='Linear-spaced';
end
        ylabel([type ' Filter']);
        xlabel(['Frequency (' xticksi 'Hz)']);
        title(['filtbankm: mode = ' w]);
    end

end