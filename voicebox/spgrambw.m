function [t,f,b]=spgrambw(s,fs,varargin)
%SPGRAMBW Draw spectrogram [T,F,B]=(s,fs,mode,bw,fmax,db,tinc,ann)
%
%  Inputs:  s         speech signal
%           FS        sample fequency (Hz)
%           MODE      optional character string specifying options (see list below)
%           BW        bandwidth resolution in Hz [default: 200]
%           FMAX      frequency range [Fmin Fstep Fmax]. If Fstep is omitted
%                     it is taken to be (Fmax-Fmin)/257, if Fmin is also omitted it is taken
%                     to be 0 (or 20Hz for mode l), if all three are omitted Fmax is taken to be FS/2.
%                     If modes m,b or e are specified then the units are in mel, bark or erb.
%                     If mode l is specified, then the units are in log10(Hz).
%           DB        either dB-range or [dB-min dB-max]
%           TINC      frame increment in seconds [0 or missing uses default=0.45/BW]
%                     or [increment offset] where sample s(n) is at offset + n/fs
%           ANN       annotation cell array: each row contains either
%                     {time 'text-string' 'font'} or {[t_start t_end] 'text-string' 'font'} where
%                     the time value is in seconds with s(n) at time offset+n/fs. The font column can
%                     omitted in which case the system font will be used. MATLAB cannot cope with
%                     so I recommend the SILDoulosIPA (serifed) or SILSophiaIPA (sans) fonts for
%                     phonetic symbols; these are now a little hard to find.
%
% Outputs:  T(NT)        time axis values (in seconds). Input sample s(n) is at time offset+n/fs.
%           F(NF)        frequency axis values in Hz or, unless mode=H, other selected frequency units
%                        according to mode: m=mel, l=log10(Hz), b=bark,e=erb-rate
%           B(NF,NT)     spectrogram values in power (or clipped dB values if 'd' option given)
%
% MODE:  'p' = show power per decade rather than power per Hz [preemphasis]
%        'P' = show power per mel/bark/erb according to y axis scaling
%        'd' = output B array is in dB rather than power
%
%        'm' = mel scale
%        'b' = bark scale
%        'e' = erb scale
%        'l' = log10 Hz frequency scale
%
%        'h' = units of the FMAX input are in Hz instead of mel/bark
%              [in this case, the Fstep parameter is used only to determine
%               the number of filters]
%        'H' = express the F output in Hz instead of mel/bark/...
%
%        'g' = draw a graph even if output arguments are present
%        'j' = jet colourmap
%        'i' = inverted colourmap (white background)
%        'c' = include a colourbar as an intensity scale
%        'w' = draw the speech waveform above the spectrogram
%        'a' = centre-align annotations rather than left-aligning them
%        't' = add time markers with annotations
%
% The BW input gives the 6dB bandwidth of the Hamming window used in the analysis.
% Equal amplitude frequency components are guaranteed to give separate peaks if they
% are this far apart. This value also determines the time resolution: the window length is
% 1.81/BW and the low-pass filter applied to amplitude modulations has a 6-dB bandwidth of
% BW/2 Hz.
%
% The units are power per Hz unless the u
% option is given in which case power per displayed unit is used
% or power per decade for the l option.

%%%% BUGS %%%%%%
% allow ANN rows to be a mixture of intervals and instants
% allow multiple ANN rows
% future options:
%       ['q' = constant q transform]
%       ['k' = add a piano keyboard to the frequency scale]
%       ['f' = label frequency axis in Hz rather than mel/bark/...]

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: spgrambw.m,v 1.11 2010/05/21 09:14:15 dmb Exp $
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

if nargin<2
    error('Usage: SPGRAMBW(s,fs,mode,bw,fmax,db,tinc)');
end
%SPGRAMBW Draw grey-scale spectrogram [T,F,B]=(s,fs,mode,bw,fmax,db,tinc)
%
% first decode the input arguments
%
ap=zeros(1,6);
j=2;
for i=1:length(varargin)
    if ischar(varargin{i})
        ap(1)=i;
    else
        ap(j)=i;
        j=j+1;
    end
end
if ap(1) && ~isempty(varargin{ap(1)})
    mode=varargin{ap(1)};
else
    mode='';  % default mode
end
if ap(2) && ~isempty(varargin{ap(2)})
    bw=varargin{ap(2)};
else
    bw=200;
end
if ap(3) && ~isempty(varargin{ap(3)})
    fmax=varargin{ap(3)};
else
    fmax=[];
end
if ap(4) && ~isempty(varargin{ap(4)})
    db=varargin{ap(4)};
else
    db=80;
end
if ap(5)>0 && ~isempty(varargin{ap(5)})
    tinc=varargin{ap(5)};
    if numel(tinc)<2
        tinc(2)=0;
    end
else
    tinc=zeros(1,2);
end
if tinc(1)<=0
    tinc(1)=1.81/(4*bw); % default frame increment
end
if ap(6)
    ann=varargin{ap(6)};
else
    ann=[];
end

% now sort out the mode flags

mdsw='  ';           % [yscale preemph]
for i=1:length(mode)
    switch mode(i)
        case {'l','m','b','e'}
            mdsw(1)=mode(i);
        case {'p','P'}
            mdsw(2)=mode(i);
    end
end
if mdsw(2)=='P'
    mdsw(2)=mdsw(1);        % preemphasis is scaling dependent
end
%
% sort out the frequency axis
%
flmin=30;                   % min frequency for l option
nfrq=257;                   % default number of frequency bins
if ~numel(fmax)             % no explicit frequency range
    switch mdsw(1)
        case 'l'
            fx=linspace(log10(flmin),log10(fs/2),nfrq);   % 20  Hz to Nyquist
        case 'm'
            fx=linspace(0,frq2mel(fs/2),nfrq);   % DC to Nyquist
        case 'b'
            fx=linspace(0,frq2bark(fs/2),nfrq);   % DC to Nyquist
        case 'e'
            fx=linspace(0,frq2erb(fs/2),nfrq);   % DC to Nyquist
        otherwise   % linear Hz scale
            fx=(0:nfrq-1)*0.5*fs/(nfrq-1);
    end
else
    if any(mode=='h')
        switch mdsw(1)
            case 'l'
                fmaxu=log10(fmax);   % 20  Hz to Nyquist
            case 'm'
                fmaxu=frq2mel(fmax);   % DC to Nyquist
            case 'b'
                fmaxu=frq2bark(fmax);   % DC to Nyquist
            case 'e'
                fmaxu=frq2erb(fmax);   % DC to Nyquist
        end
    else
        fmaxu=fmax;                 % already in the correct units
    end
    if numel(fmax)<2   % only max value specified
        if mdsw(1)=='l'
            fx=linspace(log10(flmin),fmaxu,nfrq);   % 20  Hz to fmax
        else
            fx=linspace(0,fmaxu,nfrq);   % DC to fmax
        end
    elseif numel(fmax)<3 % min and max values specified
        fx=linspace(fmaxu(1),fmaxu(2),nfrq);   % fmin to fmax
    else
        fmaxu(2)=fmax(2)*(fmaxu(3)-fmaxu(1))/(fmax(3)-fmax(1)); % scale the step size appropriately
        fx=fmaxu(1):fmaxu(2):fmaxu(3);   % fmin to fmax in steps of finc
        nfrq=length(fx);
    end
end
switch mdsw(1)          % convert the frequency range to Hz
    case 'l'
        f=10.^fx;
        frlab='log_{10}Hz';
    case 'm'
        f=mel2frq(fx);
        frlab='Mel';
    case 'b'
        f=bark2frq(fx);
        frlab='Bark';
    case 'e'
        f=erb2frq(fx);
        frlab='Erb-rate';
    otherwise
        f=fx;
        frlab='Hz';
end
if ~any(mode=='H')
    f=fx;               % output frequencies in native units unless 'H' is specified
end
%
% now calculate the spectrogram
%
winlen = fix(1.81*fs/bw);   % window length
win=0.54+0.46*cos((1-winlen:2:winlen)*pi/winlen);  % Hamming window
ninc=max(round(tinc(1)*fs),1);                 % window increment in samples
%  we nee to take account of minimum freq increment + make it exact if possible
fftlen=pow2(nextpow2(4*winlen));        % enough oversampling to get good interpolation
win=win/sqrt(sum(win.^2));       % ensure window squared sums to unity
[sf,t]=enframe(s,win,ninc);
t=tinc(2)+t/fs;                         % time axis
nfr=numel(t);                   % number of frames
b=rfft(sf,fftlen,2);
b=b.*conj(b)*2/fs;          % Power per Hz
b(:,1)=b(:,1)*0.5;   % correct for no negative zero frequency to double the power
b(:,end)=b(:,end)*0.5;   % correct for no negative nyquist frequency to double the power
fb=(0:fftlen/2)*fs/fftlen; % fft bin frequencies
dblab='Power/Hz';
switch mdsw(2)
    case {'p','l'}
        b=b.*repmat(fb*log(10),nfr,1);       % convert to power per decade
        dblab='Power/Decade';
    case 'm'
        b=b.*repmat((1+fb/700)*log(1+1000/700)/1000,nfr,1);       % convert to power per mel
        dblab='Power/Mel';
    case 'b'
        b=b.*repmat((1960+fb).^2/52547.6,nfr,1);       % convert to power per bark
        dblab='Power/Bark';
    case 'e'
        b=b.*repmat(6.23*fb.^2 + 93.39*fb + 28.52,nfr,1);       % convert to power per erb
        dblab='Power/Erb-rate';
end
%
% Now map onto the desired frequency scale
%
b=b*filtbankm(nfrq,fftlen,fs,fx(1),fx(end),['cush' mdsw(1)])';

if ~nargout || any(mode=='g') ||  any(mode=='d')
    if numel(db)<2          % find clipping limits
        plim=max(b(:))*[0.1^(0.05*db) 1];
    else
        plim=10.^(0.05*db(1:2));
    end
    bd=10*log10(min(max(b,plim(1)),plim(2)));
    if any(mode=='d')
        b=bd;       % output the dB version
    end
end
% now plot things
if ~nargout || any(mode=='g')
    imagesc(t,fx,bd');
    axis('xy');
    set(gca,'tickdir','out');
    if any(mode=='j')
        colormap('jet');
        map=colormap;
    else
        map = repmat((0:63)'/63,1,3);
    end
    if any(mode=='i')               % 'i' option = invert the colourmap
        map=map(64:-1:1,:);
    end
    colormap(map);
    if any(mode=='c')                % 'c' option = show a colourbar
        colorbar;
        cblabel([dblab ' (dB)']);
    end
    %
    % Now check if annotations or a waveform are required
    %
    dotaw=[((any(mode=='t') && size(ann,2)>1) || size(ann,2)==1) size(ann,2)>1 any(mode=='w')];
    if  any(dotaw)
        ylim=get(gca,'ylim');
        yrange = ylim(2)-ylim(1);
        zlim=ylim;
        toptaw=cumsum([0 dotaw.*[0.05 0.05 0.1]]*yrange)+ylim(2);
        zlim(2)=toptaw(4);
        set(gca,'ylim',zlim,'color',map(1,:));
        if dotaw(3)         % Plot the waveform
            smax=max(s(:));
            smin=min(s(:));
            srange=smax-smin;
            hold on
            plot(tinc(2)+(1:length(s))/fs,(s-smin)/srange*0.9*(toptaw(4)-toptaw(3))+toptaw(3),'color',map(48,:))
            hold off
        end
        if dotaw(1) || dotaw(2)
            tmk=cell2mat(ann(:,1));
            tmksel=tmk(:,1)<=t(end) & tmk(:,end)>=t(1);
            yix=1+[tmk(tmksel,1)<t(1) ones(sum(tmksel),2) tmk(tmksel,end)>t(end)]';
            tmk(:,1)=max(tmk(:,1),t(1));  % clip to axis limits
            tmk(:,end)=min(tmk(:,end),t(end));
        end
        if dotaw(1) && any(tmksel)  % draw time markers
            ymk=toptaw(1:2)*[0.8 0.4;0.2 0.6];
            switch size(tmk,2)
                case 0
                case 1      % isolated marks
                    hold on
                    plot([tmk(tmksel) tmk(tmksel)]',repmat(ymk',1,sum(tmksel)),'color',map(48,:));
                    hold off
                otherwise % draw durations

                    hold on
                    plot(tmk(tmksel,[1 1 2 2])',ymk(yix),'color',map(48,:));
                    hold off
            end
        end
        if dotaw(2) && any(tmksel) % print annotations
            if any(mode=='a')
                horal='center';
                tmk=(tmk(:,1)+tmk(:,end))*0.5;
            else
                horal='left';
                tmk=tmk(:,1);
            end
            if size(ann,2)>2
                font='Arial';
                for i=1:size(ann,1)
                    if tmksel(i)
                        if ~isempty(ann{i,3})
                            font = ann{i,3};
                        end
                        text(tmk(i),toptaw(2),ann{i,2},'color',map(48,:),'fontname',font,'VerticalAlignment','baseline','HorizontalAlignment',horal);
                    end
                end
            else
                for i=1:size(ann,1)
                    if tmksel(i)
                        text(tmk(i),toptaw(2),ann{i,2},'color',map(48,:),'VerticalAlignment','baseline','HorizontalAlignment',horal);
                    end
                end
            end
        end
    end
    xlabel(['Time (' xticksi 's)']);
    ylabel(['Frequency (' yticksi frlab ')']);
end
