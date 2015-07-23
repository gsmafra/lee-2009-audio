function [b,a]=stdspectrum(s,m,f,n,zi)
%STDSPECTRUM Generate standard acoustic/speech spectra in s- or z-domain [B,A]=(S,M,F,N,ZI)
%
%Inputs:  s  Spectrum type (either text or number - see below)
%         m  mode: char 1 specifies output type,
%                    f - frequency response (complex)
%                    m - magnitude response
%                    p - power spectrum
%                    d - decibel power spectrum
%                    t - time waveform
%                    s - s-domain filter: b(s)/a(s) [default]
%                    z - z-domain filter: b(z)/a(z)
%         f  sample frequency (modes z,t) or set of frequencies in Hz (modes f,m,p,d)
%         n  number of output samples (mode t)
%         zi initial state of filter from a previous call (mode t)
%
% Outputs:  The outputs depend on the mode selected:
%         mode = 'f', 'm', 'p' or 'd'
%             b = ouptut spectrum
%         mode = 's' or 'z'
%             b,a = numerator and denonminator of the output spectrum
%         mode = 't'
%             b = output waveform
%             a = final state of the filter - use as the zi input of a future call
%
% Spectrum types (specify either as a number or case-insensitive text abbreviation):
%   1  White      : white noise
%   2  A-Weight   : the formula for this is given in [3] and is based on
%                   the equal-loudness curves of [9]
%   3  B-Weight   : this needs to be confirmed with ANSI S1.4-1981 standard or IEC 60651
%   4  C-Weight   : the formula for this is given in [3]
%   5  LTASS-P50  : the long-term average speech spectrum that is defined by a
%                   formula on page 3 of [4] which, strangely, does not precisely
%                   match the graph shown on the same page.
%   6  LTASS-1994 : the long-term average speech spectrum that is taken from a graph in [2]
%   7  SII-intinv : The inverse spectrum of the ear's internal masking noise; this is taken
%                   from table 1 of [1]. It is inverted so that it is a bandpass rather than
%                   bandstop characteristic.
%   8  BS-468     : The weighting proposed for audio frequency noise measurement in [5] and [6].
%   9  USASI      : Noise simulating long-term programme material spectrum from [7],[8].
%                   The level is such that the power is 0dB over an infinite bandwidth
%  10  POTS       : the D spectrum from [11].
%
% References:
% [1]	Methods for the calculation of the speech intelligibility index.
%       ANSI Standard S3.5-1997 (R2007), American National Standards Institute, 1997.
% [2]	D. Byrne, H. Dillon, K. Tran, S. Arlinger, K. Wilbraham, R. Cox, B. Hayerman,
%       R. Hetu, J. Kei, C. Lui, J. Kiessling, M. N. Kotby, N. H. A. Nasser,
%       W. A. H. E. Kholy, Y. Nakanishi, H. Oyer, R. Powell, D. Stephens, R. Meredith,
%       T. Sirimanna, G. Tavartkiladze, G. I. Frolenkov, S. Westerman, and C. Ludvigsen.
%       An international comparison of long-term average speech spectra.
%       JASA, 96 (4): 2108�2120, Oct. 1994.
% [3]	CENELEC. Electroacoustics - sound level meters. Technical Report EN EN 61672-1:2003, 2003.
%       (also ANSI S1.42-2001)
% [4]	ITU-T. Artificial voices. Standard P.50, Sept. 1999.
% [5]   ITU-T. Measurement of weighted noise in sound-programme circuits.
%       Recommendation J.16, 1988.
% [6]   ITU-R. Measurement of audio-requency noise voltage level in sound broadcasting.
%       Recommendation BS.468., 1986.
% [7]   NRSC AM Reemphasis, Deemphasize, and Broadcast Audio Transmission Bandwidth Specifications,
%       EIA-549 Standard, Electronics Industries Association , July 1988.
% [8]   NRSC AM Reemphasis, Deemphasize, and Broadcast Audio Transmission Bandwidth Specifications,
%       NRSC-1-A Standard, Sept 2007, Online: http://www.nrscstandards.org/SG/NRSC-1-A.pdf 
% [9]   H. Fletcher and W. A. Munson. Loudness, its definition, measurement and calculation.
%       J. Acoust Soc Amer, 5: 82�108, Oct. 1933.
% [10]  American National Standard Specification for Sound Level Meters.
%       ANSI S1.4-1983 (R2006)/ANSI S1.4a-1985 (R2006), American National Standards Institute
% [11]	IEEE standard equipment requirements and measurement techniques for analog transmission
%       parameters for telecommunications. Standard IEEE Std 743-1995, Dec. 1995.

% Other candidates: (a) Z-weighting, (b) ISO226, (c) USASI, (d) P.48 spectra
%
% Other standards:
%    IEEE743 has several weighting filters defined
%    ITU-T 0.41 Psophometer for use on telephone-type circuits
%    Bell System Technical Reference 41009 (C-message) 
%    ISO 8041:2005 (E): Human Response to Vibration � Measuring
%    Instrumentation  
%    IEC 1260:1995, class 1 (also IEC 61260/ANSI S1.11-2004) Octave band and fractional octave band filters
%    IEC 651: Specification for Sound Level Meters

%      Copyright (C) Mike Brookes 2008
%      Version: $Id: stdspectrum.m,v 1.6 2010/09/24 15:45:49 dmb Exp $
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

% constants - could make persistent
persistent spty spz az bz fz ixsz maxj
% spty contains the name of the spectrum
% spz contains a list of the poles and zeros
%    spz(1) = gain, spz(2) = number of zeros (excluding implicit conjugates), spz(3...) zeros followed by poles
%    only one of a complex conjugate pair is given
if isempty(spz)
    spty={'White';'A-Weight';'B-Weight';'C-Weight';'LTASS-P50';'LTASS-1994';'SII-IntInv';'BS-468';'USASI';'POTS'};
    spz={[1 0];
        [7390100803.6603 4 0 0 0 0 -129.42731565506 -129.42731565506 -676.40154023295 -4636.125126885 -76618.526016858 -76618.526016858];
        [5986190155.0156 3 0 0 0 -129.42731565506 -129.42731565506 -995.88487118796 -76618.526016858 -76618.526016858];
        [5912384617.784 2 0 0 -129.42731565506 -129.42731565506 -76618.526016858 -76618.526016858];
        [1.1294790345421e+015 3 0 0 -34437.856184098 -721.94747118664+754.97798119504i -1721.704402273 -5234.2950286868 -10953.570773216+42789.342252749i];
        [19720493.192959 5 0 0 0 -22550.637578954 -11319.635610404+70239.177107659i -253.31327846696+672.10855509952i -1299.1885437118+2301.2064056419i -10646.952627978+68290.702816027i -147307.51763333];
        [6.1311266394354e+018 2 0 0 -381.08293630892 -5920.974797779 -4701.76218192+24369.279310049i 10597.854874768+39258.915154617i];
        [2.1034520039796e+024 1 0 -25903.701047817 -23615.535213635+36379.908937329i -62675.170058468 -18743.746690721+62460.156452506i];
        [72.648989380657 1 0 -2*pi*[100 320]];
        [7.8820088171767e+016 4 0 0 0 0 -452.681+1924.28i -2334+1702.73i -11264.2+8213.32i -4665.8+19828.7i];
        };
    fz=-1;
    ixsz=-1;
end
if nargin<2
    if nargout
        m='s ';     % default is Laplace transform with default unit
    else
        m='d ';
    end
end
m1=m(1);        % output format

% determine the spectrum type

if ischar(s)
    ixs=find(strcmpi(s,spty));
    if isempty(ixs)
        error('undefined spectrum type: %s',s);
    end
else
    ixs=s;
end
if ixs>size(spty,1)
    error('undefined spectrum type: %d',ixs);
end

% get s-domain function
% sb/sa is transfer function
% sz and sp are lists of zeros and poles of lengths nsp and nsz

spzi=spz{ixs};
nsz=spzi(2);
sz=spzi(3:3+nsz-1);
sz=[sz conj(sz(imag(sz)~=0))];
sp=spzi(3+nsz:end);
sp=[sp conj(sp(imag(sp)~=0))];
sb=spzi(1)*poly(sz);
sa=poly(sp);

if nargin<3 && any(m1=='fmpd')
    apz=abs([sp sz]);
    apz(apz==0)=[]; % ignore zero frequency poles/zeros
    if length(apz)==0
        apz=[100 5000];
    elseif length(apz)==1
        apz=[apz/10 apz*10];
    end
    f=logspace(log10(min(apz)*0.5/pi)-0.5,log10(max(apz)*0.5/pi)+0.5);
end
if any(m1=='fmpd')
    h=polyval(sb,2i*pi*f)./polyval(sa,2i*pi*f);
end
if any(m1=='zt') && (f~=fz || ixs~=ixsz)
    % we use an iterative method to find the best digital filter
    % we initialize the phases with either bilinear or impulse invariance
    % only using the impulase invariance if it is good (max error < 10dB)
    % we then iterate invfreqz using the s-domain magnitudes and
    % the phases of the best fit so far.
    % we use log-spaced frequencies at low frequencies and linear at high
    % we then search for various numbers of poles and zeros
    fz=f;                   % save sampling frequency for future call
    ixsz=ixs;
    if ixs==1
        bz=1;
        az=1;
        maxj=0;
    else
        warning off all % avoid lots of ill-conditioning error messages
        nflin=100;      % number of frequency samples in linear region (high freq)
        alp=1.15;        % freq ratio increment in low freq region
        f0=25*2*pi/f;  % minimum interesting frequency (25 Hz in radians)
        fx=pi/nflin/(alp-1);    % boundary between log and linear portions
        if fx<=f0 || f0>=pi
            fif=linspace(0,pi,nflin);
        elseif fx>pi
            fif=[0 logspace(log10(f0),log10(pi),ceil(log10(pi/f0)/log10(alp)))];
        else
            nlin=ceil((pi-fx)*nflin/pi);
            fif=[0 logspace(log10(f0),log10(fx),ceil(log10(fx/f0)/log10(alp))) linspace(fx+(pi-fx)/nlin,pi,nlin-1)];
        end
        h0=abs(polyval(sb,1i*fif*f)./polyval(sa,1i*fif*f));   % target magnitude spectrum
        hix=h0~=0;                                              % don't calculate dB errors if zero (e.g. at DC)
        % initialize with impulse invariance
        [bz,az]=impinvar(sb,sa,f);
        hj=freqz(bz,az,fif);
        maxj=max(abs(db(abs(hj(hix)))-db(abs(h0(hix)))));
        % or else with bilinear
        [ifb,ifa]=bilinear(sb,sa,f);
        hn=freqz(ifb,ifa,fif);
        maxi=max(abs(db(abs(hn(hix)))-db(abs(h0(hix)))));
        if maxi<maxj || maxj>10 % accept bilinear if it is better or if imp inv is bad
            maxj=maxi;
            bz=ifb;
            az=ifa;
            hj=hn;
        end
        for mm=1:length(sa)     % maximum number of poles
            for nn=1:mm         % number of zeros is always less than number of poles
                hn=hj;
                j=0;
                for i=1:30          % iterate up to 30 times (usually less)
                    h=h0.*exp(1i*angle(hn));
                    [ifb,ifa]=invfreqz(h,fif,nn,mm,[],10);
                    hn=freqz(ifb,ifa,fif);
                    maxi=max(abs(db(abs(hn(hix)))-db(abs(h0(hix)))));
                    if maxi<maxj
                        maxj=maxi;
                        bz=ifb;
                        az=ifa;
                        hj=hn;
                        j=i;
                    end
                    if i>j+5    % quit if no improvement in last five iterations
                        break
                    end
                end
            end
        end
        warning on all
    end
end
switch m1
    case 'z'
        b=bz;
        a=az;
    case 't'
        if nargin<5
            [b,a]=randfilt(bz,az,n);
        else
            [b,a]=randfilt(bz,az,n,zi);
        end
    case 'm'
        b = abs(h);
    case 'f'
        b = h;
    case 'd'
        b = db(abs(h));
    case 'p'
        b=h.*conj(h);
    case 's'
        b=sb;
        a=sa;
    otherwise
        error('Output format %s not yet implemented',m1);
end

% plot data
if ~nargout
    clf;
    if m1=='z'
        t=linspace(0,2*pi);
        rtzb=roots(b);
        rtza=roots(a);
        plot(cos(t),sin(t),':r',real(rtza),imag(rtza),'xb',real(rtzb),imag(rtzb),'ob');
        axis equal;
        axis([-1.2 1.2 -1.2 1.2]);
        title(sprintf('%s (%.0f kHz, Max err = %.1f dB)',spty{ixs},f/1000,maxj));
    else
        if any(m1=='mpd')
            semilogx(f/1000,b);
            xlabel('Frequency (kHz)');
            switch m1
                case 'd'
                    ylabel('dB');
                case 'm'
                    ylabel('Gain magnitude');
                case 'p'
                    ylabel('Power gain');
            end
        elseif m1=='f'
            subplot(212)
            semilogx(f/1000,unwrap(angle(b)));
            xlabel('Frequency (kHz)');
            ylabel('Phase (rad)');
            subplot(211)
            semilogx(f/1000,db(abs(b)));
            ylabel('Magnitude');
        elseif m1=='s'
            plot(real(sp),imag(sp),'xb',real(sz),imag(sz),'ob');
            axis equal;
            xlim=get(gca,'xlim');
            xlim(1)=min(xlim(1),-1000);
            xlim(2)=max(xlim(2),1000);
            ylim=get(gca,'ylim');
            ylim(1)=min(ylim(1),-1000);
            ylim(2)=max(ylim(2),1000);
            axis([xlim ylim]);
            hold on
            plot(xlim,[0 0],':r',[0 0],ylim,':r');
            hold off
        end
        title(spty{ixs});
    end

end