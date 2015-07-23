function [xx,ii,m,v]=psycest(iq,x,r,xp)
% Estimate multiple psychometric functions
%
% Usage: [xx,m,s,ii]=psycest(-n,p,q,xp) % initialize
%        [xx,m,s,ii]=psycest(i,x,r)    % trial result
%                    psycest(i)        % plot pdf
%                    psycest           % plot all pdfs
%
% Inputs:
%         -n        minus the number of models
%          p(:,n)   parameters for each model
%                      1  thresh
%                      2  miss
%                      3  guess
%                      4  SNR min
%                      5  SNR max
%                      6  Slope min
%                      7  slope max
%          q(:)     parameters common to all models (vector or struct)
%                      1  q.nx  number of SNR values in pdf [40]
%                      2  q.ns  number of slope values in pdf [21]
%                      3  q.nh  number of probe SNR values to evaluate [30]
%                      4  q.cs  weighting of slope relative to threshold in cost function [1]
%                      5  q.dh  minimum step size in dB for probe SNRs [0.2]
%                      6  q.sl  min slope at threshold [0.02]
%                      7  q.kp  number of std deviations of the pdfs to keep [4]
%                      8  q.hg  amount to grow expected gains in ni trials [1.3]
%                      9  q.cf  cost function: 1=variance, 2=entropy [2]
%          xp{n}(:) list of available probe SNR values
%          i        model probed
%          x        probe SNR value used
%          r        response: 0=wrong, 1=right.
%
% Outputs:
%          xx       recommended probe SNR
%          ii       recommended model to probe next
%          m(2,n)   estimated srt and slope of all models
%          v(3,n)   estimated covariance matrix entries:
%                   [var(srt) cov(srt,slope) var(slope)]'

%      Copyright (C) Mike Brookes 2009-2010
%      Version: $Id: psycest.m,v 1.2 2010/06/30 15:39:19 dmb Exp $
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bugs/Suggestions:
% (1) allow lookahead by 2 tests rather than 1
% (2)use a structure for input parameters
% (3)should only resample the pdfs when necessary
% (4)add a forgetting factor for old measurements
% (8) use quadratic interpolation to find cost function minimum
% (10) optionally output the whole pdf + its axis values
% (11) optionally output all model probe values
% (13) Should perhaps estimate and return the mean and slope compensated for the guess rate
%      i.e. if you want p, you actually test at guess+ p*(1-guess-miss)/(1-miss)
% (14) use logistic model instead of gaussian
% (15) make work for iq=0: plot only
% (16) output maximum likelihood instead of mean
% (17) remember probe snrs, models and results and recalculate the entire
%      array when changing precision
% (18) could force probe SNRs to be a multiple of minimum step size
% (19) use a non-uniform prior e.e. 1/(1+x^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% persistent variables
% sizes: nx=SNR values in pdf
%        ns=slope values in pdf (1/std dev of cumulative gaussian)
%        nh=potential test SNRs
%        ni=simultaneous models being estimated
% wq(ns*nx,ni) = prob of model; each column is a vectorized ns*nx matrix
% nr(1,9)      = parameter values: [SNR-pdf slope-pdf SNR-probe ...]
% pr(7,ni)     = input model parameters
% qr(5,ni)     = derived model parameters
% xq(nr(1),ni) = SNR values in pdf
% sq(nr(2),ni) = slope values in pdf
% mq(2,ni)     = estimated srt and slope of all models
% vq(3,ni)     = estimated covariance matrix entries:[var(srt) cov(srt,slope) var(slope)]'
% xn(1,ni)     = next probe value to use
% hn(1,ni)     = expected decrease in cost function after next probe

persistent wq xq sq nr pr qr mq vq xn hn hfact xz

% algorithm parametes that could be programmable
% dxmin=0.2; nr.dh          % minimum step size in dB
dxfact=2;           % potential trial range / threshold position in std devs
%mix=0.5;  nr.cs      % fraction of slope variance to use
% slowlim=0.02; nr.sl  % 1/max std dev of cum gaussian
%pdkeep=4;  nr.kp      % number of std deviations of the pdfs to keep
pfloor=0.001;    % pdf floor relative to uniform distribution
% hfact0=1.3; nr.hg % amount to grow expected gains in ni trials
tryexp=3;  % exponent determining trial region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iq<0  % initialization
    ni=-iq;  % number of models
    pr=repmat([0.75 0.04 0.1 -20 20 0 0.5]',1,ni);  % default parameters
    if nargin>1
        if size(x,2)>1
            pr(1:size(x,1),1:size(x,2))=x;
        else
            pr(1:size(x,1),:)=repmat(x,1,ni);
        end
    end
    nr=[40 21 30 1 0.2 0.02 4 1.3 2]';          % default parameter values
    nrf={'nx','ns','nh','cs','dh','sl','kp','hg','cf'}; % parameter field names
    numnr=length(nr);
    if nargin>2
        if isstruct(r)
            fnn=fieldnames(r);
            for i=1:length(fnn)
                mk=strcmp(fnn{i},nrf);
                if any(mk)
                    nr(mk)=r.(fnn{i});
                end
            end
        else
            nr(1:min(numnr,numel(r)))=r(:);
        end
        nr(1:3)=round(nr(1:3));
    end
    pr(6,:)=max(pr(6,:),nr(6));    % low limit of slope in prob/dB
    nsxq=nr(1)*nr(2);
    xq=(0:nr(1)-1)'*(pr(5,:)-pr(4,:))/(nr(1)-1)+repmat(pr(4,:),nr(1),1);
    sq=(0:nr(2)-1)'*(pr(7,:)-pr(6,:))/(nr(2)-1)+repmat(pr(6,:),nr(2),1);
    wq=repmat(1/nsxq,nsxq,ni);  % initialize to uniform pdf
    qr=zeros(5,ni);
    qr(1,:)=1-pr(2,:)-pr(3,:); % prob range covered by cumulative gaussian
    qr(2,:)=(pr(1,:)-pr(3,:))./qr(1,:); % cumulative gaussian at threshold
    qr(3,:)=norminv(qr(2,:)); % x position of target in std measure
    qr(4,:)=qr(1,:).*normpdf(qr(3,:));     % slope*stddev at threshold
    qr(5,:)=-norminv(min(qr(2,:),1-qr(2,:)).^tryexp);        % SNR region to try in std measure (a bit ad hoc this)
    mq=[mean(xq,1); mean(sq,1)];    % initial means
    vq=[var(xq,1,1); zeros(1,ni); var(sq,1,1)];  % initial variances
    %     hn=[1-nr(4) 0 nr(4)]*vq;  % artificially high value of cost function ensures all models are probed early on
    hn=repmat(Inf,1,ni);   % very high initial cost function
    hfact=nr(8)^(1/ni);   % growth factor to ensure that no model is neglected for too long
    xn=mq(1,:);     % start at mean value
    if nargin>=4
        if iscell(xp)
            xz=xp;
        else
            xz=repmat(num2cell(xp(:)',2),ni,1);        % else replicate for each model
        end
        for i=1:ni
            [dummy,j]=min(abs(xz{i}-mq(1,i)));      % select the closest available probe to the mean
            xn(i)=xz{i}(j);
        end
    else
        xz=cell(ni,1);          % make empty cells
    end
elseif iq>0 && nargin==3   % update pdfs with a new probe result
    nxq=nr(1);
    nsq=nr(2);
    nxh=nr(3);
    nsxq=nxq*nsq;
    thresh=pr(1,iq);        % target threshold
    guess=pr(3,iq);          % guess rate (1/choices)
    pscale=qr(1,iq);
    xtstd=qr(3,iq); % x position of target in std measure
    sfact=qr(4,iq);     % scale factor slope*Stddev
    trystd=qr(5,iq);        % SNR region to try in std measure (a bit ad hoc this)

    sqi=sq(:,iq);
    sqis=sqi/sfact;
    xqi=xq(:,iq);
    wqi=wq(:,iq);
    %
    % update probabilities with the previous test result
    %
    rz=r==0;
    wqi=wqi.*(rz+(1-2*rz)*(guess+pscale*normcdf(repmat(sqis,nxq,1)*x-reshape(sqis*xqi'-xtstd,nsxq,1)))); %  P(l | r,x)
    wqi=wqi./sum(wqi,1);        % normalize
    %     wq(:,iq)=wqi;               % save updated probabilities (not needed if we always interpolate them)

    % Calculate mean and covariance and entropy

    px=sum(reshape(wqi,nsq,nxq),1); % p(x0)
    ps=sum(reshape(wqi,nsq,nxq),2); % p(s0)
    xe=px*xqi;                            % E(x0)
    se=ps'*sqi;                            % E(s0)
    xv=px*(xqi.^2)-xe^2;               % Var(x0)
    sv=ps'*(sqi.^2)-se^2;               % Var(s0)
    sxv=wqi'*(repmat(sqi,nxq,1).*reshape(repmat(xqi',nsq,1),nsxq,1))-xe*se; % Cov(s0*x0)
    mq(:,iq)=[xe; se];     % save means
    vq(:,iq)=[xv; sxv; sv];   % save covariance matrix
    xh=(px*log(px)')*(xqi(1)-xqi(2));                     % h(x0)
    sh=(ps'*log(ps))*(sqi(1)-sqi(2));                     % h(s0)

    % now estimate the next test value
    if ~numel(xz{iq})
        dxt=max(nr(5),dxfact*trystd/(se*nxh));    % minimum step size of 0.2 dB
        xt=xe-xtstd/se+((1:nxh)-(1+nxh)/2)*dxt;  % select possible SNRs (row vector)
    else
        xzi=xz{iq};
        [xt,ixt]=min(abs(xzi-xe));   % find closest to xe
        ixt=max(1,min(1+numel(xzi)-nxh,ixt-floor((1+nxh)/2)));
        xt=xzi(ixt:min(ixt+nxh-1,numel(xzi)));
    end
    nxhp=length(xt);
    prt=guess+pscale*normcdf(repmat(sqis,nxq,1)*xt-repmat(reshape(sqis*xqi'-xtstd,nsxq,1),1,nxhp)); %  P(r | l,x)
    wqt=repmat(wqi,1,nxhp);
    pl1=prt.*wqt;  % posterior prob given success = p(l | x,r=1) unnormalized
    pl0=(wqt-pl1); % posterior prob given failure = p(l | x,r=0) unnormalized
    prh=sum(pl1,1);  % p(r | x) (row vector)
    pl1=pl1./repmat(prh,nsxq,1);  % normalize
    pl0=pl0./repmat(1-prh,nsxq,1);   % posterior prob given failure = p(l | x,r=0) normalized
    px1=squeeze(sum(reshape(pl1,nsq,nxq,[]),1)); % p(x0 | x,r=1)
    px0=squeeze(sum(reshape(pl0,nsq,nxq,[]),1)); % p(x0 | x,r=0)
    ps1=squeeze(sum(reshape(pl1,nsq,nxq,[]),2)); % p(s0 | x,r=1)
    ps0=squeeze(sum(reshape(pl0,nsq,nxq,[]),2)); % p(s0 | x,r=0)
    xet1=xqi'*px1;                            % E(x0 | x,r=1)
    xvt1=(xqi.^2)'*px1-xet1.^2;               % Var(x0 | x,r=1)
    xet0=xqi'*px0;                            % E(x0 | x,r=0)
    xvt0=(xqi.^2)'*px0-xet0.^2;               % Var(x0 | x,r=0)
    xvt=xvt1.*prh+xvt0.*(1-prh);            % E(Var(x0 | x ))
    set1=sqi'*ps1;                           % E(s0 | x,r=1)
    svt1=(sqi.^2)'*ps1-set1.^2;              % Var(s0 | x,r=1)
    set0=sqi'*ps0;                           % E(s0 | x,r=0)
    svt0=(sqi.^2)'*ps0-set0.^2;              % Var(s0 | x,r=0)
    svt=svt1.*prh+svt0.*(1-prh);            % E(Var(s0 | x ))
    xht1=sum(log(px1).*px1,1);                     % -H(x0 | x,r=1)
    xht0=sum(log(px0).*px0,1);                     % -H(x0 | x,r=0)
    xht=(xht1.*prh+xht0.*(1-prh))*(xqi(1)-xqi(2));            % h(x0 | x)
    sht1=sum(log(ps1).*ps1,1);                     % -H(s0 | x,r=1)
    sht0=sum(log(ps0).*ps0,1);                     % -H(s0 | x,r=0)
    sht=(sht1.*prh+sht0.*(1-prh))*(sqi(1)-sqi(2));            % h(s0 | x)
    switch nr(9)
        case 1
            hx=(xvt + nr(4)*svt)/(1+nr(4));               % cost function for each possible test SNR
            [hxmin,ix]=min(hx);                     % find the minimum of cost function
            hn(iq)=(xv + nr(4)*sv)/(1+nr(4))-hxmin;        % expected decrease in cost function
        case 2
            hx=(xht + nr(4)*sht)/(1+nr(4));               % cost function for each possible test SNR
            [hxmin,ix]=min(hx);                     % find the minimum of cost function
            hn(iq)=(xh + nr(4)*sh)/(1+nr(4))-hxmin;        % expected decrease in cost function
        otherwise
            error('Unrecognised cost function option');
    end
    xn(iq)=xt(ix);                              % next probe value for this model

    %%%%% degugging plots %%%%%%%
    % figure(101); plot(xt,hx); xlabel('Probe SNR'); title('Weighted Cost Function');
    % figure(102); plot(xt,xht); xlabel('Probe SNR'); title('Threshold Entropy');
    % figure(103); plot(xt,sht); xlabel('Probe SNR'); title('Slope Entropy');
    % figure(104); imagesc(xt,sqi,ps0); colorbar; hold on; plot(xt,set0); hold off; xlabel('Probe SNR'); ylabel('Slope'); title('Failure: Slope PDF and mean');
    % figure(105); imagesc(xt,sqi,ps1); colorbar; hold on; plot(xt,set1); hold off; xlabel('Probe SNR'); ylabel('Slope'); title('Success: Slope PDF and mean');
    % figure(106); imagesc(xqi,sqi,reshape(wqi,nsq,nxq)); colorbar; xlabel('Threshold SNR'); ylabel('Slope'); title('Joint pdf');
    % figure(107); imagesc(xt,xqi,px0); colorbar; hold on; plot(xt,xet0); hold off; xlabel('Probe SNR'); ylabel('Threshold'); title('Failure: Threshold PDF and mean');
    % figure(108); imagesc(xt,xqi,px1); colorbar; hold on; plot(xt,xet1); hold off; xlabel('Probe SNR'); ylabel('Threshold'); title('Success: Threshold PDF and mean');

    % rescale the pdfs if necessary

    ssd=sqrt(sv);
    xsd=sqrt(xv);
    sq2=linspace(max(nr(6),se-nr(7)*ssd),se+nr(7)*ssd,nsq)';
    xq2=linspace(xe-nr(7)*xsd,xe+nr(7)*xsd,nxq)';
    % do linear interpolation in the x directiono
    wqi=reshape(wqi,nsq,nxq);           % turn into a matrix for easy interpolation
    xqf=(xq2-xq(1,iq))/(xq(2,iq)-xq(1,iq));
    xqj=ceil(xqf);
    xqf=xqj-xqf;
    xqg=1-xqf;
    xqf((xqj<=0) | (xqj>nxq))=0;
    xqg((xqj<0) | (xqj>=nxq))=0;
    wq2=wqi(:,min(max(xqj,1),nxq)).*repmat(xqf',nsq,1)+wqi(:,min(max(xqj+1,1),nxq)).*repmat(xqg',nsq,1);
    % do linear interpolation in the s direction
    sqf=(sq2-sq(1,iq))/(sq(2,iq)-sq(1,iq));
    sqj=ceil(sqf);
    sqf=sqj-sqf;
    sqg=1-sqf;
    sqf((sqj<=0) | (sqj>nsq))=0;
    sqg((sqj<0) | (sqj>=nsq))=0;
    wq2=wq2(min(max(sqj,1),nsq),:).*repmat(sqf,1,nxq)+wq2(min(max(sqj+1,1),nsq),:).*repmat(sqg,1,nxq);
    % now normalize and apply a floor
    wq2=wq2(:);     % turn back into a vector
    wq2=wq2/sum(wq2,1);  % normalize
    wq2=max(wq2,pfloor/nsxq); %impose a floor
    wq2=wq2/sum(wq2,1);  % normalize again
    sq(:,iq)=sq2;
    xq(:,iq)=xq2;
    wq(:,iq)=wq2;
end

% now select the appropriate model to probe next

[hnmin,ii]=max(hn);         % chose model with the biggest expected decrease
hn=hn*hfact;         % increase values to ensure they all get a chance
xx=xn(ii);
m=mq;
v=vq;

if ~nargout

    % draw final 2-D distribution

    sqs=sq(:,iq);       % convert to units of prob/dB
    imagesc(xq(:,i),sqs,reshape(wq(:,i),nsq,nxq));
    hold on
    plot([xq(1,i) xq(end,i); xe xe]',[se se; [sq(1,i) sq(end,i)]]','w-');
    hold off
    colorbar;
    axis 'xy';
    xlabel(sprintf('%d%% threshold SNR (dB)',round(thresh*100)));
    ylabel('Psychometric Slope at threshold (prob/dB)');
    cblabel('Probability');
    title('Joint slope-threshold pdf');
end
