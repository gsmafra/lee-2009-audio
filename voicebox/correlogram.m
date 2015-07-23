function [y,ty]=correlogram(x,inc,nwin,nlag,m,fs)
% make correlogram,
% Inputs:
%          x   is the output of a filterbank
%              with one column per filter channel
%        inc   frame increment in samples
%       nwin   window length in samples [or window function]
%       nlag   number of lags to calculate
%          m   mode string
%              [d = subtract DC component]
%              [n = unnormalized]
%              [v = variable analysis window proportional to lag]
%              [p = output the peaks only]
%         fs   sample freq (affects only plots)
%
% Outputs:
%          y(lag,chan,frame) is correlogram
%          tf                is time of the window energy centre for each frame 
%                            x(n) is at time n

% fake input values

% x=[0 1 0 0 1 1 0 1 0 0 1; 1 2 1 2 1 2 3 2 1 2 1].';
% inc=1;
% nwin=4;
% nlag=3;
% m='';
% fs=1;

memsize=voicebox('memsize');    % set memory size to use
[nx,nc]=size(x);  % number of sampes and channels
nwl=nwin+nlag-1;
nt=pow2(nextpow2(nwl));  % transform length
nf=floor((nx-nwl+inc)/inc);  % number of frames
i1=repmat((1:nwl)',1,nc)+repmat(0:nx:nx*nc-1,nwl,1);
nb=min(nf,max(1,floor(memsize/(16*nwl*nc))));    % chunk size for calculating
nl=ceil(nf/nb);                  % number of chunks
jx0=nf-(nl-1)*nb;                % size of first chunk in frames
win=ones(nwin,1);               % window function
wincg=(1:nwin)*win.^2/(win'*win);
fwin=conj(fft(win,nt,1)); % fft of window
y=zeros(nlag,nc,nf);
% first do partial chunk
jx=jx0;
x2=zeros(nwl,nc*jx);
x2(:)=x(repmat(i1(:),1,jx)+repmat((0:jx-1)*inc,nwl*nc,1));
v=ifft(conj(fft(x2(1:nwin,:),nt,1)).*fft(x2,nt,1));
w=real(ifft(fwin(:,ones(1,nc*jx)).*fft(x2.*conj(x2),nt,1)));
w=sqrt(w(1:nlag,:).*w(ones(nlag,1),:));
if isreal(x)
    y(:,:,1:jx)=reshape(real(v(1:nlag,:))./w,nlag,nc,jx);
else
    y(:,:,1:jx)=reshape(conj(v(1:nlag,:))./w,nlag,nc,jx);
end
% now do remaining chunks
x2=zeros(nwl,nc*nb);
for il=2:nl
    ix=jx+1;            % first frame in this chunk
    jx=jx+nb;           % last frame in this chunk
    x2(:)=x(repmat(i1(:),1,nb)+repmat((ix-1:jx-1)*inc,nwl*nc,1));
    v=ifft(conj(fft(x2(1:nwin,:),nt,1)).*fft(x2,nt,1));
    w=real(ifft(fwin(:,ones(1,nc*nb)).*fft(x2.*conj(x2),nt,1)));
    w=sqrt(w(1:nlag,:).*w(ones(nlag,1),:));
    if isreal(x)
        y(:,:,ix:jx)=reshape(real(v(1:nlag,:))./w,nlag,nc,nb);
    else
        y(:,:,ix:jx)=reshape(conj(v(1:nlag,:))./w,nlag,nc,nb);
    end
end
ty=(0:nf-1)'*inc+wincg;       % calculate times of window centres


