function [v,xy,t,m]=peak2dquad(z)
%PEAK2DQUAD find quadratically-interpolated peak in a 2D array
%
%  Inputs:  Z(m,n)   is the input array
%
% Outputs:  V        is the peak value
%           XY(2,1)  is the position of the peak (in fractional subscript values)
%           T        is -1, 0, +1 for maximum, saddle point or minimum
%           M        defines the fitted quadratic: z = [x y 1]*M*[x y 1]'

%	   Copyright (C) Mike Brookes 2008
%      Version: $Id: peak2dquad.m,v 1.1 2008/11/19 16:26:00 dmb Exp $
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

persistent sz b
% first calculate the fixed matrix, b (can be stored if sz is constant)
if isempty(sz) || ~all(sz==size(z))
    sz=size(z);
    if prod(sz)<6
        error('Need at least 6 points to find peak');
    end
    xc=repmat((1:sz(1))',sz(2),1);  % x coordinate vector
    yc=reshape(repmat((1:sz(2)),sz(1),1),prod(sz),1);
    a=ones(prod(sz),6);
    a(:,1)=xc.^2;
    a(:,2)=yc.*xc;
    a(:,3)=yc.^2;
    a(:,4)=xc;
    a(:,5)=yc;
    b=(a'*a)\a';        % converts to polynomial coeficients {x^2 xy y^2 x y 1]
end

% now find the peak

c=b*z(:);
xy=zeros(2,1);
d=4*c(1)*c(3)-c(2)^2;
xy(1)=c(2)*c(5)-2*c(3)*c(4);
xy(2)=c(2)*c(4)-2*c(1)*c(5);
xy=xy/d;
w=ones(6,1);
w(1:2)=xy(1)*xy;
w(3)=xy(2)^2;
w(4:5)=xy;
v=c'*w;
t=(4*c(1)*c(3)>c(2)^2)*(2*(c(1)+c(2)>0)-1);
m=c([1 2 4; 2 3 5; 4 5 6]).*(eye(3)+1)/2;       % z = [x y 1]*m*[x y 1]'

if ~nargout
    % now plot the data

    xm=repmat((1:sz(1))',1,sz(2));  % x coordinate matrix
    ym=repmat((1:sz(2)),sz(1),1);
    pm=(c(1)*xm+c(2)*ym+c(4)).*xm+(c(3)*ym+c(5)).*ym+c(6);
    hold off
    mesh(pm','EdgeColor','k');
    hold on
    plot3(xm(:),ym(:),z(:),'rx',xy(1),xy(2),v,'bo');
    xlabel('x');
    ylabel('y');
    hold off
end


