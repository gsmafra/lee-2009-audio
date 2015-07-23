function finishat(frac,tol,fmt)
%FINISHAT print estimated finish time of a long computation (FRAC,TOL,FMT)
%
% Inputs: FRAC = fraction of total comutation that has been completed
%                As a special case, 0 is used to initialize the routine
%         TOL  = Tolerance in minutes. If the estimated time has changed by less
%                than this, then nothing will be printed. [default 1]
%         FMT  = Format string which shoul include %s at the point where the
%                estimated finish time should be inserted.
%
% Example:       finishat(0);
%                for i=1:many
%                    long computation;
%                    finishat(i/many);
%                end

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: finishat.m,v 1.2 2008/09/24 09:58:08 dmb Exp $
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

persistent oldt oldnw
if nargin<3
    fmt='Estimated finish at %s\n';
    if nargin<2
        tol=1;
    end
end
if frac<=0
    oldt=0;
    tic;
else
    nw=now;
    newt=nw+(1/frac-1)*toc/86400;
    if ~exist('oldt') || oldt==0 || (abs(newt-oldt)>tol/1440 && (nw-oldnw)>10/86400)
        oldt=newt;
        if floor(oldt)==floor(nw)
            df='HH:MM';
        else
            df='HH:MM dd-mmm-yyyy';
        end
        fprintf(2,fmt,datestr(oldt,df));    % print to std err to avoid delay
        oldnw=nw;                           %
    end
end