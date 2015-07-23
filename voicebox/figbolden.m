function figbolden(pos,pv)
% embolden the current figure
% Inputs: pos = [xmin ymin width height] gives the lower left corner position and the window size in pixels
%               [width height] leaves the lower left corner alone
%               [width] has a standard aspect ratio of 4:3
%               [-width/height] leaves the area unchanged but fixes the aspect ratio
%         pv is a cell name containing attribute-value pairs.
%            default = {'FontName' 'Ariel'; 'FontSize' 16; 'LineWidth' 2; 'MarkerSize' 8}

%      Copyright (C) Mike Brookes 2003
%      Version: $Id: figbolden.m,v 1.9 2009/09/29 21:07:35 dmb Exp $
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

ps={'Title' 'XLabel' 'YLabel' 'Children'};
if nargin<2
    pv={'FontName' 'Arial'; 'FontSize' 16; 'LineWidth' 2; 'MarkerSize' 8};
end
pp={'Symbol';'Wingdings'};      % protected fonts
if nargin<1
    pos=[];
end
scsz=get(0,'screensize');
if length(pos)
    po=get(gcf,'position');
    if length(pos)>2            % position is specified
        po(1:2)=pos(1:2);
        pos(1:2)=[];      % remove xmin,ymin
    end
    if length(pos)>1
        po(3:4)=pos(1:2);
    else
        if pos(1)>0
            po(3:4)=[1 0.75]*pos(1);
        else
            po(3:4)=[-pos(1) 1]*sqrt(-po(3)*po(4)/pos(1)); % preserve area
        end
    end
    set(gcf,'position',po);
end
hlist=get(gcf,'children');
while length(hlist)
    pl=get(hlist(1));
    %fprintf('list length = %d, handle = %f\n',length(hlist),hlist(1));
    for i=1:size(pv,1)
        if isfield(pl,pv{i,1})
            if i>1 || all(~strcmp(get(hlist(1),pv{i,1}),pp))
            set(hlist(1),pv{i,1},pv{i,2})
            end
            %fprintf('set %f %s\n',hlist(1),pv{i,1});
        end
    end
    for i=1:length(ps)
        if isfield(pl,ps{i})
            hlist=[hlist; get(hlist(1),ps{i})];
            %fprintf('add %f:%s\n',hlist(1),ps{i});
        end
    end
    hlist(1)=[];
end

