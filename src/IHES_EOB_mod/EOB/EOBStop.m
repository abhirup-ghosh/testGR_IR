function [value,isterminal,direction] = EOBStop(t,y, rmin)

%EOBStop Event function to stop EOB integration.
%
%   [value,isterminal,direction] = EOBStop(t,y, rmin)
%   Stop at minimum radius, e.g. rmin = 2.01
%   


value     = y(2)-rmin; % y(2) = r
isterminal= 1;
direction = 0;

