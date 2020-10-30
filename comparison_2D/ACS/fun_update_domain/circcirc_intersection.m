function [Bool]=circcirc_intersection(x1,y1,r1,x2,y2,r2)
%CIRCCIRC  Intersections of circles in Cartesian plane (returns 1 if two different intersection points!)
%
%  [xout,yout] = CIRCCIRC(x1,y1,r1,x2,y2,r2) finds the points
%  of intersection (if any), given two circles, each defined by center
%  and radius in x-y coordinates.  In general, two points are
%  returned.  When the circles do not intersect or are identical,
%  NaNs are returned.  When the two circles are tangent, two identical
%  points are returned.  All inputs must be scalars.
%
%  See also LINECIRC.

% Copyright 1996-2007 The MathWorks, Inc.
% Written by:  E. Brown, E. Byrns

% assert(isscalar(x1) && isscalar(y1) && isscalar(r1) && ...
%        isscalar(x2) && isscalar(y2) && isscalar(r2), ...
%     ['map:' mfilename ':mapError'], 'Inputs must be scalars')
% 
% assert(isreal([x1,y1,r1,x2,y2,r2]), ...
%     ['map:' mfilename ':mapError'], 'inputs must be real')
% 
% assert(r1 > 0 && r2 > 0, ...
%     ['map:' mfilename ':mapError'], 'radius must be positive')

% Cartesian separation of the two circle centers

r3=sqrt((x2-x1).^2+(y2-y1).^2);


if( ((r3-(r1+r2))>=-10*eps)  && r3>10*eps)  % too far apart to intersect or a single interesection point
  Bool=0;
 
elseif   ((r3<10*eps)&&(abs(r1-r2)<10*eps)) %circles identical
 Bool = 1;    
    
else %otherwise return 1 
    Bool = 1;
end
%if (r2>r3+r1);  % circle one completely inside circle two


