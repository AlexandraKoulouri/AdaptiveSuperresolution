function M = InterpolateMatrix2d(H,g,p)
% 
% Finds values in arbitrary points (x,y) on a given surface.
% The surface is defined by a triangulation of the data points
% obtained from DELAUNAY.
%
% Output:
%  f ... Interpolated values for points defined in the L-by-2 matrix p
%
% Input:
%  H ... Surface triangulation. M-by-3 "face matrix".
%  g ... Nodal cordinates for H. N-by-2 matrix.
%  f0 .. Data values for each point (row) in the g matrix
%  p ... Coordinates for the interpolation. L-by-2 matrix

% K. Karhunen, 26.04.2008

% A. Koulouri 15.05.2013


T = tsearchn(g,H,p); % find the triangle where point p_i lies (if NaN then the point is not in any triangle)
nT = length(T); %Number of p_i points

p = p';
L = [-1 -1;1 0;0 1];
f = zeros(nT,1);

M = sparse(size(p,2),size(g,1));

for ii=1:nT
  iH = T(ii);
  
  if ~isnan(iH)
    id = H(iH,:);
    gg = g(id,:)';
   %Barycentric interpolation 
    %X = inv(gg*L)*(p(:,ii)-gg(:,1));
    X = (gg*L)\(p(:,ii)-gg(:,1));
    %f(ii) = f0(id)'*[1-X(1)-X(2);X];

    M(ii,id) = M(ii,id) + [1-X(1)-X(2), X(1), X(2)]; %id index of Triangle
  else
    %nearest neighbor interpolation   if the point doesn't lie in an
    %triangle
    %dist_to_nearest = 1e9;
    nearest = 0;
    dist = sqrt((p(1,ii)-g(:,1)).^2+(p(2,ii)-g(:,2)).^2);
    Ind = find(dist<=min(dist));
    if(~isempty(Ind))
    %for j=1:size(g,1)
    %    dist = sqrt((p(1,ii)-g(j,1)).^2+(p(2,ii)-g(j,2)).^2);
    %    if(dist<dist_to_nearest)
         %   dist_to_nearest = dist;
          %  nearest = j;
    %    end
    %end
        %dist_to_nearest = dist(Ind);
        nearest = Ind;
    end
    %    
    M(ii,nearest) = M(ii,nearest) + 1;
  end  
end



