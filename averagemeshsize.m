function ams= averagemeshsize(coord,topol)
%
%  function ams= averagemeshsize(coord,topol)
%
%  Let "coord" be a matrix with n rows and 2 colomns that represents the
%  coordinates of the nodes.
%  Let "topol" be a matrix witn m rows and 3 colomns that represents the
%  nodes belonging to the triangular elements.
%  
% It generates the  average mesh size (ams) i.e. the average length of the
% triangular edges.


sumofperim=0;
for n=1:size(topol,1)
    coordloc=[coord(topol(n,1),:); ...
              coord(topol(n,2),:); ...
              coord(topol(n,3),:) ];
    perim=norm(coordloc(1,:)-coordloc(2,:),2)+...    %compute the length ot the edges using pitagora formula
          norm(coordloc(2,:)-coordloc(3,:),2)+...
          norm(coordloc(3,:)-coordloc(1,:),2);
    sumofperim=sumofperim+perim;
end

ams=sumofperim/(3*size(topol,1));    % we divide by 3 times the number of triangle because we want
                                     % the average length of the edges.

end