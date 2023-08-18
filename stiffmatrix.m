function [H,f, deltai]=stiffmatrix(coord, topol, func)
% function [H,f, deltai]=stiffmatrix(coord, topol, func)
%  
%  Let "coord" be a matrix with n rows and 2 colomns that represents the
%  coordinates of the nodes.
%  Let "topol" be a matrix witn m rows and 3 colomns that represents the
%  nodes belonging to the triangular elements.
%  Let "func" be a vector function f(x,y).
%
%  stiffmatrix  generates the stiffnes matrix H, the right hand side f related to the 
%  linear system Hu=f where u is the vector containing the value of the
%  unknow function on the nodes.
%  Moreover il generate the deltai: the areas related to every nodes (already divided by 3).

% Pattern creation for the stiffness matrix
      %set "row" and "col" vectors for the A matrix
row=ones(size(topol,1)*size(topol,2),1);
counter=0;
for n=1:size(topol,1)
    counter=counter+1;
    row(3*n-2:3*n)=row(3*n-2:3*n)*counter;
end

col=NaN(size(topol,1)*size(topol,2),1);
for n=1:size(topol,1)
    col(3*n-2:3*n)=topol(n,:);
end

      %build th pattern of H matrix

A=sparse(row,col,1);
H=A'*A; %notice that all the nonzero elements are 1 or 2 or 3 or...
H = spones(H); 



% compute all the local stiffness matrices and then assembly in the H matrix
deltas=NaN(size(topol,1),1);
for n=1:size(topol,1)
    coordloc=[coord(topol(n,1),:); ...
              coord(topol(n,2),:); ...
              coord(topol(n,3),:) ];
    delta=0.5*det([ones(3,1),coordloc]);    
    deltas(n)=delta;

    b=[coordloc(2,2)-coordloc(3,2);...
       coordloc(3,2)-coordloc(1,2);...
       coordloc(1,2)-coordloc(2,2)];
    c=[coordloc(3,1)-coordloc(2,1);...
       coordloc(1,1)-coordloc(3,1);...
       coordloc(2,1)-coordloc(1,1)];

    Hloc=0.25/delta*([b(1)*b(1), b(1)*b(2), b(1)*b(3);...
                      b(2)*b(1), b(2)*b(2), b(2)*b(3);...
                      b(3)*b(1), b(3)*b(2), b(3)*b(3)] +...
                     [c(1)*c(1), c(1)*c(2), c(1)*c(3);...
                      c(2)*c(1), c(2)*c(2), c(2)*c(3);
                      c(3)*c(1), c(3)*c(2), c(3)*c(3)]);
    for i=1:3
        row=topol(n,i);
        for j=1:3
            col=topol(n,j);
            H(row,col)=H(row,col)+Hloc(i,j);
        end
    end
end


% % we must subtract the 1 from which we started
H=H-spones(H);


% in order to compute the rhs we need to compute the area related to every node
deltai=zeros(size(coord,1),1);

for i = 1:size(coord,1)
        for j = 1:size(topol,1)
            if i == topol(j,1) || i == topol(j,2) || i == topol(j,3)
                deltai(i) = deltai(i) + deltas(j);
            end
        end
end

deltai=deltai/3;   %already divided by three

% compute rhs
f=NaN(size(H,2),1);
for i=1:length(f)
    f(i)=func(coord(i,:))*deltai(i);
end

end







