function [ D ] = MinkowskiDistance( X,Y, p )
%SUMMARY: Computes the metric distance between
%         the columns of matrices X,Y in the metric space with exponent p.
%Details: The distances are stored and returned in the matrix D.
%         D has dimension columns of X by columns of Y.
%         If the columns of one matrix have larger dimension 
%         than the ones of the other, the missing entries are paddes as
%         zeros. 
%Example: Type for example: 
%          X = magic(4); Y = magic(3);
%          D = MinkowskiDistance(X,Y,3);

[mx,nx] = size(X); 
[my,ny] = size(Y);

m = max(mx,my);
X = [X; zeros(m - mx, nx)];
Y = [Y; zeros(m - my, ny)];

D = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        D(i,j) = sum( abs(X(:,i) - Y(:,j)).^p )^(1/p);
    end
end

end

