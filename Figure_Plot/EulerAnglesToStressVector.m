function [sigma_vector_1,sigma_vector_2,sigma_vector_3] = EulerAnglesToStressVector(a,b,c)
%EULERANGLESTOSTRESSVECTOR 
%  Input: three Euler angles a, b, c
%  Output: vectors of three principal stresses
%
a=deg2rad(a);
b=deg2rad(b);
c=deg2rad(c);
% Assume an arbitrary stress tensor
StressTensor = [1 0 0;
                0 0 0;
                0 0 -1;];
%
%transformation from principal stress to geographic coordinates
R1=[cos(a)*cos(b) sin(a)*cos(b) -sin(b);
    cos(a)*sin(b)*sin(c)-sin(a)*cos(c) sin(a)*sin(b)*sin(c)+cos(a)*cos(c) cos(b)*sin(c);
    cos(a)*sin(b)*cos(c)+sin(a)*sin(c) sin(a)*sin(b)*cos(c)-cos(a)*sin(c) cos(b)*cos(c)];

Sg=R1'*StressTensor*R1;
%
[Vector,Diag_Matrix] = eig(Sg);

Value = eig(Diag_Matrix);
[~, j] = sort(Value,'descend');

sigma_vector_1 = Vector(:,j(1));
sigma_vector_2 = Vector(:,j(2));
sigma_vector_3 = Vector(:,j(3));
end

