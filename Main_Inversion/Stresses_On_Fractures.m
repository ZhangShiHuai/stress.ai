function [tau, Sn, rake] = Stresses_On_Fractures(strike,dip_angle,a,b,c,phi,S3) 
%MohrFrac_arb(S1, S2, S3,Pp,coeff friction, a, b, c, ...
%               faults=2 column matrix with strike(0-360), dip(0-90) (use
%               right-hand rule)), e.g faults=[208,60;36,62;...]
%All angle information should be input as DEGREES, stresses as MPa (output
%will be the same)
% a=trend of S1...exception when S1 is vertical a=trend SHmax-90 degrees
% b=-plunge of S1 (plunge is angle from horizontal)
% c=rake of S2, 0 if S1 or S3 is vertical, 90 if S2 is vertical

%converting fault and stress orientation info from degrees to radians 
str = strike;
dip = dip_angle;
str=deg2rad(str);
dip=deg2rad(dip);
a=deg2rad(a);
b=deg2rad(b);
c=deg2rad(c);

%defines principal stress tensor
S1 = 1; S2 = (1-S3)*phi+S3;
Pp=0; %
S=[S1 0  0
    0 S2 0
    0 0 S3];
%transformation from principal stress to geographic coordinates
R1=[cos(a)*cos(b) sin(a)*cos(b) -sin(b);
    cos(a)*sin(b)*sin(c)-sin(a)*cos(c) sin(a)*sin(b)*sin(c)+cos(a)*cos(c) cos(b)*sin(c);
    cos(a)*sin(b)*cos(c)+sin(a)*sin(c) sin(a)*sin(b)*cos(c)-cos(a)*sin(c) cos(b)*cos(c)];

Sg=R1'*S*R1;

%the following part of the code is dependent on fault orientation
%for loop makes calculations for each fault one at a time
frac_no = length(strike);
Sn = zeros(frac_no,1);
tau = zeros(frac_no,1);
rake = zeros(frac_no,1);
for k=1:frac_no
    % transformation from geographic to fault coordinate system
    R2=[cos(str(k)) sin(str(k)) 0;
         sin(str(k))*cos(dip(k)) -cos(str(k))*cos(dip(k)) -sin(dip(k));
         -sin(str(k))*sin(dip(k)) cos(str(k))*sin(dip(k)) -cos(dip(k))];

    Sf = R2*Sg*R2';
    %Solves for normal stress resolved on the fault surface
    Sn(k) = Sf(3,3);

    %Solve for rake of the slip vector
    if (Sf(3,2)>0)
        rake(k) = (atan(Sf(3,2)/(Sf(3,1))));
    elseif (Sf(3,2)<0)&&(Sf(3,1)>0)
        rake(k) = pi-atan(Sf(3,2)/(-Sf(3,1)));
    else
        rake(k) = atan((-Sf(3,2))/(-Sf(3,1)))-pi;
    end
    %Solve for shear stress resolved on the fault plane
    R3 = [cos(rake(k)) sin(rake(k)) 0;
        -sin(rake(k)) cos(rake(k)) 0;
        0 0 1];

    Sr = R3*Sf*R3';

    tau(k,1) = Sr(3,1);
    rake(k) = rad2deg(rake(k));
end

%output data
%changing signs to make output easier to understand
% rake = (sign(tau).*rake');
Sn = Sn - Pp;
tau = abs(tau);
end
