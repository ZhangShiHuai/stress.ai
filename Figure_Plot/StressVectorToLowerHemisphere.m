function [x_sigma,y_sigma] = StressVectorToLowerHemisphere(StressVector)
%STRESSVECTORTOLOWERHEMISPHERE 
%   zenithal equal-area projection
fi = rad2deg(atan(abs(StressVector(2)/StressVector(1))));
%
if (StressVector(2)>=0 && StressVector(1)>=0)
    Azimuth_Sigma = fi; 
end
%
if (StressVector(2)>=0 && StressVector(1)< 0)
    Azimuth_Sigma = 180-fi;
end
%
if (StressVector(2)< 0 && StressVector(1)< 0)
    Azimuth_Sigma = 180+fi;
end
%
if (StressVector(2)< 0 && StressVector(1)>=0)
    Azimuth_Sigma = 360-fi;
end
%
Theta_sigma = rad2deg(acos(abs(StressVector(3))));
%

% projection onto the lower hemisphere
% length on the X-Y plane
Radius_sigma = 1*sin(Theta_sigma*pi/360);

x_sigma = sqrt(2)*Radius_sigma*cos(deg2rad(Azimuth_Sigma));
y_sigma = sqrt(2)*Radius_sigma*sin(deg2rad(Azimuth_Sigma));

end

