function [varargout]=make_stereonet
% This function makes a blank lower hemisphere stereonet on a new figure and 
% optionally outputs the figure number: fignum. 
% run this before any of the other StereonetToolbox scripts, which will
% subsequently plot things on the stereonet. 
% Example:
% fignum = make_stereonet
%
% from Pollard and Fletcher 2010, Fundamentals of Structural
% Geology. 
% by Rall Walsh
% frwalsh at stanford dot edu
%  
% 
% figure number
fignum=figure; hold on % clear memory and figure, hold plot
set(gcf,'color','w')
axis equal, axis off, box off % equal scaling in x and y, no axes or box
axis ([-1 1 -1 1]) % sets scaling for x- and y-axes

plot([-2 2],[0 0],'--','Color',[0.8 0.8 0.8]) % plot x-axis
hold on
plot([0 0],[-2 2],'-','Color',[0.8 0.8 0.8]) % plot y-axis
r = 1; % radius of reference circle
TH = linspace(0,2*pi,3601); % polar angle, range 2 pi, 1/10 deg increment
[X,Y] = pol2cart(TH,r); % Cartesian coordinates of reference circle
plot(X,Y,'k','LineWidth',2) % plot reference circle

for j = 1:8 % loop to plot great circles at 10 degree increments
phid = j*(10*pi/180); % dip angle, phid
h = -r*tan(phid); rp = r/cos(phid); % x-coord of center, h, and radius, rp
X = -h + rp*cos(TH); Y = rp*sin(TH); % coordinates of pts on great circle
for q=1:length(X)   %get rid of data points outside great circle
    if ((((X(1,q)^2)+(Y(1,q)^2))^(1/2))>r);
     X(1,q)=NaN;
     Y(1,q)=NaN;
    end
end
plot(X,Y,'Color',[0.8 0.8 0.8]) % plot two sets of great circles
hold on
plot(-X,Y,'Color',[0.8 0.8 0.8])
end

for j = 1:8 % loop to plot small circles at 10 degree increments
gam = j*(10*pi/180); % cone angle, gam
k = r/cos(gam); rp = r*tan(gam); % y-coord of center, k, and radius, rp
X = rp*cos(TH); Y = k + rp*sin(TH); % coordinates of points on small circle
for q=1:length(Y)        %sort through to find points outside great circle
    if ((((X(1,q)^2)+(Y(1,q)^2))^(1/2))>r);
     X(1,q)=NaN;
     Y(1,q)=NaN;
    end
end
plot(X,Y,'--','Color',[0.8 0.8 0.8]) % plot two sets of small circles
hold on
plot(X,-Y,'--','Color',[0.8 0.8 0.8])
end
text(0,1.1,'N','fontsize',16,'horizontalalignment','center') % north 
text(1.1,0,'E','fontsize',16,'horizontalalignment','center') % north 
text(0,-1.1,'S','fontsize',16,'horizontalalignment','center') % north 
text(-1.1,0,'W','fontsize',16,'horizontalalignment','center') % north 

if nargout==1 % if outputting figure handle
    varargout{1}=fignum;
    
end
    

end
