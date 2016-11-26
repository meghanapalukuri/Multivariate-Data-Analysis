% SHADE - Band shading routine
%
% Usage: shade(x,y,y1,y2,graylevel)

function shade(x,y,y1,y2,graylevel)

if (nargin < 5)
    graylevel= 0.65;
end

cv= graylevel*ones(1,3);

% Re-sample x,y1,y2
step= (x(length(x))-x(1))/100;
xr= x(1):step:x(length(x));
y1r= interp1(x,y1,xr,'spline');
y2r= interp1(x,y2,xr,'spline');
yr= interp1(x,y,xr,'spline');

%for k= 1:length(xr),
%    ln= y1r(k):.01:y2r(k);
    %plot(xr(k)*ones(1,length(ln)),ln,'r');
%    hold on;
%end
hold on;
plot(xr,y1r,'r:','linewidth',0.5);
plot(xr,y2r,'r:','linewidth',0.5);
plot(xr,yr,'linewidth',1.5);
hold off;
