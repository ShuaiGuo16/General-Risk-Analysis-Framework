% Reinitialize geometrie values for adjoint prediction (build new discretization
% matrix). 

l.d1 = sal.d1*change(3); % 0.07 %the top rectangle
l.l1 = sal.l1*change(6); % 0.1 % 0.2 % 0.4
l.d2 = sal.d2*change(2); % 0.022 %the middle rectangle
l.l2 = sal.l2*change(5); % 0.074
l.d3 = sal.d3*change(1); % 0.03 %the bottom rectangle
l.d3u = sal.d3u*change(1); % 0.065
l.l3o = sal.l3o; % 0.06
l.l3u = sal.l3u*change(4);  % 0.1 0.164 0.228
l.l3  = l.l3o + l.l3u;

% thining
thin=1;
l.d2 = thin*l.d2;
l.d1 = thin*l.d1;
l.d3 = thin*l.d3;
l.d3u = thin*l.d3u;

% Geometrie of bottom cavitiy
% Vector for diameter

l.d3vec = [l.d3,l.d3u];
l.l3vec = [0, -l.l3o];

l.pp=spline(l.l3vec, [0, l.d3vec/2, slopeVec(ns)]);
yy = linspace(0, -l.l3o);
xx = ppval(l.pp, yy);

xx2 = [xx(1:end-1),0.5*l.d3u*ones(1, round(l.l3u/0.005))];
yy2 = [yy(1:end-1),linspace(-l.l3o,-l.l3,length(xx2) - length(xx)+1) ];

yy3 = linspace(0, -l.l3);
xx3 = spline(yy2,xx2,yy3)

l.pp = spline(yy3, xx3);




% npoints = 80;
% yy = linspace(0, -l.l3, npoints);
% xx = ppval(l.pp, yy);
% 
% intshift = round(shift*npoints); 
% xx(intshift+1) = xx(intshift+1)+l.d3u*0.03;
% xx(intshift+2) = xx(intshift+2)+l.d3u*0.07;
% xx(intshift+3:intshift+4) = xx(intshift+3:intshift+4)+l.d3u*0.12;
% xx(intshift+5) = xx(intshift+5)+l.d3u*0.07;
% xx(intshift+6) = xx(intshift+6)+l.d3u*0.03;
% 
% l.pp = spline(yy, xx);