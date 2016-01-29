% Evan Baker
% DataSurf

% This script makes a 3D surface plot of some data contained in an xlsx
% when the first row and first column are the x and y axes and the cell
% contents are the z values

function [gees,hertz,times] = DataSurf(m1,m2,m3)
close all
dims1 = size(m1); % Get the dimensions of the spreadsheet
dims2 = size(m2);
dims3 = size(m3);

% Get the axes and values out of the matr
gees1 = m1(1,2:dims1(2));
hertz1 = m1(2:dims1(1),1);
times1 = m1(2:dims1(1),2:dims1(2));

gees2 = m2(1,2:dims2(2));
hertz2 = m2(2:dims2(1),1);
times2 = m2(2:dims2(1),2:dims2(2));

gees3 = m3(1,2:dims3(2));
hertz3 = m3(2:dims3(1),1);
times3 = m3(2:dims3(1),2:dims3(2));

% Draw the surface of the data
subplot(3,2,1);surf(gees1,hertz1,times1);xlabel('Gs');ylabel('Hz');zlabel('Times');title('Trial 1');
subplot(3,2,2);contourf(gees1,hertz1,times1);xlabel('Gs');ylabel('Hz');zlabel('Times');colorbar;

subplot(3,2,3);surf(gees2,hertz2,times2);xlabel('Gs');ylabel('Hz');zlabel('Times');title('Trial 2');
subplot(3,2,4);contourf(gees2,hertz2,times2);xlabel('Gs');ylabel('Hz');zlabel('Times');colorbar;

subplot(3,2,5);surf(gees3,hertz3,times3);xlabel('Gs');ylabel('Hz');zlabel('Times');title('Average');
subplot(3,2,6);contourf(gees3,hertz3,times3);xlabel('Gs');ylabel('Hz');zlabel('Times');colorbar;
end