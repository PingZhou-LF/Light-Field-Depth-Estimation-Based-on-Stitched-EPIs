function [DLines,Chain] = Digitization(k,b,n)
%DIGITIZATION 此处显示有关此函数的摘要
%   此处显示详细说明

x  = (0:0.1:n)';
y  = k*x + b;

xd = (0:n)';
yd = k*xd + b;

% floor（y）
cy =  floor(yd);

% floor（y+0.5）
% cy =  floor(yd);

DLines=[xd,yd,cy];
Lines = [x,y];
Chain=diff(cy);
end

