% XCOLOR.m
% Generate two structures with different sets of rgb colors
% Usage:
%        xcolor;
% Inputs:
%        none
% Output:
%        structure xcoln (default MATLAB color order): blue, red, 
%        yellow, purple, green, skyblue, and brown
%        structure xcolo: blue, red, black, gray, lgray (light gray),
%        green, purple, and orange.
% 
% Christian Bustamante
% September 13, 2017

% New MATBLAB's color order (since 2014)
rgbcol = get(groot,'DefaultAxesColorOrder');
xcoln.blue    = rgbcol(1,:);
xcoln.red     = rgbcol(2,:);
xcoln.yellow  = rgbcol(3,:);
xcoln.purple  = rgbcol(4,:);
xcoln.green   = rgbcol(5,:);
xcoln.skyblue = rgbcol(6,:);
xcoln.brown   = rgbcol(7,:);
clear rgbcol

% Classic palette
xcolo.blue    = [0,0,1];
xcolo.red     = [1,0,0];
xcolo.black   = [0,0,0];
xcolo.gray    = [0.5,0.5,0.5];
xcolo.lgray   = [0.75,0.75,0.75];
xcolo.green   = [50,168,26]/256;
xcolo.purple  = [102,0,204]/256;
xcolo.orange  = [255,153,0]/256;

