
% MatLab demo script to plot Perple_X tab and ctr format files.
% JADC March 12, 2011

clear all; clf;

[x,y,a,xname,yname,zname,nvar,mvar,nrow,dnames,titl] = function_to_get_perple_x_file; %open the Perple_X file

function_for_perple_x_plots (x,y,a,xname,yname,zname,nvar,mvar,nrow,dnames,titl);

