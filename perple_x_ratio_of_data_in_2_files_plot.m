
% MatLab script to plot the ratio of properties in two Perple_X 2d tab format files.
% JADC March 30, 2011

clear all; clf;
%                                 NUMERATOR
helpdlg('The next prompts select the file and variable to be used for the numerator of the ratio','Numerator Specification');
[x,y,a,xname,yname,z1name,nvar,mvar,nrow,dnames,titl] = get_perple_x_file; %open the Perple_X file
if nvar == 1 ,errordlg(['This script is only for 3d data, I quit!']), end;
%                                 DENOMINATOR
helpdlg('The next prompts select the file and variable to be used for the denominator of the ratio','Denominator Specification');
[x,y,b,xname,yname,z2name,nvar,mvar,nrow,dnames,titl] = get_perple_x_file; %open the Perple_X file
if nvar == 1 ,errordlg(['This script is only for 3d data, I quit!']), end;
%                                 make the ratio, could check sizes first. 
a = a./b;
zname = [z1name,'/',z2name];

figure(1);

amin = min(a(:)); amax = max(a(:)); disp(['Grid data range is ',num2str(amin),' - >',num2str(amax)])

choice = questdlg('Select plot style','Plot Style','3D Surface','Auto-Contour','Contour','3D Surface');

switch choice;
    
    case '3D Surface';
        surf(x,y,a); d2 = 0;
    case 'Auto-Contour';
        [C,h]=contour(x,y,a); clabel(C,h); d2 = 1;
    case 'Contour';
        prompt = {'Minimum contour:','Maximum contour:','Contour interval:'};
        dlg = 'Contour specification';
        num_lines = 1;
        da = (amax-amin)/11;
        def = {num2str(amin+da/2),num2str(amax-da/2),num2str(da)};
        ans = inputdlg(prompt,dlg,num_lines,def);
        contours = [ans{1}:ans{3}:ans{2}];
        [C,h]= contour(x,y,a,contours);clabel(C,h,'manual'); d2 = 1;
end

if d2 == 1,
    if strcmp(titl,' ')
        titl = zname;
    else
        titl = [titl ', ' zname];
    end
end

light;shading interp;lighting phong;axis tight;xlabel(xname);ylabel(yname);zlabel(zname);title(titl);

%colorbar; %uncomment for colorbar




