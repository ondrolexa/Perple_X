
% MatLab demo script to plot Perple_X tab and ctr format files.
% JADC March 11, 2011

clear all; clf;

[x,y,a,xname,yname,zname,titl] = get_perple_x_file; %open the Perple_X file

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
    if strcmp(titl,'')
        titl = zname;
    else
        titl = [titl ', ' zname];
    end
end

light;shading interp;lighting phong;axis tight;xlabel(xname);ylabel(yname);zlabel(zname);title(titl);

%colorbar; %uncomment for colorbar





