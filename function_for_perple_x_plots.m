function [] = function_for_perple_x_plots (x,y,a,xname,yname,zname,nvar,mvar,nrow,dnames,titl)
%Generic function to make 2- and 3-d plots from Perple_X tab format files
%

figure(1);

if nvar == 1 % two cases: 1d - table -> 2d plot
    
    [kvar, ok] = listdlg('PromptString','Select the INDEPENDENT (X) variable:','ListSize',[240 500],'SelectionMode','single','ListString',dnames{1});
    if ok == 0, errordlg(['You did not choose a variable, I quit!']), end;
    
    if mvar > 2,
        [dvar, ok] = listdlg('PromptString','Select the DEPENDENT (Y) variables:','ListSize',[240 500],'ListString',dnames{1});
        if ok == 0, errordlg(['You did not choose a variable, I quit!']), end;
        [n m] = size(dvar);
        jvar = n*m;
    else   % mvar must be 2
        dvar = 2;
        jvar = 1;
    end
    
    hold all
    
    for i = 1:jvar,plot(a(kvar,1:nrow),a(dvar(i),1:nrow)),end
    
    legend(dnames{1}{dvar},'Location','Best');axis tight;xlabel(xname);title(titl);
    
elseif nvar == 2 % 2d - table -> 2/3d plot
    
    amin = min(a(:)); amax = max(a(:));
    
    choice = questdlg([zname,'range is ',num2str(amin),'->',num2str(amax) '. Select plot style:'],'Plot Style','3D Surface','Auto-Contour','Contour','3D Surface');
    
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
            c = inputdlg(prompt,dlg,num_lines,def);
            helpdlg('Carefully select contours for labeling. When done, press RETURN while the Graph window is the active window.');
            contours = [str2num(c{1}):str2num(c{3}):str2num(c{2})];
            [C,h]= contour(x,y,a,contours);clabel(C,h,'manual'); d2 = 1;
 
    end
    
    if d2 == 1,
        if strcmp(titl,' ')
            titl = zname;
        else
            titl = [deblank(titl) ', ' zname];
        end
    end
    
    light;shading interp;lighting gouraud;axis tight;xlabel(xname);ylabel(yname);zlabel(zname);title(titl);
    
    %colorbar; %uncomment for colorbar
    
end

