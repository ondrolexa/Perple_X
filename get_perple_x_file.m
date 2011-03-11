function [x,y,d,xname,yname,ok] = get_perple_x_file

clear; clf

ok = 0

while ok == 0;
    
    data_file=uigetfile('*.ctr;*.tab', 'Select a Perple_X ctr or tab file');
    enter = fopen(data_file, 'rt');
    
    fmt = fscanf(enter, '%s', 1);
    
    if strcmp(fmt,'|6.6.6')
        
        ok = 1
        
        title = fscanf(enter, '%c', 1);
        
        nvar = fscanf(enter, '%f', 1); % number of independent variables
        mvar = fscanf(enter, '%f', 1); % number of dependent variables
        
        if nvar ~= 2
            
            errordlg(['The input data is ',nvar,'-dimensional, this script is only configured for 2d']);
            
        end
        
        for i = 1:nvar
            vname(i,:) = fscanf(enter, '%8c',1);
            vmin(i)  = fscanf(enter, '%f',1);
            dv(i)    = fscanf(enter, '%f', 1);
            inc(i)   = fscanf(enter, '%f', 1);
            v(i,1:inc(i))   = vmin(i):dv(i):vmin(i)+(inc(i)-1)*dv(i);
        end
        
        dnames = textscan(enter,'%s',mvar)
        [dvar, ok] = listdlg('PromptString','Select the dependent variable','SelectionMode','single','ListString',dnames{1})
        if ok == 0, errordlg(['You did not choose a variable, I quit too!']), end;
        
        fclose(enter);
        
        a = textread(data_file, '%f','headerlines',4*(nvar+1)+1);
        b = reshape(a,mvar,inc(2),inc(1));
        d = b(dvar,1:inc(2),1:inc(1))
        
        
        x = v(1,:);
        y = v(2,:);
        xname = vname(1,:);
        yname = vname(2,:);
        
        break
        
    else
        
        choice = questdlg('Invalid file format, try again?','File error','Yes','No','Yes');
        
        switch choice
            case 'Yes'
                ok = 0
                
            case 'No'
                ok = 1
                break
                
        end
        
    end
    
    
    
end