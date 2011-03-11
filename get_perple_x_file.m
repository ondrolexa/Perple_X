function [x,y,a,xname,yname,zname,titl] = get_perple_x_file

% MatLab script to read Perple_X tab and ctr format files. 
% JADC March 11, 2011

ok = 0;

while ok == 0;
    
    data_file=uigetfile('*.ctr;*.tab', 'Select a Perple_X ctr or tab file');
    
    fid = fopen(data_file, 'rt');
    
    fmt = fscanf(fid, '%s', 1); %read revision tag
    
    if strcmp(fmt,'|6.6.6')     %valid revision
        
        ok = 1;
        
        titl  = fscanf(fid, '%c', 1); % title
        nvar  = fscanf(fid, '%f', 1); % number of independent variables
        mvar  = fscanf(fid, '%f', 1); % number of dependent variables
        
        if nvar == 2 % 2d table
            
            for i = 1:nvar                % independent variables
                vname(i,:) = fscanf(fid, '%8c', 1);
                vmin(i)    = fscanf(fid, '%f', 1);
                dv(i)      = fscanf(fid, '%f', 1);
                inc(i)     = fscanf(fid, '%f', 1);
                v(i,1:inc(i))   = vmin(i):dv(i):vmin(i)+(inc(i)-1)*dv(i);
            end
            
            if mvar == 1
                
                zname = fscanf(fid, '%c',1);
                dvar = mvar;
                
            else
                
                dnames = textscan(fid,'%s',mvar); % dependent variable names
                
                [dvar, ok] = listdlg('PromptString','Select the dependent variable:','ListSize',[200 400],'SelectionMode','single','ListString',dnames{1});
                if ok == 0, errordlg(['You did not choose a variable, I quit!']), end;
                
                zname = dnames{1}{dvar};
                
            end
            
            fclose(fid);
            
            a = textread(data_file, '%f','headerlines',4*(nvar+1)+1);
            a = reshape(a,mvar,inc(2),inc(1));
            a = reshape(a(dvar,1:inc(2),1:inc(1)),inc(2),inc(1));
            
            x = v(1,1:inc(1));
            y = v(2,1:inc(2));
            xname = vname(1,:);
            yname = vname(2,:);
            
        elseif nvar ==1 % 1d table
            
        else
            errordlg(['The input data is ',nvar,'-dimensional, this script is only configured for 2d']);
        end
        
    else
        
        choice = questdlg('Invalid file format, try again?','File error','Yes','No','Yes');
        
        switch choice;
            case 'Yes';
                ok = 0;
                
            case 'No';
                ok = 1;
                break;
        end
        
    end % end while
end

end