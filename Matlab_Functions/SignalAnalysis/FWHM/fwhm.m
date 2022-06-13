function[wdth]=fwhm(x,y,pts)


y_max=max(y);
maxpos = find(y==y_max);

    if pts<2000
        try
            data1(:,1) = x(maxpos-5:maxpos-1);
            data1(:,2) = y(maxpos-5:maxpos-1);
            data2(:,1) = x(maxpos+1:size(y,1)-2);
            data2(:,2) = y(maxpos+1:size(y,1)-2);   
        catch wdth = NaN;
        end    
    else
        try
            data1(:,1) = x(maxpos-15:maxpos-1);
            data1(:,2) = y(maxpos-15:maxpos-1);
            data2(:,1) = x(maxpos+1:size(y,1)-2);
            data2(:,2) = y(maxpos+1:size(y,1)-2);
        catch wdth = NaN;
    end
    end

    try
        x_halfmax1 = interp1(data1(:,2),data1(:,1),(y_max./2));
        x_halfmax2 = interp1(data2(:,2),data2(:,1),(y_max./2));

        wdth = x_halfmax1-x_halfmax2;
    catch 
        wdth = NaN;
    end
