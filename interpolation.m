function [A_3D ,B_3D ,ztemp]=interpolation(CoordsX_minWidth,CoordsY_minWidth,L1,...
    L2,L3,CoordsX_2D,CoordsY_2D,x0_el,y0_el,z0_el,ht)

xrange = unique(sort(CoordsX_minWidth(:)));
yrange = unique(sort(CoordsY_minWidth(:)));
temp_x = xrange;
temp_y = yrange;
i = 1;
k = 1;
A(1,1:3) = 0;
B(1,1:3) = 0;


for i = 1 :length (xrange)
    figure (4)
    hold on
    % This part is for Vertical Interpolation, so all A1,A2,A3 x points are
    % same
    A1 = temp_x(i);
    A2 = A1;
    A3 = A1;
    % Ellipse Equation in 2D
    y_ell =@(x) L2* sqrt(1-(((x-x0_el).^2)/(L1^2)))+y0_el;
    
    B3 = y_ell(A1);
    % Find a mid point between y points
    [row,col] = find( CoordsX_minWidth==A1);
    temp_1 = CoordsY_minWidth(row,col);
    % Upwards Interpolation
    B1_max =  max(temp_1,[],'all');
    B2 = (B3+B1_max)/2;
    A(end+1,:) = [A1,A2,A3];
    B(end+1,:) = [B1_max,B2,B3];
    %     plot(A(end,:),B(end,:))
    %     A_3D_2 = [A1,A2,A3];
    %     B_3D_2 = [B1_max,B2,B3];
    A_3D_2 = [A2,A3];
    B_3D_2 = [B2,B3];
    hold on
    
    % Downward Interpolation
    y_ell =@(x) -L2* sqrt(1-(((x-x0_el).^2)/(L1^2)))+y0_el;
    B3 = y_ell(A1);
    [row,col] = find( CoordsX_minWidth==A1);
    temp_1 = CoordsY_minWidth(row,col);
    B3 = y_ell(A1);
    B1_min =  min(temp_1,[],'all');
    B2 = (B3+B1_min)./2;
    A(end+1,:) = [A1,A2,A3];
    B(end+1,:) = [B1_min,B2,B3];
    %      plot(A(end,:),B(end,:),'LineWidth',1)
    %     A_3D_1 = [A1,A2,A3];
    %     B_3D_1 = [B1_min,B2,B3];
    A_3D_1 = [A2,A3];
    B_3D_1 = [B2,B3];
    
%     count = 0;
%     while true
%         hold off
%         figure(1)
%         hold on
%         % Vertical Interpolation 3D
%         [row,col] = find (CoordsX_2D== A1);
%         % Find Min and Max X in 1m width squares to limit interpolation
%         %
%         y_mid_vals = CoordsY_2D(row,col);
%         
%         if any(diff((y_mid_vals))==1)
%             %  Find the number of 1m width squares to fill fap between A_3D
%             mid_fac =unique(sort(y_mid_vals(:)))';
%             A_3D = [ A_3D_1 A1.*ones(1,length(mid_fac)) A_3D_2];
%             B_3D = [ fliplr(B_3D_1) mid_fac B_3D_2];
%             ztemp  = real( -L3*sqrt(1-(((A_3D-x0_el).^2)/ (L1^2)) - (((B_3D-y0_el).^2)/ (L2^2))))+z0_el;
%              plot3(A_3D,B_3D, ztemp,'LineWidth',2)
%             hold on
%         elseif all(diff((y_mid_vals)) < 1)
%             mid_fac =unique(sort(y_mid_vals(:)))';
%             A_3D = [ A_3D_1 A1.*ones(1,length(mid_fac)) A_3D_2];
%             B_3D = [ fliplr(B_3D_1) mid_fac B_3D_2];
%             ztemp  = real( -L3*sqrt(1-(((A_3D-x0_el).^2)/ (L1^2)) - (((B_3D-y0_el).^2)/ (L2^2))))+z0_el;
%             plot3(A_3D,B_3D,ztemp,'LineWidth',2)
%             hold on
%         elseif any(diff((y_mid_vals)) > 1)
%             mid_fac =unique(sort(y_mid_vals(:)))' ;
%             temp = length(mid_fac)/2;
%             mid_fac = [mid_fac(1:temp) NaN mid_fac(temp+1:end)];
%             A_3D = [ A_3D_1  A1.*ones(1,length(mid_fac)) A_3D_2];
%             B_3D = [ fliplr(B_3D_1) mid_fac B_3D_2];
%             ztemp  = real(- L3*sqrt(1-(((A_3D-x0_el).^2)/ (L1^2)) - (((B_3D-y0_el).^2)/ (L2^2))))+z0_el;
%             
%             plot3(A_3D,B_3D,ztemp,'LineWidth',2)
%             hold on
%         end
%         count=count+1;
%         
%         if count==1
%             hold off
%             break
%         end
%     end
end



for k = 1 :length (yrange)
    figure(4)
    hold on
    % This part is for Horizontal Interpolation, so all B1,B2,B3 y points are
    % same
    B1 = temp_y(k);
    B2 = B1;
    B3 = B1;
    % Ellipse Equation in 2D
    x_ell =@(y) L1* sqrt(1-(((y-y0_el).^2)/(L2^2)))+x0_el;
    A3 = x_ell(B1);
    % Find a mid point between x points
    [row,col] = find( CoordsY_minWidth==B1);
    temp_2 = CoordsX_minWidth(row,col);
    % Left Interpolation
    A1_max =  max(temp_2,[],'all');
    A2 = (A3+A1_max)./2;
    A(end+1,:) = [A1_max,A2,A3];
    B(end+1,:) = [B1,B2,B3];
    %     plot(A(end,:),B(end,:))
    A_3D_2 = [A1_max,A2,A3];
    B_3D_2 = [B1,B2,B3];
    
    
    hold on
    
    % Left Interpolation
    % Ellipse Equation in 2D
    x_ell =@(y) -L1* sqrt(1-(((y-y0_el).^2)/(L2^2)))+x0_el;
    A3 = x_ell(B1);
    % Find mid point between x and y
    [row,col] = find( CoordsY_minWidth==B1);
    temp_2 = CoordsX_minWidth(row,col);
    A1_min =  min(temp_2,[],'all');
    A2 = (A3+A1_min)./2;
    A(end+1,:) = [A1_min,A2,A3];
    B(end+1,:) = [B1,B2,B3];
    %    plot(A(end,:),B(end,:))
    A_3D_1 = [A1_min,A2,A3];
    B_3D_1 = [B1,B2,B3];

    count = 0;
    while true
        hold off
        figure(1)
        hold on
        % On x direction Interpolation 3D
        [row,col] = find (CoordsY_2D== B1);
        x_mid_vals = unique(CoordsX_2D(row,col));
        
        if any(diff((x_mid_vals))==1)
            mid_fac =unique(sort(x_mid_vals(:)))';
            B_3D = [ B_3D_1 B1.*ones(1,length(mid_fac)) B_3D_2];
            A_3D = [ fliplr(A_3D_1) mid_fac A_3D_2];
            ztemp  = real( -L3*sqrt(1-(((A_3D-x0_el).^2)/ (L1^2)) - (((B_3D-y0_el).^2)/ (L2^2))))+z0_el;
            plot3(A_3D,B_3D ,ztemp)
  
        elseif all(diff((x_mid_vals)) < 1)
            mid_fac =unique(sort(x_mid_vals(:)))';
            B_3D = [ B_3D_1 B1.*ones(1,length(mid_fac)) B_3D_2];
            A_3D = [ fliplr(A_3D_1) mid_fac A_3D_2];
            ztemp  = real( -L3*sqrt(1-(((A_3D-x0_el).^2)/ (L1^2)) - (((B_3D-y0_el).^2)/ (L2^2))))+z0_el;
            plot3(A_3D,B_3D ,ztemp)

          
        elseif any(diff((x_mid_vals)) > 1)
            %  Find the number of 1m width squares to fill fap between A_3D
            mid_fac =unique(sort(x_mid_vals(:)))' ;
            temp = length(mid_fac)/2;
            mid_fac = [mid_fac(1:temp) NaN mid_fac(temp+1:end)];
            B_3D = [ B_3D_1 B1.*ones(1,length(mid_fac)) B_3D_2];
            A_3D = [ fliplr(A_3D_1) mid_fac A_3D_2];
            ztemp  = real(- L3*sqrt(1-(((A_3D-x0_el).^2)/ (L1^2)) - (((B_3D-y0_el).^2)/ (L2^2))))+z0_el;
            plot3(A_3D,B_3D ,ztemp)
        end
   
 count=count+1;
        if count==1
            hold off
            break
        end
    end
    
end
A(1,:) = [];
B(1,:) = [];










end