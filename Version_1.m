clc
clear all
close all



% Coordinates for the x-axis
x1 = -15;
x2 = -10;
x3 =  10;
x4 =  15;


% Coordinates for the y-axis
y1 = -10;
y2 =  10;


% Coordinates for the z-axis
z1 = 0;
z2 = 0;
z3 = 6;
z4 = 6;

x = [ x1 x2 x3 x4];
y = [y1 y1 y2 y2];
z = [ z1 z2 z3 z4];



 [x,y] = meshgrid(x,y);
[z,~] = meshgrid(z);
figure(1)
 surf(x,y,z,'FaceColor','y','FaceAlpha',0.5)
% surf(x,y,z,'FaceColor','g','FaceAlpha',0.5)

% Inclination angle of the slope
alpha = atand( (max(z(:))-min(z(:)))/(x3-x2)); 

P = [x2 , y1 ,z1];
Q = [x2 , y2 ,z1];
R = [x3 , y1 ,z3];




V1 = P-Q;
V2 = P-R;
normal = cross(V1,V2);
syms x y z
D = [x,y,z];
% Equation of the inclined part of ground surface
draw_plane_equation = dot(normal,D-P); 
% Equation of the inclined part of ground surface in terms of z
zplane = solve(draw_plane_equation,z); 
% Store the surface eqaution
ht = matlabFunction(zplane); 

% 
% % Major and Minor Axis of Ellipsid 
L1 = 12;
L2 = 9;
L3 = L1;
% 
n1 = x2;
n2 = x3;

% Center of Ellipsoid
x0_el = 0;
y0_el = 0;
z0_el = 10;
hold on


figure(1)

t = 0:0.1:360;
xS = L1*cosd(t);
yS = L2*sind(t);

plot3(xS,yS,z0_el*ones(size(xS)))





daspect([1 1 1])
xlabel(' x [m] ') 
ylabel(' y [m] ')
zlabel(' z [m] ') 

hold off
%%
figure(2)

[sqrX,sqrY] = meshgrid(n1-5 : n2+5);
z = zeros(size(sqrX));

for k = 1:length(sqrX)
    plot([sqrX(k) sqrX(k)], [sqrY(1) sqrY(end)],'k')
    hold on
    plot([sqrX(1) sqrX(end)], [sqrY(k) sqrY(k)],'k')
    plot(sqrX,sqrY,'k')
    plot(xS,yS,'b','LineWidth',3)
end

xlabel(' x [m] ') 
ylabel(' y [m] ')
title (' Ellipsoidal Slip Surface on X-Y Plane')
%% Generate a matrix to find out the squares intersect with ellipsoid


A(:,:,1) = sqrX;
A(:,:,2) = sqrY;
Square_Matrix = zeros(2,2,2); % Generate a matrix with 4 corners of (x,y)
Square_Row_Number = length(sqrX);


% INITIALS
j =1;
count=0;
i=1;
Red_Squares =0;
Coordinates =0;
n=1;
CoordinatesX_Blue_Squares (1,1:5) = 0;
CoordinatesY_Blue_Squares (1,1:5) = 0;


figure(3)
plot(xS,yS)
hold on
while i
    Square_Matrix = A(i:i+1,j:j+1,1:2);
    j=j+1;
    
    xSquares = [Square_Matrix(2,1,1) Square_Matrix(1,1,1) Square_Matrix(1,2,1)...
        Square_Matrix(2,2,1) Square_Matrix(2,1,1)];
    ySquares = [Square_Matrix(2,1,2) Square_Matrix(1,1,2) Square_Matrix(1,2,2)...
        Square_Matrix(2,2,2) Square_Matrix(2,1,2)];
    
    xSquares_Max = max(xSquares);
    xSquares_Min = min(xSquares);
    
    ySquares_Max = max(ySquares);
    ySquares_Min = min(ySquares);
    
    [in,on] = inpolygon(xS,yS,xSquares,ySquares);
    outx = xS(in);
    outy = yS(in);
    Box_No = num2str(count);
    
    if any( outx<xSquares_Max & outx > xSquares_Min) && any(outy < ySquares_Max...
            & outy > ySquares_Min)
        fill(xSquares,ySquares,'r')
        Red_Squares(end+1) = count;
        CoordinatesX (n,1:5) = xSquares;
        CoordinatesY (n,1:5) = ySquares;
        n=n+1;
        text((xSquares_Min +xSquares_Max)/2,(ySquares_Min+ ySquares_Max)/2,Box_No)
    else
        fill(xSquares,ySquares,'b')
        text(((xSquares_Min +xSquares_Max)/2),((ySquares_Min+ ySquares_Max)/2),Box_No)
        CoordinatesX_Blue_Squares (end+1,1:5) = xSquares(:,1:5);
        
        CoordinatesY_Blue_Squares (end+1,1:5) = ySquares(:,1:5);
    end
    
    hold on
    
    
    if j==Square_Row_Number
        j=1;
        i=i+1;
    end
    
    count = count +1;
    if count == (Square_Row_Number-1)^2
        break
    end
end
xlabel(' x [m] ')
ylabel(' y [m] ')
title ( 'Squares Intersecting with the Ellipse')

CoordinatesX_Blue_Squares (1,:) = [];
CoordinatesY_Blue_Squares (1,:) = [];

%% Table
Red_Squares(1) = [];
No = 1:length(Red_Squares);
for i =1:length(CoordinatesX)
    
    formatSpec = '(%d,%d)';
    x1_y2(i) = convertCharsToStrings(sprintf(formatSpec,CoordinatesX(i,1),CoordinatesY(i,1)));
    x1_y1(i) = convertCharsToStrings(sprintf(formatSpec,CoordinatesX(i,2),CoordinatesY(i,2)));
    x2_y1(i)= convertCharsToStrings(sprintf(formatSpec,CoordinatesX(i,3),CoordinatesY(i,3)));
    x2_y2(i) = convertCharsToStrings(sprintf(formatSpec,CoordinatesX(i,4),CoordinatesY(i,4)));
    
    
end

T=table(No',Red_Squares',x1_y2',x1_y1',x2_y1',x2_y2');
T.Properties.VariableNames = {'Number','No of Effected Sqaures','X1,Y2','X1,Y1','X2,Y1','X2,Y2'};


%% Divison of Effected Squares

figure(4)

for k = 1:length(sqrX)
    plot([sqrX(k) sqrX(k)], [sqrY(1) sqrY(end)],'k')
    hold on
    plot([sqrX(1) sqrX(end)], [sqrY(k) sqrY(k)],'k')
     plot(sqrX,sqrY,'k')
end
plot(xS,yS,'b','LineWidth',1)
hold on
% daspect([1 1 1])

CoordsX(1,1:5) = 0;
CoordsY(1,1:5) = 0;

div = 4;
n = @(n) 1/(2^n);
n = n(div);


k=1;
for i= 1: length(CoordinatesX)
    
    X1_min(i,1) = min(CoordinatesX(i,:));
    X2_max(i,1) = max(CoordinatesX(i,:));
    Y1_min(i,1) = min(CoordinatesY(i,:));
    Y2_max(i,1) = max(CoordinatesY(i,:));
end


while k
    try
        X1= X1_min(k);
        X2= X2_max(k);
        Y1= Y1_min(k);
        Y2= Y2_max(k);
        
        
        
        
        
        [CoordinatesX_add_sqrs...
            ,CoordinatesY_add_sqrs,CoordinatesX_back_to_loop...
            ,CoordinatesY_back_to_loop,DifferenceX] =...
            square_produce(n,X1,X2,Y1,Y2,xS,yS);
        
        CoordinatesX_add_sqrs(1,:)= [];
        CoordinatesY_add_sqrs(1,:)= [];
        CoordinatesX_back_to_loop(1,:) = [];
        CoordinatesY_back_to_loop(1,:) = [];
        
        
        if numel(CoordinatesX_add_sqrs) >=5 && numel(CoordinatesY_add_sqrs) >= 5
            % USE FOR FoS CALCULATION %COORDS
            CoordsX=[ CoordsX ; CoordinatesX_add_sqrs];
            CoordsY=[ CoordsY ; CoordinatesY_add_sqrs];
            %         plot(CoordinatesX_add_sqrs,CoordinatesY_add_sqrs,'+')
            %         pause(0.5)
            %         hold on
        end
        
        if numel(CoordinatesX_back_to_loop) >= 5 && numel(CoordinatesY_back_to_loop) >= 5
            [valmin,locmin] = min(CoordinatesX_back_to_loop');
            [valmax,locmax] = max(CoordinatesX_back_to_loop');
            CoordinatesX =[CoordinatesX; CoordinatesX_back_to_loop];
            X1_min = [X1_min; valmin'];
            X2_max = [X2_max; valmax'];
            [valmin,locmin] = min(CoordinatesY_back_to_loop');
            [valmax,locmax] = max(CoordinatesY_back_to_loop');
            CoordinatesY =[CoordinatesY; CoordinatesX_back_to_loop];
            Y1_min = [Y1_min; valmin'];
            Y2_max = [Y2_max; valmax'];
        end
        
        k = k+1;
    catch
        break
    end
    %        if k==1360
    %            break
    %        end
    
    
    
end

hold on


for i = 1:length(CoordsX)
    plot(CoordsX(i,:),CoordsY(i,:),'r')
    hold on
end

%% Coordinates of All Squares in Ellipsoid

% Blue Squares that has no intersection with ellipsoid (inside)
CoordsX_Blue_in_Ellip(1,1:5) = 0;
CoordsY_Blue_in_Ellip(1,1:5) = 0;
for i = 1 :length(CoordinatesX_Blue_Squares(:,1))
B_x = CoordinatesX_Blue_Squares(i,:);
B_y = CoordinatesY_Blue_Squares(i,:);
xSquares = [ B_x(1)  B_x(2)  B_x(3)  B_x(4)  B_x(5)];
ySquares = [ B_y(1)  B_y(2)  B_y(3)  B_y(4)  B_y(5)];


[in,on] =inpolygon(xSquares,ySquares,xS,yS);
 if sum(in) == 5
     % Blue Squares In the Ellipsoid 2D
     CoordsX_Blue_in_Ellip(end+1,1:5) = xSquares; 
     CoordsY_Blue_in_Ellip(end+1,1:5) = ySquares;
     plot(xSquares(in),ySquares(in),'b','LineWidth',2)
     hold on
 end
end
CoordsX_Blue_in_Ellip(1,:) = [];
CoordsY_Blue_in_Ellip(1,:) = [];

% Red Squares in the Ellipsoid 2D
CoordsX(1,:) = [];
CoordsY(1,:) = [];
% Combine All 2D Coordinates
CoordsX_2D = [CoordsX(:,1:5);CoordsX_Blue_in_Ellip(:,1:5)];
CoordsY_2D = [CoordsY(:,1:5);CoordsY_Blue_in_Ellip(:,1:5)];

%%
% 
 figure(1)
 hold on
 
 
 Below_Surface_X(1,1:5) = 0;
 Below_Surface_Y(1,1:5) = 0;
 Below_Surface_Z(1,1:5) = 0;
 
 for i = 1: length( CoordsX_2D(:,1) )
     z3D  = real( -L3*sqrt(1-(((CoordsX_2D (i,:)-x0_el).^2)/ (L1^2)) - (((CoordsY_2D(i,:) -y0_el).^2)/ (L2^2))))+z0_el;
     
     F = 160*CoordsX_2D(i,:) - 240*z3D + 960; % Equation of Plane
     
     
     if all(F>=0.5 & z3D <= z4)
         plot3 (CoordsX_2D (i,:),CoordsY_2D (i,:),z3D,'b','LineWidth',2)
         Below_Surface_X(end+1,:)= CoordsX_2D (i,:);
         Below_Surface_Y(end+1,:)= CoordsY_2D (i,:);
         Below_Surface_Z(end+1,:)= z3D;
         hold on
     end
    
 end
  CoordsX_2D = Below_Surface_X(2:end,:);
     
  CoordsY_2D = Below_Surface_Y(2:end,:);
     
  CoordsZ_2D = Below_Surface_Z(2:end,:);
  
%% Rotate and Add Coordinates to 3D Slope Stability Model
 
 figure(1)
 hold on
 
 CoordsZ_2D = ones(size(CoordsX_2D)).*z0_el;
 for i = 1 :length(CoordsX_2D(:,1))    
   
     % Square Coordinates in 3D Top
     Columns_3X_Top(i,1:5)= CoordsX_2D(i,:) ;
     Columns_3Y_Top(i,1:5)= CoordsY_2D(i,:) ;
     Columns_3Z_Top(i,1:5)= CoordsZ_2D(i,:);
      plot3(Columns_3X_Top(i,1:5),Columns_3Y_Top(i,1:5),Columns_3Z_Top(i,1:5),'r','LineWidth',3)
     hold on
 end
 daspect([1 1 1])
 hold off


%% Interpolation to the edge of Ellipsoid
% Extract the red squares with 1st and 2nd mimimum width for interpolation
% 1st Minimum
choice = 1;
CoordsX_temp = CoordsX_2D(:,1:4);
CoordsY_temp = CoordsY_2D(:,1:4);
CoordsX_diff = abs(diff(abs(CoordsX_temp), 1, 2));
CoordsY_diff = abs(diff(abs(CoordsY_temp), 1, 2));
minX = min(CoordsX_diff(CoordsX_diff >0),[],1);
minY = min(CoordsY_diff(CoordsY_diff >0),[],1);
[MX,~] = find(CoordsX_diff==minX);
[MY,~] = find(CoordsY_diff==minY);
MX = unique(MX);
MY = unique(MY);
CoordsX_1stminWidth = CoordsX_2D(MX,:);
CoordsY_1stminWidth =CoordsY_2D(MY,:);
% 2nd Minimum
CoordsX_diff_temp = CoordsX_diff(:);
CoordsY_diff_temp = CoordsY_diff(:);
Sort_X = unique(sort(CoordsX_diff_temp(CoordsX_diff_temp~=0)));
Sort_Y = unique(sort(CoordsY_diff_temp(CoordsY_diff_temp~=0)));
[MX,~] = find(CoordsX_diff==Sort_X(2));
[MY,~] = find(CoordsY_diff==Sort_Y(2));
MX = unique(MX);
MY = unique(MY);
CoordsX_2ndminWidth = CoordsX_2D(MX,:);
CoordsY_2ndminWidth = CoordsY_2D(MY,:);

CoordsX_minWidth = [CoordsX_1stminWidth;CoordsX_2ndminWidth];
CoordsY_minWidth = [CoordsY_1stminWidth;CoordsY_2ndminWidth];

figure(1)
hold on
[A_3D ,B_3D ,ztemp]=interpolation(CoordsX_minWidth,CoordsY_minWidth,L1,...
       L2,L3,CoordsX_2D,CoordsY_2D, x0_el,y0_el,z0_el,ht);



%% Factor of Safety Calculation

% Find Angles
figure(1)
hold on
for i = 1: length(CoordsX_2D(:,1))
    x3D = CoordsX_2D(i,1:5);
    y3D = CoordsY_2D(i,1:5);
   
    z3D  = real( -L3*sqrt(1-(((x3D).^2)/ (L1^2)) - (((y3D).^2)/ (L2^2))))+z0_el;
    
    hold on
    
    %Column Edges
    for k = 1 :length(x3D)
   
        x_col = [x3D(k) x3D(k)];
        
        y_col = [y3D(k) y3D(k)];
        if x3D(k)<=x2
        z_plane_point = z1;
        elseif  x3D(k)>=x3
        z_plane_point = z4;
        elseif (x3D(k)>=x2 && x3D(k)<=x3)
        z_plane_point = ht(x3D(k));
        end
        z_column_top(1,k) =  z_plane_point;
        
        z_col = [z_plane_point  z3D(k)];
      
        z_height(k,1) = z_plane_point-z3D(k);
        
       
       
        hold on
        %Column Edges
        plot3(x_col,y_col,z_col,'r','LineWidth',2) 
        hold on
        % Base of Column
        fill3(x3D,y3D,z3D,[1,0,0]) 
  
    end
    
%    fill3(x3D,y3D,z_column_top,[1,0,0])
    % slope angle in x direction
    Horz_L1 = abs(diff(unique(x3D)));
    Vert = abs(diff(z3D(2:3)));
    ax1 = Vert/Horz_L1 ;
    
    Horz_L1 = abs(diff(unique(x3D)));
    Vert = abs(z3D(4)-z3D(1));
    ax2 = Vert/Horz_L1 ;
    
    ax = (ax1+ax2)/2;
    
    % slope angle in y direction
    Horz_L2 = diff(unique(y3D));
    Vert = abs(diff(z3D(1:2)));
    ay1 = Vert/Horz_L2;
    
    Horz_L2 = abs(diff(unique(y3D)));
    Vert = abs(diff(z3D(3:4)));
    ay2 = Vert/Horz_L2;
    
    ay = (ay1+ay2)/2;
    
    % Total Slope
    a_total = sqrt((ax^2)+(ay^2));
    a_total_slopedegree = atand(a_total);
    
    % Unit vector to fin slope in max direction
    e = [1;0;0]; % x direction for this model
    slope_directon = [ax ;ay;0] .* e;
    slope_directon = atand(slope_directon);
    
    % Soil Parameters to use in stress tensor
   
    g = 9.81;
    p = 2000;
    C = 13000;
    phi = 23;
    
    h = sum(z_height)/5;

     hold on
    
    FoS_loc = (tand(phi)/tand(a_total_slopedegree))+ (C/ (p*g*h*cosd(a_total_slopedegree)*sind(a_total_slopedegree)))
    if FoS_loc <= 1
    fill3(x3D,y3D,z_column_top,[0,0,1]) % Top of Column
    hold on
    end
  
end




%     nor_vec = [sind(a_total_slopedegree);0;cosd(a_total_slopedegree)]
%     t_vec = [cosd(a_total_slopedegree);0;-sind(a_total_slopedegree)]
%     sigma = [0 0 0; 0 0 0; 0 0 -p*g*h ];
%     
%     sigman = sum((sigma* nor_vec).*nor_vec)
% %     sigman = sqrt(sigman(1)^2+ sigman(2)^2+sigman(3)^2)
% %     sigmas = abs(sigma* nor_vec- sigman.*nor_vec) 
%      sigmas = sum((sigma*nor_vec).*t_vec)
% 
% %      sigmas = abs(sqrt(sigmas(1)^2+ sigmas(2)^2+sigmas(3)^2))
%     sigmas_crit = C-sigman.*tand(phi)
%     FoS_Local = sigmas_crit ./ sigmas


