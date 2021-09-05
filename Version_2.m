clc
clear
close all
% Ellipsoid Parameters

figure (1)
% % Major and Minor Axis of Ellipsid 
L1 = 8;
L2 = 6;
L3 = L1;
% 
n1 = -10;
n2 = 10;

% Center of Ellipsoid
xc = 0;
yc = 0;
zc = 8;

z_el=@(x,y) [-1 +1].*(L3*sqrt(1-((x-xc).^2/L1.^2)-((y-yc).^2/L2.^2)))+zc;
t = 0:0.1:360;

n1 = -10;
n2 = 10;

%  3D x-z plane
xS = L1*cosd(t)+xc ;
yS = L2*sind(t)+yc ;
zS = L3*sind(t)+zc ;
plot3(xS,ones(size(xS))*yc,zS,'r')
axis equal
hold on
plot3(xS,yS,ones(size(xS)).*zc,'r')
hold on
% plot3(xc*ones(size(xS)),yS,zS,'r')


% Surface Grid
[x1,y1] = meshgrid(n1:0.1:n2);
tt =linspace(n1,n2,length(x1));
zg = repmat(linspace(0,8,length(tt)),length(x1),1);
rng('default')
zg = zg + rand(size(zg));
T = griddedInterpolant(x1',y1',zg');
surf(x1,y1,zg,'EdgeColor','none');

grid on
xlabel(' x [m]')
ylabel(' y [m]')
zlabel(' z [m]')
axis equal
hold off
hold on




%% Extract Points 


t=1;
j=1;


Top= [];
Bottom= [];
X = [];
Y= [];
 




while true
    
    Range_X =  x1(t,j);
    Range_Y =  y1(t,j);
    temp = z_el(Range_X,Range_Y);
    
    if isreal(temp(1))
        if T(Range_X,Range_Y) >= temp(1)
            Top(end+1)   = T(Range_X,Range_Y);
            Bottom(end+1)= temp(1);
            X(end+1) = Range_X;
            Y(end+1) = Range_Y;
            
        end
    end
    if isreal(temp(2))
        if T(Range_X,Range_Y) >= temp(2)
            Top(end+1)   = T(Range_X,Range_Y);
            Bottom(end+1)= temp(2);
            X(end+1) = Range_X;
            Y(end+1) = Range_Y;
        end
    end
    
    
    if j == length(T.Values) && t== length(T.Values)
        break
    end
    
    if j ==length(T.Values)
        j=1;
        t=t+1;
    end
    j = j+1;
end


%%

% To seperate elements of z_el individually 
indexat = @(z_el, index) z_el(index); % to handle the indexes of the output
test1 = @(x,y,index) indexat(z_el(x,y), index); 
test  = @(x,y,i) abs(T(x,y)-test1(x,y,i)); % abs value of the difference between
% the topography and the intersecting point of the the same vertical
% distance on the slip cirle
temp  = [];

% Along X direction , find the points where y starts changing. from max
% befinning and from min beginning
for i = min(min(x1)):0.1: max(max(x1))
    
    % the minimum vertical distance at x certain point
    tempy1 = fminsearch(@(y) test(i,y,1),min(min(y1))) ;
    % Start from bottom
    if isreal(test(i,tempy1,1)) && test1(i,tempy1,1) <= T(i,tempy1) &&...% topo must be higher than the point below
            T(i,tempy1)<= test1(i,tempy1,2) % topo must be lower than the point above
        temp(end+1,:) = [i,tempy1,test1(i,tempy1,1)];
    end

    tempy2 = fminsearch(@(y) test(i,y,1),max(max(y1))) ;
    % Start from above
    if isreal(test(i,tempy2,1)) && test1(i,tempy2,1) <= T(i,tempy2) &&...
            T(i,tempy2)<= test1(i,tempy2,2)
        temp(end+1,:) = [i,tempy2,test1(i,tempy2,1)];
    end
    
     tempy3 = fminsearch(@(y) test(i,y,2),min(min(y1))) ;
    if isreal(test(i,tempy3,2)) && test1(i,tempy3,1) <= T(i,tempy3) &&...
            T(i,tempy3)<= test1(i,tempy3,2)
        temp(end+1,:) = [i,tempy3,test1(i,tempy3,2)];
    end
     
     tempy4 = fminsearch(@(y) test(i,y,2),max(max(y1))) ;
    if isreal(test(i,tempy4,2)) && test1(i,tempy4,1) <= T(i,tempy4) &&...
            T(i,tempy4)<= test1(i,tempy4,2)
        temp(end+1,:) = [i,tempy4,test1(i,tempy4,2)];
    end
end
figure(2)
plot(temp(:,1),temp(:,2))
xlabel(' x [m] ') 
ylabel(' y [m] ')
grid on


%%
temp_change=temp;
temp_new=temp(1,:);
a = temp_change(1,1:2);
k = dsearchn(temp_change(2:end,1:2),a);
temp_new(end+1,:) = temp_change(1+k,:);
temp_change(1:2,:) = [];

for i = 1:numel(temp(:,1))-4
    k = dsearchn(temp_change(:,1:2),temp_new(end,1:2));
    temp_new(end+1,:) = temp_change(k,:);
    temp_change(k,:) = [];
end
temp_new = [temp_new; temp_new(1,:)];
hold on
plot3(temp_new(:,1),temp_new(:,2),temp_new(:,3))
xS = temp_new(:,1);
yS = temp_new(:,2);
zS = temp_new(:,3);
figure(3)
plot(xS,yS,'b')
xlabel(' x [m] ') 
ylabel(' y [m] ')
grid on

%%
figure(3)


[sqrX,sqrY] = meshgrid(n1 : n2);
z = zeros(size(sqrX));

for k = 1:length(sqrX)
    plot([sqrX(k) sqrX(k)], [sqrY(1) sqrY(end)],'k')
    hold on
    plot([sqrX(1) sqrX(end)], [sqrY(k) sqrY(k)],'k')
    plot(sqrX,sqrY,'k')
  
end

hold on
plot(xS,yS,'r')
daspect([1 1 1])
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


figure(4)
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
% SQUARES 1
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

T2=table(No',Red_Squares',x1_y2',x1_y1',x2_y1',x2_y2');
T2.Properties.VariableNames = {'Number','No of Effected Sqaures','X1,Y2','X1,Y1','X2,Y1','X2,Y2'};

%% Divison of Effected Squares

figure(5)

for k = 1:length(sqrX)
    plot([sqrX(k) sqrX(k)], [sqrY(1) sqrY(end)],'k')
    hold on
    plot([sqrX(1) sqrX(end)], [sqrY(k) sqrY(k)],'k')
     plot(sqrX,sqrY,'k')
end
plot(xS,yS,'b','LineWidth',1)
hold on
daspect([1 1 1])

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





%% Factor of Safety Calculation

% Find Angles
figure(1)
hold on

 sigma_Crit = [];
    sigmass = [];
for i =100%:length(CoordsX_2D(:,1))
    Col_Top = [];
    Col_Bot = [];
    Col_Height = [];
    x3D = CoordsX_2D(i,1:5);
    y3D = CoordsY_2D(i,1:5);
    hold on
    
    %Column Edges
    for k = 1 :length(x3D)
   
        x_col = [x3D(k) x3D(k)];
        
        y_col = [y3D(k) y3D(k)];
        
        Col_Top(end+1) = T(x3D(k),y3D(k));
        Col_Bot(end+1)= min(z_el(x3D(k),y3D(k)));
        Col_Height(end+1) = Col_Top(end)- Col_Bot(end);

  
        
        %Column Edges
       plot3(x_col,y_col,[Col_Top(end) Col_Bot(end)],'r','LineWidth',2) 
        hold on
        
    end
    
    % Column Top
    fill3(x3D,y3D,Col_Top,[1,0,0])
    hold on
    % Column Bottom
    fill3(x3D,y3D,Col_Bot,[1,0,0])
   
    % slope angle in x direction
    Horz_L1 = abs(diff(unique(x3D)));
    Vert = abs(diff(Col_Bot(2:3)));
    ax1 = Vert/Horz_L1 ;
    
    Horz_L1 = abs(diff(unique(x3D)));
    Vert = abs(Col_Bot(4)-Col_Bot(1));
    ax2 = Vert/Horz_L1 ;
    ax = (ax1+ax2)/2;
    ax_ang = atand(ax);
    
    % slope angle in y direction
    Horz_L2 = diff(unique(y3D));
    Vert = abs(diff(Col_Bot(1:2)));
    ay1 = Vert/Horz_L2;
    
    
    Horz_L2 = abs(diff(unique(y3D)));
    Vert = abs(diff(Col_Bot(3:4)));
    ay2 = Vert/Horz_L2;
    ay = (ay1+ay2)/2;
    ay_ang = atand(ay);
    Length = [  Horz_L1 Horz_L2];
    % Total Slope
    a_total = sqrt((ax^2)+(ay^2));
    a_total_slopedegree = atand(a_total);
    g = 9.81; % m/s2
    p = 2000; % kg/m3
    C = 1300; % kg/m2
    phi = 23;
    
    h = sum(Col_Height)/5
   

    
        cosdip = (1+(tand(ax_ang)^2)+(tand(ay_ang)^2))^(-1/2);
        W = p*g*h.*Length;
        A = Length.*(((1-(sind(ax_ang)^2)*(sind(ay_ang)^2))^(1/2))/cosd(ax_ang)*cosd(ay_ang))
        % Soil Parameters to use in stress tensor
        sigma_Crit = (C*A+W*cosdip*tand(phi))
        sigmass = W*sind(ax_ang)

    f_loc = sum(sigma_Crit)/sum(sigmass)
    
hold on
    
    FoS_loc = (tand(phi)/tand(a_total_slopedegree))+ (C/ (p*g*h*cosd(a_total_slopedegree)*sind(a_total_slopedegree)))
    if FoS_loc <= 1
    fill3(x3D,y3D,Col_Top,[0,0,1]) % Top of Column
    hold on
    end
  
end
 FoS_loc = (tand(phi)/tand(a_total_slopedegree))+ (C/ (p*g*h*cosd(a_total_slopedegree)*sind(a_total_slopedegree)))