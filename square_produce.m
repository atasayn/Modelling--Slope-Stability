function   [CoordinatesX_add_sqrs...
    ,CoordinatesY_add_sqrs,CoordinatesX_back_to_loop...
    ,CoordinatesY_back_to_loop,DifferenceX] =...
    square_produce(n,X1,X2,Y1,Y2,xS,yS)


Pascal_DivisionX(1,1:2) = [X1 X2];
Pascal_DivisionY(1,1:2) = [Y1 Y2];
j = 1;
i = 1;
count =0;
CoordinatesX_add_sqrs(1,1:5) = 0;
CoordinatesY_add_sqrs(1,1:5) = 0;
CoordinatesX_back_to_loop(1,1:5) =0;
CoordinatesY_back_to_loop(1,1:5) =0;

% Find Mid-Point on X axis of the square each square
DifferenceX = mean([Pascal_DivisionX(1:end-1); Pascal_DivisionX(2:end)]);
Pascal_DivisionX = sort([DifferenceX  Pascal_DivisionX]);


% Find Mid-Point on Y axis of the square each square
DifferenceY = mean([Pascal_DivisionY(1:end-1); Pascal_DivisionY(2:end)]);
Pascal_DivisionY = sort([DifferenceY  Pascal_DivisionY]);

% Make a new square from the new points
[x_int_sqr,y_int_sqr]=meshgrid(Pascal_DivisionX,Pascal_DivisionY);
Square_int_Matrix = zeros(2,2,2);

B(:,:,1) = x_int_sqr;
B(:,:,2) = y_int_sqr;
Square_Row_Number = length(x_int_sqr);
if any(diff(Pascal_DivisionX) == n) || any(diff(Pascal_DivisionX) == -n)
    return
end



while i
    
    
    Square_int_Matrix = B(i:i+1,j:j+1,1:2);
    j=j+1;
    xSquares = [Square_int_Matrix(1,1,1) Square_int_Matrix(2,1,1) Square_int_Matrix(2,2,1)...
        Square_int_Matrix(1,2,1) Square_int_Matrix(1,1,1) ];
    ySquares = [Square_int_Matrix(1,1,2) Square_int_Matrix(2,1,2) Square_int_Matrix(2,2,2)...
        Square_int_Matrix(1,2,2) Square_int_Matrix(1,1,2) ];
    [in,~] = inpolygon(xSquares,ySquares,xS,yS);
    
    if numel(xSquares(in)) >= 5 && numel(ySquares(in)) >= 5
        CoordinatesX_add_sqrs(end+1,1:5) = xSquares(in);
        CoordinatesY_add_sqrs(end+1,1:5) = ySquares(in);
    end
    if numel(xSquares(in)) < 5 && numel(ySquares(in)) < 5
        CoordinatesX_back_to_loop(end+1,1:5)= xSquares;
        CoordinatesY_back_to_loop(end+1,1:5) = ySquares;
        
    end
    
    if j==Square_Row_Number
        j=1;
        i=i+1;
    end
    
    count = count +1;
    if count ==  (Square_Row_Number-1)^2
        
        break
    end
    
end


end





