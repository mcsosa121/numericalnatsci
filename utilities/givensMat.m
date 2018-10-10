function G = givensMat( a,b )
    [c,s,~] = givensRot(a,b);
    G = [c -s; s c]; 
end

