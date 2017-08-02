function validCellBool = isValidCell3D(xyzCurr,mazeMap)
[a,b,c]=size(mazeMap);
if xyzCurr(1)>=1 && xyzCurr(1)<=a && xyzCurr(2)>=1 && xyzCurr(2)<=b &&...
        xyzCurr(3)>=1 && xyzCurr(3)<=c
    validCellBool=true;
else
    validCellBool=false;
end

end

