function validCellBool = isValidCell(xyCurr,mazeMap)
[a,b]=size(mazeMap);
if xyCurr(1)>=1 && xyCurr(1)<=a && xyCurr(2)>=1 && xyCurr(2)<=b
    validCellBool=true;
else
    validCellBool=false;
end

end

