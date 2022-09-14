function bool = detectConflict(reqInfo, reserveInfo)

crossMul = (reqInfo(2)-reserveInfo(1)) * (reqInfo(1)-reserveInfo(2));

if crossMul <= 0
    bool = true;
else
    bool = false;
end

end