function a = acc(t)
    t = t/60;
    a0 = 0.4/3600; % 0.4 km/min^2 -> km/s^2
    if t>=15 && t<17
        a = a0;
    elseif t>=17 && t<19
        a = -a0;
    elseif t>=30 && t<32
        a = a0;
    elseif t>=32 && t<34
        a = -a0;
    elseif t>=45 && t<47
        a = a0;
    elseif t>=47 && t<49
        a = -a0;
    elseif t>=60 && t<62
        a = a0;
    elseif t>=62 && t<64
        a = -a0;
    else
        a= 0;
    end
end