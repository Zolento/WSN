function pt = getpt(v,d,pr)%v 频率，单位MHz;d 半径
    pt = pr-20*log10(3e8/(v*1e6*4*pi*d));
end