function w_ = auto_skew(w)
    if all(size(w)==3) %anti_skew
        w_ = [w(3,2);w(1,3);w(2,1)];
    else %skew
        w_ = [0 -w(3) w(2); w(3) 0 -w(1);-w(2) w(1) 0];
    end
end