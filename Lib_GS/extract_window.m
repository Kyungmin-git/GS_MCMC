function W = extract_window(x, t, t0, t1)
    idx = (t>=t0) & (t<=t1);
    W.x = x(idx);
    W.t = t(idx) - t(find(idx,1,'first'));
end
