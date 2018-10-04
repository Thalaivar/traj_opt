function Vw = wind_model(h, p)
    VR = p(1); hR = p(2); a = p(3);
    Vw = VR*h;
end