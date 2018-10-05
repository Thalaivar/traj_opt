function Vw = wind_model(h, p)
    VR = p(1); hR = p(2); a = p(3); choice = p(4);
    if choice == 1
        Vw = VR*h;
    else
        Vw = VR*(h/hR)^a;
    end
end