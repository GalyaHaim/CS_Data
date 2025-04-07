function x = pulse(t,width,offset)
    %mdpt = 0.5*(t(1)+t(end));
    x = (heaviside((t-offset+0.5*width)/width)-heaviside((t-offset-0.5*width)/width));
end