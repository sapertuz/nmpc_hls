function xRef_out = Quadrotor_xref(t, ~)
    x_0   = [  7, 10,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 
    x_d   = [  0,  1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 
    x_max = [ 15, 15, 15, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]; 
    x_min = [  0,  0,  0, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]; 
    N = length(t);

    xRef_out = zeros(N,12);
    
    % Line reference
%          xRef_out = x_0 + x_d.*t;
    
    % Spiral Reference
        d = 2.5;
        T = 20;
        ang = 2*pi*t/T - pi/2;
        xRef_out(:,1) = x_0(1) + cos(ang)*d; 
        xRef_out(:,2) = x_0(2) + d + sin(ang)*d;
        xRef_out(:,3) = t/8;
        xRef_out(:,6) = t*2*pi/T;

    % Ring with direction reference
%         d = 2.5;
%         T = 20;
%         ang = 2*pi*t/T - pi/2;
%         xRef_out(:,1) = x_0(1) + cos(ang)*d; 
%         xRef_out(:,2) = x_0(2) + d + sin(ang)*d;
%         xRef_out(:,6) = t*2*pi/T;

    index = xRef_out > x_max; 
    for i=1:N
        xRef_out(i,index(i,:)) = x_max(index(i,:));
    end
    index = xRef_out < x_min; 
    for i=1:N
        xRef_out(i,index(i,:)) = x_min(index(i,:));
    end