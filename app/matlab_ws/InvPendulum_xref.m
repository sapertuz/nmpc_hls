function xRef_out = InvPendulum_xref(t, ~)
    x_0   = [  0, pi, 0, 0]; 
    x_d   = [  0, 1, 0, 0]; 
    x_max = [ 1, Inf, Inf, Inf]; 
    x_min = [ -1, -Inf, -Inf, -Inf]; 
    N = length(t);

    xRef_out = zeros(N,4);
    
    index = t < 60;
    xRef_out(index,1) = .2;
    index = not(index);
    xRef_out(index,1) = -.1; 
    
    index = xRef_out > x_max; 
    for i=1:N
        xRef_out(i,index(i,:)) = x_max(index(i,:));
    end
    index = xRef_out < x_min; 
    for i=1:N
        xRef_out(i,index(i,:)) = x_min(index(i,:));
    end