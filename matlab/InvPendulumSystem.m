function y = InvPendulumSystem(~, x, u, Ts)
    global  counter
    if isempty( counter )
        counter=0; %Initializing counter
    end
    counter = counter + 1;
    
    u = u(:);
    x = x(:);

    y = InvPendStateFcn(x',u);

    if (nargin>3)
        k1 = InvPendStateFcn(x, u);
        x_temp = x + (k1*Ts/2);
        k2 = InvPendStateFcn(x_temp, u);
        x_temp = x + (k2*Ts/2);
        k3 = InvPendStateFcn(x_temp, u);
        x_temp = x + (k3*Ts);
        k4 = InvPendStateFcn(x_temp, u);
        y = x + (Ts/6)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        y = y';
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state_dot = InvPendStateFcn(state,control)
    global invPend

    ml = invPend.m*invPend.l;
    Mm = invPend.M+invPend.m;
    Mn_1 = 1/(invPend.M+invPend.m);
    h_ml = invPend.h/(invPend.m*invPend.l);
    l43 = invPend.l * 4.0/3.0;

    th = state(2);
    c = cos(th);
    s = sin(th);

    Mm_c = Mm/c;

    mls = ml*s;
    mlc = ml*c;
    u_bx1 = control(1)-invPend.b*state(2);
    state4_2 = state(4)*state(4); 

    state_dot(1) = state(2);                                                 
    state_dot(3) = state(4);                                                 
    state_dot(4) = (u_bx1 + mls*state4_2 - Mm_c*(invPend.g*s - (h_ml)*state(4))) / (mlc-Mm_c*l43); 
    state_dot(2) = (Mn_1)*(u_bx1 - mlc*state_dot(4) + mls*state4_2); 

    state_dot = state_dot(:);