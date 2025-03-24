function animateQuadrotor(time,x,u,x_ref,axis_lim)

    global drone

    l = drone.l;
    base = x(1:3)';
    q = x(4:6);
    color = [0, 0, 0];
    grosor = 2;

    for i=1:3
        s = sin(q(i));
        c = cos(q(i));
        switch i
            case 1
                Rx = [ 1,  0,  0;
                       0,  c, -s;
                       0,  s,  c];
            case 2
                Ry = [ c,  0,  s;
                       0,  1,  0;
                      -s,  0,  c];
            case 3
                Rz = [ c, -s,  0;
                       s,  c,  0;
                       0,  0,  1];
        end
    end
    A = Rz*Ry*Rx;

    prop(:,1) = base + l*A(:,2); % y
    prop(:,2) = base + l*A(:,1); % x
    prop(:,3) = base - l*A(:,2); % -y
    prop(:,4) = base - l*A(:,1); % -x
    z = A(:,3)*l/4;

    T = u(1);
    tx = u(2);
    ty = u(3);
    tz = u(4);

    clf; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,3,[1,2,4,5]) %{
    hold on; 

    view([50 40])
%     view([10 14])
    xlabel('X','FontSize',14,'FontWeight','bold')
    ylabel('Y','FontSize',14,'FontWeight','bold')
    zlabel('Z','FontSize',14,'FontWeight','bold')
    grid on
    
    if (~exist('axis_lim'))
        axis_lim = [4 10 9 17 -1 3]; 
    end
    axis(axis_lim)

    text(axis_lim(1)+0.5,axis_lim(3)+0.5,axis_lim(6)-0.5, ...
        ['t=',num2str(time),'s'],'FontSize',12)

    dibujar_linea( prop(:,1), prop(:,3), color, grosor )
    dibujar_linea( prop(:,2), prop(:,4), color, grosor )
    for i=1:4
        center = prop(:,i) + z;
        dibujar_linea( prop(:,i),center, color, grosor )
        plotCircle3D(center, A(:,3),l/3, grosor)
    end
    plot3(x_ref(:,1),x_ref(:,2),x_ref(:,3), '-r');
    title('Animation Visualization','FontSize',16)
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes_lim = 100*[-1 1 -1 1];

    subplot(2,3,3)
    plot(tz,T,'o','MarkerSize',30);
    axis(axes_lim)
    xlabel('Yaw','FontSize',12,'FontWeight','bold')
    ylabel('Thrust','FontSize',12,'FontWeight','bold')
    title('Control Input','FontSize',16)

    subplot(2,3,6)
    plot(tx,ty,'o','MarkerSize',30);
    axis(axes_lim)
    xlabel('Roll','FontSize',12,'FontWeight','bold')
    ylabel('Pitch','FontSize',12,'FontWeight','bold')
end

function dibujar_linea( p0, p1, color, grosor )
    % vectores que contienen las componentes de los puntos a dibujar
    x = [ p0(1), p1(1)];
    y = [ p0(2), p1(2)];
    z = [ p0(3), p1(3)];
    % Dibujo de las lineas
    plot3(x,y,z, 'Color', color ,'LineWidth',grosor);
end

function plotCircle3D(center,normal,radius,thick)
    center = center(:);
    center = center';
    normal = normal(:);
    normal = normal';
    theta=0:0.01:2*pi;
    v=null(normal);
    points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    plot3(points(1,:),points(2,:),points(3,:),'r-','LineWidth',thick);
end


