% PROBLEMA DE LOS TRES CUERPOS
% Luisa Zamudio <luisa.zamudio@uabc.edu.mx>
% 02/06/2021

clear, clf, clf

%% CONDICIONES DE SIMULACION 
t = 0;       % Tiempo inicial
dt = 1/365;  % Paso temporal (Dias)
tMax = 3;    % Tiempo final (Anios)
G = 1;       % Constante gravitacional universal unitaria


%% CONDICIONES INICIALES
% Masas y velocidades y posiciones iniciales de los cuerpos
% Masas en terminos de masas terrestres donde 1 MT = 6.0e24*kg
% Velocidades en terminos de km/s
% Posiciones en terminos de unidades astronomicas donde 1 UA = 150e+06*km

%% Caso tres cuerpos en vertices de un triangulo y a velocidad cero
% m0 = 4.0; m1 = 4.0; m2 = 4.0;
% 
% vx0 = 0.0; vy0 = 0.0;       x0 = 0.0;  y0 = 1.0;
% vx1 = 0.0; vy1 = 0.0;       x1 = -0.5; y1 = 0;
% vx2 = 0.0; vy2 = 0.0;       x2 = 0.5;  y2 = 0;

%% Caso especial #1 tres cuerpos 
% m0 = 4.0; m1 = 4.0; m2 = 5.0;
% 
% vx0 = 0.0; vy0 = 0.0;       x0 = 1.0;  y0 = -1.0;
% vx1 = 0.0; vy1 = 0.0;       x1 = -0.5; y1 = 1.0;
% vx2 = 0.0; vy2 = 0.0;       x2 = -2.0; y2 = -1.2;

%% Caso especial #2 tres cuerpos 
% m0 = 4.0; m1 = 4.0; m2 = 5.0;
% 
% vx0 = 0.0; vy0 = 0.0;       x0 = 1.0;  y0 = -1.0;
% vx1 = 0.0; vy1 = 0.0;       x1 = -0.5; y1 = 1.0;
% vx2 = 0.0; vy2 = 0.0;       x2 = -0.5; y2 = -3.5;

%% Caso especial #3 tres cuerpos
% m0 = 4.0; m1 = 3.0; m2 = 4.0;
%  
% vx0 = 0.0; vy0 = 0.0;       x0 = 0.0;  y0 = -0.1;
% vx1 = 0.0; vy1 = 0.0;       x1 = 2.0;  y1 = 2.0;
% vx2 = 0.0; vy2 = 0.0;       x2 = 5.0;  y2 = 0.0;

%% Caso especial #4 tres cuerpos
% m0 = 4.0; m1 = 4.0; m2 = 4;
% 
% vx0 = (1/2)*pi;   vy0 = -(1/2)*pi;      x0 = -1/2; y0 = 0.0;
% vx1 = (1/2)*pi;   vy1 = (1/2)*pi;       x1 = 1/2;  y1 = 0.0;
% vx2 = -(1/2)*pi;  vy2 = -(1/2)*pi;      x2 = 0.0;  y2 = sqrt(3)/2;

%% Solucion #1 tres cuerpos
 m0 = 4.0; m1 = 4.0; m2 = 4.0;
% Si la masa intermedia m2 disminuye las orbitas se vuelven mas amplias,
% mientras que si la masa aumenta, las orbitas se vuelven mas cerradas y 
% alargadas, y a mayor tiempo vemos el desplazamiento debido al movimiento
% de la masa intermedia

vx0 = (1/2)*pi;  vy0 = (1/2)*pi;       x0 = -1.0; y0 = 0.0;
vx1 = -(1/2)*pi; vy1 = -(1/2)*pi;      x1 = 1.0;  y1 = 0.0;
vx2 = 0.0;       vy2 = 0.0;            x2 = 0.0;  y2 = 0.0;

%% Solucion #2 tres cuerpos
% m0 = 4.0; m1 = 4.0; m2 = 4.0;
%  
% vx0 = 0.0;       vy0 = (1/2)*pi;       x0 = -1.0; y0 = 0.0;
% vx1 = 0.0;       vy1 = (1/2)*pi;       x1 = 1.0;  y1 = 0.0;
% vx2 = -(1/2)*pi; vy2 = -(1/2)*pi;      x2 = 0.0;  y2 = 0.0;

% Arreglos de posicion y velocidad de m0
X0=[x0];   Y0=[y0];
VX0=[vx0]; VY0=[vy0];

% Arreglos de posicion y velocidad de m1
X1=[x1];   Y1=[y1];
VX1=[vx1]; VY1=[vy1];

% Arreglos de posicion y velocidad de m2
X2=[x2];   Y2=[y2];
VX2=[vx2]; VY2=[vy2];

%% METODO DE RUNGE-KUTTA ORDEN 2
% Donde el metodo obtiene el valor actualizado de acuerdo con la expresion
% x(t+dt) = x(t) + f(x',t')dt
% Donde:
% x' = x(t) + 1/2f(x(t),t)dt
% t' = t + 1/2dt

while t < tMax 
    %% Calculo previo de las posiciones (x,y)
    % Posiciones previas de m0
    x0p = x0 + (1/2)*vx0*dt;  
    y0p = y0 + (1/2)*vy0*dt;  
    % Posiciones previas de m1
    x1p = x1 + (1/2)*vx1*dt;  
    y1p = y1 + (1/2)*vy1*dt;  
    % Posiciones previas de m2
    x2p = x2 + (1/2)*vx2*dt;  
    y2p = y2 + (1/2)*vy2*dt;  
    
    %% Calculo previo de las velocidades (vx,vy)
    % Velocidades previas de m0
    vx0p = vx0 - (1/2)*G*m2*((x0-x2)/(sqrt(abs(x2-x0)^2+abs(y2-y0)^2))^3)*dt - (1/2)*G*m1*((x0-x1)/(sqrt(abs(x0-x1)^2+abs(y0-y1)^2))^3)*dt;  % velocidad en x de m0
    vy0p = vy0 - (1/2)*G*m2*((y0-y2)/(sqrt(abs(x2-x0)^2+abs(y2-y0)^2))^3)*dt - (1/2)*G*m1*((y0-y1)/(sqrt(abs(x0-x1)^2+abs(y0-y1)^2))^3)*dt;  % velocidad en y de m0
    % Velocidades previas de m1
    vx1p = vx1 - (1/2)*G*m2*((x1-x2)/(sqrt(abs(x2-x1)^2+abs(y2-y1)^2))^3)*dt - (1/2)*G*m0*((x1-x0)/(sqrt(abs(x0-x1)^2+abs(y0-y1)^2))^3)*dt;  % velocidad en x de m1
    vy1p = vy1 - (1/2)*G*m2*((y1-y2)/(sqrt(abs(x2-x1)^2+abs(y2-y1)^2))^3)*dt - (1/2)*G*m0*((y1-y0)/(sqrt(abs(x0-x1)^2+abs(y0-y1)^2))^3)*dt;  % velocidad en y de m1
    % Velocidades previas de m2
    vx2p = vx2 - (1/2)*G*m0*((x2-x0)/(sqrt(abs(x2-x0)^2+abs(y2-y0)^2))^3)*dt - (1/2)*G*m1*((x2-x1)/(sqrt(abs(x2-x1)^2+abs(y2-y1)^2))^3)*dt;  % velocidad en x en m2
    vy2p = vy2 - (1/2)*G*m0*((y2-y0)/(sqrt(abs(x2-x0)^2+abs(y2-y0)^2))^3)*dt - (1/2)*G*m1*((y2-y1)/(sqrt(abs(x2-x1)^2+abs(y2-y1)^2))^3)*dt;  % velocidad en y en m2
    
    %% Calculo de las posiciones actualizadas (x,y)
    % Posiciones actualizadas de m0
    x0 = x0 + vx0p*dt;  
    y0 = y0 + vy0p*dt;  
    % Posiciones actualizadas de m1
    x1 = x1 + vx1p*dt;  
    y1 = y1 + vy1p*dt;  
    % Posiciones actualizadas de m2
    x2 = x2 + vx2p*dt;  
    y2 = y2 + vy2p*dt;  
    
    %% Calculo de las velocidades actualizadas (vx,vy)
    % Velocidades actualizadas de m0
    vx0 = vx0 - G*m2*((x0p-x2p)/(sqrt(abs(x2p-x0p)^2+abs(y2p-y0p)^2))^3)*dt - G*m1*((x0p-x1p)/(sqrt(abs(x0p-x1p)^2+abs(y0p-y1p)^2))^3)*dt;  
    vy0 = vy0 - G*m2*((y0p-y2p)/(sqrt(abs(x2p-x0p)^2+abs(y2p-y0p)^2))^3)*dt - G*m1*((y0p-y1p)/(sqrt(abs(x0p-x1p)^2+abs(y0p-y1p)^2))^3)*dt; 
    % Velocidades actualizadas de m1
    vx1 = vx1 - G*m2*((x1p-x2p)/(sqrt(abs(x2p-x1p)^2+abs(y2p-y1p)^2))^3)*dt - G*m0*((x1p-x0p)/(sqrt(abs(x0p-x1p)^2+abs(y0p-y1p)^2))^3)*dt;   
    vy1 = vy1 - G*m2*((y1p-y2p)/(sqrt(abs(x2p-x1p)^2+abs(y2p-y1p)^2))^3)*dt - G*m0*((y1p-y0p)/(sqrt(abs(x0p-x1p)^2+abs(y0p-y1p)^2))^3)*dt;    
    % Velocidades actualizadas de m2
    vx2 = vx2 - G*m0*((x2p-x0p)/(sqrt(abs(x2p-x0p)^2+abs(y2p-y0p)^2))^3)*dt - G*m1*((x2p-x1p)/(sqrt(abs(x2p-x1p)^2+abs(y2p-y1p)^2))^3)*dt;
    vy2 = vy2 - G*m0*((y2p-y0p)/(sqrt(abs(x2p-x0p)^2+abs(y2p-y0p)^2))^3)*dt - G*m1*((y2p-y1p)/(sqrt(abs(x2p-x1p)^2+abs(y2p-y1p)^2))^3)*dt;  
    
    %% Aumento de los arreglos de posicion y velocidad
    X0=[X0,x0]; Y0=[Y0,y0];
    X1=[X1,x1]; Y1=[Y1,y1];
    X2=[X2,x2]; Y2=[Y2,y2];
    VX0 = [VX0, vx0]; VY0 = [VY0, vy0];
    VX1 = [VX1, vx1]; VY1 = [VY1, vy1];
    VX2 = [VX2, vx2]; VY2 = [VY2, vy2];
  
    %% Aumento del tiempo
    t = t + dt;
 
    %% Calculo de los limites para los ejes de la figura
    xmin = min([X0 X1 X2]);
    xmax = max([X0 X1 X2]);
    ymin = min([Y0 Y1 Y2]);
    ymax = max([Y0 Y1 Y2]);
    
    %% Graficacion de las trayectorias en tiempo real
    figure(1)
    plot (X0,Y0,'-b',X1,Y1,'-r',X2,Y2,'-g','LineWidth',1.2)
    axis([xmin-1/2 xmax+1/2 ymin-1/2 ymax+1/2])
    drawnow
end

%% GRAFICACION DE LOS RESULTADOS
plot (X0,Y0,'-b',X1,Y1,'-r',X2,Y2,'-g','LineWidth',1.2)
hold on
plot(X0(end),Y0(end),'o','markerfacecolor','b','markersize',8)
hold on
plot(X1(end),Y1(end),'o','markerfacecolor','r','markersize',8)
hold on
plot(X2(end),Y2(end),'o','markerfacecolor','g','markersize',8)
hold on
plot(X0(1),Y0(1),'s','markerfacecolor','b','markersize',4)
hold on
plot(X1(1),Y1(1),'s','markerfacecolor','r','markersize',4)
hold on
plot(X2(1),Y2(1),'s','markerfacecolor','g','markersize',4)
axis([xmin-1/2 xmax+1/2 ymin-1/2 ymax+1/2])
hold off
grid on, grid minor
title('Problema de los tres cuerpos')
xlabel('distancia horizontal (UA)')
ylabel('distancia vertical (UA)')
legend({'m1','m2','m3'},'Location','northeast')
