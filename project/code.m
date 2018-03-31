% **********************************  Walking Beam Indexer Mechanism ***************************** %
%{
The Code is written based on some calculations done by hand and is not
shown or given in the code. The code is written for Kinematic Analysis of Walking Beam
Indexer Mechanism. Units - 
length - m.
velocity - m/s.
acceleration - m/s^2. 
%}

% Link lengths declared.
l1 = 0.037;
l2 = 0.087;
l3 = 0.094;
l4 = 0.1244;
l5 = 0.1276;
l7 = 0.1288;
l9 = 0.26;

% Initial angular velocities declared
w1 = -150*(2*pi/60);
w2 = 10;
w3 = 10;
al2 = 10;
al3 = 10;

% Fixed angles between link lengths declared
t4 = 180 - ((180/pi)*atan(240/198));
a = acos((l2^2 + l5^2 - (0.043)^2)/(2*l2*l5))*(180/pi);
b = acos((l5^2 + l7^2 - (0.073)^2)/(2*l5*l7))*(180/pi);
theta = [1:2:360];

% Main for loop for analysis
for t = 1:2:360

  % Displacement Analysis
  A = (l1*cosd(t) - l4*cosd(t4));
  B = (l1*sind(t) - l4*sind(t4));
  C = (l1*l4*cosd(t - t4)/l2) + ((l3^2 - l1^2 - l2^2 - l4^2)/(2*l2));
  C1 = (l1*l4*cosd(t - t4)/l3) + ((l2^2 - l1^2 - l3^2 - l4^2)/(2*l3));
  x1 = (B - sqrt(A^2 + B^2 - C1^2))/(A + C1);
  x = (B - sqrt(A^2 + B^2 - C^2))/(A + C);
  t2 = 2*atan(x)*(180/pi);
  t3 = 90 - 2*atan(x1)*(180/pi);
  if t2 <= 0
    t2 = 360 + t2;
  end
  if t3 <= 0
    t3 = 360 + t3;
  end

  % Co-ordinates of hinges
  x = l1*cosd(t) + l5*cosd(a + t2);
  y = l1*sind(t) + l5*sind(a + t2);
  x1 = l1*cosd(t);
  y1 = l1*sind(t);
  x2 = x1 + l2*cosd(t2);
  y2 = y1 + l2*sind(t2);
  x3 = x1 + l7*cosd(t2+a+b);
  y3 = y1 + l7*sind(t2+a+b);
  
  % Velocity Analysis
  J = [l2*sind(t2) l3*sind(t3); l2*cosd(t2) l3*cosd(t3)];
  W = [w2; w3];
  compX = W;
  newX = [0; 0];
  while((abs(compX(1) - newX(1)) >= 0.001)||(abs(compX(2) - newX(2)) >= 0.001))
    F = [w1*l1*sind(t)+W(1)*l2*sind(t2)+W(2)*l3*sind(t3); w1*l1*cosd(t)+W(1)*l2*cosd(t2)+W(2)*l3*cosd(t3)];
    newX = W - inv(J)*F;
    compX = W;
    W = newX;
  end
  % Velocity of output points on ternary joint.
  vx = -w1*l1*sind(t) - newX(1)*l5*sind(a + t2);
  vy = w1*l1*cosd(t) + newX(2)*l5*cosd(a + t2);
  
  % Acceleration analysis
  Al = [al2; al3];
  compA = Al;
  newA = [0; 0];
  while((abs(compA(1) - newA(1)) >= 0.001)||(abs(compA(2) - newA(2)) >= 0.001))
    F1 = [w1^2*l1*cosd(t) + l2*W(1)^2*cosd(t2) + l2*Al(1)*sind(t2) + l3*W(2)^2*cosd(t3) + l3*sind(t3)*Al(2); -w1^2*l1*sind(t) - l2*W(1)^2*sind(t2) + l2*Al(1)*cosd(t2) - l3*W(2)^2*sind(t3) + l3*cosd(t3)*Al(2)];
    newA = Al - inv(J)*F1;
    compA = Al;
    Al = newA;
  end
  % acceleration of ouput point on ternary joint
  ax = -w1^2*l1*cosd(t) - newX(1)^2*l5*cosd(a + t2) - newA(1)*l5*sind(a + t2);
  
  % Arrary Concatination and Formation
  if t == 1
    X1 = x1;
    Y1 = y1;
    X2 = x2;
    Y2 = y2;
    X3 = x3;
    Y3 = y3;
    Vx = vx;
    Vy = vy;
    Ax = ax;
    W2 = newX(1);
    W3 = newX(2);
    A2 = newA(1);
    A3 = newA(2);
    theta3 = t3;
    theta2 = t2;
    X = x;
    Y = y;
  else
    X1 = [X1 x1];
    Y1 = [Y1 y1];
    X2 = [X2 x2];
    Y2 = [Y2 y2];
    X3 = [X3 x3];
    Y3 = [Y3 y3];
    Vx = [Vx vx];
    Vy = [Vy vy];
    Ax = [Ax ax];
    W2 = [W2 newX(1)];
    W3 = [W3 newX(2)];
    A2 = [A2 newA(1)];
    A3 = [A3 newA(2)];
    theta3 = [theta3 t3];
    theta2 = [theta2 t2];
    X = [X x];
    Y = [Y y];
  end
end

% Maximum x-component velocity and acceleration points on coupler curve
%{
% Uncomment to run code
maxVx = max(Vx);
maxAx = max(Ax);
for n = 1:1:180
  if Ax(n) == maxAx
    maxn_a = n;
  end
  if Vx(n) == maxVx
    maxn_v = n;
  end
end
maxVxpt_x = X(maxn_v)
maxVxpt_y = Y(maxn_v)
maxAxpt_x = X(maxn_a)
maxAxpt_y = Y(maxn_a)
%}

theta3 = theta3 - 180;

% Maximum push distance of job (Uncomment to run code)
%Max_Push_distance = abs(max(X) - min(X))

% Co-ordinates of Hinges.
X4 = X;
Y4 = Y + 0.10;
X5 = X4 + l9;
Y5 = Y4;
X6 = X + l9;
Y6 = Y;
X7 = X2 + l9;
Y7 = Y2;
X8 = X3 + l9;
Y8 = Y3;

% Animation Code (Uncomment to see animation)
%{
for th=1:1:180
    plot(X,Y);hold on; % Uncomment to see how output point traces the coupler curve
    plot([0 X1(th)], [0 Y1(th)],'ro-');hold on;
    plot([X1(th) X2(th)], [Y1(th) Y2(th)], 'ro-'); hold on;
    plot([X2(th) -0.079], [Y2(th) 0.096], 'ro-'); hold on;
    plot([X2(th) X(th)], [Y2(th) Y(th)], 'ro-'); hold on;
    plot([X2(th) X3(th)], [Y2(th) Y3(th)], 'ro-'); hold on;
    plot([X(th) X4(th)], [Y(th) Y4(th)], 'ro-'); hold on;
    plot([X4(th) X5(th)], [Y4(th) Y5(th)], 'ro-'); hold on;
    plot([X5(th) X6(th)], [Y5(th) Y6(th)], 'ro-'); hold on;
    plot([X6(th) X7(th)], [Y6(th) Y7(th)], 'ro-'); hold on;
    plot([X7(th) -0.079+l9], [Y7(th) 0.096], 'ro-'); hold on;
    plot([X7(th) X8(th)], [Y7(th) Y8(th)], 'ro-'); hold on;
    plot([X8(th) X3(th)], [Y8(th) Y3(th)], 'ro-'); hold off;
    axis([-0.20 0.25 -0.15 0.15]);
    pause(0.0001);
end
%}

%   ************************ Results ***************************

% Displacement Analysis (Uncomment to see graph)
%plot(theta, theta2); % Coupler displacement
%plot(theta, theta3); % Rocker displacement
%plot(X, Y);  % Coupler Curve

% Velocity analysis
%plot(theta, W2); % Coupler angular velocity
%plot(theta, W3); % Rocker angular velocity
%plot(theta, Vx); % Coupler output point x-component velocity

% Acceleration Analysis
%plot(theta, A2); % Coupler angular acceleration
%plot(theta, A3); % Rocker angular acceleration
%plot(theta, Ax); % Coupler output point acceleration

% Points of maximum x-component acceleration and velocity on coupler curve.
%{
plot(X, Y);
text(maxVpt_x,maxVpt_y,"o");
text(maxApt_x,maxApt_y,"o");
%}