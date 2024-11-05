clc; clear; close all;
[mu, m0, m_dot, T, tf, r0, u0, v0, uf] = getParameters();
t = linspace(0, tf, 101);
tol = 1E-10;

solinit = bvpinit(t, @guess);
options = bvpset('RelTol',tol,'Stats','on');
sol = bvp5c(@odefun, @bcfun, solinit, options);

time = sol.x;
state = sol.y([1, 2, 3],:);
costate = sol.y([4, 5, 6],:);
control = atan2(sol.y(5, :), sol.y(6, :));
J_max = state(1, end);

figure(1);
subplot(3, 1, 1);
plot(time, state(1, :));
title("Radial Distance vs. Time");
xlabel("time");
ylabel("r");
grid on;

subplot(3, 1, 2);
plot(time, state(2, :));
title("Radial Velocity vs. Time");
xlabel("time");
ylabel("u");
grid on;

subplot(3, 1, 3);
plot(time, state(3, :));
title("Tangential Velocity vs. Time");
xlabel("time");
ylabel("v");
grid on;

figure(2);
subplot(3, 1, 1);
plot(time, costate(1, :));
title("p_1 vs. Time")
xlabel("time");
ylabel("p_1");
grid on;

subplot(3, 1, 2);
plot(time, costate(2, :));
title("p_2 vs. Time")
xlabel("time");
ylabel("p_2");
grid on;

subplot(3, 1, 3);
plot(time, costate(3, :));
title("p_3 vs. Time");
xlabel("time");
ylabel("p_3");
grid on;

for i = 1:length(control)
    if(control(i) < 0)
        control(i) = control(i) + 2*pi;
    end
end

figure(3);
plot(time, control);
title("φ vs. Time");
xlabel("time");
ylabel("φ");
grid on;

figure(4);
plot(time, control*180/pi);
title("φ (deg) vs. Time");
xlabel("time");
ylabel("φ (deg)");
grid on;

theta = linspace(0, 2*pi, 301);
x1 = state(1, 1)*cos(theta);
y1 = state(1, 1)*sin(theta);
x2 = state(1, end)*cos(theta);
y2 = state(1, end)*sin(theta);

figure(5);
plot(x1, y1, 'b');
title("Plot of Orbits");
grid on;
hold on;
plot(x2, y2, 'r');
hold on;
x_trans = state(1, :).*cos(theta);
y_trans = state(1, :).*sin(theta);
plot(x_trans, y_trans, 'g');
legend("Orbit 1", "Orbit 2", "Transfer Orbit");

[vala, valb] = size(time);
dvalue_increment = (valb-1)/20;
k = 1;
for i = 1:21
    time_inc(i) = time(k);
    control_inc(i) = control(k);
    theta_inc(i) = theta(k);
    x_trans_inc(i) = state(1, k)*cos(theta(k));
    y_trans_inc(i) = state(1, k)*sin(theta(k));
    k = k + dvalue_increment;
end

figure(6);
plot(x1, y1, 'b');
title("Plot of Orbits with Thrust Direction on Transfer");
grid on;
hold on;
plot(x2, y2, 'r');
hold on;
x_trans = state(1, :).*cos(theta);
y_trans = state(1, :).*sin(theta);
plot(x_trans, y_trans, 'g');

for i = 1:21
    a = annotation("arrow");
    a.Parent = gca;
    arrow_length = 0.3;

    a.Position = [x_trans_inc(i), y_trans_inc(i), -arrow_length*sin(control_inc(i)), arrow_length*cos(control_inc(i))];
end
legend("Orbit 1", "Orbit 2", "Transfer Orbit");

function dydt = odefun(t, y)
    [mu, m0, m_dot, T, tf, r0, u0, v0, uf] = getParameters();
    phi = atan2(y(5), y(6));

    dydt = [y(2); y(3)^2/y(1) - mu/y(1)^2 + T*sin(phi)/(m0 - abs(m_dot)*t);...
        -y(2)*y(3)/y(1) + T*cos(phi)/(m0 - abs(m_dot)*t);...
        y(3)^2/y(1)^2*y(5) - 2*mu*y(5)/y(1)^3 - y(2)*y(3)*y(6)/y(1)^2;...
        -y(4) + y(3)/y(1)*y(6); -2*y(3)*y(5)/y(1) + y(2)*y(6)/y(1)];
end

function res = bcfun(ya, yb)
    [mu, m0, m_dot, T, tf, r0, u0, v0, uf] = getParameters();

    res = [ya(1) - r0; ya(2) - u0; ya(3) - v0;...
        yb(4) - 1; yb(2) - uf; yb(3) - sqrt(mu/yb(1))];
end

function v = guess(t)
    [mu, m0, m_dot, T, tf, r0, u0, v0, uf] = getParameters();
    x0 = [r0; u0; v0];
    p0 = [1; 1; 1];
    v = [x0; p0];
end

function [mu, m0, m_dot, T, tf, r0, u0, v0, uf] = getParameters()
    mu = 1;
    m0 = 1;
    m_dot = -0.075;
    T = 0.1405;
    tf = 4;

    r0 = 1;
    u0 = 0;
    v0 = sqrt(mu/r0);
    uf = 0;
end