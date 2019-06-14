
T = 50;
Q = sqrt(10);
R = 1;

x = zeros(1,T);
y = zeros(1,T);
x_hat = zeros(1,T);
P = zeros(1,T);
J = zeros(1, T);

x(1) = 10;
y(1) = 10;
x_hat(1) = 10;
P(1) = 5;
J(1) = 1/2 - 25*(x(1)^2-1)/(1+x(1)^2)^2;

num_sig_pts = 2*height(x)+1;

for t = 2:T
    %% Actual System
    x(t) = x(t-1)/2 + 25*x(t-1)/(1+x(t-1)^2) + 8*cos(1.2*t) + Q*randn;
    y(t) = x(t)^2/20 + R*randn;
    
    %% EKF
    x_predict = x_hat(t-1)/2 + 25*x_hat(t-1)/(1+x_hat(t-1)^2) + 8*cos(1.2*t);
    P_predict = J(t-1)^2*P(t-1) + Q;

    H = x_predict/10;

    K = P_predict*H/(H^2*P_predict + R);

    x_hat(t) = x_predict + K*(y(t) - x_predict^2/20);
    P(t) = (1 - K*H)*P_predict;

    J(t) = 1/2 - 25*(x_hat(t)^2-1)/(1+x_hat(t)^2)^2;

    %% UKF
    
    
end

t_vec = 1:T;
figure(1)
plot(t_vec, x, t_vec, x_hat)