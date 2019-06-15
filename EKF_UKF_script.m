clear all;
clc;

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

%% UKF Params
x_dim = 1;
num_sig_pts = 2*x_dim+1;
sig_pts = zeros(num_sig_pts, 2);
lambda = 3 - x_dim;
%lambda = 0;
y_predict = zeros(num_sig_pts,1);
x_hat_UKF = zeros(1,T);
x_cov_hat = zeros(1,T);
k
sig_pts(1,1) = x(1);
sig_pts(1,2) = lambda/(x_dim+lambda);
for i = 2:x_dim+1
    sig_pts(i,1) = x(1) + (sqrt((x_dim + lambda)*Q));
    sig_pts(i,2) = 1/(2*(x_dim+lambda));
end
for i = x_dim+2:2*x_dim+1
    sig_pts(i,1) = x(1) - (sqrt((x_dim + lambda)*Q));
    sig_pts(i,2) = 1/(2*(x_dim+lambda));
end


for t = 2:T
    %% Actual System
    x(t) = x(t-1)/2 + 25*x(t-1)/(1+x(t-1)^2) + 8*cos(1.2*t) + Q*randn;
    y(t) = x(t)^2/20 + R*randn;
    
    %% EKF
%     x_predict = x_hat(t-1)/2 + 25*x_hat(t-1)/(1+x_hat(t-1)^2) + 8*cos(1.2*t);
%     P_predict = J(t-1)^2*P(t-1) + Q;
% 
%     H = x_predict/10;
% 
%     K = P_predict*H/(H^2*P_predict + R);
% 
%     x_hat(t) = x_predict + K*(y(t) - x_predict^2/20);
%     P(t) = (1 - K*H)*P_predict;
% 
%     J(t) = 1/2 - 25*(x_hat(t)^2-1)/(1+x_hat(t)^2)^2;

    %% UKF
    for i = 1:x_dim+1
       sig_pts(i,1) = sig_pts(i,1)/2 + 25*sig_pts(i,1)/(1+sig_pts(i,1)^2) + 8*cos(1.2*t);
    end
    x_mean_predict = dot(sig_pts(:,1),sig_pts(:,2));
    x_cov_predict = 0;
    y_hat = 0;
    S = 0;
    cross_corr = 0;
    for i = 1:2*x_dim+1
        x_cov_predict = x_cov_predict + sig_pts(i,2)*(sig_pts(i,1) - x_mean_predict)^2 + Q;
        y_predict(i) = sig_pts(i,1)^2/20;
        y_hat = y_hat + sig_pts(i,2)*y_predict(i);
    end
    for i=1:2*x_dim+1
       S = S + sig_pts(i,2)*(y_predict(i) - y_hat)^2 + R;
       cross_corr = cross_corr + sig_pts(i,2)*(sig_pts(i,1)-x_mean_predict)*(y_predict(i)-y_hat);
    end
    K = cross_corr/S;
    x_hat_UKF(t) = x_mean_predict + K*(y(t)-y_hat);
    x_cov_hat(t) = (1 - K*cross_corr)*x_cov_predict;
    
    sig_pts(1,1) = x_hat_UKF(t);
    sig_pts(1,2) = lambda/(x_dim+lambda);
    for i = 2:x_dim+1
        sig_pts(i,1) = x_hat_UKF(t) + (sqrt((x_dim + lambda)*x_cov_hat(t)));
        sig_pts(i,2) = 1/(2*(x_dim+lambda));
    end
    for i = x_dim+2:2*x_dim+1
        sig_pts(i,1) = x_hat_UKF(t) - (sqrt((x_dim + lambda)*x_cov_hat(t)));
        sig_pts(i,2) = 1/(2*(x_dim+lambda));
    end
    
end

t_vec = 1:T;
figure(1)
plot(t_vec, x, t_vec, x_hat_UKF)