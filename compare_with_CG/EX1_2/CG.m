% solve equation by CG method
% refer to book <computational methods for inverse problems> page 32
function [x_CG,num_CG,Err_CG] = CG(noise_level)
m =100; tau = 1/m;  s = (0:tau:1)'; t = 0:tau:1;
K = min(s,t).*(1-max(s,t)); 
x_real = (-6*t.^2.*(1-t).*(2-8*t+7*t.^2))';
y_real = s.^4.*(1-s).^3; 
A = tau*K*diag([1/2,ones(1,m-1),1/2]);% Rectangle integral formula
% noise = (2*rand(size(y_real))-1);save('noise','noise');
load("noise.mat");
% noise_level = 0.001; 
y=y_real.*(1+noise*noise_level);% parameters in the equation

x = zeros(m+1,1); 
g = A'*(A*x-y);
p = -g;
dlt = norm(g,2)^2;

N=10000;Err_CG=zeros(N+1,1); Err_CG(1) = norm(x-x_real,2)/norm(x_real,2);% save error
X = zeros(m+1,N+1);  % save x(t) in X
Residual =zeros(N+1,1);num=0;
for n = 1:N
    h = A'*A*p; tau = dlt/(p'*h);
    x_old = x;
    x = x_old +tau*p; % update solution
    
    X(:,n+1) = x;
    Err_CG(n+1) = norm(x-x_real,2)/norm(x_real,2);
    g = g+tau*h; % update gradient;
    dlt_old = dlt;
    dlt = norm(g,2)^2; beta = dlt/dlt_old;
    p = -g+beta*p; % update search direction;

    Residual(n) = norm(A*x - y);
    if  Residual(n)/norm(y_real - y)<1.01
        num = num+1;
        if num ==1
            x_CG = x;
            num_CG = n;
        end
    end
    if num>0 &&  n == num_CG+50
        break;
    end
    
    figure(1); 
    subplot(2,1,1),plot(t',X(:,n+1),t,x_real,'r',LineWidth=3)
    % 创建 ylabel
    ylabel({'x'});
    % 创建 xlabel
    xlabel({'n'});
    % 创建 title
    title({'CG'});
    subplot(2,1,2),plot(max(1,n-100):n+1,Err_CG(max(1,n-100):n+1))
    % 创建 title
    title({'L^2 error'});
end
