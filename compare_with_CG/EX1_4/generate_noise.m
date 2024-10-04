thet = 1.99;  % thet in (0,1)
n1 = 50; sau1 = 1/n1; m1 = 50; tau1 = 1/m1; % m1=n1
s1 = (0:sau1:1)'; t1 = 0:tau1:1;
n2 = 50; sau2 = 1/n2;m2 = 50; tau2 = 1/m2; % n2 = m2;
s2 = (0:sau2:1)'; t2 = 0:tau2:1;
K = zeros((n1+1)*(n2+1),(m1+1)*(m2+1));
for i = 1:n1+1
    for j = 1:n2+1
        KK =  sin(pi*(s1(i)-t1').^2).*sin(pi*(s2(j)-t2).^2);
        KK = diag([1/2,ones(1,m1-1),1/2])*KK*diag([1/2,ones(1,m2-1),1/2]);
       K(i+(j-1)*(n2+1),:) = reshape(KK,1,(m1+1)*(m2+1)); 
    end
end
A = tau1*tau2*K;% Rectangle integral formula
xx_real = 1+ (t1').^2.*t2.^2;%t1'.^3+t2.^3+(t1-1)'.*t2;
x_real =  reshape(xx_real,(m1+1)*(m2+1),1);
y_real = A*x_real; y_real = reshape(y_real,(n1+1)*(n2+1),1);

noise = (2*rand(size(y_real))-1); save('noise','noise');