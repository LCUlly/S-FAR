% 
m = 100; tau = 1/m; s = (0:tau:1)'; 
y_real = s.^4.*(1-s).^3; 
noise = (2*rand(size(y_real))-1);save('noise','noise');