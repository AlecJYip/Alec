clc;
clear;
close all;

%% Variables

f = 1.5e9; %Hz

c = 3e8; %m/s

lambda = c/f; %m

beta = 2*pi/lambda; %rad/m

q = 1.14;

Si = 8.5; %m

D = 18; %m

F = 0.4*D;

e = 8.85e-12;

mu = (4*pi)*(10^-7);

eta = sqrt(mu/e);

w = 2*pi*f;

%theta_rim_max = pi/2;

%theta_rim_max = pi/4;

theta_rim_max = 2*pi/180;

E_i = zeros(1,length(0:theta_rim_max/200:theta_rim_max));

E_f = zeros(1,length(0:theta_rim_max/200:theta_rim_max));

theta_rim = 0;

k = 1;

step_size = lambda/500;

range_phi = (0:step_size:2*pi);

%% Analysis
while(theta_rim<=theta_rim_max)
    theta = 0;
    n = 1;    
    range_theta = (0:step_size:theta_rim);    
    arr = 0;
    while(theta<=theta_rim)        
        phi = 0;
        p = 2*F*tan(theta/2);
        rf = F*(sec(theta/2)^2);
        zf = sqrt(p^2-rf^2);
        while(phi<=2*pi)
            
            y_hat = [sin(theta)*cos(phi) cos(theta)*cos(phi) -sin(theta)];
            
            s_i_hat = [1 0 0];
            
            H_over_I  = (cross(y_hat, s_i_hat)/magnitude(cross(y_hat, s_i_hat)))*exp(-1i*beta*rf)*(cos(theta)^q);
            
            n_hat = [-cos(theta/2) sin(theta/2) 0];
            
            J = 2*cross(n_hat, H_over_I);
            
            mag_J = magnitude(J);
                       
            func = w*mu*(1/(4*pi))*abs(mag_J*exp(1i*beta*tan(pi/2-theta)));
            
            arr = arr + (step_size^2)*func;
            phi = phi + step_size;
            
            E_f(n) = magnitude(abs(-eta*H_over_I));
            
            n = n + 1;
        end
        theta = theta + step_size;
    end
    E_i(k)= arr;
    k = k + 1;
    theta_rim = theta_rim + theta_rim_max/200;
end

magnitude_E =  eta*magnitude((cross(y_hat, s_i_hat)/magnitude(cross(y_hat, s_i_hat))));


%% Plotting
angle_rim = (0:theta_rim_max/200:theta_rim_max);
D_dBi = 20*log10(E_i/magnitude_E);
figure;plot(angle_rim*180/pi,D_dBi);hold all;xlabel('\theta [deg]');
ylabel('Directivity [dBi]');title('Figure 3');grid on;

%% Functions

function J = magnitude(Q) 

    n = 1;
    sum_1 = zeros(1, length(Q));
    while(n<=length(Q))
        
        sum_1(n) = abs(Q(n))^2;
        n = n + 1;
        
    end
    
    J = sum(sum_1);
end
