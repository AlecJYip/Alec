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

%theta_rim_max = 2*pi/180;

theta_rim_max = 2*atan(1/(4*F/D));

%theta_rim_max = 90*pi/180;

rim_step_size = lambda/1000;

numerator = zeros(1,length(0:rim_step_size:theta_rim_max));

theta_rim = 0;

k = 1;

step_size = 0.3*lambda;

range_phi = (0:step_size:2*pi);

s_i_hat = [1 0 0];

%% Analysis: Numerator

while(theta_rim<=theta_rim_max)
    theta = theta_rim;
    n = 1;    
    arr = [0 0 0];
    while(theta<=theta_rim_max)        
        
        phi = 0;
        s_i = F*(sec(theta/2)^2);
        far_field_dist = s_i*cos(theta);
        
        while(phi<=2*pi)
            
            y_hat = [sin(theta)*cos(phi) cos(theta)*cos(phi) -sin(theta)];
                       
            H_over_I  = (1/s_i)*(cross(y_hat, s_i_hat)/magnitude(cross(y_hat, s_i_hat)))...
                *exp(-1i*beta*s_i)*(cos(theta)^q);
                        
            n_hat = [-cos(theta/2) sin(theta/2) 0];
            
            J = 2*cross(n_hat, H_over_I);
            
            func = J*exp(1i*beta*far_field_dist);
            
            arr = arr + (step_size^2)*func*sin(theta)*s_i^2;
                        
            phi = phi + step_size;
                        
            n = n + 1;
        end
        theta = theta + step_size;
    end
    numerator(k)= abs(magnitude(arr))^2;
    k = k + 1;
    theta_rim = theta_rim + rim_step_size;
end

%% Analysis: Denominator

theta = 0;

arr2 = 0;

q = 0.9;

while(theta <= pi/2)
    
    arr2 = arr2 + (cos(theta)^(2*q))*sin(theta)*step_size;    
    
    theta = theta + step_size;
    
end

denominator_term_1 = 2*pi*arr2;

denominator_term_2 = 1/((4*pi)^3);

denominator = denominator_term_1*denominator_term_2;



%% Plotting
angle_rim = (0:rim_step_size:theta_rim_max);
D_dBi = 10*log10(numerator./denominator);
figure;plot(angle_rim*180/pi,D_dBi);hold all;xlabel('\theta [deg]');xlim([0 2]);
ylabel('Directivity [dBi]');title('Figure 3');grid on;

%% Functions
function J = magnitude(Q) 
    n = 1;
    sum_1 = zeros(1, length(Q));
    while(n<=length(Q))        
        sum_1(n) = abs(Q(n))^2;
        n = n + 1;        
    end    
    J = sqrt(sum(sum_1));
end
