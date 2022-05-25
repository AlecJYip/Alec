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

%theta_rim_max = 2*atan(1/(4*F/D));

%theta_rim_max = 2*pi/180; %rad

theta_rim_max = pi/2;

P = zeros(1,length(0:theta_rim_max/100:theta_rim_max));

theta_rim = 0;
k = 1;
%% Analysis
while(theta_rim<=theta_rim_max)
    theta = 0;
    range_theta = (0:0.001:theta_rim);
    range_phi = (0:0.001:2*pi);
    arr = zeros(1,length(range_theta).*length(range_phi));
    bruh = zeros(1,length(range_theta).*length(range_phi));
    f = zeros(1,length(arr));
    n = 1;
    while(theta<=theta_rim)        
        phi =0;
        while(phi<=2*pi)
            term1 = (cos(theta)^q)/sqrt((cos(theta)^2)*(sin(phi)^2)+(cos(phi)^2));
            term2 = sqrt(...
                (cos(theta)*sin(phi)*cos(theta/2))^2+...
                (cos(phi)*cos(theta/2))^2+...
                (cos(theta)*sin(phi)*sin(theta/2))^2);
            p = 2*F*tan(theta/2);
            term3 = exp(1i*beta*p*sin(theta)*cos(phi));
            func = term1*term2*term3;
            disp(abs(func))
            arr(n) = 4*pi*sin(theta)*(0.001*0.001)*abs(func)^2;
            phi = phi + 0.001;
            n = n + 1;
        end
        theta = theta + 0.001;
    end
    P(k)=sum(arr);
    k = k + 1;
    theta_rim = theta_rim + theta_rim_max/100;
end
angle_rim = (0:theta_rim_max/100:theta_rim_max);
D_dBi = 10*log10(4.*pi./P);
figure;plot(angle_rim*180/pi,D_dBi);