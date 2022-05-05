clc
clear

%CREATING THE FSE WITH COMPLEX COEFFICIENT ALGORITHM
function complex_coeff_fourier_expansion(N,t,T,signal)
t0 = linspace(0,4,numel(t));
f = 0; %fourier expansion
for k=-N:N
    fi = signal.*exp(-1i*2*pi*(1/T)*k*t); %multiply signal to exponential function
    ck = (1/T) * trapz(t0,fi); %calculate ck
%     disp(ck)
    f = f + ck*exp(1i*2*pi*(1/T)*k*t); %assigning values to the fs
end
figure
plot(t,signal,'k',LineWidth=1) %plot original signal
axis(AXIS) %sets axis
title(TITLE) %set title
xlabel(XLABEL) %set xlabel
ylabel(YLABEL) %set yalbel
grid on %set grid
hold on
plot(t,f,'r') %plot fourier series expansion
legend(YLABEL, "FSE(N="+num2str(N)+")") %set legend
end
