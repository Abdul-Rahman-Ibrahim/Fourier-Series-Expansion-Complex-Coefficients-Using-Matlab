clc
clear all

a = 9; % amplitude of square signal
b = 2; % end point of first cycle square
c = 1; % sum of this and 'b' gives the period of the square signal

%THIS IS JUST A TEST OF THE FUNCTION DEFINED BELOW
%GENERATING THE SQUARE WAVE TO TO BE REPRESENTED WITH THE FOURIER SERIES EXPANSION(FSE)
t = 0:0.01:15; %time vector with an increment of 0.01
s = zeros(size(t)); %square wave signal
for k=0:4 %generate 5 period of square signal
    for ii=1:numel(t)
        if ((t(ii)>=k*(b+c)) && (t(ii)<=k*(b+c)+b))
            s(ii) = a;
        elseif ((t(ii)>k*(b+c)+b) && (t(ii)<(k+1)*(b+c)))
            s(ii) = 0;
        end
    end
end

N = 50; %summation end limit
T = b+c; %period of the signal
f = complex_coeff_fourier_expansion(N,t,T,s)

figure()
plot(t,s,"k")
hold on
plot(t,f)

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
