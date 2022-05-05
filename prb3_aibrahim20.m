clc
clear

a = 9; % fifth digit of student number
b = 2; % sixth digit of student number
c = 1; % last digit of student number
Ns = [3,5,10,50]; %values of N

%GENERATING THE SQUARE WAVE
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

%GENERATING THE TRIANGULAR WAVE
t = 0:0.01:20; %time vector with an increment of 0.01
v = zeros(size(t)); %traingular wave vector
T = 4; %period
for k=0:4 %generate 5 period of triangular signal
    for ii=1:numel(t)
        if ((t(ii)>=k*T) && (t(ii)<=(0.5*T + k*T)))
            v(ii) = (1 - 4/T*(t(ii)-k*T));
        elseif ((t(ii)>(0.5*T + k*T)) && (t(ii)<(T + k*T)))
            v(ii) = ((4/T*(t(ii)-k*T)) - 3);
        end
    end
end

%GENERATE THE FSE OF THE SQUARE WAVE
TITLE = 'FSE of Square Wave with complex coefficient'; %title
XLABEL = 't';
YLABEL = 's(t)';
AXIS = [0,16,0,50];
t = 0:0.01:15; %time vector
t0 = linspace(-1.5,1.5,numel(t)); %limit of integration
T = b+c;
for N=Ns
     %call the FSE function
    complex_coeff_fourier_expansion(N,t,t0,T,TITLE,XLABEL,YLABEL,AXIS,s)
end

%GENERATING THE FSE OF THE SQUARE WAVE
TITLE = 'FSE of Triangular Wave with complex coefficient';
XLABEL = 't';
YLABEL = 'v(t)';
AXIS = [0,20,0,2];
t = 0:0.01:20;
t0 = linspace(0,4,numel(t));
T = 4;
for N=Ns
     %call the FSE function
    complex_coeff_fourier_expansion(N,t,t0,T,TITLE,XLABEL,YLABEL,AXIS,v)
end

%CREATING THE FSE WITH COMPLEX COEFFICIENT ALGORITHM
function complex_coeff_fourier_expansion(N,t,t0,T,TITLE,XLABEL,YLABEL,AXIS,signal)
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