%% Part 4.1- Validation of Numerical Solution

clear all;
close all;
global I_o gamma w;

% Set Initial Conditions
t(1)=0;
theta(1)=pi;
gamma=0;
I_o=1.0005;
w=1;
n=1;
dt=0.01;
dt_min=10^-6;
tolmax=1e-2;
tolmin=tolmax/100;
w_o=sqrt((I_o)^2-1);
t_final=0.995*((2*pi)/w_o);

%Euler and Improved Euler Variable Time-Step
tic
while t(n)<t_final
    if t(n)+dt>t_final
        dt=t_final-t(n);
    end
    k1=rhsLab4(theta(n),t(n));
    euler=theta(n)+dt*k1;
    k2=rhsLab4(euler,t(n));
    imp_euler=theta(n)+dt*((k1+k2)/2);
    diff=imp_euler-euler;
    if abs(diff)>tolmax
        dt=(dt/2);
        if dt<dt_min
            error('dt is too small');
        end
    else
        n=n+1;
        t(n)=t(n-1)+dt;
        theta(n)=imp_euler;
        if abs(imp_euler-euler)<tolmin
            dt=2*dt;
        end
    end
end

% Calculate Numerical Solution
for ii=1:length(theta)
    thetadot(ii)=I_o-cos(theta(ii));
end
toc

% Calculate Analytic Solution
for ii=1:length(t)
    v(ii)=((I_o)^2-1)/(I_o+cos(w_o*t(ii)+theta(1)));
end

%Graphing Analytic solution
figure;
subplot(2,1,1);
plot(t,v,'.');
xlabel('time');
ylabel('voltage');
title('Voltage Across the Josephson Device (Variable Time Step)');
hold on;

%Graphing Numeric Solution
plot(t,thetadot);
legend('Analytical','Numerical');

% Graphing dt vs. time for Numeric Solution
dt_vec=t(2:end)-t(1:end-1);
dt_vec=[dt_vec dt_vec(end)];
subplot(2,1,2);
plot(t,dt_vec);
title('dt vs. time');
xlabel('time');
ylabel('dt');

% Error Plot for Numerical Solution
figure
error3=v-thetadot;
plot(t,abs((v-thetadot)));
title('Error Plot');
xlabel('time');
ylabel('error');

% Calculating Error
error=sum(abs(v-thetadot).*dt);
error2=abs(v(end)-thetadot(end));
fprintf('integral error = %5.4f\n',error);
fprintf('end error = %5.4f\n', error2);
fprintf('tolmax = %10.10f\n', tolmax);
fprintf('tolmin = %10.10f\n', tolmin);

% Number of Steps
variable_n=n;
fprintf('The number of steps taken is: %4.2d\n',variable_n);


%% Part 4.2- Comparing Fixed and Adaptive Time Step
clear all;
global I_o gamma w;

% Initial Conditions
t(1)=0;
theta(1)=pi;
gamma=0;
I_o=1.0005;
w=1;
n=1;
dt=0.04;
tolmax=1e-2;
tolmin=tolmax/100;
w_o=sqrt((I_o)^2-1);
t_final=0.995*((2*pi)/w_o);

% Adaptive Time Step and Improved Euler Method
tic
while t<t_final
    if t+dt>t_final
        dt=t_final-t;
    end
    k1=rhsLab4(theta(n),t(n));
    euler=theta(n)+dt*k1;
    k2=rhsLab4(euler,t(n));
    imp_euler=theta(n)+dt*((k1+k2)/2);
    n=n+1;
    t(n)=t(n-1)+dt;
    theta(n)=imp_euler;
end

% Calculating Numerical Solution
for ii=1:length(theta)
    thetadot(ii)=I_o-cos(theta(ii));
end
toc

% Calculating Analytic Solution
for ii=1:length(t)
    v(ii)=((I_o)^2-1)/(I_o+cos(w_o*t(ii)+theta(1)));
end

% Graphing Analytic Solution
figure;
subplot(2,1,1);
plot(t,v,'.');
xlabel('time');
ylabel('voltage');
title('Voltage Across the Josephson Device (Constant Time Step)');
hold on;

% Graphing Numeric Solution
plot(t,thetadot);
legend('Analytical','Numerical');

% Graphing dt vs. Time
dt_vec=dt*ones(1,length(t));
subplot(2,1,2);
plot(t,dt_vec);
title('dt vs. time');
xlabel('time');
ylabel('dt');

% Calculating Error for Numerical Solution
error=sum(abs(v-thetadot).*dt);
error2=abs(v(end)-thetadot(end));
fprintf('integral error = %5.4f\n',error);
fprintf('end error = %5.4f\n', error2);
fprintf('tolmax = %10.10f\n',tolmax);
fprintf('tolmin = %10.10f\n', tolmin);

% Number of Steps
constant_n=n;
fprintf('The number of steps taken is: %4.2d\n',constant_n);

%% Part 4.2 - Varying tolmax and tolmin
clear all;
close all;

tolmin=[1e-9, 1e-8, 1e-7 1e-6, 1e-5, 1e-4,];
tolmax=[1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3];

for aa=1:length(tolmin)
    for mm=1:length(tolmax)
        %Euler and Improved Euler Variable Time-Step
        clear t theta n dt thetadot t_analytic v dt_vec error
        global I_o gamma w;
        t(1)=0;
        theta(1)=pi;
        gamma=0;
        I_o=1.0005;
        w=1;
        n=1;
        dt=0.01;
        dt_min=10^-6;
        w_o=sqrt((I_o)^2-1);
        t_final=0.995*((2*pi)/w_o);
        tic
        while t<t_final
            if t+dt>t_final
               dt=t_final-t;
            end
        k1=rhsLab4(theta(n),t(n));
        euler=theta(n)+dt*k1;
        k2=rhsLab4(euler,t(n));
        imp_euler=theta(n)+dt*((k1+k2)/2);
        diff=imp_euler-euler;
        if abs(diff)>tolmax(mm)
            dt=(dt/2);
            if dt<dt_min
                error('dt is too small');
            end
        else
            n=n+1;
            t(n)=t(n-1)+dt;
            theta(n)=imp_euler;
            if abs(imp_euler-euler)<tolmin(aa)
                dt=2*dt;
            end
        end
        
    end
    
    % Calculating Numerical Solution
    for ii=1:length(theta)
        thetadot(ii)=I_o-cos(theta(ii));
    end
    toc

    % Calculating Analytic Solution
    for ii=1:length(t)
        v(ii)=((I_o)^2-1)/(I_o+cos(w_o*t(ii)+theta(1)));
    end

    % Graphing Analytic Solution
    figure;
    subplot(2,1,1);
    plot(t,v,'.');
    xlabel('time');
    ylabel('voltage');
    str=sprintf('Voltage Across Josephson Device, tolmin=%4.4f tolmax=%4.4f',tolmin(aa), tolmax(mm));
    title(str);
    hold on;
    
    % Graphing Numerical Solution
    plot(t,thetadot);
    legend('Analytical','Numerical');

    % Graphing dt vs. Time
    dt_vec=t(2:end)-t(1:end-1);
    dt_vec=[dt_vec dt_vec(end)];
    subplot(2,1,2);
    plot(t,dt_vec);
    title('dt vs. time');
    xlabel('time');
    ylabel('dt');
    hold on

    %Error
    error=sum(abs(v-thetadot).*dt);
    disp(error);

    end
end

%% Part 5- Phase Locking

clear all;
close all;

global I_o gamma w;
t(1)=0;
theta(1)=pi;
I_o=1.01;
gamma=0.02;
w_o=sqrt((I_o)^2-1);

% When w = 0.62 and 1.1
% fortitle = [0.62 1.1];
% w_vec = [(0.62*w_o) (1.1*w_o)];

% When w is being looped
w_vec = (0.5*w_o):(0.02*w_o):(2*w_o);

n=1;
dt=0.01;
dt_min=10e-6;
tolmax=10e-5;
tolmin=tolmax/100;
t_final=2000;

%Euler and Improved Euler Variable Time-Step
for ii = 1:length(w_vec)
    n = 1;
    w = w_vec(ii);
    while t(n)<t_final
        if t(n)+dt>t_final
            dt=t_final-t(n);
        end
        k1=rhsLab4(theta(n),t(n));
        euler=theta(n)+dt*k1;
        k2=rhsLab4(euler,t(n));
        imp_euler=theta(n)+dt*((k1+k2)/2);
        diff=imp_euler-euler;
        if abs(diff)>tolmax
            dt=(dt/2);
            if dt<dt_min
                error('dt is too small');
            end
        else
            n=n+1;
            t(n)=t(n-1)+dt;
            theta(n)=imp_euler;
            if abs(imp_euler-euler)<tolmin
                dt=2*dt;
            end
        end
    end
        w_mean(ii) = (1/t_final)*(theta(end) - theta(1));
        wratio(ii) = w_mean(ii)/w;
    for jj=1:length(theta)
        thetadot(jj)=I_o-cos(theta(jj))+gamma*sin(w*t(jj));
        current(jj) = gamma*sin(w*t(jj));
    end
%  %  Graphing voltage and current vs. time for w = 0.62 and 1.1
%     figure;
%     plot(t,thetadot)
%     hold on;
%     plot(t,current);
%     str = ['Voltage and Input Current for w =', num2str(fortitle(ii)),'w_o'];
%     title(str);
%     xlabel('Time');
%     ylabel('Voltage');
end

% % Graphing wmean/w vs. w for looped w
figure;
plot(w_vec, wratio,'-o');
xlabel('w');
ylabel('w mean/w');
title('Phase-Locking');

