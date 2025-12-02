% Generalized Constantin-Lax-Majda Equation solver
%#ok<*GVMIS>
echo off, clc;
hold off;
clear all;
close all;

% Global Constants
global N; global DT; global H; global T; global F; global X; global K; global M; global MAX; global COUNT;
N = 2^14;                  % number of grid points
DT = 0.0001;               % timestep
H = 2*pi/N;               % distance between x_i and x_i+1
T = 7;                    % timescale
F = 1E-12;                % filter tolerance
X = linspace(-pi,pi-H,N); % grid points
K = [0:N/2-1, -N/2:-1];   % wavenumbers
MAX = 1E6;

% Initial Condition
eps = 0.1;
w_0 = sin(X)+eps*sin(2*X);
SolveGCLM(w_0);
%deriv_test()
%ht_test()
%calc_u_test()

function [E] = calc_u_test()
    global X; global N; global F;
    w = sin(3*X);
    u_ana = -sin(3*X)/3;
    w_c = fft(w,N);
    w_c = filter(w_c, F);
    u_num = calc_u(w_c);

    dif = abs(u_ana - u_num);
    %E = max(dif);
    E = norm(dif,2);
end

function [E] = ht_test()
    global X; global N; global F;
    w = sin(3*X);
    wht_ana = -cos(3*X);
    w_c = fft(w,N);
    w_c = filter(w_c, F);
    wht_num = ifft(ht(w_c));

    dif = abs(wht_ana - wht_num);
    %E = max(dif);
    E = norm(dif,2);
end

function [E] = deriv_test()
    global X; global N; global F;
    w = sin(3*X);
    w_c = fft(w, N);
    w_c = filter(w_c, F);
    wx_ana = 3*cos(3*X);
    wx_num_c = deriv(w_c,1);
    wx_num_c = dealias(wx_num_c);
    wx_num = ifft(wx_num_c,N);

    dif = abs(wx_ana - wx_num);
    %E = max(dif);
    E = norm(dif,2);
end


% Generalized Constantin-Lax-Majda Equation solver
function SolveGCLM(w_0)
    %a = -0.8;
    for a=-1:0.1:1
        disp("a: "+a);
        [~,~] = RK4(w_0, a); 
    end
end

% Blow-up Checker L2 norm of w_x
function b = blowup_1(w_t)
    global MAX; global N;
    n = norm(deriv(ifft(w_t,N),1),2);
    b = (n>MAX);
end

% Blow-up Checker LINF norm of u_x
function b = blowup_1(w_t)
    global MAX; global N;
    n = norm(deriv(ifft(w_t,N),1),2);
    b = (n>MAX);
end

% Derivative in Fourier Space
function d = deriv(w_c, ord)
    global K;
    % multiply by (ik)^ord
    d = ((1i.*K).^ord).*w_c;
end


% Hilbert Transform Physical
function h = ht(w_c)
    global K;
    % multiply Fourier coefficients by -i signum(k)
    h = (-1i*sign(K)).*w_c;
end

% Calculate velocity from vorticity
function u = calc_u(w_c)
    global N; global K;
    % exclude k=0 
    points = (K ~= 0);
    % calc k=0 val
    w_c(~points) = 4*log(2)*w_c(~points);
    % divide Fourier coefficients by k~=0
    w_c(points) = -w_c(points)./abs(K(points));
    % convert back to physical space
    u = ifft(w_c,N);
end

% Dealias with 2/3 rule
function w_c = dealias(w_c)
    global N; global K;
    % points above 2/3
    dealias_points = (abs(K) >= ceil(2*(N/2)/3));
    % set dealias points to zero
    w_c(dealias_points) = 0;
end

% filter |e_k| < 1E-12 to zero
function w_c = filter(w_c, threshold)
    % find points below threshold
    filter_points = (abs(w_c) <= threshold);
    % set filter points to zero
    w_c(filter_points) = 0;
end

function [f,n] = calc_RHS(w,a,t)
    global N; global F; global X; global T; global K;
    w_c = filter(fft(w,N),F); % Convert w to Fourier space
    ux_c = ht(w_c); % calculate du/dx
    wx_c = deriv(w_c,1); % calculate dw/dx
    u = calc_u(w_c); % calculate u
    w = ifft(w_c,N); %calculate w

    ux = ifft(ux_c,N); % convert du/dx to physical space
    wx = ifft(wx_c,N); % convert dw/dx to physical space

    t1 = ux.*w; % compute du/dx * w
    t2 = a*(u.*wx); % compute a*u*dw/dx

    t1 = real(ifft(filter(dealias(fft(t1,N)),F),N)); % dealias t1
    t2 = real(ifft(filter(dealias(fft(t2,N)),F),N)); % dealias t2

    f = t1 - t2;

    
    n = NaN;
    if(t ~= -1)
        mod_step = 0.5;
        if a>=0.5
            mod_step = 1.0;
        end
        if (mod(t,mod_step)==0 && t ~=T)
            disp("t: "+t);
            n = norm(wx_c,2);
            colors = {[0,0,1],[0.5,0.5,0.5]};
            figure(1);

            plot(X,real(w),'color', colors{1 + ~(t==0)});
            hold on;
        
            figure(2);
            powerSpec = filter(abs(w_c).^2,F);
            loglog(K(K>0),real(abs(powerSpec(K>0))),'color',colors{1 + ~(t==0)});
            hold on;
        
            figure(4);
            plot(X,real(ux),'color',colors{1 + ~(t==0)});
            hold on;
        
            figure(5);
            plot(X,real(wx),'color',colors{1 + ~(t==0)});
            hold on;
        end
    end
end

% Integration of the model problem using the Euler Explicit scheme;
function [w_T, b] = Euler_Explicit( w_0, a)
    global N; global F; global DT; global X; global T;

    i = 1;
    t = 0;

    f = calc_RHS(w_0,a);
    w_t = w_0;
    b = false;
    while ( t <= T && ~b)
      t = (i) * DT;
     
      w_t = w_t + f*DT;
      
      f = calc_RHS(w_t, a);

      b = (b | blowup(w_t));
      if b
        disp("BLOWUP");
      end
    
      if (mod(i-1,40000)==0)
        plot(X,w_t);
        hold on;
      end

      i = i + 1;
    end

    w_T = w_t;
end


% Integration of the model problem using the RK4 scheme;
function [w_T, b] = RK4( w_t, a)
    global DT; global X; global T; global N; global F; global K; global COUNT;
    i = 0;
    t = 0;
    b = false;
    COUNT = 1;

    norm_data = [];
    norm_times = [];
    
    for j=1:5
        figure(j); clf;
    end

    [~,n] = calc_RHS(w_t, a, t);
    norm_data(COUNT) = n;
    norm_times(COUNT) = t;
    COUNT = COUNT + 1;
    
    while ( t <= T && ~b)
      t = (i+1) * DT;
     
      [f,n] = calc_RHS(w_t, a, t);
      k1 = DT*f;

      [f,~] = calc_RHS(w_t + k1/2, a, -1);
      k2 = DT*f;
      
      [f,~] = calc_RHS(w_t + k2/2, a, -1);
      k3 = DT*f;
      
      [f,~] = calc_RHS(w_t + k3, a, -1);
      k4 = DT*f;
      
      w_t = w_t + (1/6)*( k1 + 2*k2 + 2*k3 + k4 );
      
      b = (b | blowup(w_t));

      if(~isnan(n))
        norm_data(COUNT) = n;
        norm_times(COUNT) = t;
        COUNT = COUNT + 1;
      end
      i = i + 1;
    end

    if b
        disp("Blowup");
        disp(t)
    end
    w_T = w_t;

    fg = figure(1);
    figure(fg);
    if(b)
        plot(X,real(w_t),"--",'color','r');
    else
        plot(X,real(w_t),"-",'color','r');
    end
    title("$\omega$ vs $x$: $a=$" + num2str(round(a,1)), 'Interpreter', 'latex');
    ylabel("$\omega$",'Interpreter','latex');
    xlabel("$x$",'Interpreter','latex');
    %ylim([-3,3]);
    grid on;
    filename = sprintf('../Results/w_vs_x_a%s.jpg', num2str(round(a,1)));
    exportgraphics(fg,filename,'Resolution',600)
    filename = sprintf('../Results/w_vs_x_a%s.fig', num2str(round(a,1)));
    savefig(fg,filename)

    fg = figure(2);
    figure(fg);
    w_c = filter(fft(w_t,N),F);
    powerSpec = filter(abs(w_c).^2,F);
    if(b)
        loglog(K,real(abs(powerSpec)),"--",'color','r');
    else
        loglog(K,real(abs(powerSpec)),"-",'color','r');
    end
    title("$\omega$ Fourier Spectrum: $a=$" + num2str(round(a,1)), 'Interpreter', 'latex');
    ylabel("$|\omega_k|^2$",'Interpreter','latex');
    xlabel("$k$",'Interpreter','latex');
    grid on;
    filename = sprintf('../Results/w_spectrum_a%s.jpg', num2str(round(a,1)));
    exportgraphics(fg,filename,'Resolution',600)
    filename = sprintf('../Results/w_spectrum_a%s.fig', num2str(round(a,1)));
    savefig(fg,filename)

    fg = figure(3);
    figure(fg);
    wx_c = deriv(w_c,1);
    temp = norm(wx_c,2);
    norm_data(COUNT) = temp;
    norm_times(COUNT) = t;
    val = 0;
    semilogy(norm_times(1:end-val),norm_data(1:end-val));
    title("$\|\omega_x\|$ vs Time: $a=$" + num2str(round(a,1)), 'Interpreter', 'latex');
    ylabel("$\|\omega_x\|$",'Interpreter','latex');
    xlabel("$t$",'Interpreter','latex');
    grid on;
    filename = sprintf('../Results/norm_vs_time_a%s.jpg', num2str(round(a,1)));
    exportgraphics(fg,filename,'Resolution',600)
    filename = sprintf('../Results/norm_vs_time_a%s.fig', num2str(round(a,1)));
    savefig(fg,filename)

    fg = figure(4);
    figure(fg);
    ux_c = ht(w_c);
    ux = ifft(ux_c,N);
    if(b)
        plot(X,real(ux),"--",'color','r');
    else
        plot(X,real(ux),"-",'color','r');
    end
    title("$u_x$ vs $x$: $a=$" + num2str(round(a,1)),'Interpreter', 'latex');
    ylabel("$u_x$",'Interpreter','latex');
    xlabel("x",'Interpreter','latex');
    %ylim([-2,5]);
    grid on;
    filename = sprintf('../Results/ux_vs_x_a%s.jpg',num2str(round(a,1)));
    exportgraphics(fg,filename,'Resolution',600)
    filename = sprintf('../Results/ux_vs_x_a%s.fig',num2str(round(a,1)));
    savefig(fg,filename)

    fg = figure(5);
    figure(fg);
    wx = ifft(wx_c,N);
    if(b)
        plot(X,real(wx),"--",'color','r');
    else
        plot(X,real(wx),"-",'color','r');
    end
    title("$\omega_x$ vs $x$: $a=$" + num2str(round(a,1)),'Interpreter', 'latex');
    ylabel("$\omega_x$",'Interpreter', 'latex');
    xlabel("$x$",'Interpreter', 'latex');
    %ylim([-5,5]);
    grid on;
    filename = sprintf('../Results/wx_vs_x_a%s.jpg', num2str(round(a,1)));
    exportgraphics(fg,filename,'Resolution',600)
    filename = sprintf('../Results/wx_vs_x_a%s.fig', num2str(round(a,1)));
    savefig(fg,filename)

end