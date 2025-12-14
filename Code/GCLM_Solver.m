% Generalized Constantin-Lax-Majda Equation solver
%#ok<*GVMIS>
echo off, clc;
hold off;
clear all;
close all;

% Global Constants
global N; global DT; global H; global T; global F; global X; global K2; global K; global DEALIAS_MAX;
global H2; global X2;
N = 2^14;                 % number of grid points
DT = 1E-5;                % timestep
H = 2*pi/N;               % distance between x_i and x_i+1
H2 = 2*pi/(3*N/2);
T = 7;                    % timescale
F = 1E-12;                % filter tolerance
X = linspace(-pi,pi-H,N); % grid points
X2 = linspace(-pi,pi-H2,3*N/2); % grid points
K = [0:N/2-1, -N/2:-1];   % wavenumbers
K2 =[0:(3*(N/2)/2)-1, -3*(N/2)/2:-1];
DEALIAS_MAX = 1E-10;

%SolveGCLM("Init_1");
%generateDGDTData("Init_1");
generateGraphs("Init_1");
%SolveGCLM("Init_2");
%generateDGDTData("Init_2");
generateGraphs("Init_2");


%Test functions used to verify the velocity, hilbert transform, and
%derivative were calculated correctly
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

function generateDGDTData(init)
    global X; global K; global DT;
    data_path = "../Results/Data/"+init+"/";
    
    index = 1;
    for a=1.0:-0.1:-1.0
        disp(a);
        filename = sprintf("%sa%d/l2_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data = readmatrix(filename);
        dgdt = (data(1,2:end)-data(1,1:end-1))/DT;
        disp(size(data(1,1:end-1)));
        disp(size(dgdt));
        filename = sprintf("%sa%d/l2_y_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        writematrix(dgdt, filename, 'WriteMode', 'append');
        filename = sprintf("%sa%d/l2_x_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        writematrix(data(1,1:end-1), filename, 'WriteMode', 'append');

        filename = sprintf("%sa%d/linf_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data = readmatrix(filename);
        int_vals = zeros(size(data));
        for i=2:length(data)
            int_vals(i) = trapz((0:DT:(i-1)*DT),data(1:i));
        end
        disp(size(int_vals));
        disp(size(data));
        filename = sprintf("%sa%d/linf_y_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        writematrix(data, filename, 'WriteMode', 'append');
        filename = sprintf("%sa%d/linf_x_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        writematrix(int_vals, filename, 'WriteMode', 'append');

        index = index + 1;
    end
end

% Generate graphs from data files
function generateGraphs(init)
    global X; global K; global DT;
    data_path = "../Results/Data/"+init+"/";
    image_path = "../Results/Images/"+init+"/";

    index = 1;
    alphas_l2 = zeros(12,21);
    alphas_linf = zeros(12,21);
    for i=1:12
        alphas_l2(i,1) = NaN;
        alphas_linf(i,1) = NaN;
    end
    as = (1.0:-0.1:-1.0);

    filename = sprintf("%s/blowup_times.csv",data_path);
    data = readmatrix(filename);
    fg = figure();
    clf;
    plot(as,data);
    ylim([0,8]);
    xlim([-1.1,1.1]);
    title('Under-resolved Times vs a','interpreter','latex','FontSize', 18);
    xlabel("$a$",'interpreter','latex','FontSize', 14);
    ylabel("Time",'interpreter','latex','FontSize', 14);
    grid on;
    filename = sprintf("%s/blowup_times.jpg",image_path);
    exportgraphics(fg,filename,'Resolution',600)

    for a=1.0:-0.1:-1.0
        disp(a);
        %colors = {[0,0,1],[0.5,0.5,0.5]};
        %plot(X,real(w),'color', colors{1 + ~(t==0)});

        % omega vs x plots
        filename = sprintf("%sa%d/w_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data = readmatrix(filename);
        fg = figure();
        clf;
        plot(X,data(1,:),'Color','blue');
        hold on;
        for i=5:4:(size(data, 1))
            plot(X,data(i,:),'Color',[0.5,0.5,0.5]);
        end
        plot(X,data(end,:),'Color','red');
        ylim([-3,3]);
        title("$\omega$ vs $x$ - a="+num2str(round(a,1)),"Interpreter","latex",'FontSize', 18);
        xlabel("$x$","Interpreter","latex",'FontSize', 14);
        ylabel("$\omega$","Interpreter","latex",'FontSize', 14);
        grid on;
        filename = sprintf("%sa%d/w_data_a%d_%s.jpg",image_path,index,index,num2str(round(a,1)));
        exportgraphics(fg,filename,'Resolution',600)

        % omega_x vs x plots
        filename = sprintf("%sa%d/wx_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data = readmatrix(filename);
        fg = figure();
        clf;
        plot(X,data(1,:),'Color','blue');
        hold on;
        for i=5:4:(size(data, 1))
            plot(X,data(i,:),'Color',[0.5,0.5,0.5]);
        end
        plot(X,data(end,:),'Color','red');
        ylim([-3,3]);
        title("$\omega_x$ vs $x$ - a="+num2str(round(a,1)),"Interpreter","latex",'FontSize', 18);
        xlabel("$x$","Interpreter","latex",'FontSize', 14);
        ylabel("$\omega_x$","Interpreter","latex",'FontSize', 14);
        grid on;
        filename = sprintf("%sa%d/wx_data_a%d_%s.jpg",image_path,index,index,num2str(round(a,1)));
        exportgraphics(fg,filename,'Resolution',600)

        % u vs x plots
        filename = sprintf("%sa%d/u_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data = readmatrix(filename);
        fg = figure();
        clf;
        plot(X,data(1,:),'Color','blue');
        hold on;
        for i=5:4:(size(data, 1))
            plot(X,data(i,:),'Color',[0.5,0.5,0.5]);
        end
        plot(X,data(end,:),'Color','red');
        %ylim([-3,3]);
        title("$u$ vs $x$ - a="+num2str(round(a,1)),"Interpreter","latex",'FontSize', 18);
        xlabel("$x$","Interpreter","latex",'FontSize', 14);
        ylabel("$u$","Interpreter","latex",'FontSize', 14);
        grid on;
        filename = sprintf("%sa%d/u_data_a%d_%s.jpg",image_path,index,index,num2str(round(a,1)));
        exportgraphics(fg,filename,'Resolution',600)

        % ux vs x plots
        filename = sprintf("%sa%d/ux_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data = readmatrix(filename);
        fg = figure();
        clf;
        plot(X,data(1,:),'Color','blue');
        hold on;
        for i=5:4:(size(data, 1))
            plot(X,data(i,:),'Color',[0.5,0.5,0.5]);
        end
        plot(X,data(end,:),'Color','red');
        ylim([-3,3]);
        title("$u_x$ vs $x$ - a="+num2str(round(a,1)),"Interpreter","latex",'FontSize', 18);
        xlabel("$x$","Interpreter","latex",'FontSize', 14);
        ylabel("$u_x$","Interpreter","latex",'FontSize', 14);
        grid on;
        filename = sprintf("%sa%d/ux_data_a%d_%s.jpg",image_path,index,index,num2str(round(a,1)));
        exportgraphics(fg,filename,'Resolution',600)

        % |w_c| vs K plots
        filename = sprintf("%sa%d/spec_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data = readmatrix(filename);
        fg = figure();
        clf;
        semilogy(K(K>0),data(1,:),'Color','blue');
        hold on;
        for i=5:4:(size(data, 1))
            semilogy(K(K>0),data(i,:),'Color',[0.5,0.5,0.5]);
        end
        semilogy(K(K>0),data(end,:),'Color','red');
        title("$|\hat{\omega}_k|$ vs $k$ - a="+num2str(round(a,1)),"Interpreter","latex",'FontSize', 18);
        xlabel("$k$","Interpreter","latex",'FontSize', 14);
        ylabel("$|\hat{\omega}_k|$","Interpreter","latex",'FontSize', 14);
        grid on;
        filename = sprintf("%sa%d/spec_data_a%d_%s.jpg",image_path,index,index,num2str(round(a,1)));
        exportgraphics(fg,filename,'Resolution',600)

        
        % L2 vs t plots
        filename = sprintf("%sa%d/l2_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data = readmatrix(filename);
        fg = figure();
        hold off;
        clf;
        plot(linspace(0,length(data(1,:))*DT,length(data(1,:))),data(1,:),'Color','b');
        title("$\|\omega_x\|$ vs $t$ - a="+num2str(round(a,1)),"Interpreter","latex",'FontSize', 18);
        xlabel("$t$","Interpreter","latex",'FontSize', 14);
        ylabel("$\|\omega_x\|$","Interpreter","latex",'FontSize', 14);
        grid on;
        filename = sprintf("%sa%d/l2_data_a%d_%s.jpg",image_path,index,index,num2str(round(a,1)));
        exportgraphics(fg,filename,'Resolution',600)

        % dgdt L2 vs L2 plots
        filename = sprintf("%sa%d/l2_y_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data_y = readmatrix(filename);
        filename = sprintf("%sa%d/l2_x_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data_x = readmatrix(filename);
        fg = figure();
        clf;
        loglog(data_x,data_y,'Color','blue');

        if index ~= 1
            hold on;
            x = data_x;
            y = data_y;
            x = transpose(x);
            y = transpose(y);
            s = 1250; %0.025
            for j=1:12
                %x1 = x(end*(40-j)*s:min(end,end*(41-j)*s));
                %y1 = y(end*(40-j)*s:min(end,end*(41-j)*s));
                x1 = x(end-s*j:end-s*(j-1));
                y1 = y(end-s*j:end-s*(j-1));

                ft = fittype('a*(x^b)','independent','x','dependent','y');
                startPoints = [1E-3,2];
                %lower = [1E-10,0.1];
                [f1,~] = fit(x1,y1,ft,'StartPoint',startPoints);
                plot(x1,f1(x1));
                hold on;
                alphas_l2(j,index) = f1.b;
            end
            
        end

        title("$\frac{d}{dt}\|\omega_x\|$ vs $\|\omega_x\|$ - a="+num2str(round(a,1)),"Interpreter","latex",'FontSize', 18);
        xlabel("$\|\omega_x\|$","Interpreter","latex",'FontSize', 14);
        ylabel("$\frac{d}{dt}\|\omega_x\|$","Interpreter","latex",'FontSize', 14);
        grid on;
        filename = sprintf("%sa%d/l2_dgdt_data_a%d_%s.jpg",image_path,index,index,num2str(round(a,1)));
        exportgraphics(fg,filename,'Resolution',600)

        
        %plot linf norm vs t
        filename = sprintf("%sa%d/linf_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data = readmatrix(filename);
        fg = figure();
        clf;
        plot(linspace(0,length(data_y(1,:))*DT,length(data_y(1,:))),data_y);
        title("$\|u_x\|_\infty$ vs $t$ - a="+num2str(round(a,1)),"Interpreter","latex",'FontSize', 18);
        xlabel("$t$","Interpreter","latex",'FontSize', 14);
        ylabel("$\|u_x\|_\infty$","Interpreter","latex",'FontSize', 14);
        grid on;
        filename = sprintf("%sa%d/linf_data_a%d_%s.jpg",image_path,index,index,num2str(round(a,1)));
        exportgraphics(fg,filename,'Resolution',600)

        %convert linf norm vs t to integral linf norm from 0 to t vs t
        filename = sprintf("%sa%d/linf_x_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data = readmatrix(filename);
        loglog(linspace(0,length(data(1,:))*DT,length(data(1,:))),data,'Color','blue');
        title("$\int_0^t \|u_x(\tau)\|_\infty d\tau$ vs $t$ - a="+num2str(round(a,1)),"Interpreter","latex",'FontSize', 18);
        xlabel("$t$","Interpreter","latex",'FontSize', 14);
        ylabel("$\int_0^t \|u_x(\tau)\|_\infty d\tau$","Interpreter","latex",'FontSize', 14);
        grid on;
        filename = sprintf("%sa%d/linf_int_data_a%d_%s.jpg",image_path,index,index,num2str(round(a,1)));
        exportgraphics(fg,filename,'Resolution',600)

        % dgdt LINF int vs LIF int plots
        fg = figure();
        clf;
        filename = sprintf("%sa%d/linf_x_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data_x = readmatrix(filename);
        filename = sprintf("%sa%d/linf_y_data_a%d_%s.csv",data_path,index,index,num2str(round(a,1)));
        data_y = readmatrix(filename);
        loglog(data_x,data_y,'Color','blue');

        if index ~= 1
            hold on;
            x = data_x;
            y = data_y;
            x = transpose(x);
            y = transpose(y);
            s = 1250; %0.025
            for j=1:12
                %x1 = x(end*(40-j)*s:min(end,end*(41-j)*s));
                %y1 = y(end*(40-j)*s:min(end,end*(41-j)*s));
                x1 = x(end-s*j:end-s*(j-1));
                y1 = y(end-s*j:end-s*(j-1));
    
                ft = fittype('a*(x^b)','independent','x','dependent','y');
                startPoints = [1E-3,2];
                %lower = [1E-10,0.1];
                [f1,~] = fit(x1,y1,ft,'StartPoint',startPoints);
                plot(x1,f1(x1));
                hold on;
                alphas_linf(j,index) = f1.b;
            end
            
        end

        title("$\|u_x(\tau)\|_\infty$ vs $\int_0^t\|u_x(\tau)\|_{\infty} d\tau$ - a="+num2str(round(a,1)),"Interpreter","latex",'FontSize', 18);
        xlabel("$\int_0^t\|u_x(\tau)\|_{\infty} d\tau$","Interpreter","latex",'FontSize', 14);
        ylabel("$\|u_x\|_\infty$","Interpreter","latex",'FontSize', 14);
        grid on;
        filename = sprintf("%sa%d/linf_int_dgdt_data_a%d_%s.jpg",image_path,index,index,num2str(round(a,1)));
        exportgraphics(fg,filename,'Resolution',600)
        
        close all;
        index = index + 1;
        
    end 
        

    fg = figure();
    clf;
    plot(as,alphas_l2(1,:),'color','r');
    hold on;
    for j=2:11
        plot(as,alphas_l2(j,:),'color',[0.5,0.5,0.5]);
        hold on;
    end
    plot(as,alphas_l2(12,:),'color','b');
    xlim([-1.1,1.1]);
    grid on;
    %legend('$\alpha_1$','$\alpha_2$','$\alpha_3$','interpreter','latex');
    xlabel("a",'Interpreter','latex','FontSize', 14);
    title("$\alpha$ vs $a$ - $\|\omega_x\|$",'Interpreter','latex','FontSize', 18);
    ylabel("$\alpha$",'Interpreter','latex','FontSize', 14);
    filename = sprintf("%salpha_l2_vs_a.jpg",image_path);
    exportgraphics(fg,filename,'Resolution',600)

    fg = figure();
    clf;
    plot(as,alphas_linf(1,:),'color','r');
    hold on;
    for j=2:11
        plot(as,alphas_linf(j,:),'color',[0.5,0.5,0.5]);
        hold on;
    end
    plot(as,alphas_linf(12,:),'color','b');
    xlim([-1.1,1.1]);
    grid on;
    xlabel("a",'Interpreter','latex','FontSize', 14);
    title("$\alpha$ vs $a$ - $\|u_x\|_{\infty}$",'Interpreter','latex','FontSize', 18);
    ylabel("$\alpha$",'Interpreter','latex','FontSize', 14);
    filename = sprintf("%salpha_inf_vs_a.jpg",image_path);
    exportgraphics(fg,filename,'Resolution',600)
end

% Generalized Constantin-Lax-Majda Equation solver
function SolveGCLM(init)
    global X;
    
    index = 20;
    for a=-0.9:-0.1:-1.0
        disp("a: "+a);

        % Initial conditions
        if init=="Init_1"
            w_0 = sin(X)+(0.1*sin(2*X));
        else
            w_0 = exp(sin(2*X)) - exp(cos(X));
        end
        RK4(w_0, a, index, init); 
       
        index = index + 1;
    end      
end

% L2 norm of w_x from w
function n = L2_norm(w)
    global N;
    n = norm(ifft(deriv(fft(w,N)),N),2);
end

% LINF norm of u_x from w
function n = LINF_norm(w)
    global N;
    n = norm(ifft(ht(fft(w,N)),N),Inf);
end

% Derivative in Fourier Space
function d = deriv(w_c)
    global K;
    d = (1i.*K).*w_c;
end

% Hilbert Transform in Fourier space
function h = ht(w_c)
    global K;
    h = (-1i*sign(K)).*w_c;
end

% Calculate velocity from vorticity
function w_c = calc_u(w_c)
    global N; global K;
    % exclude k=0 
    points = (K ~= 0);
    w_c(~points) = 0;
    % divide Fourier coefficients by k~=0
    w_c(points) = -sign(K(points)).*w_c(points)./K(points);
    
end

% Dealias with 2/3 rule
function w_c = dealias(w_c)
    global N; global K2; global K;
    % points above 2/3
    %dealias_points = (abs(K2) >= N/2);
    % set dealias points to zero
    %w_c(dealias_points) = 0;
    %w_c = w_c(abs(K)>=0);
    w_c = (2/3)*[w_c(1:N/2),w_c(end-N/2+1:end)];

end

% filter |w_c| <= 1E-12 to zero
function w_c = filter(w_c)
    global F;
    % find points below threshold
    filter_points = (abs(w_c) <= F);
    % set filter points to zero
    w_c(filter_points) = 0;

    
end

% Calculate u_x*w - a*u*w_x
function f = calc_RHS(w,a)
    global N; global K; global K2; global X; global X2;

    % calc values
    w_c = fft(w,N);
    ux_c = ht(w_c); % calculate du/dx in Fourier space
    wx_c = deriv(w_c); % calculate dw/dx in Fourier Space
    u_c = calc_u(w_c); % calculate u in physical space

    % pad in fourier space
    w_c2 =  (3/2)*[w_c(1:N/2),zeros(1,2^13),w_c(N/2 + 1:end)];
    ux_c2 = (3/2)*[ux_c(1:N/2),zeros(1,2^13),ux_c(N/2 + 1:end)];
    wx_c2 = (3/2)*[wx_c(1:N/2),zeros(1,2^13),wx_c(N/2 + 1:end)];
    u_c2 =  (3/2)*[u_c(1:N/2),zeros(1,2^13),u_c(N/2 + 1:end)];

    w2 = ifft(w_c2);
    ux2 = ifft(ux_c2,3*N/2); % convert du/dx to physical space
    wx2 = ifft(wx_c2,3*N/2); % convert dw/dx to physical space
    u2 = ifft(u_c2,3*N/2);

    t1 = ux2.*w2; % compute du/dx * w in physical space
    t2 = a*(u2.*wx2); % compute a*u*dw/dx in physical space

    f = t1 - t2; %f = u_x*w - a*u*w_x
    f = ifft(filter(dealias(fft(f,3*N/2))),N);
    
end

function save_data(w, a, index, init_path)
    global N; global K;
    % Save w and u data every 0.25  
    filename = sprintf('../Results/Data/%s/a%s/w_data_a%s_%s.csv', init_path, num2str(index), num2str(index), num2str(round(a,1)));
    writematrix(real(w), filename, 'WriteMode', 'append');
        
    w_c = fft(w,N);
    spec = filter(abs(w_c));
    filename = sprintf('../Results/Data/%s/a%s/spec_data_a%s_%s.csv', init_path,num2str(index), num2str(index), num2str(round(a,1)));
    writematrix(real(spec(K>0)), filename, 'WriteMode', 'append');
        
    u_c = calc_u(w_c);
    u = ifft(u_c,N);
    filename = sprintf('../Results/Data/%s/a%s/u_data_a%s_%s.csv', init_path, num2str(index), num2str(index), num2str(round(a,1)));
    writematrix(real(u), filename, 'WriteMode', 'append');
        
    ux = ifft(ht(w_c),N);
    filename = sprintf('../Results/Data/%s/a%s/ux_data_a%s_%s.csv', init_path, num2str(index), num2str(index), num2str(round(a,1)));
    writematrix(real(ux), filename, 'WriteMode', 'append');
        
    wx = ifft(deriv(w_c),N);
    filename = sprintf('../Results/Data/%s/a%s/wx_data_a%s_%s.csv', init_path, num2str(index), num2str(index), num2str(round(a,1)));
    writematrix(real(wx), filename, 'WriteMode', 'append');

end

% Integration of the model problem using the RK4 scheme;
function RK4( w_t, a, index,init_path)
    global DT; global X; global T; global N; global K; global K2;
    global DEALIAS_MAX;
    i = 1;
    t = 0;
    b = false;

    norm_data = [];
    norm_data_2 = [];
    norm_times = [];

    save_data(w_t, a, index, init_path);
    norm_data(i) = L2_norm(w_t);
    norm_data_2(i) = LINF_norm(w_t);
    norm_times(i) = t;
    disp(t);

    while ( all([(t < T),(~b)]) )

      t = i*DT;
     
      f1 = calc_RHS(w_t, a);
      k1 = DT*f1;
      f2 = calc_RHS(w_t + k1/2, a);
      k2 = DT*f2;
      f3 = calc_RHS(w_t + k2/2, a);
      k3 = DT*f3;
      f4 = calc_RHS(w_t + k3, a);
      k4 = DT*f4;

      w_t = w_t + (1/6)*( k1 + 2*k2 + 2*k3 + k4 );
  
      if(mod(t,0.25)==0)
        disp(t);
        save_data(w_t, a, index, init_path);
      end
      i = i + 1;
      norm_data(i) = L2_norm(w_t);
      norm_data_2(i) = LINF_norm(w_t);
      norm_times(i) = t;

      % underresolved check
      x = filter(abs(fft(w_t,N)));
      cut = (N/2);
      b = logical(mean(x(cut-100:cut))>DEALIAS_MAX);
    end

    if b
        disp("Under resolved");
        disp(t);
    else
        disp("DONE");
        disp(t);
    end

    save_data(w_t, a, index, init_path);
    
    filename = sprintf('../Results/Data/%s/a%s/l2_data_a%s_%s.csv', init_path, num2str(index), num2str(index), num2str(round(a,1)));
    writematrix(norm_data, filename, 'WriteMode', 'append');
    
    filename = sprintf('../Results/Data/%s/a%s/linf_data_a%s_%s.csv', init_path, num2str(index), num2str(index), num2str(round(a,1)));
    writematrix(norm_data_2, filename, 'WriteMode', 'append');

    filename = sprintf('../Results/Data/%s/blowup_times.csv', init_path);
    writematrix(t, filename, 'WriteMode', 'append');
    
end