clear all
close all
clc

addpath(genpath('TOOLBOX_DWTRed_Frame'))
addpath data
addpath functions

rng(5); % noise seed, to reproduce results

%% load data
load('Urban')
Xtrue     = Urban.Xtrue;     % ground-truth of the abundance maps
nEnd      = Urban.nEnd;      % number of materials or endmembers
materials = Urban.materials; % list of names of materials
nRow      = Urban.nRow;      % number of rows for the abundance maps
nCol      = Urban.nCol;      % number of columns for the abundance maps
S         = Urban.lib;       % library with the spectral signatures of the materials

%% Generate attenuated noisy hyperspectral data
%%% Add atenuation to satisfy the constraint <=1, attenuation is different for every pixel
eta   =  0.05;          % eta = 0 -> no attenuation
t1    = (1:nRow)/nRow;
t2    = (1:nCol)/nCol;
Latt  =  rand(2,1);
Att   =  kron(exp(-10*(t2-Latt(1)).^2),exp(-10*(t1-Latt(2)).^2)')*5;
Att   =  Att(:)';
Xatt  = (1-eta*Att).*Xtrue; % attenuated abundance maps
Y     =  S*Xatt;            % hyperspectral data

%%% Add noise
noise_snr = 1.5;                                  % signal-to-noise ratio
noise_std = mean(std(Y,[],2))*10^(-noise_snr/20); % Gaussian noise standard deviation
Y         = Y + randn(size(Y)).*noise_std;
        
%% Parameters
reg       = 0.01; % regularization parameter
time_max  = 40;   % time for running the algorithm

%% load solution
load('solution') % solution to the minimization problem (only used to plot distance to solution)

%% run PIPA
[X,obj,snr,X_Xinf,time] = PIPA(Y,S,reg,time_max,nRow,nCol,Xtrue,Xinf);

%% plot figures
color=jet(20);

%%% draw restored abondance maps
for p = 1:nEnd 
    figure  
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.4, 0.5, 0.5]);
    %%% ground-truth
    subplot(1,2,1)
    imagesc(reshape(Xtrue(p,:),nRow,nCol),[0 1]); axis image off; colorbar
    title('Ground-truth','Interpreter','latex','fontsize',15)
    %%% PIPA
    subplot(1,2,2)
    imagesc(reshape(X(p,:),nRow,nCol),[0 1]); axis image off; colorbar
    title(strcat('PIPA--',num2str(time(end),'%.2f'),'s'),'Interpreter','latex','fontsize',15)
    ax=axes('Units','Normal','Position',[.075 .075 .85 .75],'Visible','off');
    set(get(ax,'Title'),'Visible','on')
    title(materials{p},'Interpreter','latex','fontsize',20);
    h=get(ax,'Title');
end

%%% objective function
figure
plot(time,obj,'-o','Linewidth',1,'Color',color(2,:))
xlabel({'Time (s)'},'Interpreter','latex','fontsize',15)
ylabel({'Criterion'},'Interpreter','latex','fontsize',15)
xlim([0 time(end)])
ylim([0.9*min(obj),1.1*max(obj)])
title('Objective function','Interpreter','latex','fontsize',20)

%%% signal-to-noise ratio
figure
plot(time,snr,'-o','Linewidth',1,'Color',color(2,:))
xlabel({'Time (s)'},'Interpreter','latex','fontsize',15)
ylabel({'SNR'},'Interpreter','latex','fontsize',15)
xlim([0 time(end)])
ylim([1.1*min(snr),1.1*max(snr)])
title('Signal-to-noise ratio','Interpreter','latex','fontsize',20)

%%% distance to solution
figure
semilogy(time,X_Xinf,'-o','Linewidth',1,'Color',color(2,:))
xlabel({'Time (s)'},'Interpreter','latex','fontsize',15)
ylabel({'$\|x-x_{\infty}\|/\|x_{\infty}\|$'},'Interpreter','latex','fontsize',15)
xlim([0 time(end)])
ylim([10^(-3),1.1*max(X_Xinf)])
title('Distance to solution','Interpreter','latex','fontsize',20)

