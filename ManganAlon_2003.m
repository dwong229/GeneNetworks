%% Simulating Gene Networks

% Ref: Mangan + Alon PNA 2003 Structure and function of the feed-forward
% loop network motif

% mathematical model for the transciption of a variey of network structures
clear all
close all
%% initialize concentrations of transcription factors
X = 1;
Y = 1;
% Gene expression of protein Z
Z = 0;

% inducers
Sx = 1;
Sy = 0;

Alphaz = 1; % decay rate of Z.

Xactive = X*Sx;
Yactive = Y*Sy;

% TIMING
tstart = 0;
tfinal = 10;
tevent = [1 8]; % times when Sx changes
tstep = 0.02;

state0 = [Xactive Yactive Z]';

%% solve
% using ode45     
    
    % tstart to tevent(1)
    Sx = 0;
    tvec = [tstart;tevent(1)];
    statemat = [Xactive Yactive Z;Xactive Yactive Z];
    
    
    % tevent(1) to tevent(2)
    tspan = [tevent(1) tevent(2)];
    Sx = 1;
    [tvector,rxnstate] = ode45(@compute_rxnderivatives,tspan,state0);
    tvec = [tvec;tvector];
    statemat = [statemat;rxnstate];
    
    % tevent(2) to tfinal
    tspan = [tevent(2):tstep:tfinal]';
    rxnstate = [];
    zvec = statemat(end,3)*exp(-Alphaz*(tspan - tevent(2)));
    rxnstate = zeros(length(tspan),3);
    rxnstate(:,3) = zvec;
    tvec = [tvec;tspan];
    statemat = [statemat;rxnstate];
    
% plot
h1 = figure;
% subplot
subplot(2,1,1)
stairs([tstart tevent tfinal],[0 1 0 0],'LineWidth',5);
ylabel('Input Sx')
subplot(2,1,2)
plot(tvec,statemat(:,3),'LineWidth',5);
ylabel('Z')
xlabel('Time')

