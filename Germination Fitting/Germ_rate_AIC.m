% Germ_rate_AIC.m - calculates the AICc value for various model fits to the
% germination sigmodal curves from Swanson.
% Created 7/14/20 by David Elzinga
% Compute time < 1min.

clear all 
close all 
clc

% These numbers need to change when Swanson emails back.
tot_col = 1040.56;
tot_ler = 1040.56;

% Specify the raw data for Col and Ler. The raw data = the number of pollen
% grains of each accesssion that have germinated at some time point t.
global time
time = [10 20 30 40];

raw_data_col_10 = [17	75	0	0	17	8	6	2	3]';
raw_data_col_20 = [85	90	51	333	52	69	239	103	128]';
raw_data_col_30 = [249	612	695	351	624	572	299	556	612]';
raw_data_col_40 = [386	403	580	675	432	406	588	390	531]';
raw_data_col = horzcat(raw_data_col_10,raw_data_col_20,raw_data_col_30,raw_data_col_40);

raw_data_ler_10 = [0	3	4	4	2	0	1	0	0]';
raw_data_ler_20 = [55	56	13	22	42	13	36	15	31]';
raw_data_ler_30 = [236	254	286	288	275	192	400	338	388]';
raw_data_ler_40 = [100	210	300	250	254	558	401	312	362]';
raw_data_ler = horzcat(raw_data_ler_10,raw_data_ler_20,raw_data_ler_30,raw_data_ler_40);

model_time_plot = linspace(0,40,1000);

% Calculation the proportional data for Col and Ler.
global prop_data_col prop_data_ler
prop_data_col = raw_data_col ./ tot_col;
prop_data_ler = raw_data_ler ./ tot_ler;

% For each accession we will attempt to fit the following models:
%
% Gompertz_3 = Gompertz model with 3 MATLAB fitted parameters.
% Gompertz_2 = Gompertz model with 2 MATLAB fitted parameters.
% Logistic_3 = Logistic model with 3 MATLAB fitted parameters.
% Logistic_2 = Logistic model with 2 MATLAB fitted parameters.

% Specify the accession.
global accession
for accession = ["col" "ler"]
    opts = optimset('fminsearch');
    opts.Display = 'iter';
    opts.TolX = 1.e-120;
    opts.MaxFunEvals = 100000;
    opts.MaxIter = 100000;
    
    if strcmp(accession,"col")
        fig_num = 1;
        data = prop_data_col;
    else
        fig_num = 2;
        data = prop_data_ler;
    end
    
    n = length(data)*length(time);
    
    figure(fig_num)
    p(1) = scatter(10*ones(9,1),data(:,1)','bo');
    hold on 
    p(2) = scatter(20*ones(9,1),data(:,2)','bo');
    hold on 
    p(3) = scatter(30*ones(9,1),data(:,3)','bo');
    hold on 
    p(4) = scatter(40*ones(9,1),data(:,4)','bo');
    hold on
    
    
    % For the Gompertz_3 model, let MATLAB fit all 3 parameters.
    if strcmp(accession,"col")
    Gompertz_3_fit_params = fminsearchbnd(@cost_Gompertz,[0.1 10 0.25],[0 0 0],[],opts);
    else
        Gompertz_3_fit_params = fminsearchbnd(@cost_Gompertz,[0.1 10 0.25],[0 0 0],[],opts);
    end
    Gompertz_3_min_SSE = cost_Gompertz(Gompertz_3_fit_params);
    figure(fig_num)
    p(5) = plot(model_time_plot,Gompertz_3_fit_params(1) * exp(Gompertz_3_fit_params(2)/Gompertz_3_fit_params(3) * (1 - exp(-model_time_plot*Gompertz_3_fit_params(3)))));
    hold on
    Lhat = Gompertz_3_min_SSE/n;
    k = 4;
    AIC(fig_num,1) = 2*k -2*log(Lhat) + (2*k*(k+1))/(n - k - 1);
    
    % For the Gompertz 2 model, let MATLAB fit just 2 of the parameters. 
    Gompertz_2_fit_params = fminsearchbnd(@cost_Gompertz,[10 0.25],[0 0],[],opts);
    Gompertz_2_min_SSE = cost_Gompertz(Gompertz_2_fit_params);
    a = mean(data(:,end))*exp(-Gompertz_2_fit_params(1)/Gompertz_2_fit_params(2));
    figure(fig_num)
    p(6) = plot(model_time_plot,a * exp(Gompertz_2_fit_params(1)/Gompertz_2_fit_params(2) * (1 - exp(-model_time_plot*Gompertz_2_fit_params(2)))));
    hold on 
    Lhat = Gompertz_2_min_SSE/n;
    k = 3;
    AIC(fig_num,2) = 2*k -2*log(Lhat) + (2*k*(k+1))/(n - k - 1);
    
    % For the Logistic 3 model, let MATLAB fit all 3 of the parameters. 
    Logistic_3_fit_params = fminsearchbnd(@cost_Logistic,[0.5 10 0.25],[0 0 0],[],opts);
    Logistic_3_min_SSE = cost_Logistic(Logistic_3_fit_params);
    figure(fig_num)
    p(7) = plot(model_time_plot, Logistic_3_fit_params(1) ./ (1 + Logistic_3_fit_params(2)*exp(-Logistic_3_fit_params(3)*model_time_plot)));
    hold on 
    Lhat = Logistic_3_min_SSE/n;
    k = 4;
    AIC(fig_num,3) = 2*k -2*log(Lhat) + (2*k*(k+1))/(n - k - 1);
    
    % For the Logistic 2 model, let MATLAB just 3 of the parameters. 
    Logistic_2_fit_params = fminsearchbnd(@cost_Logistic,[10 0.25],[0 0],[],opts);
    Logistic_2_min_SSE = cost_Logistic(Logistic_2_fit_params);
    figure(fig_num)
    p(8) = plot(model_time_plot, mean(data(:,end)) ./ (1 + Logistic_2_fit_params(1)*exp(-Logistic_2_fit_params(2)*model_time_plot)));
    hold on 
    Lhat = Logistic_2_min_SSE/n;
    k = 3;
    AIC(fig_num,4) = 2*k -2*log(Lhat) + (2*k*(k+1))/(n - k - 1);
    
    legend([p(1) p(5) p(6) p(7) p(8)],"Data","Gompertz, k=3", "Gompertz, k=2", "Logistic, k=3", "Logistic, k=2",'Location','northwest')
    
end

















function total_cost = cost_Gompertz(x)

global time
global accession
global prop_data_col
global prop_data_ler

if strcmp(accession,"col")
    data = prop_data_col;
else 
    data = prop_data_ler;
end

% If we are fitting Gompertz 3:
if length(x) == 3
    a = x(1);
    b = x(2);
    c = x(3);
elseif length(x) == 2 % If we are fitting Gompertz 2:
    b = x(1);
    c = x(2);
    
    d = mean(data(:,end));
    a = d*exp(-b/c);
end

model = a * exp(b/c * (1 - exp(-time*c)));

total_cost = 0;       % Initialize.
for i=1:length(time)  % Loop over every entry in the time vector.
    for j = 1:1:length(data)
     this_cost = (model(i) - data(j,i))^2;  % Square-difference of model and data.
     total_cost = total_cost + this_cost; % Add in the square-difference for this datum to the total.
    end
end

end







function total_cost = cost_Logistic(x)

global time
global accession
global prop_data_col
global prop_data_ler

if strcmp(accession,"col")
    data = prop_data_col;
else 
    data = prop_data_ler;
end

% If we are fitting Logistic 3:
if length(x) == 3
    K = x(1);
    c = x(2);
    r = x(3);
elseif length(x) == 2 % If we are fitting Logistic 2:
    K = mean(data(:,end));
    c = x(1);
    r = x(2);
end

model = K ./ (1 + c*exp(-r*time));

total_cost = 0;       % Initialize.
for i=1:length(time)  % Loop over every entry in the time vector.
    for j = 1:1:length(data)
     this_cost = (model(i) - data(j,i))^2;  % Square-difference of model and data.
     total_cost = total_cost + this_cost; % Add in the square-difference for this datum to the total.
    end
end

end


