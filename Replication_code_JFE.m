%% Methodology 1 (Fernando Duarte code)
% Table 2: Unconditional and Conditional Predictive Regressions of Consumption Growth on Inflation

load('data_replication_jfe.mat','dates','inflation','cg')

Tstart=393;% start of consumption data for estimation of nominal-real covariance (Feb 1959, beginning of CG data)
Tend=1063; %final observation (Dec 2014)

g=12;
%this for loop adds second column to cg: compounded annual consumption growth
for t=Tstart:Tend-12+1 
    cg(t,2)=prod(1+cg(t:t+12-1,1))-1; 
end                                   

period=60;
C=cg(Tstart:Tend,2);
I=inflation(Tstart-2:Tend-2,1);

% One-stage uncondtional regression
disp('Unconditional predictive regression at horizon K=12 (column 4 of Table 2)')
uncond_x=[ones(570-11,1), I(end-569:end-11,1)];
[uncond_coeff,~,~,~,uncond_stats] = regress(C(end-569:end-11,1), uncond_x);
r_squared_uncond=uncond_stats(1)*100;

% Two-stage conditional regression
for t=1:length(I)-period
    clear Kw
    h=log(2)/(period);
    ts=t+period-12+1;
    for tt=1:(ts-1)
        Kw(tt,1)=exp(-abs(ts-tt-1)*h);
    end
    Kw=Kw(:,1)/sum(Kw(:,1));
    NRC(t+period,1:2)=(lscov([ones(t+period-12,1) I(1:t+period-12,1)],C(1:t+period-12,1),Kw(1:end)))';
  
end


disp('Conditional predictive regression at horizon K=12 (column 8 of Table 2)')
cond_x=[ones(570-11,1), NRC(end-569:end-11,1)+NRC(end-569:end-11,2).*I(end-569:end-11,1)];
[cond_coeff,~,~,~,cond_stats] = regress(C(end-569:end-11,1), cond_x);
r_squared_cond=cond_stats(1)*100;


row=["intercept";"slope";"R_squared (%)"];
uncond=[uncond_coeff;r_squared_uncond];
cond=[cond_coeff;r_squared_cond];
table(row,uncond,cond)

%% Methodology 1 (Lorenzo code)
% Load data from Excel file
filename = 'Inflation data.xlsx';
data = readtable(filename);

% Define dependent and independent variables
y3 = table2array(data(27:697, 2)); % Column 2: consumption growth (dependent variable)

X_core_CPI = table2array(data(27:697, 8)); % Column 8: core CPI
X_CPI_less_shelter = table2array(data(27:697, 12)); % Column 12: CPI less shelter
X_PCE = table2array(data(27:697, 16)); % Column 16: PCE

% UNCONDITIONAL REGRESSIONS
% Unconditional regression on core CPI
X_uc_core = [ones(length(X_core_CPI), 1), X_core_CPI];
[b_uc_core,~,~,~,stats_uc_core] = regress(y3, X_uc_core);
r_squared_uc_core = stats_uc_core(1)*100;
% Unconditional regression on CPI less shelter
X_uc_CPI_less_shelter = [ones(length(X_CPI_less_shelter), 1), X_CPI_less_shelter];
[b_uc_CPI_less_shelter,~,~,~,stats_uc_CPI_less_shelter] = regress(y3, X_uc_CPI_less_shelter);
r_squared_uc_CPI_less_shelter = stats_uc_CPI_less_shelter(1)*100;
% Unconditional regression on PCE
X_uc_PCE = [ones(length(X_PCE), 1), X_PCE];
[b_uc_PCE,~,~,~,stats_uc_PCE] = regress(y3, X_uc_PCE);
r_squared_uc_PCE = stats_uc_PCE(1)*100;


% CONDITIONAL REGRESSIONS
period = 60; % half-life
% Function to compute weighted regression and display results
function c=concatenate_conditional_regression(y, X, period)
    for t = 1:length(X) - period
        clear Kw
        h = log(2) / period;
        ts = t + period - 12 + 1;
        for tt = 1:(ts - 1)
            Kw(tt, 1) = exp(-abs(ts - tt - 1) * h);
        end
        Kw = Kw(:, 1) / sum(Kw(:, 1));
        NRC(t + period, 1:2) = (lscov([ones(t + period - 12, 1) X(1:t + period - 12)], y(1:t + period - 12), Kw(1:end)))';
    end
    X_cond = [ones(length(X), 1), NRC(:, 1) + NRC(:, 2) .* X];
    [b_cond,~,~,~,stats_cond] = regress(y, X_cond);
    r_squared_cond = stats_cond(1);
    c=[b_cond;r_squared_cond];
end

disp('Unconditional table')
row=["intercept";"slope";"R_squared (%)"];
core_CPI=[b_uc_core;r_squared_uc_core];
CPI_less_shelter=[b_uc_CPI_less_shelter;r_squared_uc_CPI_less_shelter];
PCE=[b_uc_PCE;r_squared_uc_PCE];
table(row,core_CPI,CPI_less_shelter,PCE)

disp('Conditional table')
row=["intercept";"slope";"R_squared (%)"];
core_CPI=concatenate_conditional_regression(y3, X_core_CPI, period);
CPI_less_shelter=concatenate_conditional_regression(y3, X_CPI_less_shelter, period);
PCE=concatenate_conditional_regression(y3, X_PCE, period);
table(row,core_CPI,CPI_less_shelter,PCE)

%% Methodology 2 (Fernando Duarte code)
% Table 3: Characterizing Inflation Beta-Sorted Portfolios and the Inflation Risk Premium
load('data_replication_jfe.mat','retHL','inflation_arma')
% retHL: 10 portfolios sorted on inflation beta (from High to Low) and H-L. DIM: 570x11
% inflation_arma: arma(1,1) innovation in inflation. DIM: 570x1
% 570 monthly observations (from July 1967 to December 2014)

% Average return (Row 3 of Table 3)
disp('Average return (Table 3, Panel B)')
mean(retHL)*1200 %from decimal to percentage --> displayed in the ouput

T=length(retHL);
for i=1:11
    Beta_pi_post(:,i)=regress(retHL(:,i),[ones(T,1) inflation_arma]);
end

disp('Estimated coefficients: intercept, slope (Table 3, Panel A)')
disp(Beta_pi_post);