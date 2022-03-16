%% STOCHASTIC MODEL DELTA %%
%
% This script file uses a stochastic model to calculate the spread of 
% COVID-19 in New Zealand (Aotearoa) as a function of vaccination rate.
%
% The model is used in Watson, L. M. (2022) Likelihood of infecting or 
% getting infected with COVID-19 as a function of vaccination status, as 
% investigated with a stochastic model for Aotearoa New Zealand for Delta 
% and Omicron variants, New Zealand Medical Journal.
%
% This script file is formatted to run over all possible vaccination rates. 
% This version tracks the number of vaccinated and unvaccinated infections

clear all; clc;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
cmap = get(gca,'ColorOrder');
figure(1); clf;
figure(2); clf;
figure(3); clf;

%% Initialization %%

R = 6; % reproduction number of delta variant
R_sym = R; % reproduction number for symptomatic individuals / clinical infections
R_asym = R/2; % reproduction number for asymptomatic individuals / subclinical infections

prob_asym = 0.33; % probability of a new case being asymptomatic / subclinical
prob_hos = 7.8/100; % probability of a new case being hospitalized

v_eff_hos = 0.87; % effectiveness of vaccine against hospitalization

% effectiveness of vaccianation from "A COVID-19 vaccination model for
% Aotearoa New Zealand" by Steyn, N., Plank, M. J., Binny, R. N., Hendy, S.
% C., Lustig, A., and Ridings, K.
e_infect = 0.7; % effectiveness of vaccination against infection
e_transmit = 0.5; % effectiveness of vaccination against transmission

% effectiveness of alert levels
% realistic scenarios from Table 3 in "A stochastic model for COVID-19
% spread and the effects of alert level 4 in Aotearoa New Zealand" by
% Plank, M. J., Binny, R. N., Hendy, S. C., Lustig, A., James, A., and
% Steyn, N. https://doi.org/10.1101/2020.04.08.20058743.
% C1 = 1; % alert level 1
% C2 = 0.72; % alert level 2
% C3 = 0.52; % alert level 3
% C4 = 0.32; % alert level 4


% distribution for time between an individual becoming infected and
% infecting another individual (generation time)
a = 5.57;
b = 4.08;
TG_distrib = makedist('Weibull','a',a,'b',b);

% recovery time (assume individuals are recovered beyond this time for computational reasons)
recovery_time = 30; % recovery time in days

% simulation time
tmax = 28; % number of days to run simulation for
t = 0:1:tmax; % time vector

% size of generation time
nT = 25;

num_realizations = 100000; % number of simulations to run

VAC_RATE = 0.5:0.1:0.9;
CONTROLS = [1];

%% Run Simulations %%

for ii = 1:length(VAC_RATE)
    
    prob_vac = VAC_RATE(ii);
    
    for jj = 1:length(CONTROLS)
        
        C = CONTROLS(jj);
                
        infect_distrib_sym = makedist('Poisson','lambda',R_sym*C); % number of secondary infections for symptomatic case
        infect_distrib_asym = makedist('Poisson','lambda',R_asym*C); % number of secondary infections for asymptomatic case

        
        for p = 1:num_realizations
            
            % Display number of realizations
            disp(strcat('Realization Number=',num2str(p),'/',num2str(num_realizations)));
            
            % initialize simulation with one symptomatic individual / clinical infection
            Ic = 1; % introduce one symptomatic case
            [mc,~] = size(Ic); % total number of symptomatic cases (including those who have recovered)
            num_infect_c = sum(Ic); % number of active symptomatic cases
            T0c = 0; % time of initial infection
            Nc = random(infect_distrib_sym); % number of individuals infected by the symptomatic case
            TGc = zeros(1,nT); % initialize generation time vector (this assumes maximum number of people infected by a single case is 15)
            TGc(1:Nc) = T0c + random(TG_distrib,1,Nc); % infection times for new cases
            
            % initialize simulation with zero asymptomatic individuals / subclinical infections
            Is = 0; % do not introduce any asymptomatic cases
            ms = 0; % total number of asymptomatic cases (including those who have recovered)
            num_infect_s = sum(Is); % number of active asymptomatic cases
            T0s = 0; % time of initial infection
            Ns = random(infect_distrib_asym); % number of individuals infected by an asymptomatic case
            TGs = zeros(1,nT); % initialize generation time vector (this assumes maximum number of people infected by a single case is 15)
            TGs(1:Ns) = T0s + random(TG_distrib,1,Ns); % infection times for new cases
            
            % number of vaccinated and unvaccianted individuals
            NUM_VAC = zeros(size(t)); % vaccinated
            NUM_UNV = zeros(size(t)); % unvaccinated
            
            NUM_VAC(1) = 0;
            NUM_UNV(1) = 1;
            
            % number of hospitalizations
            VAC_HOS = zeros(size(t)); % vaccinated
            UNV_HOS = zeros(size(t)); % unvaccinated
            
            VAC_HOS(1) = 0;
            UNV_HOS(1) = 0;
            
            for i = 2:length(t) % iterate over time steps
                
                [mc,~] = size(Ic); % number of clinical infected people
                num_infect_c(i) = sum(Ic); % total number of clinical infected people at each time step
                
                [ms,~] = size(Is); % number of subclinical infected people
                num_infect_s(i) = sum(Is); % total number of subclinical infected people at each time step
                
                N = num_infect_c(end) + num_infect_s(end); % total number of infected people
                Nc_total = num_infect_c; % running total number of clinical infections
                N_total = num_infect_c + num_infect_s; % running total number of infected people
                
                %%% iterate over clinical infections %%%
                for j = 1:mc % iterate over clinically infected individuals
                    
                    if Ic(j) == 1 % if Ic(j) == 1 case is still active, else if Ic(j) == 0, then case has recovered
                        
                        N_tmp = Nc(j); % number of people infected by current individual
                        TG_tmp = TGc(j,1:N_tmp); % time that people are infected
                        
                        for k = 1:length(TG_tmp) % iterate over number of new people who are infected by current individual
                            
                            if t(i-1) < TG_tmp(k) && t(i) > TG_tmp(k) % if generation time is in current time step
                                
                                a1 = rand(1); % random number that determines if infection is clinical or subclincal
                                a2 = rand(1); % random number that determines if exposed person is vaccinated or not
                                a3 = rand(1); % random number that determines if vaccinated individual is infected
                                a4 = rand(1); % random number that determines if vaccinated individual transmits the various or not
                                a5 = rand(1); % random number that determines if case is hospitalized
                                a6 = rand(1); % random number that determines if vaccine is effective against hospitalization
                                
                                if a2 > prob_vac % then exposed person is unvaccinated
                                    
                                    NUM_UNV(i) = NUM_UNV(i) + 1;
                                    
                                    % hospitalization
                                    if a5 < prob_hos
                                        UNV_HOS(i) = UNV_HOS(i) + 1;
                                    end
                                    
                                    if a1 > prob_asym % create new clinical infection
                                        
                                        Ic_new = [Ic; 1]; % add a new infected person
                                        T0c_new = [T0c; TG_tmp(k)]; % time that people are infected
                                        Nc_new0 = random(infect_distrib_sym); % number of people that are infected
                                        Nc_new = [Nc; Nc_new0];
                                        TGc_new0 = zeros(1,nT);
                                        TGc_new0(1:Nc_new0) = TG_tmp(k) + random(TG_distrib,1,Nc_new0);
                                        TGc_new = [TGc; TGc_new0];
                                        
                                        % update fields
                                        Ic = Ic_new;
                                        T0c = T0c_new;
                                        Nc = Nc_new;
                                        TGc = TGc_new;
                                                                                
                                    else % create new subclinical infection
                                        
                                        Is_new = [Is; 1]; % add a new infected person
                                        T0s_new = [T0s; TG_tmp(k)]; % time that people are infected
                                        Ns_new0 = random(infect_distrib_asym); % number of people that are infected
                                        Ns_new = [Ns; Ns_new0];
                                        TGs_new0 = zeros(1,nT); % time
                                        TGs_new0(1:Ns_new0) = TG_tmp(k) + random(TG_distrib,1,Ns_new0);
                                        TGs_new = [TGs; TGs_new0];
                                        
                                        % update fields
                                        Is = Is_new;
                                        T0s = T0s_new;
                                        Ns = Ns_new;
                                        TGs = TGs_new;
                                        
                                    end
                                    
                                else % exposed person is vaccinated
                                    
                                    if a3 > e_infect % person is infected then a new infection is created (which can be either clinical or subclinical)
                                        
                                        NUM_VAC(i) = NUM_VAC(i) + 1;
                                        
                                        % hospitalization
                                        if a5 < prob_hos && a6 > v_eff_hos
                                            VAC_HOS(i) = VAC_HOS(i) + 1;
                                        end
                                        
                                        if a1 > prob_asym % create new clinical infection
                                            
                                            Ic_new = [Ic; 1]; % add a new infected person
                                            T0c_new = [T0c; TG_tmp(k)]; % time that people are infected
                                            
                                            if a4 < e_transmit
                                                Nc_new0 = 0;
                                                Nc_new = [Nc; Nc_new0];
                                                TGc_new0 = zeros(1,nT);
                                                TGc_new0(1:Nc_new0) = TG_tmp(k) + random(TG_distrib,1,Nc_new0);
                                                TGc_new = [TGc; TGc_new0];
                                                
                                            else
                                                
                                                Nc_new0 = random(infect_distrib_sym); % number of people that are infected
                                                Nc_new = [Nc; Nc_new0];
                                                TGc_new0 = zeros(1,nT);
                                                TGc_new0(1:Nc_new0) = TG_tmp(k) + random(TG_distrib,1,Nc_new0);
                                                TGc_new = [TGc; TGc_new0];
                                                
                                                % update fields
                                                Ic = Ic_new;
                                                T0c = T0c_new;
                                                Nc = Nc_new;
                                                TGc = TGc_new;
                                                
                                            end
                                            
                                        else % create new subclinical infection
                                            
                                            Is_new = [Is; 1]; % add a new infected person
                                            T0s_new = [T0s; TG_tmp(k)]; % time that people are infected
                                            
                                            if a4 < e_transmit
                                                Ns_new0 = 0;
                                                
                                                Ns_new = [Ns; Ns_new0];
                                                TGs_new0 = zeros(1,nT); % time
                                                TGs_new0(1:Ns_new0) = TG_tmp(k) + random(TG_distrib,1,Ns_new0);
                                                TGs_new = [TGs; TGs_new0];
                                                
                                                % update fields
                                                Is = Is_new;
                                                T0s = T0s_new;
                                                Ns = Ns_new;
                                                TGs = TGs_new;
                                                
                                            else
                                                
                                                Ns_new0 = random(infect_distrib_asym); % number of people that are infected
                                                Ns_new = [Ns; Ns_new0];
                                                TGs_new0 = zeros(1,nT); % time
                                                TGs_new0(1:Ns_new0) = TG_tmp(k) + random(TG_distrib,1,Ns_new0);
                                                TGs_new = [TGs; TGs_new0];
                                                
                                                % update fields
                                                Is = Is_new;
                                                T0s = T0s_new;
                                                Ns = Ns_new;
                                                TGs = TGs_new;
                                                
                                            end
                                            
                                            
                                        end
                                    end
                                end
                                
                                
                                
                            end
                        end
                    end
                    
                end
            end
            
            
            %%% iterate over subclinical infections %%%
            for j = 1:ms % iterate over subclinically infected individuals
                
                if Is(j) == 1 % if Ic(j) == 1 case is still active, else if Ic(j) == 0, then case has recovered
                    
                    N_tmp = Ns(j); % number of people infected by current individual
                    TG_tmp = TGs(j,1:N_tmp); % time that people are infected
                    
                    for k = 1:length(TG_tmp) % iterate over number of new people who are infected by current individual
                        
                        if t(i-1) < TG_tmp(k) && t(i) > TG_tmp(k) % if generation time is in current time step
                            
                            a1 = rand(1); % random number that determines if infection is clinical or subclincal
                            a2 = rand(1); % random number that determines if exposed person is vaccinated or not
                            a3 = rand(1); % random number that determines if vaccinated individual is infected
                            a4 = rand(1); % random number that determines if vaccinated individual transmits the various or not
                            
                            if a2 > prob_vac % then exposed person is unvaccinated
                                
                                NUM_UNV(i) = NUM_UNV(i) + 1;
                                
                                if a1 > prob_asym % create new clinical infection
                                    
                                    Ic_new = [Ic; 1]; % add a new infected person
                                    T0c_new = [T0c; TG_tmp(k)]; % time that people are infected
                                    Nc_new0 = random(infect_distrib_sym); % number of people that are infected
                                    Nc_new = [Nc; Nc_new0];
                                    TGc_new0 = zeros(1,nT);
                                    TGc_new0(1:Nc_new0) = TG_tmp(k) + random(TG_distrib,1,Nc_new0);
                                    TGc_new = [TGc; TGc_new0];
                                    
                                    % update fields
                                    Ic = Ic_new;
                                    T0c = T0c_new;
                                    Nc = Nc_new;
                                    TGc = TGc_new;
                                    
                                else % create new subclinical infection
                                    
                                    Is_new = [Is; 1]; % add a new infected person
                                    T0s_new = [T0s; TG_tmp(k)]; % time that people are infected
                                    Ns_new0 = random(infect_distrib_asym); % number of people that are infected
                                    Ns_new = [Ns; Ns_new0];
                                    TGs_new0 = zeros(1,nT); % time
                                    TGs_new0(1:Ns_new0) = TG_tmp(k) + random(TG_distrib,1,Ns_new0);
                                    TGs_new = [TGs; TGs_new0];
                                    
                                    % update fields
                                    Is = Is_new;
                                    T0s = T0s_new;
                                    Ns = Ns_new;
                                    TGs = TGs_new;
                                    
                                end
                                
                            else % expose person is vaccinated
                                
                                if a3 > e_infect % person is infected then a new infection is created (which can be either clinical or subclinical)
                                    
                                    NUM_VAC(i) = NUM_VAC(i) + 1;
                                    
                                    if a1 > prob_asym % create new clinical infection
                                        
                                        Ic_new = [Ic; 1]; % add a new infected person
                                        T0c_new = [T0c; TG_tmp(k)]; % time that people are infected
                                        
                                        if a4 < e_transmit
                                            NC_new0 = 0;
                                            Nc_new = [Nc; Nc_new0];
                                            TGc_new0 = zeros(1,nT);
                                            TGc_new0(1:Nc_new0) = TG_tmp(k) + random(TG_distrib,1,Nc_new0);
                                            TGc_new = [TGc; TGc_new0];
                                            
                                        else
                                            
                                            Nc_new0 = random(infect_distrib_sym); % number of people that are infected
                                            Nc_new = [Nc; Nc_new0];
                                            TGc_new0 = zeros(1,nT);
                                            TGc_new0(1:Nc_new0) = TG_tmp(k) + random(TG_distrib,1,Nc_new0);
                                            TGc_new = [TGc; TGc_new0];
                                            
                                            % update fields
                                            Ic = Ic_new;
                                            T0c = T0c_new;
                                            Nc = Nc_new;
                                            TGc = TGc_new;
                                            
                                        end
                                        
                                    else % create new subclinical infection
                                        
                                        Is_new = [Is; 1]; % add a new infected person
                                        T0s_new = [T0s; TG_tmp(k)]; % time that people are infected
                                        
                                        if a4 < e_transmit
                                            Ns_new0 = 0;
                                            
                                            Ns_new = [Ns; Ns_new0];
                                            TGs_new0 = zeros(1,nT); % time
                                            TGs_new0(1:Ns_new0) = TG_tmp(k) + random(TG_distrib,1,Ns_new0);
                                            TGs_new = [TGs; TGs_new0];
                                            
                                            % update fields
                                            Is = Is_new;
                                            T0s = T0s_new;
                                            Ns = Ns_new;
                                            TGs = TGs_new;
                                            
                                        else
                                            
                                            Ns_new0 = random(infect_distrib_asym); % number of people that are infected
                                            Ns_new = [Ns; Ns_new0];
                                            TGs_new0 = zeros(1,nT); % time
                                            TGs_new0(1:Ns_new0) = TG_tmp(k) + random(TG_distrib,1,Ns_new0);
                                            TGs_new = [TGs; TGs_new0];
                                            
                                            % update fields
                                            Is = Is_new;
                                            T0s = T0s_new;
                                            Ns = Ns_new;
                                            TGs = TGs_new;
                                            
                                        end
                                        
                                        
                                    end
                                end
                            end
                            
                            
                            
                        end
                    end
                end
                
            end
            
            
            num_infect = num_infect_c+num_infect_s; % total number of infected people (clinical and subclinical)
            
            % save time series of infected people for all realizations
            clinical_infections(p,:) = num_infect_c;
            subclinical_infections(p,:) = num_infect_s;
            total_infections(p,:) = num_infect;
            
            NUMBER_VACCINATED(p,:) = NUM_VAC;
            NUMBER_UNVACCINATED(p,:) = NUM_UNV;
            
            NUMBER_VAC_HOSP(p,:) = sum(VAC_HOS);
            NUMBER_UNVAC_HOSP(p,:) = sum(UNV_HOS);
            
        end
        
        %% Save data %%
        
        filename = strcat('output_vacRate=',num2str(prob_vac*100),'%_C=',num2str(C),'_N=',num2str(num_realizations),'.mat'); % file name
        params = [R, prob_asym,prob_vac,e_infect,e_transmit,C,tmax]; % save parameters into a vector
        cd OUTPUT_HOSP
        save(filename,'total_infections','clinical_infections','subclinical_infections',...
            't','params','infect_distrib_sym','infect_distrib_asym','TG_distrib',...
            'NUMBER_VACCINATED','NUMBER_UNVACCINATED',...
            'NUMBER_VAC_HOSP','NUMBER_UNVAC_HOSP');
        cd ..
        
    end
end





