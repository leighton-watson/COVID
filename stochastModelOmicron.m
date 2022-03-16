%% STOCHASTIC MODEL OMICRON %%
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
% This version tracks the number of unvaccinated, vaccinated, and boosted
% infections. This version does not model hospitalization rates. 

% Parameters are taken from COVID-19 Modeling Aotearoa
% A preliminary assessment of the potential impact of the Omicron variant
% of SARS-CoV-2 in Aotearoa New Zealand
% Vattiato G, Maclaren O, Lustig A, Binny RN, Hendy SC, Plank MJ (2022)

clear all; clc;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
cmap = get(gca,'ColorOrder');

%% Initialization %%

% Reproduction number
R = 2.6; % reproduction number of omicron variant in the presence of red traffic light setting
Rdist = makedist('Poisson','lambda',R);

% Generation time
Gdist = makedist('norm',3.3,1.3);

% Vaccine effectiveness
% Assume two doses are 15+ weeks after injection
vac2_infect = 0.14;
vac2_transmit = 0.03;
% Assume two doses are 2-5+ weeks after injection
vac3_infect = 0.58;
vac3_transmit = 0.26;

% Portion of population vaccinated
dose0 = 0.1; % unvaccinated
dose2 = 0.4; % vaccinated
dose3 = 0.5; % boosted

% Simulation time
tmax = 31; % number of days to run simulation for
t = 0:1:tmax; % time vector
nt = length(t);
num_realizations = 100000; % number of simulations to run

%% Run Simulations %%

tic
for p = 1:num_realizations
    
    % Display number of realizations
    disp(strcat('Realization Number=',num2str(p),'/',num2str(num_realizations)));
    
    % counters for number of unvaccinated/vaccinated/boosted infected by
    % unvaccinated/vaccinated/boosted
    
    B2B = 0; % boosted infects boosted
    B2V = 0; % boosted infects vaccinated
    B2U = 0; % boosted infects unvaccinated
    
    V2B = 0; % vaccinated infected boosted
    V2V = 0; % vaccinated infects vaccinated
    V2U = 0; % vaccinated infects unvaccinated
    
    U2B = 0; % unvaccinated infects boosted
    U2V = 0; % unvaccinated infects vaccinated
    U2U = 0; % unvaccinated infects unvaccinated
    
    % number of boosted, vaccinated, and unvaccinated infected individuals
    NUM_BOS = zeros(1,nt+20);
    NUM_VAC = zeros(1,nt+20);
    NUM_UNV = zeros(1,nt+20);
    
    % initialize simulation infected individuals
    NUM_UNV(1) = 1;
    NUM_VAC(1) = 0;
    NUM_BOS(1) = 0;
    
    for i = 1:length(t) % iterate over time steps
        
        %%% iterate over unvaccinated infections %%%
        for j = 1:NUM_UNV(i)
            N = random(Rdist); % number of new infections
            TG = abs(random(Gdist,1,N)); % generation time of new infections
            
            for k = 1:N % iterate over number of new infections
                
                a1 = rand(1); % is new infection boosted, vaccinated, or unvaccinated
                a2 = rand(1); % is vaccinated individual is infected or not
                
                if a1 < dose3 % potential new case is boosted
                    if a2 > vac3_infect % boosted case is infected
                        tidx = t(i) + ceil(TG(k)); % time of infection
                        U2B = U2B + 1; % unvaccinated case infects boosted
                        NUM_BOS(tidx) = NUM_BOS(tidx) + 1; % add new boosted infection
                    end
                elseif a1 < dose3+dose2 % potential new case is vaccianted
                    if a2 > vac2_infect
                        tidx = t(i) + ceil(TG(k)); % time of infection
                        U2V = U2V + 1; % unvaccinated case infects vaccinated
                        NUM_VAC(tidx) = NUM_VAC(tidx) + 1; % add new vaccinated infection
                    end
                else % new cases is unvaccinated
                    tidx = t(i) + ceil(TG(k)); % time of infection
                    U2U = U2U + 1; % unvaccinated case infects unvaccinated
                    NUM_UNV(tidx) = NUM_UNV(tidx) + 1; % add new unvaccinated infection
                end
            end
        end
        
        %%% iterate over vaccinated infections
        for j = 1:NUM_VAC(i)
            
            a3 = rand(1); % vaccine effectiveness against transmission
            
            if a3 > vac2_transmit % virus is transmitted
                
                N = random(Rdist); % number of new infections
                TG = abs(random(Gdist,1,N)); % generation time of new infections
                
                for k = 1:N % iterate over number of new infections
                    
                    a1 = rand(1); % is new infection boosted, vaccinated, or unvaccinated
                    a2 = rand(1); % is vaccinated individual is infected or not
                    
                    if a1 < dose3 % potential new case is boosted
                        if a2 > vac3_infect % boosted case is infected
                            tidx = t(i) + ceil(TG(k)); % time of infection
                            V2B = V2B + 1; % unvaccinated case infects boosted
                            NUM_BOS(tidx) = NUM_BOS(tidx) + 1; % add new boosted infection
                        end
                    elseif a1 < dose3+dose2 % potential new case is vaccianted
                        if a2 > vac2_infect
                            tidx = t(i) + ceil(TG(k)); % time of infection
                            V2V = V2V + 1; % unvaccinated case infects vaccinated
                            NUM_VAC(tidx) = NUM_VAC(tidx) + 1; % add new vaccinated infection
                        end
                    else % new cases is unvaccinated
                        tidx = t(i) + ceil(TG(k)); % time of infection
                        V2U = V2U + 1; % unvaccinated case infects unvaccinated
                        NUM_UNV(tidx) = NUM_UNV(tidx) + 1; % add new unvaccinated infection
                    end
                    
                    
                end
            end
        end
        
        %%% iterate over boosted infections
        for j = 1:NUM_BOS(i)
            
            a3 = rand(1); % vaccine effectiveness against transmission
            
            if a3 > vac3_transmit % virus is transmitted
                
                N = random(Rdist); % number of new infections
                TG = abs(random(Gdist,1,N)); % generation time of new infections
                
                for k = 1:N % iterate over number of new infections
                    
                    a1 = rand(1); % is new infection boosted, vaccinated, or unvaccinated
                    a2 = rand(1); % is vaccinated individual is infected or not
                    
                    if a1 < dose3 % potential new case is boosted
                        if a2 > vac3_infect % boosted case is infected
                            tidx = t(i) + ceil(TG(k)); % time of infection
                            B2B = B2B + 1; % unvaccinated case infects boosted
                            NUM_BOS(tidx) = NUM_BOS(tidx) + 1; % add new boosted infection
                        end
                    elseif a1 < dose3+dose2 % potential new case is vaccianted
                        if a2 > vac2_infect
                            tidx = t(i) + ceil(TG(k)); % time of infection
                            B2V = B2V + 1; % unvaccinated case infects vaccinated
                            NUM_VAC(tidx) = NUM_VAC(tidx) + 1; % add new vaccinated infection
                        end
                    else % new cases is unvaccinated
                        tidx = t(i) + ceil(TG(k)); % time of infection
                        B2U = B2U + 1; % unvaccinated case infects unvaccinated
                        NUM_UNV(tidx) = NUM_UNV(tidx) + 1; % add new unvaccinated infection
                    end
                    
                    
                end
            end
        end
        
    end
    
    % save outputs across realizations
    NUM_UNV_SAVE(p,:) = NUM_UNV;
    NUM_VAC_SAVE(p,:) = NUM_VAC;
    NUM_BOS_SAVE(p,:) = NUM_BOS;
        
    U2U_SAVE(p) = U2U;
    U2V_SAVE(p) = U2V;
    U2B_SAVE(p) = U2B;
    
    V2U_SAVE(p) = V2U;
    V2V_SAVE(p) = V2V;
    V2B_SAVE(p) = V2B;
    
    B2U_SAVE(p) = B2U;
    B2V_SAVE(p) = B2V;
    B2B_SAVE(p) = B2B;
    
    
    
end
toc

% save outputs

filename = "OMICRON_SEED_" + "_dose0=" + num2str(dose0) + "_dose2=" + num2str(dose2) + ...
    "_dose3=" + num2str(dose3) + "_unv=" + num2str(NUM_UNV(1)) + ...
    "_vac=" + num2str(NUM_VAC(1)) + "_bos=" + num2str(NUM_BOS(1)) + ".mat";
save(filename);





