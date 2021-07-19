function [] = Figure()
load PAR.mat
EnsembleSize= size(par_WT,2);

Pert = 1; 
NoSteps = Pert*400; %number of steps
varpert = {'ind'}; %define perturbation variable --> ind
varpert1 = sym(varpert);

p={'rU','kcatL','kcatG','KmL','KmG','betaL','betaG','a1','a2','alphax','ind','m','eL','eG'};
par1 = par_KO;
par2 = par_WT;
par3 = par_2xCra;

psym = sym(p);

syms rU kcatL kcatG KmL KmG betaL betaG a1 a2 alphax ind m eL eG

rL = kcatL * eL * m/(m + KmL);
rG = kcatG * eG * m/(m + KmG);

mue = rL/alphax;
%KO
TR_L1 = betaL;
TR_G1 = betaG * ind;
%base strain
TR_L2 = betaL * m^a1;
TR_G2 = betaG * ind;        
%2xCra
TR_L3 = betaL * m^a1;
TR_G3 = betaG * ind * m^a2;

dilL = eL*mue;
dilG = eG*mue;
dilm = m * mue;
%define rate vector r, stoichiometric matrix S and mass balance F = S*r
%KO
r1 = [rU;rL;rG;TR_L1;TR_G1;dilL;dilG;dilm];
S1 = [1,-1,-1,0,0,0,0,-1; 0, 0,0,1,0,-1,0,0; 0,0,0,0,1,0,-1,0];
F1 = S1*r1;

%base strain
r2 = [rU;rL;rG;TR_L2;TR_G2;dilL;dilG;dilm];
S2 = [1,-1,-1,0,0,0,0,-1; 0, 0,0,1,0,-1,0,0; 0,0,0,0,1,0,-1,0];
F2 = S2*r2;
%2xCra
r3 = [rU;rL;rG;TR_L3;TR_G3;dilL;dilG;dilm];
S3 = [1,-1,-1,0,0,0,0,-1; 0, 0,0,1,0,-1,0,0; 0,0,0,0,1,0,-1,0];
F3 = S3*r3;

%set up derivatives as matlab function for higher computation speed
DFDX1 = matlabFunction(jacobian(F1,[m, eL, eG]),'vars',{psym});
DFDP1 = matlabFunction(jacobian(F1,varpert1),'vars',{psym});    
%
DFDX2 = matlabFunction(jacobian(F2,[m, eL, eG]),'vars',{psym});
DFDP2 = matlabFunction(jacobian(F2,varpert1),'vars',{psym});    
%
DFDX3 = matlabFunction(jacobian(F3,[m, eL, eG]),'vars',{psym});
DFDP3 = matlabFunction(jacobian(F3,varpert1),'vars',{psym});

%initialize result matrices
indx1 = nan(EnsembleSize, NoSteps+1);           
indx2 = nan(EnsembleSize, NoSteps+1);           
indx3 = nan(EnsembleSize, NoSteps+1); 

%glycerol flux
rGx1 = nan(EnsembleSize, NoSteps+1);           
rGx2 = nan(EnsembleSize, NoSteps+1);           
rGx3 = nan(EnsembleSize, NoSteps+1);

%lower glycolysis flux
rLx1 = nan(EnsembleSize, NoSteps+1);           
rLx2 = nan(EnsembleSize, NoSteps+1);           
rLx3 = nan(EnsembleSize, NoSteps+1); 

%metabolite concentrations
c1 = nan(EnsembleSize, NoSteps+1);           
c2 = nan(EnsembleSize, NoSteps+1);           
c3 = nan(EnsembleSize, NoSteps+1); 

parfor k = 1:EnsembleSize
    disp(k)
    parx1 = par1(:,k);
    parx2 = par2(:,k);
    parx3 = par3(:,k);
  
    Varini = parx1(find(strcmp(p,varpert))); %initial state of pert variable
    Varmax = Pert; %max value of the pert variable
    InitialState1 = [parx1(end-2) parx1(end-1) parx1(end)];
    InitialState2 = [parx2(end-2) parx2(end-1) parx2(end)];
    InitialState3 = [parx3(end-2) parx3(end-1) parx3(end)];
    
    options = odeset('Events',@event_function,'RelTol', 1e-5, 'AbsTol', 1e-5);
    TimeIn=clock;
    [t1, conc1] = ode15s(@continuation_method,0:1/NoSteps:1,InitialState1,options,psym,parx1,Varini,Varmax,DFDX1,DFDP1,TimeIn,varpert);
    TimeIn=clock;
    [t2, conc2] = ode15s(@continuation_method,0:1/NoSteps:1,InitialState2,options,psym,parx2,Varini,Varmax,DFDX2,DFDP2,TimeIn,varpert);
    TimeIn=clock;
    [t3, conc3] = ode15s(@continuation_method,0:1/NoSteps:1,InitialState3,options,psym,parx3,Varini,Varmax,DFDX3,DFDP3,TimeIn,varpert);
                 
    %check steady state at tend
    [indx1n,rL1,rG1,concx1] = getvars(t1,conc1,Varini,Varmax,NoSteps,F1,psym,parx1,rL,rG,DFDX1,varpert);
    indx1(k,:) = indx1n;
    rGx1(k,:) = rG1;
    rLx1(k,:) = rL1;
    c1(k,:) = concx1;
    
    
    [indx2n,rL2,rG2,concx2] = getvars(t2,conc2,Varini,Varmax,NoSteps,F2,psym,parx2,rL,rG,DFDX2,varpert);
    indx2(k,:) = indx2n;  
    rGx2(k,:) = rG2;
    rLx2(k,:) = rL2;
    c2(k,:) = concx2;


    [indx3n,rL3,rG3,concx3] = getvars(t3,conc3,Varini,Varmax,NoSteps,F3,psym,parx3,rL,rG,DFDX3,varpert);
    indx3(k,:) = indx3n;
    rGx3(k,:) = rG3;
    rLx3(k,:) = rL3;
    c3(k,:) = concx3;
end


%Robustness Analysis plot
figure(1)
subplot(1,2,2)
rob1 = sum(~isnan(indx1))/EnsembleSize;
rob2 = sum(~isnan(indx2))/EnsembleSize;
rob3 = sum(~isnan(indx3))/EnsembleSize;
plot(rob3)
plot(0:1:NoSteps,[rob1;rob2;rob3])
legend('KO','base strain','2xCra')
xticks([0 NoSteps/2 NoSteps])
xticklabels([0 50 100])
yticks([0 0.5 1])
xlabel('Induction, %');
ylabel('Fraction of stable models, %');
ylim([0 1])

%Cumulative sum plot
%Resuls for KO strain model
for k = 1:EnsembleSize
    inx1 = min(find(isnan(rGx1(k,:))));
    if isempty(inx1)
       fluxx1(k,1) = rGx1(k,NoSteps+1);
    else
       fluxx1(k,1) = rGx1(k,inx1-1);
    end
end

%Results for base strain model
for k = 1:EnsembleSize
    inx2 = min(find(isnan(rGx2(k,:))));
    if isempty(inx2)
       fluxx2(k,1) = rGx2(k,NoSteps+1);
    else
       fluxx2(k,1) = rGx2(k,inx2-1);
    end
end

%Results for 2x Cra model
for k = 1:EnsembleSize
    inx3 = min(find(isnan(rGx3(k,:))));
    if isempty(inx3)
       fluxx3(k,1) = rGx3(k,NoSteps+1);
    else
       fluxx3(k,1) = rGx3(k,inx3-1);
    end
end


flux_steps = 1:0.01:40.87;

for k = 1:length(flux_steps)
    cumsumx1(k,1) = sum(fluxx1 >= flux_steps(k));   
    cumsumx2(k,1) = sum(fluxx2 >= flux_steps(k));    
    cumsumx3(k,1) = sum(fluxx3 >= flux_steps(k));
end

%
subplot(1,2,1)
xxx = flux_steps*0.002*60;
plot(xxx*2,cumsumx1/EnsembleSize)
hold on
plot(xxx*2,cumsumx2/EnsembleSize)
hold on
plot(xxx*2,cumsumx3/EnsembleSize)
xlim([0 10])
xlabel('glycerol flux')
ylabel('fraction of models, %')
legend('KO model','Base Strain','2x Cra')
yticks([0 0.5 1])
end
function dx= continuation_method(t,x,psym,par,Varini,Varfinal,DFDX,DFDP,~,varpert)
    par(end-2:end) = x;
    U = Varini + t.*(Varfinal-Varini);
    par(find(ismember(psym,varpert))) = U;
    XJac = DFDX(par');
    DFDP1 = DFDP(par');
    dx = -XJac\DFDP1*(Varfinal-Varini);
end

function [value,isterminal,direction] = event_function(t,x,psym,par,Varini,Varfinal,DFDX,~,TimeIn,varpert)
    par(end-2:end)=x;
    U = Varini + t.*(Varfinal-Varini);
    par(find(ismember(psym,varpert))) = U;
    XJac = DFDX(par');
    value = max(real(eig(XJac))); % when value = 0, an event is triggered
    if max(real(eig(XJac)))> -1E-05 || any(par<0) %stability criteria
        value = 0;
    end
    isterminal = 1; % terminate after the first event
    direction = 0;  % get all the zeros
end
