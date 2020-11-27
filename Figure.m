load PAR.mat
EnsembleSize= size(par_WT,2);

Pert = 1; 
NoSteps = Pert*400; %number of steps
varpert = {'ind'}; %define perturbation variable, here it is ind
varpert1 = sym(varpert);

p={'rU','kcatL','kcatG','KmL','KmG','betaL','betaG','a1','a2','alpha','ind','m','eL','eG'};
par1 = par_KO;
par2 = par_WT;
par3 = par_2xCra;

psym = sym(p);

syms rU kcatL kcatG KmL KmG betaL betaG a1 a2 alpha ind m eL eG

rL = kcatL * eL * m/(m + KmL);
rG = kcatG * eG * m/(m + KmG);
mue = rL/alpha;
%KO
TR_L1 = betaL;
TR_G1 = betaG * ind;
%base strain
TR_L2 = betaL * m^a1;
TR_G2 = betaG * ind;        
%2xCra
TR_L3 = betaL * m^a1;
TR_G3 = betaG * ind * m^a2;

%%
dilL = eL*mue;
dilG = eG*mue;

%define rate vector r, stoichiometric matrix S and mass balance F = S*r
%KO
r1 = [rU;rL;rG;TR_L1;TR_G1;dilL;dilG];
S1 = [1,-1,-1,0,0,0,0; 0, 0,0,1,0,-1,0; 0,0,0,0,1,0,-1];
F1 = S1*r1;

%base strain
r2 = [rU;rL;rG;TR_L2;TR_G2;dilL;dilG];
S2 = [1,-1,-1,0,0,0,0; 0, 0,0,1,0,-1,0; 0,0,0,0,1,0,-1];
F2 = S2*r2;
%2xCra
r3 = [rU;rL;rG;TR_L3;TR_G3;dilL;dilG];
S3 = [1,-1,-1,0,0,0,0; 0, 0,0,1,0,-1,0; 0,0,0,0,1,0,-1];
F3 = S3*r3;

%set up derivatives as matlab function for enhanced speed
DFDX1 = matlabFunction(jacobian(F1,[m, eL, eG]),'vars',{psym});
DFDP1 = matlabFunction(jacobian(F1,varpert1),'vars',{psym});    
%
DFDX2 = matlabFunction(jacobian(F2,[m, eL, eG]),'vars',{psym});
DFDP2 = matlabFunction(jacobian(F2,varpert1),'vars',{psym});    
%
DFDX3 = matlabFunction(jacobian(F3,[m, eL, eG]),'vars',{psym});
DFDP3 = matlabFunction(jacobian(F3,varpert1),'vars',{psym});

indx1 = nan(EnsembleSize, NoSteps+1);           
indx2 = nan(EnsembleSize, NoSteps+1);           
indx3 = nan(EnsembleSize, NoSteps+1); 

rGx1 = nan(EnsembleSize, NoSteps+1);           
rGx2 = nan(EnsembleSize, NoSteps+1);           
rGx3 = nan(EnsembleSize, NoSteps+1); 


for k = 1:EnsembleSize
    disp(k)
    parx1 = par1(:,k);
    parx2 = par2(:,k);
    parx3 = par3(:,k);
  
    Varini = parx1(find(strcmp(p,varpert))); %initial state of pert variable
    Varmax = Pert; %max value of the pert variable
    InitialState1 = [parx1(end-2) parx1(end-1) parx1(end)];
    InitialState2 = [parx2(end-2) parx2(end-1) parx2(end)];
    InitialState3 = [parx3(end-2) parx3(end-1) parx3(end)];
    
    options = odeset('Events',@event_function,'RelTol', 1e-6, 'AbsTol', 1e-6);
    TimeIn=clock;
    [t1, conc1] = ode15s(@dxFunc,0:1/NoSteps:1,InitialState1,options,psym,parx1,Varini,Varmax,DFDX1,DFDP1,TimeIn,varpert);
    TimeIn=clock;
    [t2, conc2] = ode15s(@dxFunc,0:1/NoSteps:1,InitialState2,options,psym,parx2,Varini,Varmax,DFDX2,DFDP2,TimeIn,varpert);
    TimeIn=clock;
    [t3, conc3] = ode15s(@dxFunc,0:1/NoSteps:1,InitialState3,options,psym,parx3,Varini,Varmax,DFDX3,DFDP3,TimeIn,varpert);
                
    %check steady state at tend
    [indx1n] = getvars(t1,conc1,Varini,Varmax,NoSteps,F1,psym,parx1,rL,rG,DFDX1,varpert);
    indx1(k,:) = indx1n;
    
    [indx2n] = getvars(t2,conc2,Varini,Varmax,NoSteps,F2,psym,parx2,rL,rG,DFDX2,varpert);
    indx2(k,:) = indx2n;  

    [indx3n] = getvars(t3,conc3,Varini,Varmax,NoSteps,F3,psym,parx3,rL,rG,DFDX3,varpert);
    indx3(k,:) = indx3n;

end

figure(1)
rob1 = sum(~isnan(indx1))/EnsembleSize;
rob2 = sum(~isnan(indx2))/EnsembleSize;
rob3 = sum(~isnan(indx3))/EnsembleSize;
plot(rob3)
plot(0:1:NoSteps,[rob1;rob2;rob3])
legend('KO','base strain','2xCra')
xticks([0 NoSteps/2 NoSteps])
xticklabels([0 50 100])
xlabel('Induction, %');
ylabel('Robustness, %');
ylim([0 1])


function dx= dxFunc(t,x,psym,par,Varini,Varfinal,DFDX,DFDP,~,varpert)
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
    if max(real(eig(XJac)))> -1E-05 || any(par<0)
        value = 0;
    end
    isterminal = 1; % terminate after the first event
    direction = 0;  % get all the zeros
end
