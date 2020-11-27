%this file generates the .mat-file PAR.mat. It contains steady state
%solutions of three models (base strain model, knock-out model, 2xCra model)

ParSize = 1000; %adjust number of Parameter sets
tspan = 0:1:3000;
par_end = zeros(14,ParSize);
opts = odeset('RelTol',1e-07,'AbsTol',1E-7);

for i1 = 1:ParSize
    disp(i1)
    value = 0; %repeat when stability criteria are not met
    while value == 0
        [p,par] = Sample(1,1); %samples one parameter set from case gg = 1
        y0 = par(end-2:end); %initial conditions
        for af = 0:2 %loop through different regulations
            
            %calculate constraint parameters
            par(strcmp('kcatL',p)) = par(strcmp('rU',p))./(par(strcmp('eL',p))/(1 + par(strcmp('KmL',p))));
            par(strcmp('betaL',p)) = par(strcmp('eL',p)) *par(strcmp('rU',p))/5747;
      
            [~,y]  =  ode23s(@(t,c) odemodel(t,c,p,par,af),tspan,y0,opts);
            
            %Parameters;
            rU = par(find(strcmp(p,'rU')),1);
            kcatL = par(find(strcmp(p,'kcatL')),1);
            kcatG = par(find(strcmp(p,'kcatG')),1);
            KmL = par(find(strcmp(p,'KmL')),1);
            KmG = par(find(strcmp(p,'KmG')),1);
            betaL = par(find(strcmp(p,'betaL')),1);
            betaG = par(find(strcmp(p,'betaG')),1);
            a1 = par(find(strcmp(p,'a1')),1);
            a2 = par(find(strcmp(p,'a2')),1);
            alpha = par(find(strcmp(p,'alpha')),1);
            ind = par(find(strcmp(p,'ind')),1);

            %Variables at t(end)
            m  =  y(end,1);
            eL =  y(end,2);
            eG =  y(end,3);
            %Definition of the growth rate
            mue  =  (kcatL * eL * m/(m + KmL))/alpha;
            
            %Mass balances
            dmdt = rU - kcatL * eL * m/(m + KmL) - kcatG * eG * m/(m + KmG);
            
            if af == 0 %KO
            deLdt = betaL - eL*mue;
            deGdt = betaG * ind - eG*mue;
            
            elseif af == 1 %base strain model
            deLdt = betaL * m^a1  - eL*mue;
            deGdt = betaG * ind - eG*mue;    
            
            elseif af == 2 % 2x Cra model
            deLdt = betaL * m^a1 - eL*mue;
            deGdt = betaG * ind * m^a2 - eG*mue;
            end
       
            %check if solutions are in steady state
            F = [dmdt;deLdt;deGdt];
            if af == 0
                SS_KO = sum(abs(F));
            elseif af == 1
                SS_WT = sum(abs(F));
            elseif af == 2
                SS_2xCra = sum(abs(F));
            end
            
            %calculation of the Jacobian Matrix
            if af == 0
            J = [ (eG*kcatG*m)/(KmG + m)^2 - (eL*kcatL)/(KmL + m) - (eG*kcatG)/(KmG + m) + (eL*kcatL*m)/(KmL + m)^2,              -(kcatL*m)/(KmL + m),            -(kcatG*m)/(KmG + m)
                               (eL^2*kcatL*m)/(alpha*(KmL + m)^2) - (eL^2*kcatL)/(alpha*(KmL + m)), -(2*eL*kcatL*m)/(alpha*(KmL + m)),                               0
                             (eG*eL*kcatL*m)/(alpha*(KmL + m)^2) - (eG*eL*kcatL)/(alpha*(KmL + m)),   -(eG*kcatL*m)/(alpha*(KmL + m)), -(eL*kcatL*m)/(alpha*(KmL + m))];

            elseif af == 1
            J = [ (eG*kcatG*m)/(KmG + m)^2 - (eL*kcatL)/(KmL + m) - (eG*kcatG)/(KmG + m) + (eL*kcatL*m)/(KmL + m)^2,              -(kcatL*m)/(KmL + m),            -(kcatG*m)/(KmG + m)
                    a1*betaL*m^(a1 - 1) - (eL^2*kcatL)/(alpha*(KmL + m)) + (eL^2*kcatL*m)/(alpha*(KmL + m)^2), -(2*eL*kcatL*m)/(alpha*(KmL + m)),                               0
                             (eG*eL*kcatL*m)/(alpha*(KmL + m)^2) - (eG*eL*kcatL)/(alpha*(KmL + m)),   -(eG*kcatL*m)/(alpha*(KmL + m)), -(eL*kcatL*m)/(alpha*(KmL + m))];

            elseif af == 2
            J = [ (eG*kcatG*m)/(KmG + m)^2 - (eL*kcatL)/(KmL + m) - (eG*kcatG)/(KmG + m) + (eL*kcatL*m)/(KmL + m)^2,              -(kcatL*m)/(KmL + m),            -(kcatG*m)/(KmG + m)
                    a1*betaL*m^(a1 - 1) - (eL^2*kcatL)/(alpha*(KmL + m)) + (eL^2*kcatL*m)/(alpha*(KmL + m)^2), -(2*eL*kcatL*m)/(alpha*(KmL + m)),                               0
                    a2*betaG*ind*m^(a2 - 1) - (eG*eL*kcatL)/(alpha*(KmL + m)) + (eG*eL*kcatL*m)/(alpha*(KmL + m)^2),   -(eG*kcatL*m)/(alpha*(KmL + m)), -(eL*kcatL*m)/(alpha*(KmL + m))];
   
            end
            
            %calculate real parts of the Jacobian
            if af == 0
                val = real(eig(J));
                ew_KO = max(val(1:2));
                sol_KO = [m;eL;eG];
            elseif af == 1
                val = real(eig(J));
                ew_WT = max(val(1:2));
                sol_WT = [m;eL;eG];
            elseif af == 2
                val = real(eig(J));
                ew_2xCra = max(val(1:2));
                sol_2xCra = [m;eL;eG];
            end
        end %af
        
        %check for steady state and eigenvalues (stability)
        if  SS_KO < 1E-08 && SS_WT < 1E-08 && SS_2xCra < 1E-08...
            && ew_KO < -1E-05 && ew_WT < -1E-05 && ew_2xCra < -1E-05
            
            value = 1;
            par_end(:,i1) = par;    
            sol1(:,i1) = sol_KO;
            sol2(:,i1) = sol_2xCra;
            sol3(:,i1)  = sol_WT;       
        end
   end
end

%steady state solution vectors
par = par_end;
%KO
par(end-2:end,:) = sol1;
par_KO = par;
%WT
par(end-2:end,:) = sol3;
par_WT = par;
%2xCra
par(end-2:end,:) = sol2;
par_2xCra = par;

save('PAR.mat','par_WT','par_KO','par_2xCra')



function dcdt  =  odemodel(~,c,p,par,af,~)

par(end-2:end,1) = c;

%get kinetic parameters
rU = par(find(strcmp(p,'rU')),1);
kcatL = par(find(strcmp(p,'kcatL')),1);
kcatG = par(find(strcmp(p,'kcatG')),1);
KmL = par(find(strcmp(p,'KmL')),1);
KmG = par(find(strcmp(p,'KmG')),1);
betaL = par(find(strcmp(p,'betaL')),1);
betaG = par(find(strcmp(p,'betaG')),1);
a1 = par(find(strcmp(p,'a1')),1);
a2 = par(find(strcmp(p,'a2')),1);
alpha = par(find(strcmp(p,'alpha')),1);
ind = par(find(strcmp(p,'ind')),1);

m = par(find(strcmp(p,'m')),1);
eL = par(find(strcmp(p,'eL')),1);
eG = par(find(strcmp(p,'eG')),1);

rL = kcatL * eL * m/(m + KmL);
rG = kcatG * eG * m/(m + KmG);

dcdt(1,1) = rU - rL - rG;


mue = rL/alpha;

if af == 0 %KO
   dcdt(2,1) = betaL - eL*mue;
   dcdt(3,1) = betaG * ind - eG*mue; 
elseif af == 1 %base strain
   dcdt(2,1) = betaL * m^a1  - eL*mue;
   dcdt(3,1) = betaG * ind - eG*mue;

elseif af == 2 %2x Cra
   dcdt(2,1) = betaL * m^a1 - eL*mue;
   dcdt(3,1) = betaG * ind * m^a2 - eG*mue;
end

end

