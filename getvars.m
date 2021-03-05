function [indx,rLx,rGx,concx] = getvars(t,conc,Varini,Varmax,NoSteps,F1,psym,parx1,rL,rG,DFDX,varpert)

indx = nan(1, NoSteps+1);           

rGx = NaN(1,NoSteps+1);
rLx = NaN(1,NoSteps+1);
concx  = NaN(1,NoSteps+1);
concx(1,1:length(conc(:,1))) = conc(:,1);

for j = 1:length(t)
    U1 = Varini + t(j).*(Varmax-Varini);
    parx1(find(ismember(psym,varpert))) = U1;
    parx1(end-2:end) = conc(j,:);
    indx(1,j) = double(U1/Varmax)*100;
    rLx(1,j) = double(subs(rL,psym,parx1'));
    rGx(1,j) = double(subs(rG,psym,parx1'));
end

end

