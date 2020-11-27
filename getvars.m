function [indx] = getvars(t,conc,Varini,Varmax,NoSteps,F1,psym,parx1,rL,rG,DFDX,varpert)

indx = nan(1, NoSteps+1);           


for j = 1:length(t)
    U1 = Varini + t(j).*(Varmax-Varini);
    parx1(find(ismember(psym,varpert))) = U1;
    parx1(end-2:end) = conc(j,:);
    indx(1,j) = double(U1/Varmax)*100;
end

end

