function [p,par]=Sample(ParSize,gg)

p={'rU','kcatL','kcatG','KmL','KmG','betaL','betaG','a1','a2','alpha','ind','m','eL','eG'};

%parameters average flux
if gg == 1 %called by create_SS_solutions_par -> for figure4b
    par_lb=[66.67, 16080/3.33, 1705/3.33, 0.1,  0.1 , 0.00027/3.33, 0.0085,    1, 1, 5747,1, 1, 0.023,1]';
    par_ub=[66.67, 16080*3.33, 1705*3.33, 10, 10, 0.00027*3.33, 0.0085,       2, 2, 5747,1, 1, 0.023,1]';
end

par=zeros(length(par_lb),ParSize);
%random sampling
for i=1:ParSize
    par(:,i)= 10.^((log10(par_lb)-log10(par_ub)).*rand(length(par_lb), 1)+log10(par_ub));
end

par(end) = 0; %eG
par(end-3) = 0; %ind
end