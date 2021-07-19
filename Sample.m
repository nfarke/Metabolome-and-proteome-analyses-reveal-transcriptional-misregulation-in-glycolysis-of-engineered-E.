function [p,par]=Sample(ParSize,gg)

p={'rU','kcatL','kcatG','KmL','KmG','betaL','betaG','a1','a2','alphax','ind','m','eL','eG'};

%parameters average flux
if gg == 1 %called by create_SS_solutions_par -> for figure4b
    par_lb=[40.87, 16080/3.33, 1705/3.33, 0.01,  0.01 , 0.000238/3.33,       0.0017,    1, 1, 4086,1, 1, 0.0238,1]';
    par_ub=[40.87, 16080*3,    1705*3,      10,     10, 0.000238*3,       0.0017,       2, 2, 4086,1, 1, 0.0238,1]';
end

par=zeros(length(par_lb),ParSize);
%random sampling
for i=1:ParSize
    par(:,i)= 10.^((log10(par_lb)-log10(par_ub)).*rand(length(par_lb), 1)+log10(par_ub));
end

par(end) = 0; %eG
par(end-3) = 0; %ind
end