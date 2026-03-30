%function [r]=solve_n_springs(ics,times)

close 
clear all

ics=0:0.5:5.0;

times=linspace(0,10,100)';

start_time = times(1);

num_cells = length(ics);

eqns=zeros(num_cells);

k=18.2816648;

eqns(1,end) = 1;
eqns(2,end-1:end) = [1/k,1];

for row=3:num_cells
    eqns(row,:) = [1/k*eqns(row-1,2:end),0] + 2*eqns(row-1,:)  - eqns(row-2,:);
end

p = [1/k*eqns(end,:),0] + [0, eqns(end,:)] - [0, eqns(end-1,:)];
p_sym = poly2sym(p);

L = roots(p)';

A=zeros(num_cells);
L_pow = zeros(num_cells);

for row = 1:num_cells
    L_pow(num_cells+1-row,:) = L.^(row-1);
    A(row,:) = eqns(row,:)*L_pow;
end

b=-(1:num_cells) +ics;

coeficients = A\b';

r=zeros(num_cells,length(times));

for row = 1:num_cells
    r(row,:) = ((eqns(row,:)*(L_pow)).*coeficients')*exp((times-start_time)*L)' + row; 
end

sep = ((eqns(end,:)*(L_pow)).*coeficients')*exp((times-start_time)*L)' +11 -...
      (((eqns(1,:)*(L_pow)).*coeficients')*exp((times-start_time)*L)' +1) ...
      

plot(times,r');

figure 

plot(times,sep, 'w');%r(end,:)-r(1,:))
hold on 
plot([0,10],[9,9],'w:')   % vs. 'k:'
plot([1,1],[4.5,10.5],'w:')

sep_1 = ((eqns(end,:)*(L_pow)).*coeficients')*exp(L)' +11 -...
      (((eqns(1,:)*(L_pow)).*coeficients')*exp(L)' +1)