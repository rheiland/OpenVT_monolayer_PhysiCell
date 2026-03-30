function relaxation_exact

% find the value of k/eta such that w(1) = 9
k_eta = fzero(@(k_eta) exact_width(k_eta, 1) - 9, 18); % yields 18.2816647214633

t = linspace(0,5,51)'; % temporal resolution
w = exact_width(k_eta, t); % in cell diameters

close all
plot(t, w);
hold on 
xline(1, 'k:')
yline(9, 'k:')
ylim([5 11])

writetable(table(t, w), 'relaxation_exact.csv');

end


function w = exact_width(k_eta, t)

N = 11; % number of cells
x = linspace(0,(N-1)/2,N)'; % initial cell positions in cell diameters

eqns = zeros(N);
eqns(1,end) = 1;
eqns(2,end-1:end) = [1/k_eta, 1];
for i = 3:N
    eqns(i,:) = [eqns(i-1,2:end)/k_eta, 0] + 2*eqns(i-1,:) - eqns(i-2,:);
end

P = [eqns(end,:)/k_eta, 0] + [0, eqns(end,:)] - [0, eqns(end-1,:)];
lambda = roots(P)';

A = zeros(N);
lambda_pow = zeros(N);
for i = 1:N
    lambda_pow(N+1-i,:) = lambda.^(i-1);
    A(i,:) = eqns(i,:) * lambda_pow;
end

b = x - (1:N)';

coeff = A \ b;

w = N - 1 + (eqns(end,:)-eqns(1,:)) * lambda_pow .* coeff' * exp(t*lambda)';
w = w';

end