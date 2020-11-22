% Excitatory neurons Inhibitory neurons
%Ne = 800; 
%Ni = 200;
division = 0.2;
ms = 2000;

%re = rand(Ne,1); 
%ri = rand(Ni,1);
%a = [0.02*ones(Ne,1); 0.02+0.08*ri]; %time scale of rec var, smaller, slower recovery
%b = [0.2*ones(Ne,1); 0.25-0.05*ri]; %sensitivity of u
%c = [-65+15*re.^2; -65*ones(Ni,1)]; %after spike reset value
%d = [8-6*re.^2; 2*ones(Ni,1)]; %after spike reset of u

% Create neurons
neuron = rand(1000, 1);
a = zeros(1000);
b = zeros(1000);
c = zeros(1000);
d = zeros(1000);


% Assign parameters to neurons
for i = 1:1000
    if neuron(i) > division
        a(i) = 0.02; % time scale of rec var, smaller, slower recovery
        b(i) = 0.2; % sensitivity of u
        c(i) = -65 + 15 * rand(1)^2; % after spike reset value
        d(i) = 8 - 6 * rand(1)^2; % after spike reset of u
    else
        a(i) = 0.02 + 0.08 * rand(1); % time scale of rec var, smaller, slower recovery
        b(i) = 0.25 - 0.05 * rand(1); % sensitivity of u
        c(i) = -65; % after spike reset value
        d(i) = 2; % after spike reset of u
    end    
end

%S = [0.5*rand(Ne+Ni,Ne), -rand(Ne+Ni,Ni)];

% Add neuron connections
S = zeros(1000, 1000);
for i = 1:1000
    for j = 1:1000
       if abs(i - j) <= 25 || abs(i - j) > 975 
           S(i, j) = 1;
       else
           S(i, j) = 0;
       end
    end
end
S = S - eye(size(S));

v = -65 * ones(1000,1); % Initial values of v
u = b .* v; % Initial values of u
firings = []; % spike timings

for t=1:ms % simulation of 1000 ms
    I = rand(1000, 1); 
    I(neuron > division) = 5.5 .* I(neuron > division);
    I(neuron < division) = 2 .* I(neuron < division);
    
    % I = [5*randn(Ne,1); 2*randn(Ni,1)]; % thalamic input
    
    fired = find(v>=30); % indices of spikes
    
    firings = [firings; t+0*fired,fired];
    v(fired) = c(fired);
    u(fired) = u(fired) + d(fired);
    I = I + sum(S(:,fired),2);
    v = v + 0.5*(0.04*v.^2 + 5*v + 140 - u + I); % step 0.5 ms
    v = v + 0.5*(0.04*v.^2 + 5*v + 140-u + I); % for numerical
    u = u + a.*(b.*v - u); % stability

end

%h = animatedline('Marker','.');
x = firings(:,1);
y = firings(:,2);
plot(x,y,'.')

%z = (zeros(length(x),1));

% % for k = 1:length(x)
% %     xlim([0 1000])
% %     ylim([0 1000])
% %     scatter(x(k),y(k),'.');
% %     %refreshdata(hplot,'caller')
% %     %drawnow
% %     pause(.008)
% % end
% for k = 1:length(x)
%     figure;
%     axes('position',[0 1000 0 1000]);
%     plot1 = scatter(x(k),y(k),'.');
%     xlim([0 1000]);
%     ylim([0 1000]);
%     pause(.008)
%     %set(gca,'Color','none');
%     %set(gca,'CLim',[0, 1E-4]);
% end
% %drawnow
% %plot(firings(:,1),firings(:,2))
