clear
%%%%%%define derivative%%%%%%%%% 
dudt = @(u,v,a1,b) a1/(1 + v^b) - u; 
dvdt = @(u,v,a2,g) a2/(1 + u^g) - v; 
%%%%%%%%%%%Define Initial Condition%%%%%%%%%%%%%
y0 = [0 0];
u_0 = y0(1); 
v_0 = y0(2); 
a1 = 2; 
a2 = 4; 
b = 2; 
g = 2; 

tspan = [0 12];

h = 0.1; 

t_0 = tspan(1); 
t_f = tspan(2); 
t1 = t_0:h:t_f;
nStep = t_f/h;

%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nStep
    u(i) = u_0 + h*dudt(u_0,v_0,a1,b);
    v(i) = v_0 + h*dvdt(u_0,v_0,a2,g);
    u_0 = u(i); 
    v_0 = v(i); 
end

u = [y0(1) u]; 
v = [y0(2) v]; 
figure(1)
subplot(1,2,1)
plot(t1,u);
hold on 
plot(t1,v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,y] = ode45(@(t,y) funToggleSwitch(t,y,a1,a2,b,g),tspan,y0);

subplot(1,2,2)
plot(t,y(:,1)); 
hold on 
plot(t,y(:,2));