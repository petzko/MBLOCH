clear; clc;
barrier =3 ;
well = 10;

bias = 2.0;



structure = [barrier well barrier well barrier well barrier well barrier ];
[N,M] = size(structure);
L = sum(structure);
dx = 0.01;
barrier_N = round(barrier/dx);
well_N = round(well/dx);
V = [];

for i = 1:M
   if mod(i,2) == 1
       V = [V ;ones(barrier_N,1)];
   else
        V = [V ;zeros(well_N,1)];
   end
end

% plot(V)
length(V);
voltages = linspace(0,bias,length(V))';
x = linspace(0,L,length(V));
biasedV = V - voltages;
plot(x,biasedV,'LineWidth',3.0)
axis([-10 L+10 (min(biasedV)-0.3) (max(biasedV)+0.3)]);


