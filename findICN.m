function ICN = findICN(x, y, phi, le)

a=le/3; %length OA
b=le; %length OB
h=le*2/3; %length OH
%le=1*150; %length of links
n=100; %length of end-effector


alpha   = y + n*sin(phi);
epsilon = y - n*sin(phi);
beta    = x + n*cos(phi) - a;
zeta    = x - n*cos(phi) - b;
rho     = sqrt(x*x + y*y);
gamma   = (rho*rho + n*n + 2*a*x - 3*a*a + 2*a*n*cos(phi) + 2*n*(x*cos(phi) + y*sin(phi))) / (2*le);
eta     = (rho*rho + n*n - 2*b*x + b*b + 2*b*n*cos(phi) - 2*n*(x*cos(phi) + y*sin(phi))) / (2*le);   
delta   = sqrt( alpha*alpha + beta*beta);
iota    = sqrt( epsilon*epsilon + zeta*zeta);

theta1  = 2*n*pi + acos(gamma/delta) + atan(alpha/beta);
theta2  = 2*n*pi + acos(eta/iota) + atan(epsilon/zeta);
rho     = sqrt(x*x + y*y);

O = [0 ; 0];
A = [-a ; 0];
B = [-b ; 0];
H = [h ; 0];
G = [h + le*cos(theta2) ; -le*sin(theta2)];
C = [-b + le*cos(theta1) ; -le*sin(theta1)];
E = [x ; y];
D = [x + n*cos(phi) ; y + n*sin(phi)];
F = [x - n*cos(phi) ; y - n*sin(phi)];

bCap  = (B - O)/le;
aCap  = (A - O)/le;
hCap  = (H - O)/le;

u1Cap = (C - B)/le;
v1Cap = (D - C)/le;
w1Cap = (E - D)/le;

u2Cap = (G - H)/le;
v2Cap = (F - G)/le;
w2Cap = (E - F)/le;

u3Cap = (E - A)/rho;

Em = [0 -1; 1 0];

M = [transpose(v1Cap) -n*transpose(v1Cap)*Em*w1Cap ; transpose(v2Cap) -n*transpose(v2Cap)*Em*w2Cap ; transpose(u3Cap) 0];
N = [le*transpose(v1Cap)*Em*u1Cap 0 0 ; 0 le*transpose(v2Cap)*Em*u2Cap 0 ; 0 0 1];

% Ji = inv(N)*M;
% Jf = inv(M)*N;
% 
% JiNorma = normalize(Ji);
% JfNorma = normalize(Jf);
% 
% JiEigen=eig(JiNorma);
% JiEigenMax=max(JiEigen);
% JiEigenMin=min(JiEigen);
% 
% JfEigen=eig(JfNorma);
% JfEigenMax=max(JfEigen);
% JfEigenMin=min(JfEigen);
% 
% %ICN = Inverse condition number
% 
% JiICN = JiEigenMax / JiEigenMin;
% JfICN = JfEigenMax / JfEigenMin;

if det(M)==0 | det(N)==0
    %fprintf("here");
    Micn=0;
    Nicn=0;
else
    Micn=1/cond(M);
    Nicn=1/cond(N);
end

ICN = [Micn ; Nicn];
end
