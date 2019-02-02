close all;clear all;clc;
%n = 10;L=1;H=1;K=.5;s=10; % n must be even number
n = input('Number of Divisons= '); %must be even
L = input('Length= ');
H = input('Hieght= ');
K = input('K= ');
s = input('Source= ');
error = input('error= ');
deltaX = L/(n-1);
deltaY = H/(n-1);
leftB = input('Left Boundary= ');
%rightB = input('right');
upB = input('Upper Boundary= ');
downB = input('Down Boundary= ');
T = zeros(n*n,1);
A = zeros(n*n,n*n);
size(A);
j = 1;
for i=1:n*n %row
    %for j = 1:n*n %column 
    if i==1
        A(i,i+1) = -(K*deltaY)/deltaX;
        A(i,i+n) = -(K*deltaX)/deltaY;
        A(i,i) = -(A(i,i+1)+ A(i,i+n)+(-2*(K*deltaY)/deltaX)+(-2*(K*deltaX)/deltaY));
    elseif i~=n*n
        if i<n
            A(i,i-1) =  -(K*deltaY)/deltaX; %wets
            A(i,i+1) = A(i,i-1); % east
            A(i,i+n) = -(K*deltaX)/deltaY; % south
            A(i,i) =  -(A(i,i-1)+A(i,i+1)+A(i,i+n)+(-2*(K*deltaX)/deltaY));
        elseif i==n
            A(i,i+n) = -(K*deltaX)/deltaY;
            A(i,i-1) =  -(K*deltaY)/deltaX;
            A(i,i) = -(A(i,i-1)+A(i,i+n)+(-(2*K*deltaX)/deltaY));
        elseif i>n
            if mod(i,n) == 0 %The last column
                A(i,i+n) = -(K*deltaX)/deltaY;
                A(i,i-1) =  -(K*deltaY)/deltaX;
                A(i,i-n) = A(i,i+n);
                A(i,i) = -(A(i,i+n)+A(i,i-1)+A(i,i-n));
                j = j + 1;
                j;
            elseif i==(j*n + 1) %The first column
                if i==(n*n - n + 1)
                    A(i,i-n) = -(K*deltaX)/deltaY;
                    A(i,i+1) =  -(K*deltaY)/deltaX;
                    A(i,i) = -(A(i,i-n)+A(i,i+1)+-(2*((K*deltaY)/deltaX))+(-2*(K*deltaX)/deltaY));
                else    
                    A(i,i+n) = -(K*deltaX)/deltaY;
                    A(i,i+1) =  -(K*deltaY)/deltaX;
                    A(i,i-n) = A(i,i+n);
                    %a(i,i+1) = a(i,i-1);
                    A(i,i)  = -(A(i,i+n)+A(i,i+1)+A(i,i-n)+(-2*(K*deltaY)/deltaX));
                    'insideif';
                    i;
                end       
            elseif i > (n*n)-n+1
                A(i,i-n) = -(K*deltaX)/deltaY;
                A(i,i-1) =  -(K*deltaY)/deltaX;
                A(i,i+1) = A(i,i-1);
                A(i,i) = -(A(i,i-n)+A(i,i-1)+A(i,i+1)+(-2*((K*deltaX)/deltaY)));
                
            else
                A(i,i+n) = -(K*deltaX)/deltaY;
                A(i,i-1) =  -(K*deltaY)/deltaX;
                A(i,i-n) = A(i,i+n);
                A(i,i+1) = A(i,i-1);
                A(i,i)  = -(A(i,i+n)+A(i,i+1)+A(i,i-n)+A(i,i+1));
                
            end
        end
    else
        A(i,i-n) = -(K*deltaX)/deltaY;
        A(i,i-1) =  -(K*deltaY)/deltaX;
        A(i,i) = -(A(i,i-n)+A(i,i-1)-(2*((K*deltaX)/deltaY))-(2*((K*deltaY)/deltaX)));
    end
    
end
b = zeros(n*n,1);
S = s*deltaX*deltaY;
b(1,1) = ((2*K*deltaY/deltaX)*leftB) + ((2*K*deltaX/deltaY)*upB) +S;
b(2:n) = ((2*K*deltaX/deltaY)*upB) +S;
for c = (n+1) : (n*n - n)
    b(c) = S;
end
c = 1;
while c<=(n-1)
    b(c*n+1) = S + (2*(K*deltaY)/deltaX)*leftB;
    c = c+1;
end
for c = (n*n - (n/2)) : (n*n)
    b(c) = S + ((2*K*deltaX)/deltaY) * downB;
end
b(n*n-(n-1)) =  S + ((2*K*deltaX)/deltaY) * downB +((2*K*deltaY)/deltaX) * leftB ;
x = linspace(0,L,n);
y = linspace(0,H,n);
[X,Y] = meshgrid(x,y);
%T = A\b;
C = A;
B = b./diag(A);
for i = 1:n*n
    C(i,i) = 0;
end
er=1;
Tinit = zeros(n*n,1);
f = zeros(n,n);
it = 0;
err = 1;
while err > error
    T = ((-C*Tinit)./diag(A)) + B;
    %error = max
    %j
    R = sum(abs(T-Tinit));
    F = sum(abs(diag(A).*T));
    err = R/F
    t = reshape(T,n,n);
    %t = t';
    %er = max(abs((A*T)- b))
    Tinit = T;
    %T';
    it = it +1;
    t(:,end) = t(:,end-1);
    errors(it) = err;
end
figure(1);
title('Tempartue Distribuation')
t = reshape(T,n,n);
surface(Y,X,t);
figure(2);
plot(1:it,errors);
semilogy(errors)
xlim = ([-1,it]);
xlabel('iterations')
ylabel('error')
title('Iteartions VS Error')

            
            
            
            



