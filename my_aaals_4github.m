function [A B ZZ_a ZZ_b] = my_aaals (X, r, scaling) 
% AAALS: ALS algorithm for AA 

   siz = size (X)
 if( length(siz) == 3)
%  supposed to be [time, lat lon];
   n = siz(1); la = siz(2); lo = siz(3);
   X = reshape(X, n, la*lo);
  end
       
if (nargin == 1)
  r = 4; scaling = 0;
elseif (nargin == 2)
  scaling = 0;
end
%=====================================================

[n,p] = size(X); %n - cases; p - variables
disp('Number of archetypes'), r
%pr = p*r; 
%nz = pr - p; %required number of zeros

% Mean centered data matrix with unit length columns

 stdX = std(X); 
 if (scaling == 0)
  Z = X;
 elseif (scaling == 1)
   Z=sqrt(1/(n-1))*(X-ones(n,1)*mean(X))./(ones(n,1)*stdX); scaling=1;
%  Z=sqrt(1/(n-1)) * X ./ (ones(n,1)*stdX); scaling=1;

%  Z=sqrt(1/(n-1))*(X-ones(n,1)*mean(X)) / mean(stdX); scaling=1; 
%   Z=sqrt(1/(n-1))*X / mean(stdX); scaling=1;
  end
%   Z = X;

% Initializing A and B
% -------------------- of NT
 A = rand(n,r); A = diag(1./sum(A,2))*A;
 B = rand(n,r); B = B*diag(1./sum(B,1));
 Y = Z' * B;

% MCMC Intial AAs
%----------------
  mm = 400;  Y = iniAA (Z, r, mm); Y = Y';
  A = Z*Y*inv(Y'*Y); A = diag(1./sum(A,2))*A;
  B = A*inv(A'*A)  ; B = B*diag(1./sum(B,1));
  II1 = A<0; A(II1) = rand;  A = diag(1./sum(A,2))*A;
  II1 = B<0; B(II1) = rand;  B = B*diag(1./sum(B,1));

% Objective function
fX = norm(X,'fro');
f0 = norm((eye(n)-A*B')*X,'fro')/fX;
disp('Initial value of objective function:' ), f0
Frec = f0;

continue_ok = 'Y'; max_iter = 0;

tic % Start counting CPU time
%timecheck = cputime;  

while continue_ok == 'Y' 
    max_iter = max_iter +1;

% find A row by row
W = Y;
WW = W'*W;

tspa =20.; tspan = [0; tspa];
tspan = (0:.1:tspa); 

nve = ones (1, r); pve = ones (1, n);
% can add this property to options to force >= 0 solution:  'NonNegative',nve
           options=odeset('RelTol',1e-10,'AbsTol',1e-10,'Jacobian',@FJac);
% options = [];

for i = 1:n
    v = Z(i,:)';
    u0 = A(i,:)';
%    tspan = [0; tspa];
    [t,u] = ode15s(@myfun, tspan, u0, options, W, WW, v);
%     [t,u] = ode15s(@myfun, tspan, u0, @myJ, W, WW, v);

%     disp ('length t'), t', length(t)
%     plot (u(:, 4)); hold on; plot (u0);
%     plot(u(:,1)); plot(u(:,2)); plot(u(:,3)), plot(u(:,4)); plot(u(5,:));
%     u0'
%     [u(1,1) u(1,2) u(1,3) u(1,4) u(1,5)]
%     [u(2,1) u(2,2) u(2,3) u(2,4) u(2,5)]
%      [u(length(t),1) u(length(t),2) u(length(t),3) u(length(t),4) u(length(t),5)]
%     return; 

    A(i,:) = u(length(t),:);
end
%A

% update Y with new A
Y = Z'*A/(A'*A);

% find B column by column
W = Z';
WW = W'*W;
for j = 1:r
    v = Y(:,j);
    u0 = B(:,j);
%    tspan = [0; tspa];
    [t,uu] = ode15s(@myfun, tspan, u0, options, W, WW, v);
%     [t,u] = ode15s(@myfun, tspan, u0, @myJ, W, WW, v);
    B(:,j) = uu(length(t),:);

%     plot (u(:, 100)); hold on; 
%     plot(u(:,150)); plot(u(:,300)); plot(u(:,389)), plot(u(:,405)); plot(u(:,599));
%     [u0(100)  u0(150)  u0(300)  u0(389)  u0(405)]
%     [u(1,100) u(1,150) u(1,300) u(1,389) u(1,405)]
%     [u(2,100) u(2,150) u(2,300) u(2,389) u(2,405)]
%     return;

end
%B

% update Y with new B
Y = Z'*B;
    
    % Objective function
    f = norm(Z - A*Y','fro')/fX;
    Frec = [Frec f];
    
    % Terminate?
    diff = f0 - f ;
    if (diff <= .00001)||(diff<0)
        continue_ok = 'N';
    end
%    [A, B]
    f0 = f;
end 

     figure; 
     subplot(2,2,1)
     for ii = 1:r
       plot(u(:,ii)); hold on
     end
       hold off
     subplot(2,2,2)
     for ii=1:n
      plot(uu(:,ii)); hold on
     end

      hold off

toc % End counting CPU time

fprintf('The fit is %12.8f \n',Frec)
disp('Number of iterations:' )
length(Frec)
disp('Value of objective function at minimum:' )%, f0
fprintf('The fit is %12.8f \n',f)

max_iter

figure;
plotind = 1:max_iter+1;
 plot(plotind, Frec,'ro')

% archetypes

% I added this to (?) improve B

%   B =  A * inv(A' * A);

if (scaling == 1)
  Z = sqrt(n-1) * Z .* (ones(n,1)*stdX);
% Z = sqrt(n-1) * Z * mean(stdX);
end

 % compute the Z'*A - but normalize A first
  A_n = A;

  for i=1:r
    xox=A(:,i);
    A_n(:,i)= A(:,i)/sum(xox);
  end

ZZ_b = Z'*B; ZZ_a = Z'*A_n;
if (length(siz) == 3)
  ZZ_a = reshape (ZZ_a', r, siz(2), siz(3));
  ZZ_b = reshape (ZZ_b', r, siz(2), siz(3));
end  

end
% archetypes:

%==================================================
    function dudt = myfun(t,u,W,WW,v)
      dudt =   - (diag(u) - u*u')*(WW*u - W'*v);
    end
   
    function dfdy = FJac (t, u, W, WW, v)
     pp = length(u);
     ww1 = WW*u; v1 = W'*v;
     A1 = (diag(u) - 2*u*u')*WW; A2 = diag(ww1-v1);
     A3 = u'*(ww1 - v1)*eye(pp); A4 = u*v1';
     dfdy = - (A1 + A2 - A3 + A4);
%      dfdy = WW;
    end
 
    function d2ud2t = myJ(t,WW)
      d2ud2t = WW ;
    end
