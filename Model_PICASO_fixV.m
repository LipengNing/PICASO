function [U2,Diff,Fcoef,s_est]=Model_PICASO_fixV(s,u,b)
%input: s the normalized dMRI signal from one voxel
% u: the gradient direction vectors
% b: a vector of b-values
id=find(b>1);
q=repmat(sqrt(b(id)),1,3).*u(id,:);
%[T,~,~]=Fit_1T(s,u,b);% estimate one-tensor model to initialize gradient 
T=EstTensor(q,s(id));
[V,~,~]=svd(T);
v=V(:,1);
options = optimset('Display','off','MaxIter',5000,'MaxFunEvals',5000);  
%

X0=         [.1      .1         .5     .5]';  
param.lb =  [0;      0;         0;     0]; 
param.ub =  [1;     1;         3;    3]; 

h_fn = model_BlochTorrey(u,b,v); 
%%

[Fcoef, ~,~] = lsqnonlin(@fun,X0,param.lb,param.ub,options); 
v_perp=null(v');
V=[v v_perp];

U2=V*diag([Fcoef(1)*Fcoef(3) Fcoef(2)*Fcoef(4) Fcoef(2)*Fcoef(4)])*V';

Diff=V*diag([Fcoef(3) Fcoef(4) Fcoef(4)])*V';


s_est=h_fn(Fcoef);

%figure(1);clf;plot(h_fn(Fcoef));hold on;plot(s,'r');
function fd = fun(X) 
fd = [((h_fn(X))-double((s)))
    ];
end

end


function h_fn = model_BlochTorrey(u,b,v)
  h_fn = @(X) BlochTorrey_signal(X,u,b,v);
end


function s=BlochTorrey_signal(X,u,b,v)

v_perp=null(v');
V=[v v_perp];

U2D=V*diag([X(1) X(2) X(2)])*V';
Diff=V*diag([X(3) X(4) X(4)])*V';

s=sum((u*U2D).*u,2)+(1-sum((u*U2D).*u,2)).*exp(-b.*sum((u*Diff).*u,2));
end

