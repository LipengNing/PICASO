function ComputePICASO(Input,Mask,OutputName)
%Input: Input-- dMRI file in nrrd format
%       Mask-- brain mask in nrrd format
%       Output-- the name of outputfiles
% The code can be suitably modified in the case of multiple data files and
% different file format. For example, the following code is an example for
% the case of multiple input files:
%Input1=[Input '-b1000-unring.nrrd'];
%Input2=[Input '-b3000-unring.nrrd'];
%[S1, u1, b1, voxel, Sm, s0] = nhdr_diff_multiB(Input1);
%[S2, u2, b2, voxel, Sm, s0] = nhdr_diff_multiB(Input2);
%[nx,ny,nz,N]=size(S1);
%S1=reshape(S1,nx,ny,nz,N);
%S2=reshape(S2,nx,ny,nz,N);
%b=[b1(1:30);b2(1:30)];
%u=[u1(1:N,:);u2(1:30,:)];

[S, u, b, ~, ~, ~] = nhdr_diff_multiB(Input);
[nx,ny,nz,N]=size(S);


mask=nrrdLoad(Mask);
b=b/1000;
u=u(1:N,:);

Disturb_per=zeros(nx,ny,nz);
Disturb_par=zeros(nx,ny,nz);
Diff_per=zeros(nx,ny,nz);
Diff_par=zeros(nx,ny,nz);
S=reshape(S,nx*ny*nz,N);
parfor i=1:nx*ny*nz
    if(mask(i))
        s=double(S(i,:))';
        if(mean(s)>0.01)
            s(s>1)=1;
            s(s<0)=0;
            [~,~,coef]=Model_PICASO_fixV(s,u,b);
            Disturb_per(i)=coef(2)*coef(4);
            Disturb_par(i)=coef(1)*coef(3);
            Diff_per(i)=coef(4);
            Diff_par(i)=coef(3);
        end
    end
end
output=[OutputName '.mat'];
save(output,'Disturb_per','Disturb_par','Diff_per','Diff_par');