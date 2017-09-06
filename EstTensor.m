function D=EstTensor(q,E)
%compute the DTI tensor using q value and normalized dMRI signals E
E=E(:);
E(E>1)=1;
E=E(E>1e-4);
q=q(E>1e-4,:);
s=-log(E);
A=[q(:,1).^2  q(:,2).^2  q(:,3).^2  2*q(:,1).*q(:,2)  2*q(:,1).*q(:,3)  2*q(:,2).*q(:,3)];
d=pinv(A'*A)*A'*s;d=real(d);
D=[d(1) d(4) d(5);
   d(4) d(2) d(6);
   d(5) d(6) d(3)];
[U,e]=eig(D);
e=diag(e);
e(e<5e-6)=5e-6;
D=U*diag(e)*U';
end

