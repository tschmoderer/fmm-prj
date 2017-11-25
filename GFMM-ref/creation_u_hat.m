function [u_hatp,u_hatm]=creation_u_hat(M,N,Fm,Fp,u);

u_hatp=10000*ones(M,N);
u_hatm=10000*ones(M,N);

for(k=1:length(Fp))
u_hatp(Fp(k,1),Fp(k,2))=u(Fp(k,1),Fp(k,2));
end

for(k=1:length(Fm))
u_hatm(Fm(k,1),Fm(k,2))=u(Fm(k,1),Fm(k,2));
end
