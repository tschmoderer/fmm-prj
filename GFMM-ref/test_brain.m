clear all
close all
clc

nu=imread('left_sag_081.tif');
nu=nu(:,:,1);
[M N]=size(nu);
nu=double(nu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%construction of the field theta%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h=1;
deltat=0.5;
n_lig=[1:M];
n_col=[1:N];
x=h*n_lig';
y=h*n_col;
x0=x*ones(1,N);
y0=ones(M,1)*y;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %       Creation of theta      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fp=[];
Fm=[];



  
 
 for(k=1:M)
    for(l=1:N)
        if((k-70)^2/55^2+(l-95)^2/75^2<=1 & (k-68)^2+(l-80)^2>=5^2)
            theta(k,l)=1.0;
            char_theta(k,l)=1.0;
        else
            theta(k,l)=-1.0;
            char_theta(k,l)=0.0;

        end

    end
 end

figure;imagesc(nu);colormap(gray);
hold on;
contour(1:N,1:M,theta,...
    [0 0],'r','LineWidth',2);
axis off;

Itermax=1500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Construction of F- denoted by Fm%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ind_lig,ind_col]=find(theta==1);
for(k=1:length(ind_lig))
    
    if(ind_lig(k)~=M&(theta(ind_lig(k)+1,ind_col(k))==-1))
        Fm(end+1,:)=[ind_lig(k)+1,ind_col(k)];
        
    end
    if(ind_lig(k)~=1&(theta(ind_lig(k)-1,ind_col(k))==-1))
        Fm(end+1,:)=[ind_lig(k)-1,ind_col(k)];
        
    end
    if(ind_col(k)~=N&(theta(ind_lig(k),ind_col(k)+1)==-1))
        Fm(end+1,:)=[ind_lig(k),ind_col(k)+1];
        
    end
    if(ind_col(k)~=1&(theta(ind_lig(k),ind_col(k)-1)==-1))
        Fm(end+1,:)=[ind_lig(k),ind_col(k)-1];
        
    end
end
Fm=unique(Fm,'rows');
clear ind_lig ind_col


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Construction of F+ denoted by Fp%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ind_lig,ind_col]=find(theta==-1);

for(k=1:length(ind_lig))
    
    if(ind_lig(k)~=M&(theta(ind_lig(k)+1,ind_col(k))==1))
        Fp(end+1,:)=[ind_lig(k)+1,ind_col(k)];
        
    end
    if(ind_lig(k)~=1&(theta(ind_lig(k)-1,ind_col(k))==1))
        Fp(end+1,:)=[ind_lig(k)-1,ind_col(k)];
        
    end
    if(ind_col(k)~=N&(theta(ind_lig(k),ind_col(k)+1)==1))
        Fp(end+1,:)=[ind_lig(k),ind_col(k)+1];
        
    end
    if(ind_col(k)~=1&(theta(ind_lig(k),ind_col(k)-1)==1))
        Fp(end+1,:)=[ind_lig(k),ind_col(k)-1];
        
    end
end
Fp=unique(Fp,'rows');
clear ind_lig ind_col

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Initialization of the fronts%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u=zeros(M,N);
tn_1=0.0;
for(l=1:Itermax)
    
    
    c2=sum(sum((1-char_theta).*nu))/sum(sum((1-char_theta)));
    c1=sum(sum((char_theta).*nu))/sum(sum((char_theta)));
   
    
    for(i=1:M)
        for(j=1:N)
        c(i,j)=(((nu(i,j)-c2)^2-(nu(i,j)-c1)^2));  
        end
    end
    

    
    
    c_hat=c;
    
    for(i=2:M-1)
        for(j=2:N-1)
            
            if((c(i-1,j)*c(i,j)<0 & abs(c(i,j)) <=abs(c(i-1,j)))|(c(i+1,j)*c(i,j)<0 & abs(c(i,j)) <=abs(c(i+1,j))) |...
                    (c(i,j-1)*c(i,j)<0 & abs(c(i,j)) <=abs(c(i,j-1)))|(c(i,j+1)*c(i,j)<0 & abs(c(i,j)) <=abs(c(i,j+1))))
                
                c_hat(i,j)=0.0;
                
            else
                c_hat(i,j)=c(i,j);
                
            end
        end
    end
    
    
   
    
    [u_hatp,u_hatm]=creation_u_hat(M,N,Fm,Fp,u);
    
    u_tilde=creation_u_tilde(M,N,Fm,Fp,c_hat,h,u_hatp,u_hatm);
    tn_tilde=min(u_tilde(:,3));
    tn_hat=min(tn_tilde,tn_1+deltat);
    tn=max(tn_1,tn_hat);
    if(tn==tn_1+deltat&tn<tn_tilde)
       tn_1=tn; 
    else
    NAm=[];
    NAp=[];
    p=find(u_tilde(1:length(Fm),3)==tn_tilde);
    q=find(u_tilde(length(Fm)+1:length(Fm)+length(Fp),3)==tn_tilde);
    
    if(length(p)>=1)
        for(k=1:length(p))
            NAm(k,1)=u_tilde(p(k),1);
            NAm(k,2)=u_tilde(p(k),2);
            
        end
    end
    clear k
    
    if(length(q)>=1)
        for(k=1:length(q))
            NAp(k,1)=u_tilde(length(Fm)+q(k),1);
            NAp(k,2)=u_tilde(length(Fm)+q(k),2);
            
        end
    end
    
    clear p;
    clear q;
    
    
    theta=creation_theta(NAm,NAp,M,N,theta);
    for(ko=1:M)
        for(lu=1:N)
            if(theta(ko,lu)==1.0)
                char_theta(ko,lu)=1.0;
            else
                char_theta(ko,lu)=0.0;
                
            end
            
        end
    end
    if(rem(l,100)==0)
        figure;imagesc(nu),colormap(gray);
        hold on; contour(1:N,1:M,theta,[0,0],'r','LineWidth',2);
        axis off
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%Construction of Fnewm and Fnewp%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fnewm=[];
    Fnewp=[];
    [ind_lig,ind_col]=find(theta==1);
    
    for(k=1:length(ind_lig))
        
        if(ind_lig(k)~=M&(theta(ind_lig(k)+1,ind_col(k))==-1))
            Fnewm(end+1,:)=[ind_lig(k)+1,ind_col(k)];
            
        end
        if(ind_lig(k)~=1&(theta(ind_lig(k)-1,ind_col(k))==-1))
            Fnewm(end+1,:)=[ind_lig(k)-1,ind_col(k)];
            
        end
        if(ind_col(k)~=N&(theta(ind_lig(k),ind_col(k)+1)==-1))
            Fnewm(end+1,:)=[ind_lig(k),ind_col(k)+1];
            
        end
        if(ind_col(k)~=1&(theta(ind_lig(k),ind_col(k)-1)==-1))
            Fnewm(end+1,:)=[ind_lig(k),ind_col(k)-1];
            
        end
    end
    Fnewm=unique(Fnewm,'rows');
    
    clear ind_lig ind_col
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%Construction of F+ denoted by Fp%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [ind_lig,ind_col]=find(theta==-1);
    
    for(k=1:length(ind_lig))
        
        if(ind_lig(k)~=M&(theta(ind_lig(k)+1,ind_col(k))==1))
            Fnewp(end+1,:)=[ind_lig(k)+1,ind_col(k)];
            
        end
        if(ind_lig(k)~=1&(theta(ind_lig(k)-1,ind_col(k))==1))
            Fnewp(end+1,:)=[ind_lig(k)-1,ind_col(k)];
            
        end
        if(ind_col(k)~=N&(theta(ind_lig(k),ind_col(k)+1)==1))
            Fnewp(end+1,:)=[ind_lig(k),ind_col(k)+1];
            
        end
        if(ind_col(k)~=1&(theta(ind_lig(k),ind_col(k)-1)==1))
            Fnewp(end+1,:)=[ind_lig(k),ind_col(k)-1];
            
        end
    end
    Fnewp=unique(Fnewp,'rows');
    clear ind_lig ind_col;
    
    if(size(NAm,1)>=1)
        for(k=1:size(NAm,1))
            u(NAm(k,1),NAm(k,2))=tn;
        end
    end
    
    if(size(NAp,1)>=1)
        for(k=1:size(NAp,1))
            u(NAp(k,1),NAp(k,2))=tn;
        end
    end
    
    NA=[NAm;NAp];
    F=[Fm;Fp];
    VNA=[];
    
    for(k=1:size(NA,1))
        if(NA(k,1)~=1)
            VNA(end+1,:)=[NA(k,1)-1 NA(k,2)];
        end
        if(NA(k,1)~=M)
            VNA(end+1,:)=[NA(k,1)+1, NA(k,2)];   
        end
        if(NA(k,2)~=1)
            VNA(end+1,:)=[NA(k,1), NA(k,2)-1];   
        end
        if(NA(k,2)~=N)
            VNA(end+1,:)=[NA(k,1), NA(k,2)+1];   
        end
        
    end
    VNA=unique(VNA,'rows');
    
    [C,I]=setdiff(VNA,F,'rows');
    
    if(length(I)>=1)
        for(k=1:length(I))
            u(VNA(I(k),1),VNA(I(k),2))=tn;
        end
    end
    clear Fm Fp;
    
    Fm=Fnewm;
    Fp=Fnewp;
    tn_1=tn;
    clear Fnewm Fnewp clear VNA NAp NAm F;
    clear u_tilde;
    clear c_hat;
    clear u_hatp u_hatm;
end
    
    
end




