function [u_tilde]=creation_u_tilde(M,N,Fm,Fp,c_hat,h,u_hatp,u_hatm);
u_tilde=[];
u_tilde(:,1)=Fm(:,1);
u_tilde(:,2)=Fm(:,2);

for(k=1:length(Fm))
    
    if(c_hat(Fm(k,1),Fm(k,2))<=0)
        
        u_tilde(k,3)=10000;
        
    else
        
        if(Fm(k,2)~=N)
            V1=u_hatp(Fm(k,1),Fm(k,2)+1);
        else
            V1=10000;    
        end
        if(Fm(k,2)~=1)
            V2=u_hatp(Fm(k,1),Fm(k,2)-1);
        else
            V2=10000;
        end
        if(Fm(k,1)~=M)
            V3=u_hatp(Fm(k,1)+1,Fm(k,2));
        else
            V3=10000;
        end
        if(Fm(k,1)~=1)
            V4=u_hatp(Fm(k,1)-1,Fm(k,2));
        else
            V4=10000;
        end
        
        minic=min(V1,V2);
        minil=min(V3,V4);
        
        if(minic==10000&minil~=10000)
            
            u_tilde(k,3)=minil+h/abs(c_hat(Fm(k,1),Fm(k,2)));       
            
        elseif (minic~=10000&minil==10000)
            u_tilde(k,3)=minic+h/abs(c_hat(Fm(k,1),Fm(k,2)));
            
        elseif(minic~=10000&minil~=10000)
            %disp('cas3'); 
            if(V1+V2+V3+V4<2000)
            V1
            V2
            V3
            V4
            break
            end
            discri=4.0*(minil+minic)*(minil+minic)-8.0*(minil*minil+minic*minic-(h*h)/(c_hat(Fm(k,1),Fm(k,2))*c_hat(Fm(k,1),Fm(k,2))));
            
            if(discri<0)
                u_tilde(k,3)=min(minic,minil)+h/abs(c_hat(Fm(k,1),Fm(k,2)));      
                
            else
                
                if(max((2*(minil+minic)+sqrt(discri))/4,(2*(minil+minic)-sqrt(discri))/4)>=max(minil,minic))
                    
                    u_tilde(k,3)=max((2*(minil+minic)+sqrt(discri))/4,(2*(minil+minic)-sqrt(discri))/4);
                    
                  
                    
                else
                    
                    u_tilde(k,3)=min(minic,minil)+h/abs(c_hat(Fm(k,1),Fm(k,2)));   
                    
                end
                
            end
            
        else disp('erreur Fm');
            
            
        end
    end
    
    
    
end

clear k;
u_tilde(length(Fm)+1:length(Fm)+length(Fp),1)=Fp(:,1);
u_tilde(length(Fm)+1:length(Fm)+length(Fp),2)=Fp(:,2);

for(k=1:length(Fp))
    
    if(c_hat(Fp(k,1),Fp(k,2))>=0)
        
        u_tilde(length(Fm)+k,3)=10000;
        
    else
        if(Fp(k,2)~=N)
            V1=u_hatm(Fp(k,1),Fp(k,2)+1);
        else
            V1=10000;    
        end
        if(Fp(k,2)~=1)
            V2=u_hatm(Fp(k,1),Fp(k,2)-1);
        else
            V2=10000;
        end
        if(Fp(k,1)~=M)
            V3=u_hatm(Fp(k,1)+1,Fp(k,2));
        else
            V3=10000;
        end
        if(Fp(k,1)~=1)    
            V4=u_hatm(Fp(k,1)-1,Fp(k,2));
        else
            V4=10000;
        end
        
        minic=min(V1,V2);
        minil=min(V3,V4);
        
        if(minic==10000&minil~=10000)
            
            u_tilde(length(Fm)+k,3)=minil+h/abs(c_hat(Fp(k,1),Fp(k,2)));       
            
        elseif (minic~=10000&minil==10000)
            
            u_tilde(length(Fm)+k,3)=minic+h/abs(c_hat(Fp(k,1),Fp(k,2)));
            
        elseif(minic~=10000&minil~=10000)
            
            %disp('cas6')
            if(V1+V2+V3+V4<2000)
            V1
            V2
            V3
            V4
            break
            end
            discri=4.0*(minil+minic)*(minil+minic)-8.0*(minil*minil+minic*minic-(h*h)/(c_hat(Fp(k,1),Fp(k,2))*c_hat(Fp(k,1),Fp(k,2))));
            
            if(discri<0)    
                u_tilde(length(Fm)+k,3)=min(minic,minil)+h/abs(c_hat(Fp(k,1),Fp(k,2)));   
            else 
                
                
                if(max((2*(minil+minic)+sqrt(discri))/4,(2*(minil+minic)-sqrt(discri))/4)>=max(minil,minic))
                    
                    u_tilde(length(Fm)+k,3)=max((2*(minil+minic)+sqrt(discri))/4,(2*(minil+minic)-sqrt(discri))/4);
                    
                  
                    
                else
                    
                    u_tilde(length(Fm)+k,3)=min(minic,minil)+h/abs(c_hat(Fp(k,1),Fp(k,2)));     
                    
                end
            end
            
            
        else disp('erreur Fp');
            
            
            
        end
    end
    
    
    
end