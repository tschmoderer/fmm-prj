clearvars;
clc;

N = 100;

x_min = -5;
x_max = 5;
y_min = -5;
y_max = 5;

dx = (x_max-x_min)/N;
dy = (y_max-y_min)/N;

u     = 1000*ones(N,N);
Theta = -1*ones(N,N);
Z = zeros(N,N);

c = ones(N,N);

t = 0;

%initial front
I = [];

for i=1:N
    for j=1:N
        x = x_min+j*dx;
        y = y_min+i*dy;
        if(f(x,y)<=0)
            u(i,j) = 0;
            Z(i,j) = 2;
            Theta(i,j) = 1;
            I = [I; i j];
        end
    end
end

NB = narrow(Theta);
Ut = utiles(Theta);


for n=1:200
   NB = narrow(Theta);
   tu  = zeros(1,length(NB(:,1)));
   for k=1:length(NB)
       i1 = NB(k,1);
       i2 = NB(k,2);
       
       %voisins
       u_i = u(i1,i2);
       if i1==1; u_N = 0; else u_N = u(i1-1,i2); end;%nord
       if i1==N; u_S = 0; else u_S = u(i1+1,i2); end;%sud
       if i2==1; u_O = 0; else u_O = u(i1,i2-1); end;%ouest
       if i2==N; u_E = 0; else u_E = u(i1,i2+1); end;%est
       
       %cas impossible
       if u_i<u_N && u_i<u_S && u_i<u_O && u_i<u_E
           error('Imposible de résoudre');
       else
           if u_i<u_N && u_i<u_S %le premier max vaut 0
               if u_i-u_O>u_i-u_E
                   tu(k) = u_O+dx/c(i1,i2);
               else
                   tu(k) = u_E+dx/c(i1,i2);
               end
           elseif u_i<u_O && u_i<u_E %le deuxieme max vaut 0
               if u_i-u_N>u_i-u_S
                   tu(k) = u_N+dx/c(i1,i2);
               else
                   tu(k) = u_S+dx/c(i1,i2);
               end
           else %aucun max ne vaut 0
               if u_i-u_N>u_i-u_S
                   if u_i-u_O>u_i-u_E
                       D    = 4*(u_N+u_O)-8*(u_N^2+u_O^2-(dx^2)/(c(i1,i2)^2));
                       tu(k) = real(0.5*(u_N+u_O)+0.25*sqrt(D));
                   else
                       D    = 4*(u_N+u_E)-8*(u_N^2+u_E^2-(dx^2)/(c(i1,i2)^2));
                       tu(k) = real(0.5*(u_N+u_E)+0.25*sqrt(D));
                   end
               else
                   if u_i-u_O>u_i-u_E
                       D    = 4*(u_S+u_O)-8*(u_S^2+u_O^2-(dx^2)/(c(i1,i2)^2));
                       tu(k) = real(0.5*(u_S+u_O)+0.25*sqrt(D));
                   else
                       D    = 4*(u_S+u_E)-8*(u_S^2+u_E^2-(dx^2)/(c(i1,i2)^2));
                       tu(k) = real(0.5*(u_S+u_E)+0.25*sqrt(D));
                   end
               end
           end
       end     
   end
   
   m  = min(tu(tu~=0));
   na = find(tu==m);
   NA = NB(na,:);
   
   for k=1:length(NA)
       Theta(NA(k,1),NA(k,2)) = 1;
   end
   
   Ut = utiles(Theta);
   for k=1:length(NA)
       if not(isempty(find(Ut(:,1)==NA(k,1) & Ut(:,2)==NA(k,2))))
           u(NA(k,1),NA(k,2)) = t;
       end
   end
   
    %affichage
    [XX YY] = meshgrid(x_min:dx:x_max-dx,y_min:dx:y_max-dx);
    surf(XX,YY,u);
    xlabel('x')
    ylabel('y')
    zlabel('T')
    title(['time : ' num2str(t)]);
    view(0,90)
    drawnow
    pause(0.5)
    
    t = t+m
end

