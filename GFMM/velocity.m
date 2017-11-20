function c = velocity(theta,N)
    c = ones(N+3);
    c(find(theta < 0)) = -1;
   % c(floor(N/2)-5:floor(N/2)+5,floor(N/2)-5:floor(N/2)+5) = -1;
   % c(:,1:floor(N/2)) = -1;
   % c = -ones(N+3);
    s = size(c);
    for i = 2:s(1)-1
        for j = 2:s(2)-1
            if (c(i,j)*c(i+1,j) < 0 && abs(c(i,j)) <= abs(c(i+1,j))) || (c(i,j)*c(i-1,j) < 0 && abs(c(i,j)) <= abs(c(i-1,j))) || (c(i,j)*c(i,j+1) < 0 && abs(c(i,j)) <= abs(c(i,j+1))) ||  (c(i,j)*c(i,j-1) < 0 && abs(c(i,j)) <= abs(c(i,j-1)))   
                c(i,j) = 0;
            end  
        end
    end
    
end

