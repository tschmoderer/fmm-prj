function theta=creation_theta(NAm,NAp,M,N,theta);

if(size(NAm,1)>=1)
    for(k=1:size(NAm,1))
        theta(NAm(k,1),NAm(k,2))=1;
    end
end


if(size(NAp,1)>=1)
    for(k=1:size(NAp,1))
        theta(NAp(k,1),NAp(k,2))=-1;
    end
end