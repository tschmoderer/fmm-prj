function V = voisins(i,j)


V = [];
for k=-1:1:1
    for l=-1:1:1
        if i+k>0 && j+l>0 && (i+k~=i || j+l~=j)
            V = [V; i+k j+l]
        end
    end
end
end
