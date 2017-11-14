function U=utiles(Theta)
NB = narrow(Theta);
U = [];
for k=1:length(NB)
   v = voisins(NB(k,1),NB(k,2));
   for l=1:length(v)
      if(Theta(v(l,1),v(l,2))==1)
          U = [U; v(l,1) v(l,2)];
      end
   end
end
U = unique(U,'rows');
end