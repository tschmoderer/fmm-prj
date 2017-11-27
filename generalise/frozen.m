function F=frozen(Theta)

[lig col] = find(Theta==1);
n = size(Theta);
F(:,1) = lig;
F(:,2) = col;

end